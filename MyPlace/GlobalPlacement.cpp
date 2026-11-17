// ============================================================================
// 全局布局实现（任务5《全局布局设计》）
//
// 目标函数：f(x,y) = Σ WL(e;x,y) + λ·D(x,y)
//   WL：WA-LSE 平滑线长（公式(3)(4)(5)(6)，数值稳定实现）
//   D ：静电类比密度惩罚（公式(7)~(11)、(20)(21)）
// 求解：Nesterov 加速梯度法 + 密度溢出率驱动的停止判定
// ============================================================================

#include "GlobalPlacement.h"
#include "Helper.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
constexpr double kPi = 3.14159265358979323846;

// 密度梯度符号：文档式(9) 为 ∇D = -q·ξ（ξ 为电场）。
// 若实测发现密度溢出率随迭代上升，把该常数取反即可。
constexpr double kDensityGradientSign = -1.0;
}

// ============================================================================
// 基础工具
// ============================================================================

bool GlobalPlacement::IsValidFinite(double val)
{
    return std::isfinite(val);
}

// 软边界惩罚梯度：单元越界时把它拉回芯片内部（安全保护）
double GlobalPlacement::BoundaryPenaltyGrad(double pos, double ll, double ur, double k) const
{
    if (pos < ll)
        return k * (ll - pos);
    if (pos > ur)
        return -k * (pos - ur);
    return 0.0;
}

double GlobalPlacement::BinCoordX(double x) const
{
    if (!(binWidth > 0.0))
        return 0.0;
    // 密度网格的第 i 个采样点位于第 i 个 bin 的中心，因此连续坐标为 (x-x0)/binWidth
    return (x - static_cast<double>(chipLL.x)) / binWidth;
}

double GlobalPlacement::BinCoordY(double y) const
{
    if (!(binHeight > 0.0))
        return 0.0;
    return (y - static_cast<double>(chipLL.y)) / binHeight;
}

double GlobalPlacement::SampleField(const std::vector<double>& field, double bx, double by) const
{
    return bilinearSample(field, binM, bx, by);
}

// 公式(6)：γ = 8·w_b·10^(20/9·(τ-0.1)-1)，w_b 为网格宽度，τ 为密度溢出率。
// 密度溢出大时 γ 很大（线长模型非常平滑），随着布局均匀（τ→0）γ 迅速减小、
// 线长项逐渐恢复精度。
double GlobalPlacement::WirelengthSmoothingGamma() const
{
    const double tau = clampValue(lastTau, 0.01, 1.0);
    const double wb = std::max(binWidth, 1e-6);
    const double exponent = (20.0 / 9.0) * (tau - 0.1) - 1.0;
    const double gamma = gammaBinScale * wb * std::pow(10.0, exponent);
    if (!IsValidFinite(gamma) || gamma < 1e-3)
        return 1e-3;
    return gamma;
}

// ============================================================================
// 1. 初始化阶段
// ============================================================================

// 由 SiteRow 求出 core 区域边界；并向内侧收一个极小余量，
// 保证写回的 float 坐标经"左下角 = 中心 - 半宽"换算后不会因舍入而越界。
void GlobalPlacement::SetChipBoundary(PlaceData* db)
{
    double xmin = std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();

    for (const SiteRow& row : db->SiteRows)
    {
        xmin = std::min(xmin, static_cast<double>(row.start.x));
        xmax = std::max(xmax, static_cast<double>(row.end.x));
        ymin = std::min(ymin, static_cast<double>(row.bottom));
        ymax = std::max(ymax, static_cast<double>(row.bottom + row.height));
    }

    // 数据里没有 SiteRow 时退化为所有单元的包围盒
    if (!(xmax > xmin) || !(ymax > ymin))
    {
        for (Module* m : db->Nodes)
        {
            if (m->isFixed) continue;
            double lx, ly, rx, ry;
            getModuleBBox(m, lx, ly, rx, ry);
            xmin = std::min(xmin, lx); xmax = std::max(xmax, rx);
            ymin = std::min(ymin, ly); ymax = std::max(ymax, ry);
        }
        for (Module* m : db->Terminals)
        {
            double lx, ly, rx, ry;
            getModuleBBox(m, lx, ly, rx, ry);
            xmin = std::min(xmin, lx); xmax = std::max(xmax, rx);
            ymin = std::min(ymin, ly); ymax = std::max(ymax, ry);
        }
    }

    const double width = xmax - xmin;
    const double height = ymax - ymin;
    const double margin = 1e-4 * std::max(width, height) + 1e-3;

    chipLL.x = static_cast<float>(xmin + margin);
    chipLL.y = static_cast<float>(ymin + margin);
    chipUR.x = static_cast<float>(xmax - margin);
    chipUR.y = static_cast<float>(ymax - margin);

    chipArea = (static_cast<double>(chipUR.x) - static_cast<double>(chipLL.x))
             * (static_cast<double>(chipUR.y) - static_cast<double>(chipLL.y));
}

// 1-1 填充单元初始化
//   filler 总面积 = targetDensity·chipArea - 标准单元面积 - 宏单元面积；
//   单个 filler 面积取"中间 80% 单元"的平均面积（正方形 filler）；
//   随机撒在 core 内，并与原始节点合并成 NodesAndFillers。
void GlobalPlacement::FillerInit(PlaceData* db)
{
    fillerNodes.clear();
    nodesAndFillers.clear();

    stdCellArea = 0.0;
    macroArea = 0.0;
    std::vector<double> cellAreas;
    cellAreas.reserve(db->Nodes.size());

    for (Module* m : db->Nodes)
    {
        if (m->isFixed) continue;
        const double area = std::max(0.0, static_cast<double>(m->area));
        if (!IsValidFinite(area) || area <= 0.0) continue;
        if (m->isMacro)
        {
            macroArea += area;
        }
        else
        {
            stdCellArea += area;
            cellAreas.push_back(area);
        }
    }

    // filler 只填充"可放置"的空白：先在芯片面积中扣除固定终端与不可放置区域，
    // 再按目标密度补齐标准单元/宏单元剩下的部分（文档公式在终端面积为 0 时与之等价）。
    const double placeableArea = std::max(0.0, chipArea - terminalCoverage - darkCoverage);
    const double fillerTotalArea = targetBinDensity * placeableArea - stdCellArea - macroArea;

    if (fillerTotalArea > 1e-9 && !cellAreas.empty())
    {
        // 取面积位于中间 80% 的单元求平均，避免个别超大/超小单元拉偏 filler 尺寸
        std::sort(cellAreas.begin(), cellAreas.end());
        std::size_t lo = cellAreas.size() / 10;
        std::size_t hi = cellAreas.size() - lo;
        if (hi <= lo) { lo = 0; hi = cellAreas.size(); }
        double sum = 0.0;
        for (std::size_t i = lo; i < hi; ++i) sum += cellAreas[i];
        const double avgArea = sum / static_cast<double>(hi - lo);

        double oneFillerArea = avgArea;
        long long fillerCount = 0;
        if (oneFillerArea > 1e-12)
            fillerCount = static_cast<long long>(fillerTotalArea / oneFillerArea);

        // 上限保护：超大电路上 filler 数量过多会显著拖慢迭代
        const long long kMaxFillerCount = 600000;
        if (fillerCount > kMaxFillerCount)
        {
            fillerCount = kMaxFillerCount;
            oneFillerArea = fillerTotalArea / static_cast<double>(fillerCount);
        }

        if (fillerCount > 0)
        {
            const double side = std::sqrt(oneFillerArea);
            const double usableW = static_cast<double>(chipUR.x) - static_cast<double>(chipLL.x) - side;
            const double usableH = static_cast<double>(chipUR.y) - static_cast<double>(chipLL.y) - side;

            std::mt19937 rng(42);  // 固定种子，保证结果可复现
            std::uniform_real_distribution<double> distX(0.0, std::max(usableW, 0.0));
            std::uniform_real_distribution<double> distY(0.0, std::max(usableH, 0.0));

            fillerNodes.reserve(static_cast<std::size_t>(fillerCount));
            for (long long i = 0; i < fillerCount; ++i)
            {
                Module* filler = new Module();
                filler->name = "filler";
                filler->isFiller = true;
                filler->isFixed = false;
                filler->isMacro = false;
                filler->width = static_cast<float>(side);
                filler->height = static_cast<float>(side);
                filler->area = static_cast<float>(oneFillerArea);
                filler->center.x = static_cast<float>(static_cast<double>(chipLL.x) + side * 0.5 + distX(rng));
                filler->center.y = static_cast<float>(static_cast<double>(chipLL.y) + side * 0.5 + distY(rng));
                fillerNodes.push_back(filler);
            }
        }
    }

    nodesAndFillers.reserve(db->Nodes.size() + fillerNodes.size());
    for (Module* m : db->Nodes) nodesAndFillers.push_back(m);
    for (Module* f : fillerNodes) nodesAndFillers.push_back(f);
}

// 1-2 密度 Bin 网格初始化
//   网格数由"平均单元面积 + 目标密度"决定（每个 bin 约容纳 8×8 个平均单元，
//   并限制在 [16,64] 内以保证 DCT 求电场的开销可控）；
//   终端/不可放置区域的密度按文档乘 targetDensity。
void GlobalPlacement::BinInit(PlaceData* db)
{
    double totalArea = 0.0;
    std::size_t count = 0;
    for (Module* m : db->Nodes)
    {
        if (m->isFixed) continue;
        const double area = static_cast<double>(m->area);
        if (area > 0.0 && IsValidFinite(area)) { totalArea += area; ++count; }
    }
    const double avgCellArea = (count > 0) ? (totalArea / static_cast<double>(count)) : 1.0;

    const double ideal = std::sqrt(std::max(1.0, chipArea * targetBinDensity / (64.0 * avgCellArea)));
    binM = clampIndex(static_cast<int>(std::lround(ideal)), 16, 64);

    binWidth = (static_cast<double>(chipUR.x) - static_cast<double>(chipLL.x)) / binM;
    binHeight = (static_cast<double>(chipUR.y) - static_cast<double>(chipLL.y)) / binM;
    binArea = binWidth * binHeight;

    bins.assign(static_cast<std::size_t>(binM), std::vector<BinGrid>(static_cast<std::size_t>(binM)));
    for (int i = 0; i < binM; ++i)
    {
        for (int j = 0; j < binM; ++j)
        {
            BinGrid& b = bins[i][j];
            b.x0 = static_cast<double>(chipLL.x) + i * binWidth;
            b.y0 = static_cast<double>(chipLL.y) + j * binHeight;
            b.x1 = b.x0 + binWidth;
            b.y1 = b.y0 + binHeight;
            b.area = binArea;
            b.nodeDensity = 0.0;
            b.fillerDensity = 0.0;
            b.terminalDensity = 0.0;
            b.darkDensity = 0.0;
        }
    }

    // 固定终端面积（文档：targetDensity × 重叠面积）
    for (Module* term : db->Terminals)
    {
        double lx, ly, rx, ry;
        getModuleBBox(term, lx, ly, rx, ry);
        const int i0 = clampIndex(static_cast<int>(std::floor((lx - chipLL.x) / binWidth)), 0, binM - 1);
        const int i1 = clampIndex(static_cast<int>(std::floor((rx - chipLL.x) / binWidth)), 0, binM - 1);
        const int j0 = clampIndex(static_cast<int>(std::floor((ly - chipLL.y) / binHeight)), 0, binM - 1);
        const int j1 = clampIndex(static_cast<int>(std::floor((ry - chipLL.y) / binHeight)), 0, binM - 1);
        for (int i = i0; i <= i1; ++i)
            for (int j = j0; j <= j1; ++j)
            {
                const double ov = rectOverlapArea(lx, ly, rx, ry, bins[i][j].x0, bins[i][j].y0,
                                                  bins[i][j].x1, bins[i][j].y1);
                if (ov > 0.0) bins[i][j].terminalDensity += targetBinDensity * ov;
            }
    }

    // 不可放置区域（SiteRow 覆盖范围之外）面积：先统计每个 bin 被行覆盖的面积
    std::vector<double> rowCoverage(static_cast<std::size_t>(binM) * binM, 0.0);
    for (const SiteRow& row : db->SiteRows)
    {
        const double rx0 = static_cast<double>(row.start.x);
        const double rx1 = static_cast<double>(row.end.x);
        const double ry0 = static_cast<double>(row.bottom);
        const double ry1 = static_cast<double>(row.bottom + row.height);
        if (!(rx1 > rx0) || !(ry1 > ry0)) continue;

        const int i0 = clampIndex(static_cast<int>(std::floor((rx0 - chipLL.x) / binWidth)), 0, binM - 1);
        const int i1 = clampIndex(static_cast<int>(std::floor((rx1 - chipLL.x) / binWidth)), 0, binM - 1);
        const int j0 = clampIndex(static_cast<int>(std::floor((ry0 - chipLL.y) / binHeight)), 0, binM - 1);
        const int j1 = clampIndex(static_cast<int>(std::floor((ry1 - chipLL.y) / binHeight)), 0, binM - 1);
        for (int i = i0; i <= i1; ++i)
            for (int j = j0; j <= j1; ++j)
            {
                const double ov = rectOverlapArea(rx0, ry0, rx1, ry1, bins[i][j].x0, bins[i][j].y0,
                                                  bins[i][j].x1, bins[i][j].y1);
                rowCoverage[static_cast<std::size_t>(i) * binM + j] += ov;
            }
    }
    for (int i = 0; i < binM; ++i)
        for (int j = 0; j < binM; ++j)
        {
            BinGrid& b = bins[i][j];
            const double dark = std::max(0.0, b.area - rowCoverage[static_cast<std::size_t>(i) * binM + j]);
            b.darkDensity = targetBinDensity * dark;
        }

    // 固定终端/不可放置区域实际占用的面积（除以 targetDensity 还原原始覆盖面积），
    // 供 FillerInit 计算"还能填多少 filler"使用。
    terminalCoverage = 0.0;
    darkCoverage = 0.0;
    for (const auto& column : bins)
        for (const BinGrid& b : column)
        {
            terminalCoverage += b.terminalDensity;
            darkCoverage += b.darkDensity;
        }
    if (targetBinDensity > 0.0)
    {
        terminalCoverage /= targetBinDensity;
        darkCoverage /= targetBinDensity;
    }

    // Bookshelf 的 terminal 尺寸经常被人为放大：本数据里 63 个 terminal 的覆盖面积
    // 就达到芯片面积的 43%，若把它们全部当成密度障碍，目标密度根本无法达到
    // （单元只能挤进剩下的 57%）。此时只把 terminal 当作固定引脚参与线长，
    // 不计入密度预算（与 ePlace 处理 ISPD05 这类 terminal 的做法一致）。
    if (terminalCoverage > 0.15 * chipArea)
    {
        std::printf("注意：terminal 覆盖面积占芯片 %.1f%%，远超正常可放置比例，"
                    "本算法只把 terminal 当作固定引脚，不计入密度预算\n",
                    100.0 * terminalCoverage / chipArea);
        for (auto& column : bins)
            for (BinGrid& b : column)
                b.terminalDensity = 0.0;
        terminalCoverage = 0.0;
    }

    binDensityDirty = true;
}

// 1-3 梯度向量初始化
//   线长梯度长度 = 原始节点数；密度/总梯度长度 = 节点 + filler 数。
void GlobalPlacement::GradientVectorInitialization(PlaceData* db)
{
    const std::size_t nodeCount = db->Nodes.size();
    const std::size_t allCount = nodesAndFillers.size();

    gradWireX.assign(nodeCount, 0.0);
    gradWireY.assign(nodeCount, 0.0);
    gradDenX.assign(allCount, 0.0);
    gradDenY.assign(allCount, 0.0);
    gradTotalX.assign(allCount, 0.0);
    gradTotalY.assign(allCount, 0.0);
    vx.assign(allCount, 0.0);
    vy.assign(allCount, 0.0);

    pinPosX.clear();
    pinPosY.clear();
    pinNodes.clear();
    pinExp.clear();
}

// 建立 Net::netPins → 节点下标（db->Nodes 下标）的 CSR 索引；固定单元的 pin 记为 -1。
// 这样每次迭代计算线长梯度时不必再做 Module* → 下标的查找。
void GlobalPlacement::BuildNetPinIndex(PlaceData* db)
{
    std::unordered_map<const Module*, std::int32_t> nodeIndex;
    nodeIndex.reserve(db->Nodes.size() * 2 + 16);
    for (std::size_t i = 0; i < db->Nodes.size(); ++i)
        nodeIndex.emplace(db->Nodes[i], static_cast<std::int32_t>(i));

    netPinStart.assign(db->Nets.size() + 1, 0);
    std::size_t total = 0;
    for (std::size_t n = 0; n < db->Nets.size(); ++n)
    {
        netPinStart[n] = total;
        total += db->Nets[n]->netPins.size();
    }
    netPinStart[db->Nets.size()] = total;

    netPinIndex.assign(total, -1);
    for (std::size_t n = 0; n < db->Nets.size(); ++n)
    {
        const Net* net = db->Nets[n];
        for (std::size_t k = 0; k < net->netPins.size(); ++k)
        {
            const Pin* pin = net->netPins[k];
            if (!pin || !pin->module || pin->module->isFixed) continue;
            auto it = nodeIndex.find(pin->module);
            if (it != nodeIndex.end())
                netPinIndex[netPinStart[n] + k] = it->second;
        }
    }
}

// 1-5 惩罚系数 lambda 初始化
//   λ = Σ|线长梯度| / Σ|密度梯度|，把两项梯度拉到同一量级。
double GlobalPlacement::PenaltyFactorInit(PlaceData* db)
{
    UpdateBinDensities(db);
    lastTau = ComputeTauOverflow(db);
    GetTotalGradient(db, 1.0);

    double numerator = 0.0;
    for (double v : gradWireX) if (IsValidFinite(v)) numerator += std::fabs(v);
    for (double v : gradWireY) if (IsValidFinite(v)) numerator += std::fabs(v);

    double denominator = 0.0;
    for (double v : gradDenX) if (IsValidFinite(v)) denominator += std::fabs(v);
    for (double v : gradDenY) if (IsValidFinite(v)) denominator += std::fabs(v);

    if (!(denominator > 1e-12) || !(numerator > 1e-12))
        return 1.0;

    const double lambda = numerator / denominator;
    if (!IsValidFinite(lambda) || lambda <= 0.0)
        return 1.0;
    return lambda;
}

// ============================================================================
// 2. 梯度计算
// ============================================================================

// 1-4-1 更新 bin 内的单元密度(nodeDensity)与填充密度(fillerDensity)
void GlobalPlacement::UpdateBinDensities(PlaceData* /*db*/)
{
    if (!binDensityDirty)
        return;

    for (auto& column : bins)
    {
        for (BinGrid& b : column)
        {
            b.nodeDensity = 0.0;
            b.fillerDensity = 0.0;
        }
    }

    const double invBinW = 1.0 / binWidth;
    const double invBinH = 1.0 / binHeight;

    for (Module* mod : nodesAndFillers)
    {
        if (mod->isFixed) continue;

        double lx, ly, rx, ry;
        getModuleBBox(mod, lx, ly, rx, ry);
        if (rx <= chipLL.x || lx >= chipUR.x || ry <= chipLL.y || ly >= chipUR.y) continue;

        const int i0 = clampIndex(static_cast<int>(std::floor((lx - chipLL.x) * invBinW)), 0, binM - 1);
        const int i1 = clampIndex(static_cast<int>(std::floor((rx - chipLL.x) * invBinW)), 0, binM - 1);
        const int j0 = clampIndex(static_cast<int>(std::floor((ly - chipLL.y) * invBinH)), 0, binM - 1);
        const int j1 = clampIndex(static_cast<int>(std::floor((ry - chipLL.y) * invBinH)), 0, binM - 1);

        // 文档 1-4-2-1 的 LengthScale：单元尺寸小于网格时按尺寸比例缩放面积贡献。
        // 注意：该缩放不满足面积守恒（ePlace 原文式(3)直接用重叠面积），
        // 若一直生效会让小尺寸标准单元的密度贡献被大幅削弱、溢出率恒为 0，
        // 因此默认关闭，仅作可选项保留。
        double scale = 1.0;
        if (useLengthScaling)
        {
            const double mw = rx - lx;
            const double mh = ry - ly;
            const double sx = (mw < binWidth) ? (mw / binWidth) : 1.0;
            const double sy = (mh < binHeight) ? (mh / binHeight) : 1.0;
            scale = sx * sy;
        }

        const double factor = (mod->isFiller || mod->isMacro) ? (targetBinDensity * scale) : scale;

        for (int i = i0; i <= i1; ++i)
        {
            for (int j = j0; j <= j1; ++j)
            {
                BinGrid& b = bins[i][j];
                const double ov = rectOverlapArea(lx, ly, rx, ry, b.x0, b.y0, b.x1, b.y1);
                if (ov <= 0.0) continue;
                if (mod->isFiller)
                    b.fillerDensity += factor * ov;
                else
                    b.nodeDensity += factor * ov;
            }
        }
    }

    binDensityDirty = false;
}

// 1-4-2-3 计算密度分布的 DCT 系数 a_{j,k}（公式(20)）
void GlobalPlacement::ComputeDensitySpectrum()
{
    const int m = binM;
    const std::size_t total = static_cast<std::size_t>(m) * static_cast<std::size_t>(m);

    densityGrid.assign(total, 0.0);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j)
        {
            const BinGrid& b = bins[i][j];
            densityGrid[static_cast<std::size_t>(i) * m + j] = safeDivide(b.CellArea(), b.area);
        }

    // 查表：cos / sin(π·j·(x+0.5)/m)
    // 采样点取 bin 中心（半整数坐标）而不是整数坐标，这样 cos 基是正交的
    // （DCT-II）。若直接取 cos(πjx/m)，基在 0..m-1 上不正交，
    // 连均匀密度都会算出非零电场，从而产生把单元整体推向芯片某一侧的系统性假力。
    cosTable.resize(total);
    sinTable.resize(total);
    for (int j = 0; j < m; ++j)
    {
        for (int x = 0; x < m; ++x)
        {
            const double angle = kPi * static_cast<double>(j) * (static_cast<double>(x) + 0.5)
                               / static_cast<double>(m);
            const std::size_t idx = static_cast<std::size_t>(j) * m + x;
            cosTable[idx] = std::cos(angle);
            sinTable[idx] = std::sin(angle);
        }
    }

    // a_{j,k} = (1/m²)·Σ_x Σ_y ρ(x,y)·cos(w_j x)·cos(w_k y)
    // 利用可分离性拆成两次一维求和：O(m³) 而不是朴素实现的 O(m⁴)
    dctCoeff.assign(total, 0.0);
    std::vector<double> partial(static_cast<std::size_t>(m), 0.0);
    const double invM2 = 1.0 / (static_cast<double>(m) * static_cast<double>(m));

    for (int k = 0; k < m; ++k)
    {
        const double* cosK = &cosTable[static_cast<std::size_t>(k) * m];
        for (int i = 0; i < m; ++i)
        {
            const double* row = &densityGrid[static_cast<std::size_t>(i) * m];
            double sum = 0.0;
            for (int j = 0; j < m; ++j) sum += row[j] * cosK[j];
            partial[i] = sum;
        }
        for (int j = 0; j < m; ++j)
        {
            const double* cosJ = &cosTable[static_cast<std::size_t>(j) * m];
            double sum = 0.0;
            for (int i = 0; i < m; ++i) sum += cosJ[i] * partial[i];
            dctCoeff[static_cast<std::size_t>(j) * m + k] = sum * invM2;
        }
    }
}

// 1-4-2-5 由 DCT 系数计算电场 ξx、ξy（公式(10)）
//   ξx(i,j) = ΣΣ a[jj][kk]·w_jj/(w_jj²+w_kk²)·sin(w_jj·i)·cos(w_kk·j)
//   ξy(i,j) = ΣΣ a[jj][kk]·w_kk/(w_jj²+w_kk²)·cos(w_jj·i)·sin(w_kk·j)
// 同样利用可分离性把 O(m⁴) 降到 O(m³)。
void GlobalPlacement::ComputeElectricField()
{
    const int m = binM;
    const std::size_t total = static_cast<std::size_t>(m) * static_cast<std::size_t>(m);

    fieldX.assign(total, 0.0);
    fieldY.assign(total, 0.0);
    if (total == 0) return;

    std::vector<double> w(static_cast<std::size_t>(m), 0.0);
    for (int j = 0; j < m; ++j) w[j] = kPi * static_cast<double>(j) / static_cast<double>(m);

    // ---------- ξx ----------
    // 第一步：T[jj][y] = Σ_kk a[jj][kk]·w_jj/(w_jj²+w_kk²)·cos(w_kk·y)
    fieldWork.assign(total, 0.0);
    for (int jj = 0; jj < m; ++jj)
    {
        if (w[jj] == 0.0) continue;  // 直流分量对电势梯度无贡献
        for (int y = 0; y < m; ++y)
        {
            double sum = 0.0;
            for (int kk = 0; kk < m; ++kk)
            {
                const double denom = w[jj] * w[jj] + w[kk] * w[kk];
                if (denom < 1e-12) continue;
                sum += dctCoeff[static_cast<std::size_t>(jj) * m + kk] * (w[jj] / denom)
                     * cosTable[static_cast<std::size_t>(kk) * m + y];
            }
            fieldWork[static_cast<std::size_t>(jj) * m + y] = sum;
        }
    }
    // 第二步：ξx(i,y) = Σ_jj sin(w_jj·i)·T[jj][y]
    for (int i = 0; i < m; ++i)
    {
        for (int y = 0; y < m; ++y)
        {
            double sum = 0.0;
            for (int jj = 0; jj < m; ++jj)
                sum += sinTable[static_cast<std::size_t>(jj) * m + i] * fieldWork[static_cast<std::size_t>(jj) * m + y];
            fieldX[static_cast<std::size_t>(i) * m + y] = sum;
        }
    }

    // ---------- ξy ----------
    fieldWork.assign(total, 0.0);
    for (int kk = 0; kk < m; ++kk)
    {
        if (w[kk] == 0.0) continue;
        for (int x = 0; x < m; ++x)
        {
            double sum = 0.0;
            for (int jj = 0; jj < m; ++jj)
            {
                const double denom = w[jj] * w[jj] + w[kk] * w[kk];
                if (denom < 1e-12) continue;
                sum += dctCoeff[static_cast<std::size_t>(jj) * m + kk] * (w[kk] / denom)
                     * cosTable[static_cast<std::size_t>(jj) * m + x];
            }
            fieldWork[static_cast<std::size_t>(kk) * m + x] = sum;
        }
    }
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            double sum = 0.0;
            for (int kk = 0; kk < m; ++kk)
                sum += sinTable[static_cast<std::size_t>(kk) * m + j] * fieldWork[static_cast<std::size_t>(kk) * m + i];
            fieldY[static_cast<std::size_t>(i) * m + j] = sum;
        }
    }
}

// 1-4-2 更新密度梯度 ∇D（公式(9)：∂D/∂x_i = -q_i·ξx，∂D/∂y_i = -q_i·ξy）
void GlobalPlacement::ComputeDensityGradient(PlaceData* /*db*/)
{
    std::fill(gradDenX.begin(), gradDenX.end(), 0.0);
    std::fill(gradDenY.begin(), gradDenY.end(), 0.0);

    ComputeDensitySpectrum();
    ComputeElectricField();

    for (std::size_t idx = 0; idx < nodesAndFillers.size(); ++idx)
    {
        Module* mod = nodesAndFillers[idx];
        if (mod->isFixed) continue;

        const double cx = static_cast<double>(mod->center.x);
        const double cy = static_cast<double>(mod->center.y);
        const double bx = BinCoordX(cx);
        const double by = BinCoordY(cy);

        const double ex = SampleField(fieldX, bx, by);
        const double ey = SampleField(fieldY, bx, by);

        // 电荷量 q 取该单元在密度统计中的有效面积（宏单元/filler 与密度贡献保持一致）
        const double area = static_cast<double>(mod->area);
        const double q = (mod->isMacro || mod->isFiller) ? (targetBinDensity * area) : area;

        gradDenX[idx] = kDensityGradientSign * q * ex;
        gradDenY[idx] = kDensityGradientSign * q * ey;

        if (!IsValidFinite(gradDenX[idx])) gradDenX[idx] = 0.0;
        if (!IsValidFinite(gradDenY[idx])) gradDenY[idx] = 0.0;
    }
}

// 1-4-3 更新线长梯度 ∇WL（WA-LSE 模型，公式(3)(4)(5)）
//   WL_e(x) = γ·[ ln Σ exp(x_i/γ) + ln Σ exp(-x_i/γ) ]
//   ∂WL_e/∂x_i = exp(x_i/γ)/Σexp(x_j/γ) - exp(-x_i/γ)/Σexp(-x_j/γ)
// 用"网络 max/min 坐标"做平移保证指数不溢出（与公式(4)(5)归一化后等价）。
void GlobalPlacement::ComputeWirelengthGradient(PlaceData* db)
{
    std::fill(gradWireX.begin(), gradWireX.end(), 0.0);
    std::fill(gradWireY.begin(), gradWireY.end(), 0.0);
    wireGradMax = 0.0;

    const double gamma = WirelengthSmoothingGamma();
    const double invGamma = 1.0 / gamma;

    for (std::size_t n = 0; n < db->Nets.size(); ++n)
    {
        Net* net = db->Nets[n];
        if (!net) continue;
        const std::vector<Pin*>& pins = net->netPins;
        const std::size_t degree = pins.size();
        if (degree < 2) continue;

        pinPosX.resize(degree);
        pinPosY.resize(degree);
        pinNodes.resize(degree);
        pinExp.resize(degree * 4);

        bool valid = true;
        std::size_t movableCount = 0;
        double xmin = std::numeric_limits<double>::max();
        double xmax = -std::numeric_limits<double>::max();
        double ymin = std::numeric_limits<double>::max();
        double ymax = -std::numeric_limits<double>::max();

        for (std::size_t k = 0; k < degree; ++k)
        {
            const Pin* pin = pins[k];
            if (!pin || !pin->module) { valid = false; break; }
            const double px = static_cast<double>(pin->module->center.x) + static_cast<double>(pin->offset.x);
            const double py = static_cast<double>(pin->module->center.y) + static_cast<double>(pin->offset.y);
            if (!IsValidFinite(px) || !IsValidFinite(py)) { valid = false; break; }

            pinPosX[k] = px;
            pinPosY[k] = py;
            pinNodes[k] = netPinIndex[netPinStart[n] + k];
            if (pinNodes[k] >= 0) ++movableCount;

            xmin = std::min(xmin, px); xmax = std::max(xmax, px);
            ymin = std::min(ymin, py); ymax = std::max(ymax, py);
        }
        // 全部引脚都固定的网络不产生梯度
        if (!valid || movableCount == 0) continue;

        double sumPosX = 0.0, sumNegX = 0.0, sumPosY = 0.0, sumNegY = 0.0;
        for (std::size_t k = 0; k < degree; ++k)
        {
            const double ePosX = safeExp((pinPosX[k] - xmax) * invGamma);
            const double eNegX = safeExp((xmin - pinPosX[k]) * invGamma);
            const double ePosY = safeExp((pinPosY[k] - ymax) * invGamma);
            const double eNegY = safeExp((ymin - pinPosY[k]) * invGamma);
            pinExp[4 * k + 0] = ePosX;
            pinExp[4 * k + 1] = eNegX;
            pinExp[4 * k + 2] = ePosY;
            pinExp[4 * k + 3] = eNegY;
            sumPosX += ePosX; sumNegX += eNegX;
            sumPosY += ePosY; sumNegY += eNegY;
        }
        if (!(sumPosX > 0.0) || !(sumNegX > 0.0) || !(sumPosY > 0.0) || !(sumNegY > 0.0)) continue;

        const double weight = (std::isfinite(net->weight) && net->weight > 0.0f)
                            ? static_cast<double>(net->weight) : 1.0;
        const double invSumPosX = weight / sumPosX;
        const double invSumNegX = weight / sumNegX;
        const double invSumPosY = weight / sumPosY;
        const double invSumNegY = weight / sumNegY;

        for (std::size_t k = 0; k < degree; ++k)
        {
            const std::int32_t nodeIdx = pinNodes[k];
            if (nodeIdx < 0) continue;  // 固定引脚只参与包围盒，不产生梯度
            const std::size_t idx = static_cast<std::size_t>(nodeIdx);
            gradWireX[idx] += pinExp[4 * k + 0] * invSumPosX - pinExp[4 * k + 1] * invSumNegX;
            gradWireY[idx] += pinExp[4 * k + 2] * invSumPosY - pinExp[4 * k + 3] * invSumNegY;
        }
    }

    for (std::size_t i = 0; i < gradWireX.size(); ++i)
    {
        if (!IsValidFinite(gradWireX[i])) gradWireX[i] = 0.0;
        if (!IsValidFinite(gradWireY[i])) gradWireY[i] = 0.0;
        wireGradMax = std::max(wireGradMax, std::fabs(gradWireX[i]));
        wireGradMax = std::max(wireGradMax, std::fabs(gradWireY[i]));
    }
}

// 1-4-4 更新总梯度 ∇f = ∇WL + λ·∇D（公式(2)）
void GlobalPlacement::GetTotalGradient(PlaceData* db, double lambda)
{
    UpdateBinDensities(db);
    ComputeDensityGradient(db);
    ComputeWirelengthGradient(db);

    const std::size_t nodeCount = gradWireX.size();

    for (std::size_t i = 0; i < nodesAndFillers.size(); ++i)
    {
        Module* mod = nodesAndFillers[i];
        double gx = lambda * gradDenX[i];
        double gy = lambda * gradDenY[i];

        if (i < nodeCount)
        {
            gx += gradWireX[i];
            gy += gradWireY[i];
        }

        if (!mod->isFixed)
        {
            const double cx = static_cast<double>(mod->center.x);
            const double cy = static_cast<double>(mod->center.y);
            gx += BoundaryPenaltyGrad(cx, static_cast<double>(chipLL.x), static_cast<double>(chipUR.x), boundK);
            gy += BoundaryPenaltyGrad(cy, static_cast<double>(chipLL.y), static_cast<double>(chipUR.y), boundK);
        }

        if (!IsValidFinite(gx)) gx = 0.0;
        if (!IsValidFinite(gy)) gy = 0.0;
        gradTotalX[i] = gx;
        gradTotalY[i] = gy;
    }
}

// ============================================================================
// 3. 位置更新与参数自适应
// ============================================================================

// Nesterov 动量更新：
//   v_k = β·v_{k-1} + (1-β)·∇f
//   x_{k+1} = x_k - α·v_k
// 学习率 α 按"本次迭代允许的单步位移 stepLimit"归一化（α = stepLimit / 典型动量），
// 这样不同规模、不同梯度量级的电路都能用同一套超参数稳定收敛。
void GlobalPlacement::NesterovUpdate(PlaceData* db, double beta, double stepLimit)
{
    double maxMomentum = 0.0;
    for (std::size_t i = 0; i < nodesAndFillers.size(); ++i)
    {
        if (nodesAndFillers[i]->isFixed) continue;
        vx[i] = beta * vx[i] + (1.0 - beta) * gradTotalX[i];
        vy[i] = beta * vy[i] + (1.0 - beta) * gradTotalY[i];
        if (!IsValidFinite(vx[i])) vx[i] = 0.0;
        if (!IsValidFinite(vy[i])) vy[i] = 0.0;
        maxMomentum = std::max(maxMomentum, std::fabs(vx[i]));
        maxMomentum = std::max(maxMomentum, std::fabs(vy[i]));
    }
    (void)db;

    if (!(maxMomentum > 0.0))
        return;

    // 步长尺度用 |v| 的中位数（抽样估计）而不是最大值：
    // 少数处于极端密集区的单元梯度可能比普通单元大几个数量级，
    // 用最大值归一化会把所有普通单元的位移压到几乎为 0（实测会完全停滞）。
    double scale = maxMomentum;
    {
        momentumSample.clear();
        const std::size_t total = nodesAndFillers.size();
        const std::size_t stride = std::max<std::size_t>(1, total / 8192);
        momentumSample.reserve(total / stride + 1);
        for (std::size_t i = 0; i < total; i += stride)
        {
            if (nodesAndFillers[i]->isFixed) continue;
            momentumSample.push_back(std::max(std::fabs(vx[i]), std::fabs(vy[i])));
        }
        if (!momentumSample.empty())
        {
            const std::size_t mid = momentumSample.size() / 2;
            std::nth_element(momentumSample.begin(), momentumSample.begin() + mid, momentumSample.end());
            const double median = momentumSample[mid];
            if (median > 0.0) scale = median;
        }
    }
    // 下界：线长梯度的量级是 O(1)（WA 权重是概率），拿它做参考可以避免
    // "整体梯度都接近 0（例如密度已经均匀、只剩数值噪声）时把噪声放大成满步长"。
    scale = std::max(scale, 1e-4 * wireGradMax);
    if (!(scale > 0.0))
        return;

    const double alpha = stepLimit / scale;
    const double maxDelta = 4.0 * stepLimit;  // 单次位移上限，避免极端单元数值爆炸
    const double centerX = 0.5 * (static_cast<double>(chipLL.x) + static_cast<double>(chipUR.x));
    const double centerY = 0.5 * (static_cast<double>(chipLL.y) + static_cast<double>(chipUR.y));

    for (std::size_t i = 0; i < nodesAndFillers.size(); ++i)
    {
        Module* mod = nodesAndFillers[i];
        if (mod->isFixed) continue;

        const double dx = clampValue(-alpha * vx[i], -maxDelta, maxDelta);
        const double dy = clampValue(-alpha * vy[i], -maxDelta, maxDelta);
        double nx = static_cast<double>(mod->center.x) + dx;
        double ny = static_cast<double>(mod->center.y) + dy;

        const double hw = static_cast<double>(mod->width) * 0.5;
        const double hh = static_cast<double>(mod->height) * 0.5;
        const double loX = static_cast<double>(chipLL.x) + hw;
        const double hiX = static_cast<double>(chipUR.x) - hw;
        const double loY = static_cast<double>(chipLL.y) + hh;
        const double hiY = static_cast<double>(chipUR.y) - hh;

        nx = (loX <= hiX) ? clampValue(nx, loX, hiX) : centerX;
        ny = (loY <= hiY) ? clampValue(ny, loY, hiY) : centerY;

        mod->center.x = static_cast<float>(nx);
        mod->center.y = static_cast<float>(ny);
    }

    binDensityDirty = true;
}

// 随机散布：仅当输入 .pl 的密度溢出率过高（单元全部堆在原点/角点）时作为兜底起点。
// 这不是本模块的核心算法，只是避免全局布局从一个"所有单元重合在一个点"的
// 退化初始解出发——那种状态下芯片边界处的电场恒为 0（DCT 解满足零通量边界条件），
// 单元得不到任何密度方向的推力，迭代会停滞。
void GlobalPlacement::ScatterInit(PlaceData* db)
{
    std::mt19937 rng(20261017u);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (Module* m : db->Nodes)
    {
        if (m->isFixed) continue;
        const double hw = static_cast<double>(m->width) * 0.5;
        const double hh = static_cast<double>(m->height) * 0.5;
        const double loX = static_cast<double>(chipLL.x) + hw;
        const double hiX = static_cast<double>(chipUR.x) - hw;
        const double loY = static_cast<double>(chipLL.y) + hh;
        const double hiY = static_cast<double>(chipUR.y) - hh;
        const double x = (loX < hiX) ? (loX + (hiX - loX) * unit(rng)) : 0.5 * (loX + hiX);
        const double y = (loY < hiY) ? (loY + (hiY - loY) * unit(rng)) : 0.5 * (loY + hiY);
        m->center.x = static_cast<float>(x);
        m->center.y = static_cast<float>(y);
    }
    binDensityDirty = true;
}

// 把所有可移动单元投影回芯片内部（filler 一并处理）
void GlobalPlacement::ClampModulesToChip(PlaceData* db)
{
    (void)db;
    const double centerX = 0.5 * (static_cast<double>(chipLL.x) + static_cast<double>(chipUR.x));
    const double centerY = 0.5 * (static_cast<double>(chipLL.y) + static_cast<double>(chipUR.y));

    bool changed = false;
    for (Module* mod : nodesAndFillers)
    {
        if (mod->isFixed) continue;
        const double hw = static_cast<double>(mod->width) * 0.5;
        const double hh = static_cast<double>(mod->height) * 0.5;
        const double loX = static_cast<double>(chipLL.x) + hw;
        const double hiX = static_cast<double>(chipUR.x) - hw;
        const double loY = static_cast<double>(chipLL.y) + hh;
        const double hiY = static_cast<double>(chipUR.y) - hh;

        const double cx = static_cast<double>(mod->center.x);
        const double cy = static_cast<double>(mod->center.y);
        const double nx = (loX <= hiX) ? clampValue(cx, loX, hiX) : centerX;
        const double ny = (loY <= hiY) ? clampValue(cy, loY, hiY) : centerY;
        if (nx != cx || ny != cy)
        {
            mod->center.x = static_cast<float>(nx);
            mod->center.y = static_cast<float>(ny);
            changed = true;
        }
    }
    if (changed) binDensityDirty = true;
}

// 2-3 更新 lambda（公式(37)：λ_k = μ_k·λ_{k-1}）
double GlobalPlacement::UpdateLambda(double oldLambda, double deltaHPWL)
{
    // μ = 1.1^(1 - ΔHPWL/ΔHPWL_ref)
    // ΔHPWL 取最近若干次迭代的窗口值，避免单步噪声让 λ 来回振荡；
    // ΔHPWL_ref 按电路规模放大（文档给的是 3.5e5，只适用于小规模算例）；
    // μ 限制在 [0.95, 1.05]，保证 λ 平滑变化。
    const double referenceDelta = std::max(3.5e5, 1e-3 * initHPWLRef);
    double mu = std::pow(1.1, -(deltaHPWL / referenceDelta) + 1.0);
    mu = clampValue(mu, 0.95, 1.05);

    double lambda = oldLambda * mu;
    if (!IsValidFinite(lambda) || lambda <= 0.0)
        lambda = oldLambda;

    return clampValue(lambda, 1e-9, 1e12);
}

// ============================================================================
// 4. 评估
// ============================================================================

// 3.3 输出线长：HPWL = Σ_e (max x - min x) + (max y - min y)
// 引脚位置 = 单元中心 + 引脚偏移（与 placement_report 的统计口径一致）
double GlobalPlacement::ComputeHPWL(PlaceData* db)
{
    double total = 0.0;
    for (const Net* net : db->Nets)
    {
        if (!net || net->netPins.size() < 2) continue;

        double xmin = std::numeric_limits<double>::max();
        double xmax = -std::numeric_limits<double>::max();
        double ymin = std::numeric_limits<double>::max();
        double ymax = -std::numeric_limits<double>::max();
        bool valid = true;

        for (const Pin* pin : net->netPins)
        {
            if (!pin || !pin->module) { valid = false; break; }
            const double px = static_cast<double>(pin->module->center.x) + static_cast<double>(pin->offset.x);
            const double py = static_cast<double>(pin->module->center.y) + static_cast<double>(pin->offset.y);
            if (!IsValidFinite(px) || !IsValidFinite(py)) { valid = false; break; }
            xmin = std::min(xmin, px); xmax = std::max(xmax, px);
            ymin = std::min(ymin, py); ymax = std::max(ymax, py);
        }
        if (!valid) continue;
        total += (xmax - xmin) + (ymax - ymin);
    }
    return total;
}

// 停止判定用的密度溢出率：
//   τ = 总溢出面积 / 归一化总单元面积
//   总溢出面积 = Σ_b max(0, b 内单元面积 - targetDensity·b 面积)
//   归一化总单元面积 = 标准单元总面积 + 宏单元总面积·targetDensity
double GlobalPlacement::ComputeTauOverflow(PlaceData* /*db*/)
{
    double overflowArea = 0.0;
    for (const auto& column : bins)
    {
        for (const BinGrid& b : column)
        {
            const double limit = targetBinDensity * b.area;
            const double cellArea = b.CellArea();
            if (cellArea > limit)
                overflowArea += (cellArea - limit);
        }
    }

    const double normalizedArea = stdCellArea + targetBinDensity * macroArea;
    if (!(normalizedArea > 1e-12))
        return 0.0;

    const double tau = overflowArea / normalizedArea;
    return IsValidFinite(tau) ? std::max(0.0, tau) : 1.0;
}

// ============================================================================
// 5. 主流程
// ============================================================================

void GlobalPlacement::run(PlaceData* db)
{
    if (!db)
        return;

    std::printf("==== ePlace 静电类比全局布局开始 ====\n");

    SetChipBoundary(db);
    if (!(static_cast<double>(chipUR.x) > static_cast<double>(chipLL.x))
        || !(static_cast<double>(chipUR.y) > static_cast<double>(chipLL.y)))
    {
        std::printf("错误：无法确定有效的布局区域（SiteRow 为空且单元包围盒无效）\n");
        return;
    }

    // BinInit 先于 FillerInit：filler 总量需要知道终端/不可放置区域已占用的面积
    BinInit(db);
    FillerInit(db);
    GradientVectorInitialization(db);
    BuildNetPinIndex(db);

    // 初始 .pl 可能把单元堆在芯片外（例如全部放在原点），先投影回 core 内部，
    // 否则这些单元不会计入任何 bin 的密度，密度梯度会低估初始的严重堆叠。
    ClampModulesToChip(db);

    // 退化初始解兜底：初始 .pl 若把所有单元都放在同一个点，
    // 芯片边界处电场恒为 0（DCT 解满足零通量边界条件），迭代会停滞。
    UpdateBinDensities(db);
    lastTau = ComputeTauOverflow(db);
    if (lastTau > 2.0)
    {
        std::printf("初始布局全部堆叠(tau=%.4f)，先做一次随机散布作为全局布局起点\n", lastTau);
        ScatterInit(db);
        ClampModulesToChip(db);
        UpdateBinDensities(db);
        lastTau = ComputeTauOverflow(db);
        std::printf("随机散布后 tau = %.4f\n", lastTau);
    }

    const double initHPWL = ComputeHPWL(db);
    initHPWLRef = initHPWL;
    double lambda = PenaltyFactorInit(db);
    const double lambdaInit = lambda;
    double tau = lastTau;

    std::printf("布局区域: (%.1f, %.1f) - (%.1f, %.1f)\n",
                chipLL.x, chipLL.y, chipUR.x, chipUR.y);
    std::printf("单元: %zu 个（标准单元面积 %.1f，宏单元面积 %.1f），filler: %zu 个\n",
                db->Nodes.size(), stdCellArea, macroArea, fillerNodes.size());
    std::printf("密度网格: %d x %d（bin 尺寸 %.2f x %.2f）\n", binM, binM, binWidth, binHeight);
    std::printf("初始 HPWL = %.3f，初始 tau = %.4f，lambda = %.6e\n", initHPWL, tau, lambda);

    iterLogFile.open("placement_iter.log");
    if (iterLogFile.is_open())
        iterLogFile << "iter,hpwl,tau,lambda\n";

    const double chipShortSide = std::min(static_cast<double>(chipUR.x) - static_cast<double>(chipLL.x),
                                          static_cast<double>(chipUR.y) - static_cast<double>(chipLL.y));
    int iter = 0;
    int doneIters = 0;

    for (; iter < maxIter; ++iter)
    {
        // 后期逐步收紧步长上限，减少收敛阶段的抖动；
        // 同时把单步位移限制在半个网格以内——步长大于网格尺寸时，
        // 单元会"跳过"密度网格，迭代变成数值不稳定（实测会导致密度堆积）。
        const double decay = std::max(0.2, 1.0 - 1.8 * static_cast<double>(iter)
                                            / static_cast<double>(std::max(1, maxIter)));
        const double stepLimit = std::min(moveFraction * chipShortSide * decay,
                                          0.5 * std::min(binWidth, binHeight));

        GetTotalGradient(db, lambda);
        NesterovUpdate(db, betaMomentum, stepLimit);
        ClampModulesToChip(db);

        const double currHPWL = ComputeHPWL(db);

        UpdateBinDensities(db);
        tau = ComputeTauOverflow(db);

        hpwlHistory.push_back(currHPWL);
        if (hpwlHistory.size() > 10) hpwlHistory.erase(hpwlHistory.begin());
        lambda = UpdateLambda(lambda, currHPWL - hpwlHistory.front());
        lambda = clampValue(lambda, lambdaInit * 1e-3, lambdaInit * 1e6);

        std::printf("iter %4d | HPWL %.1f | tau %.4f | lambda %.6e\n", iter, currHPWL, tau, lambda);
        if (iterLogFile.is_open())
            iterLogFile << iter << "," << currHPWL << "," << tau << "," << lambda << "\n";
        doneIters = iter + 1;

        if (tau < targetTau)
        {
            std::printf(">>> 密度溢出率 tau = %.4f 低于目标 %.4f，迭代终止\n", tau, targetTau);
            break;
        }
    }

    if (iterLogFile.is_open())
        iterLogFile.close();

    const double finalHPWL = ComputeHPWL(db);
    std::printf("==== 全局布局完成：迭代 %d 次，HPWL %.3f -> %.3f（变化 %+.2f%%），tau = %.4f ====\n",
                doneIters, initHPWL, finalHPWL,
                (initHPWL > 0.0) ? (100.0 * (finalHPWL - initHPWL) / initHPWL) : 0.0, tau);

    for (Module* filler : fillerNodes)
        delete filler;
    fillerNodes.clear();
    nodesAndFillers.clear();
}
