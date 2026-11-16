#include "GlobalPlacement.h"
#include "Helper.h"
#include <cmath>
#include <cstdio>
#include <random>
#include <algorithm>

static constexpr double PI = 3.14159265358979323846;

void GlobalPlacement::FillerInit(PlaceData* db)
{
    double chipW = static_cast<double>(chipUR.x) - static_cast<double>(chipLL.x);
    double chipH = static_cast<double>(chipUR.y) - static_cast<double>(chipLL.y);
    double chipArea = chipW * chipH;

    double movableArea = 0.0;
    std::vector<double> cellAreas;
    for(auto* m : db->Nodes)
    {
        if(!m->isFixed && !m->isMacro)
        {
            movableArea += m->area;
            cellAreas.push_back(m->area);
        }
    }

    double fillerTotalArea = targetBinDensity * chipArea - movableArea;
    if(fillerTotalArea < 1e-9)
        return;

    double avgCellArea = 0.0;
    for(auto a : cellAreas) avgCellArea += a;
    if(!cellAreas.empty()) avgCellArea /= double(cellAreas.size());

    double fillerW = std::sqrt(avgCellArea);
    double fillerH = std::sqrt(avgCellArea);
    double oneFillerArea = fillerW * fillerH;
    int fillerCnt = static_cast<int>(fillerTotalArea / oneFillerArea);

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> distX(static_cast<double>(chipLL.x), static_cast<double>(chipUR.x) - fillerW);
    std::uniform_real_distribution<double> distY(static_cast<double>(chipLL.y), static_cast<double>(chipUR.y) - fillerH);

    for(int i = 0; i < fillerCnt; ++i)
    {
        Module* f = new Module();
        f->name = "filler_" + std::to_string(i);
        f->isFiller = true;
        f->isFixed = false;
        f->isMacro = false;
        f->area = oneFillerArea;
        f->width = static_cast<float>(fillerW);
        f->height = static_cast<float>(fillerH);
        f->center.x = static_cast<float>(distX(rng));
        f->center.y = static_cast<float>(distY(rng));
        fillerNodes.push_back(f);
    }

    nodesAndFillers.clear();
    for(auto* m : db->Nodes) nodesAndFillers.push_back(m);
    for(auto* f : fillerNodes) nodesAndFillers.push_back(f);
}

void GlobalPlacement::BinInit(PlaceData* db)
{
    bins.assign(binM, std::vector<BinGrid>(binM));
    double chipW = static_cast<double>(chipUR.x) - static_cast<double>(chipLL.x);
    double chipH = static_cast<double>(chipUR.y) - static_cast<double>(chipLL.y);
    double binW = chipW / binM;
    double binH = chipH / binM;

    for(int i = 0; i < binM; ++i)
    {
        for(int j = 0; j < binM; ++j)
        {
            BinGrid& b = bins[i][j];
            b.x0 = static_cast<double>(chipLL.x) + i * binW;
            b.y0 = static_cast<double>(chipLL.y) + j * binH;
            b.x1 = b.x0 + binW;
            b.y1 = b.y0 + binH;
            b.nodeDensity = 0.0;
            b.fillerDensity = 0.0;
            b.terminalDensity = 0.0;
            b.darkDensity = 0.0;
        }
    }
    (void)db;
}

void GlobalPlacement::GradientVectorInitialization(PlaceData* db)
{
    size_t origNodeCnt = db->Nodes.size();
    size_t allCnt = nodesAndFillers.size();

    gradWireX.assign(origNodeCnt, 0.0);
    gradWireY.assign(origNodeCnt, 0.0);

    gradDenX.assign(allCnt, 0.0);
    gradDenY.assign(allCnt, 0.0);

    gradTotalX.assign(allCnt, 0.0);
    gradTotalY.assign(allCnt, 0.0);

    vx.assign(allCnt, 0.0);
    vy.assign(allCnt, 0.0);
    (void)db;
}

double GlobalPlacement::PenaltyFactorInit(PlaceData* db)
{
    GetTotalGradient(db, 1.0);
    double sumWL = 0.0, sumDen = 0.0;
    for(auto v : gradWireX) sumWL += std::fabs(v);
    for(auto v : gradWireY) sumWL += std::fabs(v);
    for(auto v : gradDenX) sumDen += std::fabs(v);
    for(auto v : gradDenY) sumDen += std::fabs(v);

    if(sumDen < 1e-12)
        return 1e-3;
    return sumWL / sumDen;
}

void GlobalPlacement::UpdateBinDensities(PlaceData* /*db*/)
{
    for(auto& row : bins)
    {
        for(auto& b : row)
        {
            b.nodeDensity = 0.0;
            b.fillerDensity = 0.0;
        }
    }

    for(size_t idx = 0; idx < nodesAndFillers.size(); idx++)
    {
        Module* m = nodesAndFillers[idx];
        if(m->isFixed) continue;
        double lx, ly, rx, ry;
        getModuleBBox(m, lx, ly, rx, ry);

        for(int i = 0; i < binM; i++)
        {
            for(int j = 0; j < binM; j++)
            {
                BinGrid& b = bins[i][j];
                double ov = rectOverlapArea(lx, ly, rx, ry, b.x0, b.y0, b.x1, b.y1);
                if(ov < 1e-12) continue;
                if(m->isFiller)
                    b.fillerDensity += ov;
                else
                    b.nodeDensity += ov;
            }
        }
    }
}

void GlobalPlacement::ComputeDensityGradient(PlaceData* /*db*/)
{
    gradDenX.assign(gradDenX.size(), 0.0);
    gradDenY.assign(gradDenY.size(), 0.0);
    double chipW = static_cast<double>(chipUR.x) - static_cast<double>(chipLL.x);
    double chipH = static_cast<double>(chipUR.y) - static_cast<double>(chipLL.y);
    double binW = chipW / binM;
    double binH = chipH / binM;
    double binArea = binW * binH;

    // 手写简化DCT‑电场（ePlace静电类比，无FFTW）
    std::vector<std::vector<double>> rho(binM, std::vector<double>(binM,0.0));
    for(int i=0;i<binM;i++)
    {
        for(int j=0;j<binM;j++)
        {
            double totalArea = bins[i][j].nodeDensity + bins[i][j].fillerDensity;
            rho[i][j] = totalArea / binArea;
        }
    }

    std::vector<std::vector<double>> a(binM, std::vector<double>(binM,0.0));
    for(int j=0;j<binM;j++)
    {
        for(int k=0;k<binM;k++)
        {
            double sum = 0.0;
            for(int x=0;x<binM;x++)
            {
                for(int y=0;y<binM;y++)
                {
                    double wj = (j * PI) / binM;
                    double wk = (k * PI) / binM;
                    sum += rho[x][y] * std::cos(wj * x) * std::cos(wk * y);
                }
            }
            a[j][k] = sum / (double)(binM*binM);
        }
    }

    for(size_t idx = 0; idx < nodesAndFillers.size(); idx++)
    {
        Module* m = nodesAndFillers[idx];
        if(m->isFixed) continue;
        double cx = static_cast<double>(m->center.x);
        double cy = static_cast<double>(m->center.y);
        int i = static_cast<int>((cx - static_cast<double>(chipLL.x))/binW);
        int j = static_cast<int>((cy - static_cast<double>(chipLL.y))/binH);
        i = std::max(0, std::min(binM-1,i));
        j = std::max(0, std::min(binM-1,j));

        double xi_x = 0.0, xi_y = 0.0;
        for(int jj=0;jj<binM;jj++)
        {
            for(int kk=0;kk<binM;kk++)
            {
                double wj = (jj*PI)/binM;
                double wk = (kk*PI)/binM;
                double denom = wj*wj + wk*wk;
                if(denom < 1e-12) continue;
                xi_x += a[jj][kk] * wj / denom
                        * std::sin(wj * i) * std::cos(wk * j);
                xi_y += a[jj][kk] * wk / denom
                        * std::cos(wj * i) * std::sin(wk * j);
            }
        }
        double q = static_cast<double>(m->area);
        gradDenX[idx] = - q * xi_x;
        gradDenY[idx] = - q * xi_y;
    }
}

void GlobalPlacement::ComputeWirelengthGradient(PlaceData* db)
{
    gradWireX.assign(gradWireX.size(),0.0);
    gradWireY.assign(gradWireY.size(),0.0);
    const double gamma = 200.0;

    for(auto* net : db->Nets)
    {
        if(net->netPins.size() < 2) continue;
        double xmin = 1e18, xmax = -1e18;
        double ymin = 1e18, ymax = -1e18;
        for(auto* pin : net->netPins)
        {
            double cx = static_cast<double>(pin->module->center.x);
            double cy = static_cast<double>(pin->module->center.y);
            xmin = std::min(xmin, cx);
            xmax = std::max(xmax, cx);
            ymin = std::min(ymin, cy);
            ymax = std::max(ymax, cy);
        }
        for(size_t idx = 0; idx < db->Nodes.size(); idx++)
        {
            Module* m = db->Nodes[idx];
            if(m->isFixed) continue;
            double cx = static_cast<double>(m->center.x);
            double cy = static_cast<double>(m->center.y);
            double dx = 0.0, dy = 0.0;
            if(std::fabs(cx - xmin) < gamma) dx -= 1.0/gamma;
            if(std::fabs(cx - xmax) < gamma) dx += 1.0/gamma;
            if(std::fabs(cy - ymin) < gamma) dy -= 1.0/gamma;
            if(std::fabs(cy - ymax) < gamma) dy += 1.0/gamma;
            gradWireX[idx] += dx;
            gradWireY[idx] += dy;
        }
    }
}

void GlobalPlacement::GetTotalGradient(PlaceData* db, double lambda)
{
    UpdateBinDensities(db);
    ComputeDensityGradient(db);
    ComputeWirelengthGradient(db);

    for(size_t i = 0; i < nodesAndFillers.size(); i++)
    {
        gradTotalX[i] = lambda * gradDenX[i];
        gradTotalY[i] = lambda * gradDenY[i];
        if(i < db->Nodes.size())
        {
            gradTotalX[i] += gradWireX[i];
            gradTotalY[i] += gradWireY[i];
        }
    }
}

void GlobalPlacement::NesterovUpdate(PlaceData* /*db*/, double alpha, double beta)
{
    for(size_t i = 0; i < nodesAndFillers.size(); i++)
    {
        Module* m = nodesAndFillers[i];
        if(m->isFixed) continue;

        // Nesterov前瞻动量更新
        vx[i] = beta * vx[i] + (1.0 - beta) * gradTotalX[i];
        vy[i] = beta * vy[i] + (1.0 - beta) * gradTotalY[i];

        double dx = - alpha * vx[i];
        double dy = - alpha * vy[i];

        double newX = static_cast<double>(m->center.x) + dx;
        double newY = static_cast<double>(m->center.y) + dy;

        // 钳位到芯片边界
        newX = std::max(static_cast<double>(chipLL.x), std::min(static_cast<double>(chipUR.x), newX));
        newY = std::max(static_cast<double>(chipLL.y), std::min(static_cast<double>(chipUR.y), newY));

        m->center.x = static_cast<float>(newX);
        m->center.y = static_cast<float>(newY);
    }
}

double GlobalPlacement::UpdateLambda(double oldLambda, double deltaHPWL)
{
    const double REF_HPWL = 3.5e5;
    double mu = std::pow(1.1, -(deltaHPWL / REF_HPWL) + 1.0);
    mu = std::max(0.75, std::min(1.1, mu));
    return oldLambda * mu;
}

double GlobalPlacement::ComputeHPWL(PlaceData* db)
{
    double total = 0.0;
    for(auto* net : db->Nets)
    {
        if(net->netPins.size() < 2) continue;
        double xmin = 1e18, xmax = -1e18;
        double ymin = 1e18, ymax = -1e18;
        for(auto* pin : net->netPins)
        {
            double cx = static_cast<double>(pin->module->center.x);
            double cy = static_cast<double>(pin->module->center.y);
            xmin = std::min(xmin, cx);
            xmax = std::max(xmax, cx);
            ymin = std::min(ymin, cy);
            ymax = std::max(ymax, cy);
        }
        total += (xmax - xmin) + (ymax - ymin);
    }
    return total;
}

double GlobalPlacement::ComputeTauOverflow(PlaceData* /*db*/)
{
    double chipW = static_cast<double>(chipUR.x) - static_cast<double>(chipLL.x);
    double chipH = static_cast<double>(chipUR.y) - static_cast<double>(chipLL.y);
    double binW = chipW / binM;
    double binH = chipH / binM;
    double binArea = binW * binH;
    double overflowArea = 0.0;
    double totalMovableArea = 0.0;

    for(auto* m : nodesAndFillers)
    {
        if(!m->isFixed) totalMovableArea += static_cast<double>(m->area);
    }
    if(totalMovableArea < 1e-12) return 1.0;

    for(int i = 0; i < binM; i++)
    {
        for(int j = 0; j < binM; j++)
        {
            double rho = (bins[i][j].nodeDensity + bins[i][j].fillerDensity) / binArea;
            if(rho > targetBinDensity)
            {
                overflowArea += (rho - targetBinDensity) * binArea;
            }
        }
    }
    return overflowArea / totalMovableArea;
}

void GlobalPlacement::run(PlaceData* db)
{
    printf("==== ePlace Nesterov全局布局开始 ====\n");

    // 获取芯片边界从SiteRow，用double临时变量规避float/double std::min/max冲突
    double cx_ll = 1e18, cy_ll = 1e18;
    double cx_ur = -1e18, cy_ur = -1e18;
    for(auto& row : db->SiteRows)
    {
        cx_ll  = std::min(cx_ll,  static_cast<double>(row.start.x));
        cy_ll  = std::min(cy_ll,  static_cast<double>(row.bottom));
        cx_ur  = std::max(cx_ur,  static_cast<double>(row.end.x));
        cy_ur  = std::max(cy_ur,  static_cast<double>(row.bottom + row.height));
    }
    chipLL.x = static_cast<float>(cx_ll);
    chipLL.y = static_cast<float>(cy_ll);
    chipUR.x = static_cast<float>(cx_ur);
    chipUR.y = static_cast<float>(cy_ur);

    FillerInit(db);
    BinInit(db);
    GradientVectorInitialization(db);
    double lambda = PenaltyFactorInit(db);

    double prevHPWL = ComputeHPWL(db);
    double tau = 1.0;
    const double stepAlpha = 1e-4;

    for(int iter = 0; iter < maxIter; ++iter)
    {
        GetTotalGradient(db, lambda);
        NesterovUpdate(db, stepAlpha, betaMomentum);

        double currHPWL = ComputeHPWL(db);
        double deltaHPWL = currHPWL - prevHPWL;
        lambda = UpdateLambda(lambda, deltaHPWL);
        prevHPWL = currHPWL;

        tau = ComputeTauOverflow(db);
        printf("iter:%-4d | HPWL:%.2f | tau:%.4f | lambda:%.6e\n",
               iter, currHPWL, tau, lambda);

        if(tau < targetTau)
        {
            printf(">>> 密度溢出tau达到目标，迭代终止\n");
            break;
        }
    }

    // 释放filler内存
    for(auto* f : fillerNodes) delete f;
    fillerNodes.clear();
    nodesAndFillers.clear();
    printf("==== ePlace全局布局完成 ====\n");
}
