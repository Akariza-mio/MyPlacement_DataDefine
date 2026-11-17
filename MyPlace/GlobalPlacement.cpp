#include "GlobalPlacement.h"
#include "Helper.h"
#include <cmath>
#include <cstdio>
#include <random>
#include <algorithm>
#include <limits>

static constexpr double PI = 3.14159265358979323846;

bool GlobalPlacement::IsValidFinite(double val)
{
    return std::isfinite(val);
}

double GlobalPlacement::BoundaryPenaltyGrad(double pos, double ll, double ur, double k)
{
    double g = 0.0;
    if(pos < ll)
    {
        g = k * (ll - pos);
    }
    if(pos > ur)
    {
        g = -k * (pos - ur);
    }
    return g;
}

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
            movableArea += static_cast<double>(m->area);
            cellAreas.push_back(static_cast<double>(m->area));
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
        f->area = static_cast<float>(oneFillerArea);
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

            //统计固定终端密度
            for(auto* term : db->Terminals)
            {
                if(!term->isFixed) continue;
                double lx,ly,rx,ry;
                getModuleBBox(term,lx,ly,rx,ry);
                double ov = rectOverlapArea(lx,ly,rx,ry,b.x0,b.y0,b.x1,b.y1);
                b.terminalDensity += ov;
            }

            //统计不可放置dark区域
            double binTotalArea = binW*binH;
            double rowOverlap = 0.0;
            for(auto& row : db->SiteRows)
            {
                double rlx = static_cast<double>(row.start.x);
                double rly = static_cast<double>(row.bottom);
                double rrx = static_cast<double>(row.end.x);
                double rry = static_cast<double>(row.bottom + row.height);
                rowOverlap += rectOverlapArea(b.x0,b.y0,b.x1,b.y1, rlx,rly,rrx,rry);
            }
            b.darkDensity = binTotalArea - rowOverlap;
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
    for(auto v : gradWireX) if(IsValidFinite(v)) sumWL += std::fabs(v);
    for(auto v : gradWireY) if(IsValidFinite(v)) sumWL += std::fabs(v);
    for(auto v : gradDenX) if(IsValidFinite(v)) sumDen += std::fabs(v);
    for(auto v : gradDenY) if(IsValidFinite(v)) sumDen += std::fabs(v);

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
    double chipW = static_cast<double>(chipUR.x) - static_cast<double>(chipLL.x);
    double chipH = static_cast<double>(chipUR.y) - static_cast<double>(chipLL.y);
    double binW = chipW / binM;
    double binH = chipH / binM;

    for(size_t idx = 0; idx < nodesAndFillers.size(); idx++)
    {
        Module* m = nodesAndFillers[idx];
        if(m->isFixed) continue;
        double lx, ly, rx, ry;
        getModuleBBox(m, lx, ly, rx, ry);
        double mw = rx-lx;
        double mh = ry-ly;

        for(int i = 0; i < binM; i++)
        {
            for(int j = 0; j < binM; j++)
            {
                BinGrid& b = bins[i][j];
                double ov = rectOverlapArea(lx, ly, rx, ry, b.x0, b.y0, b.x1, b.y1);
                if(ov < 1e-12) continue;
                //LengthScale小单元缩放
                double scaleX = (mw < binW) ? mw / binW : 1.0;
                double scaleY = (mh < binH) ? mh / binH : 1.0;
                double scale = std::max(scaleX, scaleY);
                double eff_ov = ov * scale;

                if(m->isFiller)
                    b.fillerDensity += eff_ov;
                else
                    b.nodeDensity += eff_ov;
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

    std::vector<std::vector<double>> rho(binM, std::vector<double>(binM,0.0));
    for(int i=0;i<binM;i++)
    {
        for(int j=0;j<binM;j++)
        {
            double totalArea = bins[i][j].nodeDensity + bins[i][j].fillerDensity
                              + bins[i][j].terminalDensity + bins[i][j].darkDensity;
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

        if(!IsValidFinite(gradDenX[idx])) gradDenX[idx] = 0.0;
        if(!IsValidFinite(gradDenY[idx])) gradDenY[idx] = 0.0;
    }
}

//轻量化采样WA‑LSE平滑线长梯度：大引脚网络均匀采样加速
void GlobalPlacement::ComputeWirelengthGradientWA(PlaceData* db)
{
    for(size_t i = 0; i < gradWireX.size(); i++)
    {
        gradWireX[i] = 0.0;
        gradWireY[i] = 0.0;
    }
    const double gamma = gammaWA;

    for(auto* net : db->Nets)
    {
        auto& pins = net->netPins;
        int np = static_cast<int>(pins.size());
        if(np <= 1) continue;

        std::vector<std::pair<double,double>> pinPos;
        std::vector<size_t> pinNodeIdx;
        pinPos.reserve(np);
        pinNodeIdx.reserve(np);

        for(auto* pin : pins)
        {
            Module* m = pin->module;
            if(m->isFixed) continue;
            //查找节点在db‑Nodes中的下标
            size_t idx = std::string::npos;
            for(size_t k=0;k<db->Nodes.size();k++)
            {
                if(db->Nodes[k]==m) { idx = k; break; }
            }
            if(idx == std::string::npos) continue;
            double x = static_cast<double>(m->center.x);
            double y = static_cast<double>(m->center.y);
            pinPos.emplace_back(x,y);
            pinNodeIdx.push_back(idx);
        }
        int validN = static_cast<int>(pinPos.size());
        if(validN <=1) continue;

        std::vector<std::pair<double,double>> usePins;
        std::vector<size_t> useIdx;

        if(validN > samplePinThresh)
        {
            //大网均匀采样
            int keep = samplePinThresh;
            double step = (validN-1)/(double)(keep-1);
            for(int s=0;s<keep;s++)
            {
                int sel = static_cast<int>(std::round(s*step));
                usePins.push_back(pinPos[sel]);
                useIdx.push_back(pinNodeIdx[sel]);
            }
        }
        else
        {
            usePins.swap(pinPos);
            useIdx.swap(pinNodeIdx);
        }

        int K = static_cast<int>(usePins.size());
        double sumExpX=0,sumExpNegX=0;
        double sumExpY=0,sumExpNegY=0;
        for(auto &p : usePins)
        {
            double x = p.first;
            double y = p.second;
            sumExpX     += std::exp(x/gamma);
            sumExpNegX  += std::exp(-x/gamma);
            sumExpY     += std::exp(y/gamma);
            sumExpNegY  += std::exp(-y/gamma);
        }
        if(std::fabs(sumExpX)<1e-12 || std::fabs(sumExpNegX)<1e-12) continue;
        if(std::fabs(sumExpY)<1e-12 || std::fabs(sumExpNegY)<1e-12) continue;

        for(int s=0;s<K;s++)
        {
            size_t mid = useIdx[s];
            double x = usePins[s].first;
            double y = usePins[s].second;
            double gx = std::exp(x/gamma)/sumExpX  - std::exp(-x/gamma)/sumExpNegX;
            double gy = std::exp(y/gamma)/sumExpY  - std::exp(-y/gamma)/sumExpNegY;
            gradWireX[mid] += gx;
            gradWireY[mid] += gy;
        }
    }
    //数值保护
    for(size_t i=0;i<gradWireX.size();i++)
    {
        if(!IsValidFinite(gradWireX[i])) gradWireX[i]=0.0;
        if(!IsValidFinite(gradWireY[i])) gradWireY[i]=0.0;
    }
}

void GlobalPlacement::GetTotalGradient(PlaceData* db, double lambda)
{
    UpdateBinDensities(db);
    ComputeDensityGradient(db);
    ComputeWirelengthGradientWA(db);

    for(size_t i = 0; i < nodesAndFillers.size(); i++)
    {
        gradTotalX[i] = lambda * gradDenX[i];
        gradTotalY[i] = lambda * gradDenY[i];
        if(i < db->Nodes.size())
        {
            gradTotalX[i] += gradWireX[i];
            gradTotalY[i] += gradWireY[i];
        }
        //叠加软边界惩罚梯度
        Module* m = nodesAndFillers[i];
        if(!m->isFixed)
        {
            double cx = static_cast<double>(m->center.x);
            double cy = static_cast<double>(m->center.y);
            gradTotalX[i] += BoundaryPenaltyGrad(cx, static_cast<double>(chipLL.x), static_cast<double>(chipUR.x), boundK);
            gradTotalY[i] += BoundaryPenaltyGrad(cy, static_cast<double>(chipLL.y), static_cast<double>(chipUR.y), boundK);
        }
        if(!IsValidFinite(gradTotalX[i])) gradTotalX[i]=0.0;
        if(!IsValidFinite(gradTotalY[i])) gradTotalY[i]=0.0;
    }
}

void GlobalPlacement::NesterovUpdate(PlaceData* /*db*/, double alpha, double beta)
{
    for(size_t i = 0; i < nodesAndFillers.size(); i++)
    {
        Module* m = nodesAndFillers[i];
        if(m->isFixed) continue;

        vx[i] = beta * vx[i] + (1.0 - beta) * gradTotalX[i];
        vy[i] = beta * vy[i] + (1.0 - beta) * gradTotalY[i];

        double dx = - alpha * vx[i];
        double dy = - alpha * vy[i];

        double newX = static_cast<double>(m->center.x) + dx;
        double newY = static_cast<double>(m->center.y) + dy;

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
            double rho = (bins[i][j].nodeDensity + bins[i][j].fillerDensity
                          + bins[i][j].terminalDensity + bins[i][j].darkDensity) / binArea;
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
    iterLogFile.open("placement_iter.log");
    if(iterLogFile.is_open())
    {
        iterLogFile << "iter,HPWL,tau,lambda\n";
    }

    //读取芯片边界来自SiteRow
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
        if(iterLogFile.is_open())
        {
            iterLogFile << iter << "," << currHPWL << "," << tau << "," << lambda << "\n";
        }

        if(tau < targetTau)
        {
            printf(">>> 密度溢出tau达到目标，迭代终止\n");
            break;
        }
    }

    iterLogFile.close();
    //释放filler内存
    for(auto* f : fillerNodes) delete f;
    fillerNodes.clear();
    nodesAndFillers.clear();
    printf("==== ePlace全局布局完成 ====\n");
}
