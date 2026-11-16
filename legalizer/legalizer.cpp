#include "legalizer.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace std;

static const double EPS = 1e-3;

namespace
{

struct FreeIvl
{
    double x0, x1;
};

//行自由区间 = [start.x, end.x] − 宏障碍(与行带面积相交者, 裁剪到行内, 合并)
vector<FreeIvl> computeFreeIvls(const SiteRow &r, const vector<Module *> &terms)
{
    vector<pair<double, double>> obs;
    for (Module *t : terms)
    {
        if ((double)t->ll.y < (double)r.start.y + r.height - EPS &&
            (double)t->ll.y + t->height > (double)r.start.y + EPS)
        {
            double x0 = max((double)t->ll.x, (double)r.start.x);
            double x1 = min((double)t->ll.x + t->width, (double)r.end.x);
            if (x1 > x0 + EPS) obs.push_back({x0, x1});
        }
    }
    sort(obs.begin(), obs.end());
    vector<FreeIvl> out;
    double cur = r.start.x;
    for (const auto &o : obs)
    {
        if (o.first > cur + EPS)
            out.push_back({cur, min(o.first, (double)r.end.x)});
        cur = max(cur, o.second);
    }
    if ((double)r.end.x > cur + EPS)
        out.push_back({cur, (double)r.end.x});
    return out;
}

//把 x 吸附到行的 site 网格上: x' = start + round((x-start)/step)*step
inline double snapToGrid(double x, double start, double step)
{
    return start + floor((x - start) / step + 0.5) * step;
}

//在该行找距 tx 最近的可行位置; pv 为该行已放置单元(x 有序)
//缝隙 = 自由区间内 [prevEnd, 下一单元起点] 及尾缝隙 [最后单元终点, 区间终点]
//返回位置已吸附到 site 网格(约束 2), 并保证仍在缝隙内(约束 3/4)
bool bestSlot(const SiteRow &row, const vector<FreeIvl> &fiv, const vector<Module *> &pv,
              double width, double tx, double &outX)
{
    double best = 1e30, bx = 0;
    for (const FreeIvl &iv : fiv)
    {
        if (iv.x1 - iv.x0 < width - EPS) continue;
        //二分定位本区间内的已放置单元(它们是 pv 的连续切片)
        auto b = lower_bound(pv.begin(), pv.end(), iv.x0,
                             [](const Module *m, double v) { return m->ll.x < v; });
        auto e = upper_bound(b, pv.end(), iv.x1,
                             [](double v, const Module *m) { return v < m->ll.x; });
        double prevEnd = iv.x0;
        for (auto it = b; it != e; ++it)
        {
            double gap0 = prevEnd, gap1 = (*it)->ll.x;
            if (gap1 - gap0 >= width - EPS)
            {
                double cand = min(max(tx, gap0), gap1 - width);
                cand = snapToGrid(cand, row.start.x, row.step); //site 对齐
                if (cand < gap0 - EPS) cand += row.step;        //吸附越出左界则右移一格
                if (cand + width > gap1 + EPS) cand -= row.step; //越出右界则左移一格
                if (cand >= gap0 - EPS && cand + width <= gap1 + EPS &&
                    fabs(cand - tx) < best)
                {
                    best = fabs(cand - tx);
                    bx = cand;
                }
            }
            prevEnd = (double)(*it)->ll.x + (*it)->width;
            //缝隙左端已越过目标且已有候选, 更远的缝隙只会更差
            if (gap0 > tx + best + EPS && best < 1e30) break;
        }
        if (iv.x1 - prevEnd >= width - EPS) //尾缝隙
        {
            double cand = min(max(tx, prevEnd), iv.x1 - width);
            cand = snapToGrid(cand, row.start.x, row.step);
            if (cand < prevEnd - EPS) cand += row.step;
            if (cand + width > iv.x1 + EPS) cand -= row.step;
            if (cand >= prevEnd - EPS && cand + width <= iv.x1 + EPS &&
                fabs(cand - tx) < best)
            {
                best = fabs(cand - tx);
                bx = cand;
            }
        }
    }
    if (best < 1e30) { outX = bx; return true; }
    return false;
}

} // namespace

Legalizer::Stats Legalizer::legalize(PlaceData *db, bool verbose)
{
    Stats S;
    if (db->SiteRows.empty())
    {
        cerr << "[Legalizer] 无 SiteRows" << endl;
        return S;
    }
    const size_t NR = db->SiteRows.size();
    const double rowH = db->SiteRows[0].height;
    const double y0 = db->SiteRows[0].start.y;
    const double rowW = db->SiteRows[0].end.x - db->SiteRows[0].start.x;

    //---- 1. 每行自由区间 ----
    vector<vector<FreeIvl>> freeI(NR);
    for (size_t i = 0; i < NR; i++)
        freeI[i] = computeFreeIvls(db->SiteRows[i], db->Terminals);

    //---- 2. 每行已放置单元(x 有序), 支持二分插入 ----
    vector<vector<Module *>> placed(NR);

    //---- 3. 单元按初始 x 升序(Abacus 经典处理顺序) ----
    vector<Module *> cells = db->Nodes;
    sort(cells.begin(), cells.end(),
         [](const Module *a, const Module *b) { return a->ll.x < b->ll.x; });

    //溢出堆放区(核心区上方), 正常不会触发(利用率仅 33%)
    double minX = 1e30, maxY = -1e30;
    for (const SiteRow &r : db->SiteRows)
    {
        minX = min(minX, (double)r.start.x);
        maxY = max(maxY, (double)r.start.y + r.height);
    }
    size_t overCol = 0, overRow = 0;

    for (Module *m : cells)
    {
        const double tx = m->ll.x, ty = m->ll.y;
        int r0 = (int)llround((ty - y0) / rowH);
        r0 = max(0, min((int)NR - 1, r0));

        double bestCost = 1e30;
        int bestRow = -1;
        double bestX = 0;

        //---- 行搜索: 从 r0 上下交替扩展; 双向 y 代价都不可能更优即停 ----
        for (int k = 0;; k++)
        {
            bool alive = false;
            if (r0 + k < (int)NR)
            {
                double yCost = fabs((double)db->SiteRows[r0 + k].start.y - ty);
                if (yCost < bestCost)
                {
                    alive = true;
                    double x;
                    if (bestSlot(db->SiteRows[r0 + k], freeI[r0 + k], placed[r0 + k], m->width, tx, x))
                    {
                        double cost = yCost + fabs(x - tx);
                        if (cost < bestCost) { bestCost = cost; bestRow = r0 + k; bestX = x; }
                    }
                }
            }
            if (k > 0 && r0 - k >= 0)
            {
                double yCost = fabs((double)db->SiteRows[r0 - k].start.y - ty);
                if (yCost < bestCost)
                {
                    alive = true;
                    double x;
                    if (bestSlot(db->SiteRows[r0 - k], freeI[r0 - k], placed[r0 - k], m->width, tx, x))
                    {
                        double cost = yCost + fabs(x - tx);
                        if (cost < bestCost) { bestCost = cost; bestRow = r0 - k; bestX = x; }
                    }
                }
            }
            if (!alive) break; //两个方向都已越界或 y 代价不可能更优
        }

        if (bestRow >= 0)
        {
            double dx = bestX - tx;
            double dy = (double)db->SiteRows[bestRow].start.y - ty;
            m->ll.x = (float)bestX;
            m->ll.y = db->SiteRows[bestRow].start.y;
            m->center.x = m->ll.x + m->width / 2.0f;
            m->center.y = m->ll.y + m->height / 2.0f;
            //二分插入保持该行有序
            auto &pv = placed[bestRow];
            auto it = upper_bound(pv.begin(), pv.end(), bestX,
                                  [](double v, const Module *mod) { return v < mod->ll.x; });
            pv.insert(it, m);
            S.legalized++;
            S.sumDisp += fabs(dx) + fabs(dy);
            S.sumSqDisp += dx * dx + dy * dy;
            S.maxDisp = max(S.maxDisp, fabs(dx) + fabs(dy));
        }
        else
        {
            //溢出: 堆放在核心区上方
            m->ll.x = (float)(minX + (double)(overCol % 40) * (rowW / 40.0));
            m->ll.y = (float)(maxY + 100.0 + (double)overRow * rowH);
            if (++overCol % 40 == 0) overRow++;
            m->center.x = m->ll.x + m->width / 2.0f;
            m->center.y = m->ll.y + m->height / 2.0f;
            S.overflow++;
        }
    }

    if (verbose)
    {
        printf("---------- 合法化 [Tetris] ----------\n");
        printf("放置单元 : %zu / %zu (溢出 %zu)\n", S.legalized, cells.size(), S.overflow);
        printf("位移统计 : 总 %.0f, 平均 %.2f, 最大 %.0f (曼哈顿)\n",
               S.sumDisp, S.sumDisp / max<size_t>(1, S.legalized), S.maxDisp);
        printf("平方位移: 总 %.3e (均方根 %.2f)\n",
               S.sumSqDisp, sqrt(S.sumSqDisp / max<size_t>(1, S.legalized)));
        printf("-------------------------------------\n");
    }
    return S;
}
