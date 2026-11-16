#include "legalcheck.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

using namespace std;

static const double EPS = 1e-3; //坐标均为整数值，浮点比较留千分之一容差

//闭区间 [x0, x1]
struct Ivl
{
    double x0, x1;
};

//把区间序列排序并合并为互不重叠的序列
static vector<Ivl> mergeIvls(vector<Ivl> v)
{
    sort(v.begin(), v.end(), [](const Ivl &a, const Ivl &b) { return a.x0 < b.x0; });
    vector<Ivl> out;
    for (const Ivl &iv : v)
    {
        if (!out.empty() && iv.x0 <= out.back().x1 + EPS)
            out.back().x1 = max(out.back().x1, iv.x1);
        else
            out.push_back(iv);
    }
    return out;
}

//cell [cx0,cx1] 与合并后障碍区间的总交长
static double overlapLen(const vector<Ivl> &obs, double cx0, double cx1)
{
    double s = 0;
    for (const Ivl &iv : obs)
    {
        if (iv.x0 >= cx1) break; //obs 已按 x0 排序，可提前终止
        s += max(0.0, min(cx1, iv.x1) - max(cx0, iv.x0));
    }
    return s;
}

LegalReport LegalChecker::check(PlaceData *db, const string &tag, bool verbose)
{
    LegalReport R;
    R.checked = db->Nodes.size();
    if (db->SiteRows.empty())
    {
        cerr << "[LegalChecker] 无 SiteRows" << endl;
        return R;
    }

    //---- 行索引: 底边 y -> 行 ----
    map<double, const SiteRow *> rowMap;
    for (const SiteRow &r : db->SiteRows)
        rowMap[(double)r.start.y] = &r;

    //---- 约束 1-3: 逐单元判定 ----
    for (Module *m : db->Nodes)
    {
        //lower_bound + 容差查找，避免浮点 map 精确匹配问题
        auto it = rowMap.lower_bound((double)m->ll.y - EPS);
        if (it == rowMap.end() || fabs(it->first - m->ll.y) > EPS ||
            fabs(m->height - it->second->height) > EPS)
        {
            R.offRowCnt++;
            continue; //不落行的单元不参与后续行内检查
        }
        const SiteRow &r = *it->second;
        if (m->ll.x < r.start.x - EPS || m->ll.x + m->width > r.end.x + EPS)
            R.outOfRowCnt++;
        else
        {
            double k = (m->ll.x - r.start.x) / r.step;
            if (fabs(k - floor(k + 0.5)) > EPS)
                R.misalignCnt++;
        }
    }

    //---- 按行分组（约束 4 的载体）----
    map<double, vector<Module *>> byRow;
    for (Module *m : db->Nodes)
        byRow[(double)m->ll.y].push_back(m);

    //---- 固定宏在各行带上的障碍区间（裁剪到行 x 范围并合并）----
    map<double, vector<Ivl>> obstacles;
    for (const auto &kv : byRow)
    {
        const double rb = kv.first;            //行底边
        const SiteRow *row = nullptr;
        {
            auto it = rowMap.lower_bound(rb - EPS);
            if (it == rowMap.end() || fabs(it->first - rb) > EPS) continue;
            row = it->second;
        }
        vector<Ivl> ivs;
        for (Module *t : db->Terminals)
        {
            //宏与行带 (rb, rb+rowH) 有面积相交才构成该行障碍
            if (t->ll.y < rb + row->height - EPS && t->ll.y + t->height > rb + EPS)
            {
                Ivl iv;
                iv.x0 = max((double)t->ll.x, (double)row->start.x);
                iv.x1 = min((double)t->ll.x + t->width, (double)row->end.x);
                if (iv.x1 > iv.x0 + EPS) ivs.push_back(iv);
            }
        }
        obstacles[rb] = mergeIvls(ivs);
    }

    //---- 约束 4a: 可动-可动（组内按 x 排序，一次扫描）----
    for (auto &kv : byRow)
    {
        vector<Module *> &v = kv.second;
        sort(v.begin(), v.end(),
             [](const Module *a, const Module *b) { return a->ll.x < b->ll.x; });
        double sumW = 0, unionLen = 0, segStart = 0, segEnd = 0;
        bool open = false;
        double maxEnd = -1e30; //已扫过单元的最右边界
        for (Module *m : v)
        {
            double x0 = m->ll.x, x1 = (double)m->ll.x + m->width;
            if (x0 < maxEnd - EPS)
                R.mmOverlapCells++; //起点落在前人右边界的左侧 -> 参与重叠
            maxEnd = max(maxEnd, x1);
            sumW += m->width;
            //区间并集长度（v 已按 x0 有序）
            if (!open || x0 > segEnd + EPS)
            {
                if (open) unionLen += segEnd - segStart;
                segStart = x0;
                segEnd = x1;
                open = true;
            }
            else
                segEnd = max(segEnd, x1);
        }
        if (open) unionLen += segEnd - segStart;
        R.mmOverlapArea += (sumW - unionLen) * v[0]->height;
    }

    //---- 约束 4b: 可动-固定 ----
    for (const auto &kv : byRow)
    {
        auto oit = obstacles.find(kv.first);
        if (oit == obstacles.end() || oit->second.empty()) continue;
        const vector<Ivl> &obs = oit->second;
        for (Module *m : kv.second)
        {
            double ov = overlapLen(obs, m->ll.x, m->ll.x + m->width);
            if (ov > EPS)
            {
                R.mfOverlapCells++;
                R.mfOverlapPairs++; //合并后的障碍区间互不重叠，交 1 段 = 1 对
                R.mfOverlapArea += ov * m->height;
                R.mfViolators.push_back(m);
            }
        }
    }

    if (verbose)
    {
        printf("---------- 合法化检查 [%s] ----------\n", tag.c_str());
        printf("检查单元   : %zu\n", R.checked);
        printf("行覆盖违规 : %zu   site对齐违规: %zu   越界违规: %zu\n",
               R.offRowCnt, R.misalignCnt, R.outOfRowCnt);
        printf("可动-可动  : 重叠单元 %zu, 重叠面积 %.0f\n", R.mmOverlapCells, R.mmOverlapArea);
        printf("可动-固定  : 重叠单元 %zu, 重叠对 %zu, 重叠面积 %.0f\n",
               R.mfOverlapCells, R.mfOverlapPairs, R.mfOverlapArea);
        printf("结论       : %s\n", R.legal() ? "LEGAL (合法)" : "ILLEGAL (存在违规)");
        printf("-------------------------------------\n");
    }
    return R;
}
