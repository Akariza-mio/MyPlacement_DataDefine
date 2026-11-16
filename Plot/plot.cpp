#define cimg_display 0 //纯输出位图，不启用 X11/GDI 显示后端，免除外链依赖
#include "CImg.h"

#include "plot.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <ostream>
#include <sstream>
#include <vector>

using namespace cimg_library;
using namespace std;

//==================== 1. 布局结果分析 ====================
void PlacementReport::analyze(PlaceData *db, const string &title)
{
    //---- 面积统计 ----
    double movableArea = 0, fixedArea = 0;
    for (Module *m : db->Nodes)
        movableArea += (double)m->width * m->height;
    for (Module *m : db->Terminals)
        fixedArea += (double)m->width * m->height;

    //---- 核心区包围盒（SiteRows 并集的外沿）----
    double minX = 1e30, minY = 1e30, maxX = -1e30, maxY = -1e30;
    for (const SiteRow &r : db->SiteRows)
    {
        minX = min(minX, (double)r.start.x);
        minY = min(minY, (double)r.start.y);
        maxX = max(maxX, (double)r.end.x);
        maxY = max(maxY, (double)r.start.y + r.height);
    }
    double coreW = maxX - minX, coreH = maxY - minY;
    double rowArea = 0;
    for (const SiteRow &r : db->SiteRows)
        rowArea += ((double)r.end.x - r.start.x) * r.height;

    //---- HPWL（半周长线长模型）: pin 绝对坐标 = ll + offset（仅朝向 N 成立）----
    double hpwl = 0;
    size_t deg1 = 0;
    for (Net *n : db->Nets)
    {
        if (n->netPins.empty()) { deg1++; continue; }
        double x0 = 1e30, y0 = 1e30, x1 = -1e30, y1 = -1e30;
        for (Pin *p : n->netPins)
        {
            double px = p->module->ll.x + p->offset.x;
            double py = p->module->ll.y + p->offset.y;
            x0 = min(x0, px); y0 = min(y0, py);
            x1 = max(x1, px); y1 = max(y1, py);
        }
        hpwl += (x1 - x0) + (y1 - y0);
    }

    //---- 20x20 网格密度（可动单元按中心归属；固定单元不计入）----
    const int NB = 20;
    double bw = coreW / NB, bh = coreH / NB, binArea = bw * bh;
    vector<vector<double>> dens(NB, vector<double>(NB, 0.0));
    for (Module *m : db->Nodes)
    {
        int bx = (int)((m->center.x - minX) / bw);
        int by = (int)((m->center.y - minY) / bh);
        bx = max(0, min(NB - 1, bx));
        by = max(0, min(NB - 1, by));
        dens[by][bx] += (double)m->width * m->height;
    }
    double dmax = 0, dsum = 0;
    int dmaxI = 0, dmaxJ = 0;
    for (int i = 0; i < NB; i++)
        for (int j = 0; j < NB; j++)
        {
            double d = dens[i][j] / binArea;
            dsum += d;
            if (d > dmax) { dmax = d; dmaxI = i; dmaxJ = j; }
        }

    printf("================ 布局结果分析 %s ================\n", title.c_str());
    printf("可动单元  : %zu 个, 总面积 %.0f\n", db->Nodes.size(), movableArea);
    printf("固定单元  : %zu 个, 总面积 %.0f\n", db->Terminals.size(), fixedArea);
    printf("核心区    : [%.0f, %.0f] x [%.0f, %.0f], %d 行, 行总面积 %.0f\n",
           minX, maxX, minY, maxY, (int)db->SiteRows.size(), rowArea);
    printf("利用率    : %.2f%% (可动面积/行面积), %.2f%% (可动面积/核心区面积)\n",
           100.0 * movableArea / rowArea, 100.0 * movableArea / (coreW * coreH));
    printf("HPWL      : 总 %.0f (单位与 BookShelf 坐标一致), 平均每网 %.2f, 空网 %zu\n",
           hpwl, hpwl / max<size_t>(1, db->Nets.size() - deg1), deg1);
    printf("网格密度  : 平均 %.1f%%, 最大 %.1f%% (bin[%d,%d])\n",
           100.0 * dsum / (NB * NB), 100.0 * dmax, dmaxJ, dmaxI);
    printf("=====================================================\n");
}

//==================== 2. 演示用初始散布（legalizer 最小雏形）====================
//策略：按宽度降序排序后逐行贪心装填（first-fit 的退化形式）
//注意：这是为可视化生成测试输入的临时实现，正式合法化算法应考虑位移最小化
void PlacementReport::demoSpreadPlacement(PlaceData *db)
{
    if (db->SiteRows.empty())
    {
        cerr << "[demoSpreadPlacement] 无 SiteRows，跳过" << endl;
        return;
    }
    const double rowH = db->SiteRows[0].height;

    vector<Module *> cells = db->Nodes;
    sort(cells.begin(), cells.end(),
         [](const Module *a, const Module *b) { return a->width > b->width; });

    //溢出堆放区：核心区上方
    double minX = 1e30, maxY = -1e30;
    for (const SiteRow &r : db->SiteRows)
    {
        minX = min(minX, (double)r.start.x);
        maxY = max(maxY, (double)r.start.y + r.height);
    }
    size_t placed = 0, overflow = 0, overCol = 0, overRow = 0;

    size_t ri = 0;
    double cursor = db->SiteRows[0].start.x;
    for (Module *m : cells)
    {
        //多行高度单元或装不下的单元 → 溢出区（adaptec1 中不会触发）
        bool fits = (m->height <= rowH + 1e-6);
        if (fits)
        {
            while (ri < db->SiteRows.size())
            {
                const SiteRow &r = db->SiteRows[ri];
                if (cursor + m->width <= r.end.x + 1e-6) break;
                ri++;
                if (ri < db->SiteRows.size()) cursor = db->SiteRows[ri].start.x;
            }
            fits = (ri < db->SiteRows.size());
        }
        if (fits)
        {
            const SiteRow &r = db->SiteRows[ri];
            m->ll.x = cursor;
            m->ll.y = r.start.y;
            cursor += m->width;
            placed++;
        }
        else
        {
            m->ll.x = minX + (double)(overCol % 50) * 250.0;
            m->ll.y = maxY + 100.0 + (double)(overRow) * rowH;
            if (++overCol % 50 == 0) overRow++;
            overflow++;
        }
        m->center.x = m->ll.x + m->width / 2.0f;
        m->center.y = m->ll.y + m->height / 2.0f;
    }
    printf("[demoSpreadPlacement] 按行贪心装填完成: 装入 %zu, 溢出 %zu\n", placed, overflow);
}

//==================== 2b. GP 代理输入 + HPWL 辅助 ====================
//全局布局由他人模块产出; 开发合法化器时用均匀随机散布模拟其输出
void PlacementReport::simulateGlobalPlacement(PlaceData *db, unsigned seed)
{
    if (db->SiteRows.empty()) return;
    double cx0 = 1e30, cy0 = 1e30, cx1 = -1e30, cy1 = -1e30;
    for (const SiteRow &r : db->SiteRows)
    {
        cx0 = min(cx0, (double)r.start.x); cy0 = min(cy0, (double)r.start.y);
        cx1 = max(cx1, (double)r.end.x);   cy1 = max(cy1, (double)r.start.y + r.height);
    }
    //确定性线性同余发生器, 保证结果可复现且跨平台一致
    unsigned long long s = seed;
    s = s * 2862933555777941757ULL + 3037000493ULL;
    auto rnd = [&]()
    {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return (double)((s >> 33) & 0x7FFFFFFFULL) / 2147483647.0;
    };
    for (Module *m : db->Nodes)
    {
        m->ll.x = (float)(cx0 + rnd() * (cx1 - cx0 - m->width));
        m->ll.y = (float)(cy0 + rnd() * (cy1 - cy0 - m->height));
        m->center.x = m->ll.x + m->width / 2.0f;
        m->center.y = m->ll.y + m->height / 2.0f;
    }
    printf("[simulateGlobalPlacement] 均匀随机散布完成 (seed=%u)\n", seed);
}

double PlacementReport::hpwl(PlaceData *db)
{
    double hp = 0;
    for (Net *n : db->Nets)
    {
        if (n->netPins.empty()) continue;
        double x0 = 1e30, y0 = 1e30, x1 = -1e30, y1 = -1e30;
        for (Pin *p : n->netPins)
        {
            double px = p->module->ll.x + p->offset.x;
            double py = p->module->ll.y + p->offset.y;
            x0 = min(x0, px); y0 = min(y0, py);
            x1 = max(x1, px); y1 = max(y1, py);
        }
        hp += (x1 - x0) + (y1 - y0);
    }
    return hp;
}

//==================== 3. CImg 位图输出 ====================
bool PlacementReport::plotImage(PlaceData *db, const string &bmpPath, int imgSize)
{
    //---- 包围盒：全部单元 + 核心区 ----
    double bx0 = 1e30, by0 = 1e30, bx1 = -1e30, by1 = -1e30;
    auto expand = [&](double x, double y, double w, double h)
    {
        bx0 = min(bx0, x); by0 = min(by0, y);
        bx1 = max(bx1, x + w); by1 = max(by1, y + h);
    };
    for (Module *m : db->Nodes) expand(m->ll.x, m->ll.y, m->width, m->height);
    for (Module *m : db->Terminals) expand(m->ll.x, m->ll.y, m->width, m->height);
    for (const SiteRow &r : db->SiteRows)
        expand(r.start.x, r.start.y, r.end.x - r.start.x, r.height);
    const double pad = (bx1 - bx0 + by1 - by0) * 0.01;
    bx0 -= pad; by0 -= pad; bx1 += pad; by1 += pad;
    const double sceneW = bx1 - bx0, sceneH = by1 - by0;

    const int M = 48; //像素边距
    const double scale = (double)(imgSize - 2 * M) / max(sceneW, sceneH);
    auto X = [&](double x) { return (int)((x - bx0) * scale) + M; };
    auto Y = [&](double y) { return imgSize - (int)((y - by0) * scale) - M; }; //图像 y 轴向下，需翻转

    CImg<unsigned char> img(imgSize, imgSize, 1, 3);
    img.fill(255);

    //---- 核心区底色（浅灰）----
    double cx0 = 1e30, cy0 = 1e30, cx1 = -1e30, cy1 = -1e30;
    for (const SiteRow &r : db->SiteRows)
    {
        cx0 = min(cx0, (double)r.start.x); cy0 = min(cy0, (double)r.start.y);
        cx1 = max(cx1, (double)r.end.x);   cy1 = max(cy1, (double)r.start.y + r.height);
    }
    unsigned char rowbg[] = {233, 233, 229};
    unsigned char movC[] = {70, 110, 205};
    unsigned char fixC[] = {205, 62, 62};
    unsigned char black[] = {25, 25, 25};
    unsigned char white[] = {255, 255, 255};
    if (!db->SiteRows.empty())
        img.draw_rectangle(X(cx0), Y(cy1), X(cx1), Y(cy0), rowbg);

    //---- 可动单元（蓝）----
    for (Module *m : db->Nodes)
        img.draw_rectangle(X(m->ll.x), Y(m->ll.y + m->height),
                           X(m->ll.x + m->width), Y(m->ll.y), movC);

    //---- 固定单元（红，后画置顶）----
    for (Module *m : db->Terminals)
        img.draw_rectangle(X(m->ll.x), Y(m->ll.y + m->height),
                           X(m->ll.x + m->width), Y(m->ll.y), fixC);

    //---- 核心区边框 ----
    if (!db->SiteRows.empty())
    {
        img.draw_line(X(cx0), Y(cy0), X(cx1), Y(cy0), black);
        img.draw_line(X(cx1), Y(cy0), X(cx1), Y(cy1), black);
        img.draw_line(X(cx1), Y(cy1), X(cx0), Y(cy1), black);
        img.draw_line(X(cx0), Y(cy1), X(cx0), Y(cy0), black);
    }

    //---- 标题与图例 ----
    char buf[256];
    snprintf(buf, sizeof(buf), "adaptec1 | movable %zu | fixed %zu | legalized (Tetris)",
             db->Nodes.size(), db->Terminals.size());
    img.draw_text(12, 10, buf, black, white, 1.0f, 13);
    int ly = imgSize - 26;
    img.draw_rectangle(12, ly, 26, ly + 14, fixC);
    img.draw_text(32, ly + 1, "fixed cell", black, white, 1.0f, 13);
    img.draw_rectangle(160, ly, 174, ly + 14, movC);
    img.draw_text(180, ly + 1, "movable cell", black, white, 1.0f, 13);
    img.draw_rectangle(330, ly, 344, ly + 14, rowbg);
    img.draw_text(350, ly + 1, "core rows", black, white, 1.0f, 13);

    try
    {
        img.save(bmpPath.c_str());
    }
    catch (CImgException &e)
    {
        cerr << "[plotImage] 保存失败: " << e.what() << endl;
        return false;
    }
    cout << "[plotImage] 位图已输出: " << bmpPath << " (" << imgSize << "x" << imgSize << ")" << endl;
    return true;
}

//==================== 3b. 违规视图（合法化检查结果可视化）====================
bool PlacementReport::plotViolations(PlaceData *db, const vector<Module *> &violators,
                                     const string &bmpPath, int imgSize)
{
    if (db->SiteRows.empty()) return false;

    //违规单元集合（O(1) 查找；注意 Module::idx 是解析序号，不能当数组下标）
    set<const Module *> violSet(violators.begin(), violators.end());

    //核心区即画布（违规均发生在核心区内）
    double cx0 = 1e30, cy0 = 1e30, cx1 = -1e30, cy1 = -1e30;
    for (const SiteRow &r : db->SiteRows)
    {
        cx0 = min(cx0, (double)r.start.x); cy0 = min(cy0, (double)r.start.y);
        cx1 = max(cx1, (double)r.end.x);   cy1 = max(cy1, (double)r.start.y + r.height);
    }
    const double coreW = cx1 - cx0, coreH = cy1 - cy0;
    const int M = 48;
    const double scale = (double)(imgSize - 2 * M) / max(coreW, coreH);
    auto X = [&](double x) { return (int)((x - cx0) * scale) + M; };
    auto Y = [&](double y) { return imgSize - (int)((y - cy0) * scale) - M; };

    CImg<unsigned char> img(imgSize, imgSize, 1, 3);
    img.fill(255);
    unsigned char rowbg[] = {233, 233, 229};
    unsigned char movC[] = {150, 180, 230};   //正常可动单元: 淡蓝
    unsigned char violC[] = {239, 159, 39};   //违规单元: 橙黄警示
    unsigned char fixC[] = {205, 62, 62};     //固定宏: 红
    unsigned char black[] = {25, 25, 25};
    unsigned char white[] = {255, 255, 255};

    img.draw_rectangle(X(cx0), Y(cy1), X(cx1), Y(cy0), rowbg);

    //正常可动单元（淡蓝，先画）
    for (Module *m : db->Nodes)
        if (!violSet.count(m))
            img.draw_rectangle(X(m->ll.x), Y(m->ll.y + m->height),
                               X(m->ll.x + m->width), Y(m->ll.y), movC);
    //违规单元（橙黄，置顶突出）
    for (Module *m : violators)
        img.draw_rectangle(X(m->ll.x), Y(m->ll.y + m->height),
                           X(m->ll.x + m->width), Y(m->ll.y), violC);
    //固定宏（红色覆盖其上，违规区域呈现"红上有橙"的叠压关系）
    for (Module *m : db->Terminals)
        img.draw_rectangle(X(m->ll.x), Y(m->ll.y + m->height),
                           X(m->ll.x + m->width), Y(m->ll.y), fixC);

    char buf[256];
    snprintf(buf, sizeof(buf), "violations: %zu movable cells overlap fixed macros",
             violators.size());
    img.draw_text(12, 10, buf, black, white, 1.0f, 13);
    int ly = imgSize - 26;
    img.draw_rectangle(12, ly, 26, ly + 14, violC);
    img.draw_text(32, ly + 1, "violator", black, white, 1.0f, 13);
    img.draw_rectangle(160, ly, 174, ly + 14, movC);
    img.draw_text(180, ly + 1, "legal cell", black, white, 1.0f, 13);
    img.draw_rectangle(300, ly, 314, ly + 14, fixC);
    img.draw_text(320, ly + 1, "fixed macro", black, white, 1.0f, 13);

    try { img.save(bmpPath.c_str()); }
    catch (CImgException &e)
    {
        cerr << "[plotViolations] 保存失败: " << e.what() << endl;
        return false;
    }
    cout << "[plotViolations] 违规视图已输出: " << bmpPath << endl;
    return true;
}

//==================== 4. GDSII 输出 ====================
//GDSII 流格式：记录 = [2字节大端长度][记录类型][数据类型][载荷]
//本输出器图层约定: L1=核心区轮廓  L2=布置行  L10=固定单元  L11=可动单元
//坐标单位: 数据库单位 1nm（UNITS 记录声明 user unit=1um, 1dbu=1nm）
namespace
{
inline void putU16be(ostream &f, unsigned v)
{
    f.put((char)((v >> 8) & 0xFF));
    f.put((char)(v & 0xFF));
}
inline void putI32be(ostream &f, long long v)
{
    unsigned u = (unsigned)v;
    f.put((char)((u >> 24) & 0xFF));
    f.put((char)((u >> 16) & 0xFF));
    f.put((char)((u >> 8) & 0xFF));
    f.put((char)(u & 0xFF));
}
inline void putRec(ostream &f, unsigned rt, unsigned dt) //无数据记录
{
    putU16be(f, 4);
    f.put((char)rt);
    f.put((char)dt);
}
inline void putRecU16(ostream &f, unsigned rt, unsigned dt, unsigned v)
{
    putU16be(f, 6);
    f.put((char)rt);
    f.put((char)dt);
    putU16be(f, v);
}
inline void putStrRec(ostream &f, unsigned rt, const string &s) //ASCII 记录，长度须为偶数
{
    string t = s;
    if (t.size() % 2) t += ' ';
    putU16be(f, (unsigned)(4 + t.size()));
    f.put((char)rt);
    f.put((char)0x06);
    f.write(t.data(), t.size());
}
//GDSII REAL8: 1字节(符号+16的指数，偏置64) + 7字节尾数(大端，表示 0.xxx 十六进制小数)
inline void putReal8(ostream &f, double v)
{
    unsigned char b[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    if (v != 0)
    {
        bool neg = (v < 0);
        if (neg) v = -v;
        int e = 0;
        while (v >= 1.0) { v /= 16.0; e++; }
        while (v < 1.0 / 16.0) { v *= 16.0; e--; }
        unsigned long long mant = (unsigned long long)llround(v * 72057594037927936.0); //16^14 = 2^56，7字节尾数满量程（0.xxxxxxxxxxxxxx 共14个十六进制小数位）
        b[0] = (unsigned char)((neg ? 0x80 : 0x00) | ((e + 64) & 0x7F));
        for (int i = 7; i >= 1; i--)
        {
            b[i] = (unsigned char)(mant & 0xFF);
            mant >>= 8;
        }
    }
    f.write((char *)b, 8);
}
inline void putTimestamps(ostream &f, unsigned rt) //BGNLIB/BGNSTR: 12 个 u16 时间戳
{
    putU16be(f, 28);
    f.put((char)rt);
    f.put((char)0x02);
    const unsigned ts[6] = {2026, 11, 9, 14, 40, 0};
    for (int k = 0; k < 2; k++)
        for (int i = 0; i < 6; i++)
            putU16be(f, ts[i]);
}
inline void putRect(ostream &f, int layer, long long x0, long long y0, long long x1, long long y1)
{
    putRec(f, 0x08, 0x00);                 //BOUNDARY
    putRecU16(f, 0x0D, 0x02, layer);       //LAYER
    putRecU16(f, 0x0E, 0x02, 0);           //DATATYPE
    putU16be(f, 44);                       //XY: 5 个闭合点
    f.put((char)0x10);
    f.put((char)0x03);
    putI32be(f, x0); putI32be(f, y0);
    putI32be(f, x1); putI32be(f, y0);
    putI32be(f, x1); putI32be(f, y1);
    putI32be(f, x0); putI32be(f, y1);
    putI32be(f, x0); putI32be(f, y0);
    putRec(f, 0x11, 0x00);                 //ENDEL
}
inline long long nm(double v) { return llround(v * 1000.0); } //um -> nm(DBU)
} // namespace

bool PlacementReport::writeGDS(PlaceData *db, const string &gdsPath)
{
    ofstream f(gdsPath, ios::binary);
    if (!f.is_open())
    {
        cerr << "[writeGDS] 无法创建文件: " << gdsPath << endl;
        return false;
    }
    putRecU16(f, 0x00, 0x02, 600); //HEADER, 版本 6.00
    putTimestamps(f, 0x01);        //BGNLIB
    putStrRec(f, 0x02, "PLACEMENT_LIB"); //LIBNAME
    putU16be(f, 20);               //UNITS: user unit=1um, 1 dbu=1nm
    f.put((char)0x03);
    f.put((char)0x05);
    putReal8(f, 0.001); //user units per DBU
    putReal8(f, 1e-9);  //meters per DBU
    putTimestamps(f, 0x05);       //BGNSTR
    putStrRec(f, 0x06, "TOP");    //STRNAME

    //L1: 核心区轮廓
    if (!db->SiteRows.empty())
    {
        double cx0 = 1e30, cy0 = 1e30, cx1 = -1e30, cy1 = -1e30;
        for (const SiteRow &r : db->SiteRows)
        {
            cx0 = min(cx0, (double)r.start.x); cy0 = min(cy0, (double)r.start.y);
            cx1 = max(cx1, (double)r.end.x);   cy1 = max(cy1, (double)r.start.y + r.height);
        }
        putRect(f, 1, nm(cx0), nm(cy0), nm(cx1), nm(cy1));
        //L2: 布置行
        for (const SiteRow &r : db->SiteRows)
            putRect(f, 2, nm(r.start.x), nm(r.start.y), nm(r.end.x), nm(r.start.y + r.height));
    }
    //L10: 固定单元
    for (Module *m : db->Terminals)
        putRect(f, 10, nm(m->ll.x), nm(m->ll.y), nm(m->ll.x + m->width), nm(m->ll.y + m->height));
    //L11: 可动单元
    for (Module *m : db->Nodes)
        putRect(f, 11, nm(m->ll.x), nm(m->ll.y), nm(m->ll.x + m->width), nm(m->ll.y + m->height));

    putRec(f, 0x07, 0x00); //ENDSTR
    putRec(f, 0x04, 0x00); //ENDLIB

    const double mb = (double)((long long)f.tellp()) / 1048576.0;
    f.close();
    printf("[writeGDS] GDSII 已输出: %s (%.1f MB, 图层 L1=core L2=rows L10=fixed L11=movable)\n",
           gdsPath.c_str(), mb);
    return true;
}
