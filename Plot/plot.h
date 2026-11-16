#pragma once
#include "placedata.h"
#include <string>
#include <vector>

//布局结果分析与可视化
//1. analyze          : HPWL 线长 / 核心区利用率 / 网格密度 统计
//2. demoSpreadPlacement : 演示用初始散布（按行贪心装填，legalizer 最小雏形，非正式算法）
//3. plotImage        : CImg 位图输出（BMP，无需外部图像库）
//4. writeGDS         : GDSII 流格式输出（二进制，KLayout 可直接打开）
class PlacementReport
{
public:
    static void analyze(PlaceData *db, const std::string &title = "");

    //输入: db 中可动单元当前坐标被改写为按行贪心装填的结果
    static void demoSpreadPlacement(PlaceData *db);

    //GP 代理输入: 均匀随机散布在核心区内(模拟全局布局输出, 允许不落行/重叠/压宏)
    static void simulateGlobalPlacement(PlaceData *db, unsigned seed = 42);

    //HPWL 线长(半周长模型, 当前坐标)
    static double hpwl(PlaceData *db);

    static bool plotImage(PlaceData *db, const std::string &bmpPath, int imgSize = 1600);

    //违规视图: violators 中的单元以警示色高亮（合法化检查结果可视化）
    static bool plotViolations(PlaceData *db, const std::vector<Module *> &violators,
                               const std::string &bmpPath, int imgSize = 1600);

    static bool writeGDS(PlaceData *db, const std::string &gdsPath);
};
