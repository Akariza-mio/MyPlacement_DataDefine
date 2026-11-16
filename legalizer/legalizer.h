#pragma once
#include "placedata.h"
#include <string>

//Tetris 式布局合法化器
//输入: db->Nodes 当前坐标视为全局布局输出(允许不落行/不对齐/重叠/压宏)
//输出: 每个可动单元获得合法坐标(落行/对齐/不越界/不重叠/避开固定宏)
//算法三要素:
//  1. 行自由区间: 行 [起点,终点] 减去与该行带相交的固定宏(裁剪+合并)
//  2. 最近可行槽位: 行内已放置单元把自由区间切成缝隙, 每个缝隙取
//     clamp(目标x, 缝隙左端, 缝隙右端-宽), 全部缝隙取最近
//  3. 行搜索: 从目标行上下交替扩展, 代价 = |Δy|+|Δx|, 双向 y 代价
//     均超过当前最优即早停(行距随偏移单调增)
//处理顺序: 单元按初始 x 升序(Abacus 经典顺序)
class Legalizer
{
public:
    struct Stats
    {
        size_t legalized = 0; //成功放置数
        size_t overflow = 0;  //放不下(核心区满)的单元数
        double sumDisp = 0;   //曼哈顿位移和 Σ(|dx|+|dy|)
        double sumSqDisp = 0; //欧氏平方位移和 Σ(dx²+dy²)
        double maxDisp = 0;   //最大曼哈顿位移
    };

    static Stats legalize(PlaceData *db, bool verbose = true);
};
