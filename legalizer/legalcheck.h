#pragma once
#include "placedata.h"
#include <string>
#include <vector>

//布局合法化检查结果
//约束 1 行覆盖: ll.y == 某行底边且高度 == 行高（多行高单元暂判为不落行）
//约束 2 site 对齐: (ll.x - 行起点x) / step 为整数
//约束 3 不越界: [ll.x, ll.x+width] ⊆ [start.x, end.x]
//约束 4 无重叠: 可动-可动互不重叠; 可动不压固定宏
struct LegalReport
{
    size_t checked = 0;        //参与检查的可动单元数
    size_t offRowCnt = 0;      //违反约束 1
    size_t misalignCnt = 0;    //违反约束 2
    size_t outOfRowCnt = 0;    //违反约束 3
    size_t mmOverlapCells = 0; //参与可动-可动重叠的单元数
    double mmOverlapArea = 0;  //可动-可动重叠面积
    size_t mfOverlapCells = 0; //与固定宏重叠的可动单元数
    size_t mfOverlapPairs = 0; //可动-固定重叠对数
    double mfOverlapArea = 0;  //可动-固定重叠面积
    std::vector<Module *> mfViolators; //与固定宏重叠的单元清单（供可视化）

    bool legal() const
    {
        return offRowCnt == 0 && misalignCnt == 0 && outOfRowCnt == 0 &&
               mmOverlapCells == 0 && mfOverlapCells == 0;
    }
};

class LegalChecker
{
public:
    //对 db->Nodes 当前坐标做合法化检查；tag 用于输出标识
    //复杂度 O(n log n)：按行分组 + 组内排序扫描
    static LegalReport check(PlaceData *db, const std::string &tag, bool verbose = true);
};
