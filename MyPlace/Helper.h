#pragma once
#include "objects.h"

// 计算两个轴对齐矩形重叠面积
inline double rectOverlapArea(double ax0, double ay0, double ax1, double ay1,
    double bx0, double by0, double bx1, double by1)
{
    const double x0 = (ax0 > bx0) ? ax0 : bx0;
    const double y0 = (ay0 > by0) ? ay0 : by0;
    const double x1 = (ax1 < bx1) ? ax1 : bx1;
    const double y1 = (ay1 < by1) ? ay1 : by1;
    if (x1 <= x0 || y1 <= y0)
        return 0.0;
    return (x1 - x0) * (y1 - y0);
}

// 获取module包围盒：left bottom right top
inline void getModuleBBox(const Module* m, double& lx, double& ly, double& rx, double& ry)
{
    const double hw = m->width * 0.5;
    const double hh = m->height * 0.5;
    lx = m->center.x - hw;
    ly = m->center.y - hh;
    rx = m->center.x + hw;
    ry = m->center.y + hh;
}
