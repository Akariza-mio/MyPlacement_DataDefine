#pragma once
// ============================================================================
// Helper.h：布局流程通用的几何/数学小工具（header-only）
// ============================================================================

#include "objects.h"

#include <algorithm>
#include <cmath>

// ---------------------------------------------------------------------------
// 几何：轴对齐矩形
// ---------------------------------------------------------------------------

// 计算两个轴对齐矩形 [ax0,ax1]×[ay0,ay1] 与 [bx0,bx1]×[by0,by1] 的重叠面积。
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

// 获取 module 的包围盒：left bottom right top
inline void getModuleBBox(const Module* m, double& lx, double& ly, double& rx, double& ry)
{
    const double hw = static_cast<double>(m->width) * 0.5;
    const double hh = static_cast<double>(m->height) * 0.5;
    lx = static_cast<double>(m->center.x) - hw;
    ly = static_cast<double>(m->center.y) - hh;
    rx = static_cast<double>(m->center.x) + hw;
    ry = static_cast<double>(m->center.y) + hh;
}

// 两个一维区间的重叠长度（负值截断为 0）
inline double overlapLength1D(double a0, double a1, double b0, double b1)
{
    const double lo = (a0 > b0) ? a0 : b0;
    const double hi = (a1 < b1) ? a1 : b1;
    return (hi > lo) ? (hi - lo) : 0.0;
}

// ---------------------------------------------------------------------------
// 数值：截断、安全除法、安全指数
// ---------------------------------------------------------------------------

inline double clampValue(double v, double lo, double hi)
{
    return (v < lo) ? lo : ((v > hi) ? hi : v);
}

inline int clampIndex(int v, int lo, int hi)
{
    return (v < lo) ? lo : ((v > hi) ? hi : v);
}

// 分母接近 0 时返回 fallback，避免 inf/nan 污染整条数据流
inline double safeDivide(double numerator, double denominator, double fallback = 0.0)
{
    if (std::fabs(denominator) > 1e-30 && std::isfinite(numerator) && std::isfinite(denominator))
        return numerator / denominator;
    return fallback;
}

// 带上下截断的 exp：坐标量级远大于平滑参数 γ 时也不会溢出
inline double safeExp(double z)
{
    return std::exp(clampValue(z, -700.0, 700.0));
}

// ---------------------------------------------------------------------------
// 网格场：双线性插值
// ---------------------------------------------------------------------------
// grid 为 m×m 的场，按 grid[i*m+j] 存放（i 为 x 方向下标，j 为 y 方向下标）；
// (bx,by) 为连续的网格坐标，网格点 (i,j) 位于 (i,j)。
inline double bilinearSample(const std::vector<double>& grid, int m, double bx, double by)
{
    if (m <= 0 || grid.size() < static_cast<std::size_t>(m) * static_cast<std::size_t>(m))
        return 0.0;
    const double cx = clampValue(bx, 0.0, static_cast<double>(m - 1));
    const double cy = clampValue(by, 0.0, static_cast<double>(m - 1));
    const int i0 = static_cast<int>(std::floor(cx));
    const int j0 = static_cast<int>(std::floor(cy));
    const int i1 = (i0 + 1 < m) ? (i0 + 1) : i0;
    const int j1 = (j0 + 1 < m) ? (j0 + 1) : j0;
    const double fx = cx - i0;
    const double fy = cy - j0;

    const double v00 = grid[static_cast<std::size_t>(i0) * m + j0];
    const double v10 = grid[static_cast<std::size_t>(i1) * m + j0];
    const double v01 = grid[static_cast<std::size_t>(i0) * m + j1];
    const double v11 = grid[static_cast<std::size_t>(i1) * m + j1];
    const double bottom = v00 + (v10 - v00) * fx;
    const double top = v01 + (v11 - v01) * fx;
    return bottom + (top - bottom) * fy;
}
