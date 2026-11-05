#pragma once

#include "myplace.h"
#include <cstddef>
#include <string>
#include <vector>

namespace initial_place {

// 在布局前保存一份位置和属性，布局后才能知道固定单元有没有被改动。
// 这里只保存指针，不拥有单元；使用快照期间，原单元必须仍然有效。
struct ModuleSnapshot {
    const Module* module = nullptr;
    std::string name;
    POS_2D center;
    float width = 0;
    float height = 0;
    float orientation = 0;
    bool fixed = false;
};

using PlacementSnapshot = std::vector<ModuleSnapshot>;

struct PlacementReport {
    std::size_t totalModules = 0;
    std::size_t movableModules = 0;
    std::size_t fixedModules = 0;
    std::size_t outOfBounds = 0;
    std::size_t fixedChanged = 0; // 布局前固定的单元，其位置或保存的属性发生变化。
    std::size_t invalidModules = 0; // 每个单元只计一次：几何/名称无效或属性被改动。
    double hpwl = 0;
};

// 空/重复单元指针会抛出 std::invalid_argument；快照不修改数据库。
PlacementSnapshot capturePlacement(const PlaceData& db);
// 快照必须恰好覆盖当前全部单元。集合错误和连通关系错误抛出 invalid_argument。
// 几何问题计入 invalidModules；中心非有限时 hpwl 为 NaN，调用方必须拒绝保存。
// 此报告不要求初始布局无重叠，也不检查 SiteRow 对齐或执行合法化。
PlacementReport inspectPlacement(const PlaceData& db, const CoreRegion& region,
                                 const PlacementSnapshot& before);
// 所有 pin 按 center+offset 求各网包围框；空/单 pin 网贡献为 0。
// 同时检查数据库正反向归属；坏指针关系或非有限坐标/偏移抛出 invalid_argument。
double calculateHpwl(const PlaceData& db);

} // namespace initial_place
