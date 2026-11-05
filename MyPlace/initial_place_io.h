#pragma once

#include "myplace.h"
#include <string>

namespace initial_place {
inline constexpr const char* kParserRevision = "7312fdf17e44c96b75a8ad1f408de0478b6ab847";

// 负责一次完整读入，并拥有读取器 new 出来的单元、网络和管脚。
// db 可交给布局器使用，但不要往它的指针列表中插入外部对象。
struct LoadedCircuit {
    PlaceData db{};
    CoreRegion region{};
    LoadedCircuit() = default;
    ~LoadedCircuit();
    LoadedCircuit(const LoadedCircuit&) = delete;
    LoadedCircuit& operator=(const LoadedCircuit&) = delete;

    // directory 和 benchmark 都采用 UTF-8；Windows 的旧读取器路径在内部转换。
    // 单元坐标、尺寸和管脚偏移沿用数据库的 float 精度，不保留任意十进制原文精度。
    // 固定坐标还须能经过“左下角 -> float中心 -> 导出左下角”精确往返，否则拒绝。
    // 成功或失败后都不能再次 load；重新运行请新建 LoadedCircuit。
    // 这保证左下角转中心的操作只执行一次。失败时数据库保持为空。
    void load(const std::string& directory, const std::string& benchmark);

private:
    bool attempted_ = false;
    bool ownsDatabase_ = false;
    void clear() noexcept;
};
}
