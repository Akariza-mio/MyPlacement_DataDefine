#include "myplace.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
// 将一个方向的中心允许区间转换为 float 能准确表示的安全端点。
// lower：中心坐标允许的最小值；upper：中心坐标允许的最大值，二者均为 double。
// moduleName：单元名称，仅在报错时帮助定位出问题的单元。
// 返回 pair：first 是安全下限，second 是安全上限；它们不是 x/y 两个坐标。
// center 存在 float 中，端点向区间内部取整，避免 double 转 float 后轻微越界。
std::pair<float, float> centerRange(double lower, double upper,
                                   const std::string& moduleName)
{
    float first = static_cast<float>(lower); // 下限转换成 float 后的值。
    float last = static_cast<float>(upper);  // 上限转换成 float 后的值。
    // nextafter 将浮点数向指定方向移动一个可表示的步长，不是简单地加减 1。
    if (static_cast<double>(first) < lower) {
        first = std::nextafter(first, std::numeric_limits<float>::infinity());
    }
    if (static_cast<double>(last) > upper) {
        last = std::nextafter(last, -std::numeric_limits<float>::infinity());
    }
    if (first > last) {
        throw std::invalid_argument("No representable center for module: " + moduleName);
    }
    return std::make_pair(first, last);
}

// 暂存一个单元及其中心允许范围，所有检查通过后再统一写回位置。
struct PlacementRange {
    Module* module;             // 指向原单元；通过它修改 center 会更新原数据。
    std::pair<float, float> x;  // 横坐标区间：[x.first, x.second]。
    std::pair<float, float> y;  // 纵坐标区间：[y.first, y.second]。
};
}

// region 是区域边界，seed 是随机种子；具体参数约定见 myplace.h 的声明。
void MyPlacer::randomPlace(const CoreRegion& region, unsigned int seed)
{
    if (db == nullptr) {
        throw std::invalid_argument("randomPlace requires a PlaceData database");
    }
    // isfinite 排除 NaN（无效数字）和正负无穷大，避免后续范围计算失效。
    if (!std::isfinite(region.xmin) || !std::isfinite(region.ymin)
        || !std::isfinite(region.xmax) || !std::isfinite(region.ymax)
        || region.xmin >= region.xmax || region.ymin >= region.ymax) {
        throw std::invalid_argument("CoreRegion must have finite bounds and positive size");
    }

    const double coreWidth = static_cast<double>(region.xmax) - region.xmin;   // 区域总宽度。
    const double coreHeight = static_cast<double>(region.ymax) - region.ymin; // 区域总高度。
    std::vector<PlacementRange> ranges; // 收集每个可移动单元及其安全坐标范围。
    std::unordered_set<const Module*> seen; // 同一个单元地址只能出现一次，包括固定单元。

    // 先检查并计算每个可移动单元的范围；这里不检查单元之间的重叠。
    // group 依次指向 Nodes、Terminals 两个列表；是否移动始终由 isFixed 决定。
    for (const std::vector<Module*>* group : {&db->Nodes, &db->Terminals}) {
        // module 是当前单元的地址；module->width 表示读取这个单元的宽度。
        for (Module* module : *group) {
            if (module == nullptr) {
                throw std::invalid_argument("Null module in PlaceData");
            }
            if (!seen.insert(module).second) {
                throw std::invalid_argument("Duplicate module pointer in PlaceData: " + module->name);
            }
            if (module->isFixed) {
                continue;
            }
            if (!std::isfinite(module->width) || !std::isfinite(module->height)
                || module->width <= 0 || module->height <= 0) {
                throw std::invalid_argument("Invalid dimensions for module: " + module->name);
            }
            if (module->width > coreWidth || module->height > coreHeight) {
                throw std::invalid_argument("Module does not fit CoreRegion: " + module->name);
            }

            // 中心到左右边界各留半个宽度，到上下边界各留半个高度。
            const auto x = centerRange(region.xmin + module->width / 2.0,
                                       region.xmax - module->width / 2.0, module->name);
            const auto y = centerRange(region.ymin + module->height / 2.0,
                                       region.ymax - module->height / 2.0, module->name);
            ranges.push_back(PlacementRange{module, x, y}); // 记录范围，暂不改变位置。
        }
    }

    // 生成器只初始化一次。原地更新中心，尺寸、方向等属性保持不变。
    std::mt19937 rng(seed); // rng：伪随机数生成器；seed 决定这次序列的起点。
    for (const PlacementRange& range : ranges) {
        // range 是当前单元的范围记录；randomX/randomY 分别规定横、纵坐标取值范围。
        std::uniform_real_distribution<double> randomX(range.x.first, range.x.second);
        std::uniform_real_distribution<double> randomY(range.y.first, range.y.second);
        // randomX(rng) 从生成器取数；转换为项目坐标所用的 float 后写回原单元。
        range.module->center.x = static_cast<float>(randomX(rng));
        range.module->center.y = static_cast<float>(randomY(rng));
    }
}

// filename 是输出路径；这里保存现有位置，不调用 randomPlace，也不改变任何单元。
void MyPlacer::writePl(const std::string& filename) const
{
    if (db == nullptr) {
        throw std::invalid_argument("writePl requires a PlaceData database");
    }
    if (filename.empty() || filename.find('\0') != std::string::npos) {
        throw std::invalid_argument("Output filename must be nonempty and contain no null character");
    }

    // 第一遍只检查数据。ofstream 打开已有文件时会清空旧内容，所以先完成检查。
    // names 用于防止同名单元被重复写出；它不改变 Nodes、Terminals 的原有顺序。
    std::unordered_set<std::string> names;
    for (const std::vector<Module*>* group : {&db->Nodes, &db->Terminals}) {
        for (const Module* module : *group) {
            // const Module*：本函数通过这个指针只能读取单元，不能修改单元属性。
            if (module == nullptr) {
                throw std::invalid_argument("Null module in PlaceData");
            }
            // .pl 的名称占一个文本字段，不能含空白或换行，也不能以注释符 # 开头。
            const std::string& name = module->name;
            if (name.empty() || name.front() == '#'
                || name.find_first_of(" \t\r\n\f\v") != std::string::npos
                || name.find('\0') != std::string::npos) {
                throw std::invalid_argument("Invalid module name in .pl output");
            }
            if (!names.insert(name).second) {
                throw std::invalid_argument("Duplicate module name in .pl output: " + name);
            }
            // 不检查重叠或区域边界；导出只保存当前位置。零尺寸的固定 IO 也可以保存。
            if (!std::isfinite(module->center.x) || !std::isfinite(module->center.y)
                || !std::isfinite(module->width) || !std::isfinite(module->height)
                || module->width < 0 || module->height < 0) {
                throw std::invalid_argument("Invalid position or dimensions for module: " + name);
            }
            // 现有模板只有 float orientation，没有其他方向的编码约定。
            // 本阶段明确约定 0 表示 N；不能把未知方向悄悄写成 N，导致数据失真。
            if (module->orientation != 0.0f) {
                throw std::invalid_argument("Unsupported orientation for module: " + name
                                            + " (only 0 = N is currently defined)");
            }
        }
    }

    // output 是文件输出流。trunc 表示重新写入时覆盖已有内容；父目录需事先存在。
    std::ofstream output(filename, std::ios::out | std::ios::trunc);
    if (!output.is_open()) {
        throw std::runtime_error("Cannot open .pl output file: " + filename);
    }
    // 固定使用小数点，避免系统地区设置把坐标写成带逗号的小数。
    output.imbue(std::locale::classic());
    // 坐标转换用 double 计算，保留足够有效数字，减少导出后重新读取的精度损失。
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << "UCLA pl 1.0\n\n"; // Bookshelf .pl 的文件头。

    // 第二遍输出：先 Nodes，再 Terminals；是否固定只看 isFixed，不按列表名称推断。
    for (const std::vector<Module*>* group : {&db->Nodes, &db->Terminals}) {
        for (const Module* module : *group) {
            // 内存中保存的是中心，.pl 要求左下角；只计算局部变量，不回写 center。
            const double leftX = static_cast<double>(module->center.x) - module->width / 2.0;
            const double bottomY = static_cast<double>(module->center.y) - module->height / 2.0;

            // 一行格式：单元名 左下角x 左下角y : 方向 [固定标记]。
            output << module->name << ' ' << leftX << ' ' << bottomY << " : N";
            if (module->isFixed) {
                output << " /FIXED"; // 固定单元同样输出，但保留固定标记和原位置。
            }
            output << '\n';
        }
    }

    // close 会把缓冲区内容写入文件；磁盘空间不足等错误可能到这一步才出现。
    output.close();
    if (!output) {
        throw std::runtime_error("Failed to write or close .pl output file: " + filename);
    }
}
