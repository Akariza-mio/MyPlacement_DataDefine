#include "placement_report.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace initial_place {
namespace {

// 指针先排序，再检查重复；排序后的数组也可用来快速查找归属。
// 使用 std::less 比较不同对象的地址，不依赖对象在内存中的排列方式。
template<class T>
void sortAndCheckPointers(std::vector<const T*>& pointers, const char* source) {
    for (const T* pointer : pointers) {
        if (pointer == nullptr) {
            throw std::invalid_argument(std::string("Null pointer in ") + source);
        }
    }
    std::sort(pointers.begin(), pointers.end(), std::less<const T*>());
    if (std::adjacent_find(pointers.begin(), pointers.end()) != pointers.end()) {
        throw std::invalid_argument(std::string("Duplicate pointer in ") + source);
    }
}

template<class T>
std::vector<const T*> checkedPointers(const std::vector<T*>& source,
                                      const char* description) {
    std::vector<const T*> result(source.begin(), source.end());
    sortAndCheckPointers(result, description);
    return result;
}

std::vector<const Module*> collectModules(const PlaceData& db) {
    std::vector<const Module*> modules;
    modules.reserve(db.Nodes.size() + db.Terminals.size());
    modules.insert(modules.end(), db.Nodes.begin(), db.Nodes.end());
    modules.insert(modules.end(), db.Terminals.begin(), db.Terminals.end());
    sortAndCheckPointers(modules, "PlaceData::Nodes/Terminals");
    return modules;
}

template<class T>
bool containsPointer(const std::vector<const T*>& sorted, const T* pointer) {
    return std::binary_search(sorted.begin(), sorted.end(), pointer,
                              std::less<const T*>());
}

// 不能直接使用 Pin::idx：手工图可能没有设置 idx，外部算法也可能重新编号。
std::size_t findPin(const std::vector<const Pin*>& pins, const Pin* pin,
                    const char* source) {
    if (pin == nullptr) {
        throw std::invalid_argument(std::string("Null pin in ") + source);
    }
    const auto position = std::lower_bound(pins.begin(), pins.end(), pin,
                                            std::less<const Pin*>());
    if (position == pins.end() || *position != pin) {
        throw std::invalid_argument(std::string("Pin is absent from PlaceData::Pins: ")
                                    + source);
    }
    return static_cast<std::size_t>(position - pins.begin());
}

bool finiteCenter(const Module& module) {
    return std::isfinite(module.center.x) && std::isfinite(module.center.y);
}

bool validGeometry(const Module& module) {
    if (!finiteCenter(module) || !std::isfinite(module.width)
        || !std::isfinite(module.height) || !std::isfinite(module.orientation)) {
        return false;
    }
    // 固定 IO 可以是一个零面积点；可移动单元必须有正宽度和正高度。
    return module.isFixed ? module.width >= 0 && module.height >= 0
                          : module.width > 0 && module.height > 0;
}

bool validName(const std::string& name) {
    return !name.empty() && name.front() != '#'
        && name.find_first_of(" \t\r\n\f\v") == std::string::npos
        && name.find('\0') == std::string::npos;
}

void sortUniqueNets(std::vector<const Net*>& nets) {
    std::sort(nets.begin(), nets.end(), std::less<const Net*>());
    nets.erase(std::unique(nets.begin(), nets.end()), nets.end());
}

double inspectHpwl(const PlaceData& db, const std::vector<const Module*>& modules,
                   bool allowInvalidModuleCenter) {
    const auto nets = checkedPointers(db.Nets, "PlaceData::Nets");
    const auto pins = checkedPointers(db.Pins, "PlaceData::Pins");
    bool coordinatesValid = true;
    for (const Module* module : modules) {
        if (!finiteCenter(*module)) {
            coordinatesValid = false;
            if (!allowInvalidModuleCenter) {
                throw std::invalid_argument("Nonfinite module center in HPWL: " + module->name);
            }
        }
    }

    // 先确认 pin 指向本数据库里的对象，再读取这些对象的坐标。
    for (const Pin* pin : pins) {
        if (pin->module == nullptr || !containsPointer(modules, pin->module)) {
            throw std::invalid_argument("Pin module is absent from PlaceData::Nodes/Terminals");
        }
        if (pin->net == nullptr || !containsPointer(nets, pin->net)) {
            throw std::invalid_argument("Pin net is absent from PlaceData::Nets");
        }
        if (!std::isfinite(pin->offset.x) || !std::isfinite(pin->offset.y)) {
            throw std::invalid_argument("Nonfinite pin offset for module: " + pin->module->name);
        }
    }

    // 每个 pin 只用一个字节记录两件事：是否在所属 net 和 module 中各出现一次。
    // 不为每个网或每个单元新建大集合，额外内存随总 pin 数线性增长。
    constexpr unsigned char seenInNet = 1;
    constexpr unsigned char seenInModule = 2;
    std::vector<unsigned char> seen(pins.size(), 0);
    double hpwl = 0;
    for (const Net* net : nets) {
        double xmin = std::numeric_limits<double>::infinity();
        double ymin = std::numeric_limits<double>::infinity();
        double xmax = -std::numeric_limits<double>::infinity();
        double ymax = -std::numeric_limits<double>::infinity();
        for (const Pin* pin : net->netPins) {
            const std::size_t index = findPin(pins, pin, "Net::netPins");
            if (pin->net != net) {
                throw std::invalid_argument("Pin::net does not match its containing Net::netPins");
            }
            if ((seen[index] & seenInNet) != 0) {
                throw std::invalid_argument("Duplicate pin pointer in Net::netPins");
            }
            seen[index] |= seenInNet;

            // 管脚位置是单元中心加管脚偏移，不是单元左下角，也不能只用中心代替。
            const double x = static_cast<double>(pin->module->center.x) + pin->offset.x;
            const double y = static_cast<double>(pin->module->center.y) + pin->offset.y;
            if (std::isfinite(x) && std::isfinite(y)) {
                xmin = std::min(xmin, x); xmax = std::max(xmax, x);
                ymin = std::min(ymin, y); ymax = std::max(ymax, y);
            }
        }
        // 空网和单管脚网都没有需要连接的两个点，贡献为 0。
        if (net->netPins.size() > 1 && coordinatesValid) {
            hpwl += (xmax - xmin) + (ymax - ymin);
        }
    }

    // 复用两个临时数组核对单元的反向关联，容量只随最大的单元度数增长。
    std::vector<const Net*> expectedNets;
    std::vector<const Net*> actualNets;
    for (const Module* module : modules) {
        expectedNets.clear();
        actualNets.clear();
        for (const Pin* pin : module->modulePins) {
            const std::size_t index = findPin(pins, pin, "Module::modulePins");
            if (pin->module != module) {
                throw std::invalid_argument("Pin::module does not match Module::modulePins: " + module->name);
            }
            if ((seen[index] & seenInModule) != 0) {
                throw std::invalid_argument("Duplicate pin pointer in Module::modulePins: " + module->name);
            }
            seen[index] |= seenInModule;
            expectedNets.push_back(pin->net);
        }
        for (const Net* net : module->nets) {
            if (net == nullptr) {
                throw std::invalid_argument("Null net in Module::nets: " + module->name);
            }
            actualNets.push_back(net);
        }
        // 同一个单元可在一条网中有多个不同 pin，因此这里只对网引用去重。
        // 真实 pin 仍全部保留，HPWL 的包围框已经逐 pin 计算过。
        sortUniqueNets(expectedNets);
        sortUniqueNets(actualNets);
        if (expectedNets != actualNets) {
            throw std::invalid_argument("Module::nets does not match Module::modulePins: " + module->name);
        }
    }
    for (unsigned char flags : seen) {
        if (flags != (seenInNet | seenInModule)) {
            throw std::invalid_argument("A PlaceData::Pins entry is missing from its net or module");
        }
    }
    return coordinatesValid ? hpwl : std::numeric_limits<double>::quiet_NaN();
}

} // namespace

PlacementSnapshot capturePlacement(const PlaceData& db) {
    // 先排除空指针和重复单元，防止保存一份有歧义的快照。
    const auto modules = collectModules(db);
    PlacementSnapshot snapshot;
    snapshot.reserve(modules.size());
    for (const Module* module : modules) {
        snapshot.push_back(ModuleSnapshot{module, module->name, module->center,
                                         module->width, module->height,
                                         module->orientation, module->isFixed});
    }
    return snapshot;
}

PlacementReport inspectPlacement(const PlaceData& db, const CoreRegion& region,
                                 const PlacementSnapshot& before) {
    if (!std::isfinite(region.xmin) || !std::isfinite(region.ymin)
        || !std::isfinite(region.xmax) || !std::isfinite(region.ymax)
        || region.xmin >= region.xmax || region.ymin >= region.ymax) {
        throw std::invalid_argument("Placement report requires finite CoreRegion bounds and positive size");
    }
    const auto modules = collectModules(db);
    if (before.size() != modules.size()) {
        throw std::invalid_argument("PlacementSnapshot count does not match PlaceData module count");
    }
    std::vector<const ModuleSnapshot*> previous;
    previous.reserve(before.size());
    for (const ModuleSnapshot& item : before) {
        if (item.module == nullptr) {
            throw std::invalid_argument("Null module pointer in PlacementSnapshot");
        }
        previous.push_back(&item);
    }
    std::sort(previous.begin(), previous.end(), [](const ModuleSnapshot* a, const ModuleSnapshot* b) {
        return std::less<const Module*>()(a->module, b->module);
    });
    for (std::size_t i = 1; i < previous.size(); ++i) {
        if (previous[i - 1]->module == previous[i]->module) {
            throw std::invalid_argument("Duplicate module pointer in PlacementSnapshot");
        }
    }
    for (std::size_t i = 0; i < previous.size(); ++i) {
        if (previous[i]->module != modules[i]) {
            throw std::invalid_argument("PlacementSnapshot module set does not match PlaceData");
        }
    }

    PlacementReport report;
    report.totalModules = modules.size();
    for (std::size_t i = 0; i < modules.size(); ++i) {
        const Module& module = *modules[i];
        const ModuleSnapshot& old = *previous[i];
        if (module.isFixed) ++report.fixedModules;
        else ++report.movableModules;

        const bool metadataChanged = module.name != old.name || module.width != old.width
            || module.height != old.height || module.orientation != old.orientation
            || module.isFixed != old.fixed;
        const bool geometryValid = validGeometry(module);
        if (!geometryValid || !validName(module.name) || metadataChanged) {
            ++report.invalidModules;
        }
        // 固定单元要求完全相等，不用容差掩盖细小移动；固定单元可位于 core 外。
        if (old.fixed && (metadataChanged || module.center.x != old.center.x
                         || module.center.y != old.center.y)) {
            ++report.fixedChanged;
        }
        if (!module.isFixed && geometryValid) {
            // 半宽和四条边都用 double 计算，避免 float 加法把微小越界舍入掉。
            const double left = static_cast<double>(module.center.x) - module.width / 2.0;
            const double right = static_cast<double>(module.center.x) + module.width / 2.0;
            const double bottom = static_cast<double>(module.center.y) - module.height / 2.0;
            const double top = static_cast<double>(module.center.y) + module.height / 2.0;
            if (left < region.xmin || right > region.xmax
                || bottom < region.ymin || top > region.ymax) {
                ++report.outOfBounds;
            }
        }
    }
    report.hpwl = inspectHpwl(db, modules, true);
    return report;
}

double calculateHpwl(const PlaceData& db) {
    return inspectHpwl(db, collectModules(db), false);
}

} // namespace initial_place
