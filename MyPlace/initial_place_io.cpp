#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

#include "initial_place_io.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace initial_place {
namespace {
using Tokens = std::vector<std::string>;

// 旧读取器默认文件完全正确。先逐行检查，避免它使用没读到的宽高等数值。
class CheckedFile {
public:
    explicit CheckedFile(const std::string& filename)
        : filename_(filename), input_(std::filesystem::u8path(filename)) {
        if (!input_) throw std::runtime_error("Cannot open input file: " + filename);
    }
    [[noreturn]] void fail(const std::string& reason) const {
        throw std::runtime_error(filename_ + ":" + std::to_string(line_) + ": " + reason);
    }
    bool next(Tokens& tokens) {
        std::string line;
        while (std::getline(input_, line)) {
            ++line_;
            if (line.empty() || line.front() == '#') continue;
            tokens.clear();
            std::istringstream fields(line);
            for (std::string word; fields >> word;) tokens.push_back(std::move(word));
            // 不改写原文件；明确拒绝旧读取器不能安全处理的写法。
            if (tokens.empty()) fail("Whitespace-only lines are unsupported by the reused parser");
            if (tokens[0].front() == '#') fail("Indented comments are unsupported; put # in column 1");
            if (tokens[0] == "UCLA" && line.front() != 'U') fail("The UCLA header must start in column 1");
            return true;
        }
        if (input_.bad()) fail("Error while reading input file");
        return false;
    }
    void header(const char* kind) {
        Tokens words;
        if (!next(words) || words != Tokens{"UCLA", kind, "1.0"}) fail("Expected UCLA " + std::string(kind) + " 1.0 header");
    }
private:
    std::string filename_;
    std::ifstream input_;
    std::size_t line_ = 0;
};

void shape(CheckedFile& file, const Tokens& words, std::size_t count, const char* description) {
    if (words.size() != count) file.fail(std::string("Invalid fields in ") + description);
}

int integer(CheckedFile& file, const std::string& text, const char* description) {
    int value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() || value < 0)
        file.fail(std::string("Expected a nonnegative integer for ") + description);
    return value;
}

template<class Real>
Real real(CheckedFile& file, const std::string& text, const char* description) {
    Real value = 0;
    // 与队友代码使用相同的流式数值读取规则；先确认每个数都能成功提取。
    // 这样极小数值或特殊写法也不会让旧读取器中断后留下未初始化字段。
    std::istringstream field(text);
    if (!(field >> value) || field.peek() != std::char_traits<char>::eof() || !std::isfinite(value))
        file.fail(std::string("Expected a finite number for ") + description);
    return value;
}

int countHeader(CheckedFile& file, const Tokens& words, const char* description) {
    shape(file, words, 3, description);
    if (words[1] != ":") file.fail("Expected ':' in count header");
    return integer(file, words[2], description);
}

void safeName(CheckedFile& file, const std::string& name) {
    if (name.empty() || name.front() == 'U' || name.front() == '#' || name.find('\0') != std::string::npos
        || name == "NumNodes" || name == "NumTerminals" || name == "NumNets"
        || name == "NumPins" || name == "NetDegree")
        file.fail("Unsupported name (reserved words and U-prefixed names are unsafe for the reused parser): " + name);
}

struct NodeInfo {
    float width = 0;
    float height = 0;
    bool fixed = false;
    bool positionSeen = false;
};
using NodeMap = std::unordered_map<std::string, NodeInfo>;

struct Bounds {
    double xmin = std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    void include(double left, double bottom, double right, double top) {
        xmin = std::min(xmin, left); ymin = std::min(ymin, bottom);
        xmax = std::max(xmax, right); ymax = std::max(ymax, top);
    }
};

struct Summary {
    int modules = 0;
    int fixed = 0;
    int nets = 0;
    int pins = 0;
    int rows = 0;
    Bounds bounds;
};

NodeMap inspectNodes(const std::string& filename, Summary& summary) {
    CheckedFile file(filename); file.header("nodes");
    NodeMap nodes;
    bool totalSeen = false, fixedSeen = false, recordsStarted = false;
    std::size_t fixedCount = 0;
    Tokens words;
    while (file.next(words)) {
        if (words[0] == "NumNodes" || words[0] == "NumTerminals") {
            if (recordsStarted) file.fail("Node count headers must precede node records");
            const bool total = words[0] == "NumNodes";
            bool& seen = total ? totalSeen : fixedSeen;
            if (seen) file.fail("Duplicate node count header");
            seen = true;
            (total ? summary.modules : summary.fixed) = countHeader(file, words, words[0].c_str());
            continue;
        }
        if (!totalSeen || !fixedSeen) file.fail("Missing NumNodes or NumTerminals before node records");
        recordsStarted = true;
        if (words.size() != 3 && words.size() != 4) file.fail("Expected: node width height [terminal]");
        safeName(file, words[0]);
        NodeInfo node;
        node.width = real<float>(file, words[1], "node width");
        node.height = real<float>(file, words[2], "node height");
        node.fixed = words.size() == 4;
        if (node.fixed && words[3] != "terminal") file.fail("Only the terminal node classification is supported");
        if (node.width < 0 || node.height < 0 || (!node.fixed && (node.width == 0 || node.height == 0)))
            file.fail("Movable nodes need positive sizes; fixed IO sizes may be zero");
        const double area = static_cast<double>(node.width) * node.height;
        const float storedArea = static_cast<float>(area);
        if (!std::isfinite(storedArea) || (area > 0 && storedArea == 0)) file.fail("Node area cannot be represented as float");
        if (!nodes.emplace(words[0], node).second) file.fail("Duplicate node name: " + words[0]);
        fixedCount += node.fixed;
        if (nodes.size() > static_cast<std::size_t>(summary.modules)) file.fail("More node records than NumNodes");
    }
    if (!totalSeen || !fixedSeen || nodes.size() != static_cast<std::size_t>(summary.modules)
        || fixedCount != static_cast<std::size_t>(summary.fixed)) file.fail("Node or terminal count mismatch");
    return nodes;
}

void inspectPositions(const std::string& filename, NodeMap& nodes) {
    CheckedFile file(filename); file.header("pl");
    Tokens words;
    std::size_t count = 0;
    while (file.next(words)) {
        if (words.size() != 5 && words.size() != 6) file.fail("Expected: node x y : N [/FIXED]");
        safeName(file, words[0]);
        auto found = nodes.find(words[0]);
        if (found == nodes.end()) file.fail("Unknown node in .pl: " + words[0]);
        NodeInfo& node = found->second;
        if (node.positionSeen) file.fail("Duplicate node position: " + words[0]);
        if (words[3] != ":") file.fail("Expected ':' before node orientation");
        if (words[4] != "N") file.fail("Unsupported orientation for " + words[0] + "; only N is supported");
        const bool fixed = words.size() == 6;
        if (fixed && words[5] != "/FIXED") file.fail("Unsupported placement flag: " + words[5]);
        if (fixed != node.fixed) file.fail("Fixed flag disagrees with .nodes classification: " + words[0]);
        const float x = real<float>(file, words[1], "placement x");
        const float y = real<float>(file, words[2], "placement y");
        const float centerX = static_cast<float>(static_cast<double>(x) + node.width / 2.0);
        const float centerY = static_cast<float>(static_cast<double>(y) + node.height / 2.0);
        if (!std::isfinite(centerX) || !std::isfinite(centerY))
            file.fail("Node center is outside the finite float range: " + words[0]);
        // 固定位置不仅要在算法中保持中心不变，保存时也必须回到原左下角。
        // 坐标过大、单元过小时，加半宽/半高可能被 float 舍掉；提前拒绝。
        if (fixed && (static_cast<double>(centerX) - node.width / 2.0 != static_cast<double>(x)
                   || static_cast<double>(centerY) - node.height / 2.0 != static_cast<double>(y)))
            file.fail("Fixed node cannot preserve its lower-left coordinates through float center conversion: " + words[0]);
        node.positionSeen = true;
        ++count;
    }
    if (count != nodes.size()) file.fail("Missing node positions in .pl");
}

void inspectRows(const std::string& filename, Summary& summary) {
    CheckedFile file(filename); file.header("scl");
    Tokens words;
    bool countSeen = false, inRow = false;
    int actualRows = 0;
    unsigned fieldsSeen = 0;
    double bottom = 0, height = 0, width = 0, spacing = 0, origin = 0;
    int sites = 0;
    while (file.next(words)) {
        const std::string& tag = words[0];
        if (tag == "NumRows") {
            if (countSeen || inRow || actualRows != 0) file.fail("Duplicate or misplaced NumRows");
            summary.rows = countHeader(file, words, "NumRows"); countSeen = true;
        } else if (tag == "CoreRow") {
            if (!countSeen || inRow) file.fail("Missing NumRows or nested CoreRow");
            shape(file, words, 2, "CoreRow");
            if (words[1] != "Horizontal") file.fail("Only horizontal placement rows are supported");
            inRow = true; fieldsSeen = 0;
        } else if (tag == "End") {
            shape(file, words, 1, "End");
            if (!inRow || fieldsSeen != 127) file.fail("Incomplete row; all seven row fields are required");
            if (width != spacing) file.fail("Sitewidth must equal Sitespacing; other site layouts are unsupported");
            const double right = origin + static_cast<double>(sites) * spacing;
            const double top = bottom + height;
            if (!std::isfinite(right) || !std::isfinite(top) || right <= origin || top <= bottom
                || !std::isfinite(static_cast<float>(origin)) || !std::isfinite(static_cast<float>(right))
                || !std::isfinite(static_cast<float>(bottom)) || !std::isfinite(static_cast<float>(top)))
                file.fail("Row bounds are not finite, positive, representable coordinates");
            summary.bounds.include(origin, bottom, right, top);
            inRow = false;
            if (++actualRows > summary.rows) file.fail("More rows than NumRows");
        } else {
            if (!inRow) file.fail("Unexpected field outside CoreRow: " + tag);
            unsigned bit = 0;
            if (tag == "Coordinate") bit = 1;
            else if (tag == "Height") bit = 2;
            else if (tag == "Sitewidth") bit = 4;
            else if (tag == "Sitespacing") bit = 8;
            else if (tag == "Siteorient") bit = 16;
            else if (tag == "Sitesymmetry") bit = 32;
            else if (tag == "SubrowOrigin") bit = 64;
            else file.fail("Unsupported row field: " + tag);
            if (fieldsSeen & bit) file.fail("Duplicate row field or multiple row segments: " + tag);
            if (tag == "SubrowOrigin") {
                shape(file, words, 6, "SubrowOrigin");
                if ((fieldsSeen & 15) != 15 || words[1] != ":" || words[3] != "NumSites" || words[4] != ":")
                    file.fail("SubrowOrigin requires earlier coordinate, height, site width and spacing fields");
                origin = real<double>(file, words[2], "subrow origin");
                // 读取器把起点保存在 float，也要确保这个存储步骤不会溢出。
                (void)real<float>(file, words[2], "subrow origin");
                sites = integer(file, words[5], "NumSites");
                if (sites == 0) file.fail("NumSites must be positive");
            } else {
                shape(file, words, 3, "row field");
                if (words[1] != ":") file.fail("Expected ':' in row field");
                if (tag == "Coordinate") bottom = real<double>(file, words[2], "row coordinate");
                else if (tag == "Height") height = real<double>(file, words[2], "row height");
                else if (tag == "Sitewidth") width = real<double>(file, words[2], "site width");
                else if (tag == "Sitespacing") spacing = real<double>(file, words[2], "site spacing");
                // 目前只承诺课程数据中已经验证的行格式，遇到扩展明确停止。
                else if (words[2] != "1") file.fail("Only the verified Siteorient/Sitesymmetry value 1 is supported");
                if ((tag == "Height" && height <= 0) || (tag == "Sitewidth" && width <= 0)
                    || (tag == "Sitespacing" && spacing <= 0)) file.fail("Row height, site width and spacing must be positive");
            }
            fieldsSeen |= bit;
        }
    }
    if (!countSeen || inRow || actualRows != summary.rows || actualRows == 0) file.fail("Missing/incomplete rows or NumRows mismatch");
}

void inspectNets(const std::string& filename, const NodeMap& nodes, Summary& summary) {
    CheckedFile file(filename); file.header("nets");
    Tokens words;
    bool netCountSeen = false, pinCountSeen = false, recordsStarted = false;
    std::size_t nets = 0, pins = 0;
    int remaining = 0;
    std::unordered_set<std::string> names;
    while (file.next(words)) {
        const std::string& tag = words[0];
        if (tag == "NumNets" || tag == "NumPins") {
            if (recordsStarted) file.fail("Net count headers must precede connection records");
            const bool netHeader = tag == "NumNets";
            bool& seen = netHeader ? netCountSeen : pinCountSeen;
            if (seen) file.fail("Duplicate net/pin count header");
            seen = true;
            (netHeader ? summary.nets : summary.pins) = countHeader(file, words, tag.c_str());
        } else if (tag == "NetDegree") {
            if (!netCountSeen || !pinCountSeen || remaining != 0) file.fail("Incomplete previous net or missing net/pin counts");
            shape(file, words, 4, "NetDegree");
            if (words[1] != ":") file.fail("Expected ':' after NetDegree");
            remaining = integer(file, words[2], "NetDegree");
            safeName(file, words[3]);
            if (!names.insert(words[3]).second) file.fail("Duplicate net name: " + words[3]);
            recordsStarted = true;
            if (++nets > static_cast<std::size_t>(summary.nets)) file.fail("More networks than NumNets");
            if (static_cast<std::size_t>(remaining) > static_cast<std::size_t>(summary.pins) - pins)
                file.fail("NetDegree exceeds the remaining NumPins total");
        } else {
            if (!recordsStarted || remaining <= 0) file.fail("Pin record outside a declared NetDegree");
            shape(file, words, 5, "pin");
            safeName(file, tag);
            if (nodes.find(tag) == nodes.end()) file.fail("Unknown module in net: " + tag);
            if ((words[1] != "I" && words[1] != "O" && words[1] != "B") || words[2] != ":")
                file.fail("Expected pin direction I/O/B and ':' before offsets");
            (void)real<float>(file, words[3], "pin x offset");
            (void)real<float>(file, words[4], "pin y offset");
            --remaining;
            if (++pins > static_cast<std::size_t>(summary.pins)) file.fail("More pins than NumPins");
        }
    }
    if (!netCountSeen || !pinCountSeen || remaining != 0 || nets != static_cast<std::size_t>(summary.nets)
        || pins != static_cast<std::size_t>(summary.pins)) file.fail("Network/pin count or NetDegree mismatch");
}

Summary inspect(const std::string& prefix) {
    // 只存核对名称所需的小索引，不保存另一份百万管脚数据库。
    // 函数结束后索引被释放，再让队友的读取器建立正式对象。
    Summary summary;
    auto nodes = inspectNodes(prefix + ".nodes", summary);
    inspectPositions(prefix + ".pl", nodes);
    inspectRows(prefix + ".scl", summary);
    inspectNets(prefix + ".nets", nodes, summary);
    return summary;
}

std::string parserPath(const std::string& utf8) {
#ifdef _WIN32
    // 正式接口统一用 UTF-8；原读取器的窄字符 ifstream 在 Windows 使用系统代码页。
    // 只在最后调用原读取器时转换路径，不允许把不能表示的字悄悄换成问号。
    if (utf8.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Input path is too long");
    const int size = static_cast<int>(utf8.size());
    const int wideSize = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, utf8.data(), size, nullptr, 0);
    if (!wideSize) throw std::invalid_argument("Input path is not valid UTF-8");
    std::wstring wide(static_cast<std::size_t>(wideSize), L'\0');
    if (!MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, utf8.data(), size, wide.data(), wideSize))
        throw std::runtime_error("Cannot decode input path as UTF-8");
    if (GetACP() == CP_UTF8) return utf8;
    BOOL substituted = FALSE;
    const int narrowSize = WideCharToMultiByte(CP_ACP, WC_NO_BEST_FIT_CHARS, wide.data(), wideSize,
                                               nullptr, 0, nullptr, &substituted);
    if (!narrowSize || substituted) throw std::runtime_error("Input path cannot be represented in the reused parser's Windows code page");
    std::string narrow(static_cast<std::size_t>(narrowSize), '\0');
    if (!WideCharToMultiByte(CP_ACP, WC_NO_BEST_FIT_CHARS, wide.data(), wideSize,
                            narrow.data(), narrowSize, nullptr, &substituted) || substituted)
        throw std::runtime_error("Cannot convert input path for the reused parser");
    return narrow;
#else
    return utf8;
#endif
}

float inward(double value, bool lower) {
    float converted = static_cast<float>(value);
    if (!std::isfinite(value) || !std::isfinite(converted)) throw std::runtime_error("Core bounds exceed finite float coordinates");
    if ((lower && static_cast<double>(converted) < value) || (!lower && static_cast<double>(converted) > value))
        converted = std::nextafter(converted, lower ? std::numeric_limits<float>::infinity()
                                                   : -std::numeric_limits<float>::infinity());
    return converted;
}

void arrange(PlaceData& db, CoreRegion& region, const Summary& summary) {
    if (db.moduleCount != summary.modules || db.MacroCount != summary.fixed || db.netCount != summary.nets
        || db.pinCount != summary.pins || db.Nodes.size() + db.Terminals.size() != static_cast<std::size_t>(summary.modules)
        || db.Terminals.size() != static_cast<std::size_t>(summary.fixed)
        || db.moduleMap.size() != static_cast<std::size_t>(summary.modules)
        || db.Nets.size() != static_cast<std::size_t>(summary.nets)
        || db.Pins.size() != static_cast<std::size_t>(summary.pins)
        || db.SiteRows.size() != static_cast<std::size_t>(summary.rows))
        throw std::runtime_error("Reused parser output does not match validated input counts");

    Bounds rows;
    for (const SiteRow& row : db.SiteRows) {
        if (!std::isfinite(row.start.x) || !std::isfinite(row.end.x) || !std::isfinite(row.bottom)
            || !std::isfinite(row.height) || row.height <= 0 || row.end.x <= row.start.x)
            throw std::runtime_error("Reused parser produced an invalid placement row");
        // end.y 仍在行底部；上边界必须加上行高。
        rows.include(row.start.x, row.bottom, row.end.x, row.bottom + row.height);
    }
    // 读取器的行端点是 float，可能向外舍入；与预检的精确边界取交集。
    region = {inward(std::max(rows.xmin, summary.bounds.xmin), true),
              inward(std::max(rows.ymin, summary.bounds.ymin), true),
              inward(std::min(rows.xmax, summary.bounds.xmax), false),
              inward(std::min(rows.ymax, summary.bounds.ymax), false)};
    if (region.xmin >= region.xmax || region.ymin >= region.ymax)
        throw std::runtime_error("Placement region has no representable positive size");

    for (const auto* group : {&db.Nodes, &db.Terminals}) {
        for (Module* module : *group) {
            if (!module || db.moduleMap.at(module->name) != module)
                throw std::runtime_error("Reused parser produced invalid module references");
            const bool terminal = group == &db.Terminals;
            if (module->isFixed != terminal) throw std::runtime_error("Reused parser fixed classification mismatch");
            module->center.x = static_cast<float>(static_cast<double>(module->center.x) + module->width / 2.0);
            module->center.y = static_cast<float>(static_cast<double>(module->center.y) + module->height / 2.0);
            module->area = static_cast<float>(static_cast<double>(module->width) * module->height);
            module->orientation = 0; // 已经逐条确认文件方向为 N。
            if (!std::isfinite(module->center.x) || !std::isfinite(module->center.y) || !std::isfinite(module->area))
                throw std::runtime_error("Invalid converted coordinates or area: " + module->name);
            if (!module->isFixed) {
                const double xlo = region.xmin + module->width / 2.0;
                const double xhi = region.xmax - module->width / 2.0;
                const double ylo = region.ymin + module->height / 2.0;
                const double yhi = region.ymax - module->height / 2.0;
                if (inward(xlo, true) > inward(xhi, false) || inward(ylo, true) > inward(yhi, false))
                    throw std::runtime_error("Movable module cannot fit placement region: " + module->name);
            }
        }
    }
    for (const Pin* pin : db.Pins) {
        if (!pin || !pin->module || !pin->net || !std::isfinite(pin->offset.x) || !std::isfinite(pin->offset.y))
            throw std::runtime_error("Reused parser produced an invalid pin");
    }
}
}

void LoadedCircuit::clear() noexcept {
    if (ownsDatabase_) {
        // 指针列表本身不会 delete 对象；由这个持有者集中释放一次。
        for (Pin* pin : db.Pins) delete pin;
        for (Net* net : db.Nets) delete net;
        for (Module* module : db.Nodes) delete module;
        for (Module* module : db.Terminals) delete module;
        db = PlaceData{};
        ownsDatabase_ = false;
    }
    region = CoreRegion{};
}

LoadedCircuit::~LoadedCircuit() { clear(); }

void LoadedCircuit::load(const std::string& directory, const std::string& benchmark) {
    if (attempted_) throw std::logic_error("LoadedCircuit supports only one load attempt; create a new object");
    attempted_ = true;
    if (!db.Nodes.empty() || !db.Terminals.empty() || !db.Pins.empty() || !db.Nets.empty()
        || !db.SiteRows.empty() || !db.moduleMap.empty() || db.moduleCount || db.MacroCount || db.netCount || db.pinCount)
        throw std::logic_error("LoadedCircuit must have an empty database before load");
    try {
        if (directory.empty() || benchmark.empty() || benchmark == "." || benchmark == ".."
            || benchmark.find_first_of("/\\\r\n") != std::string::npos || benchmark.find('\0') != std::string::npos)
            throw std::invalid_argument("Provide a directory and a plain benchmark name without path separators");
        if (directory.find('\0') != std::string::npos) throw std::invalid_argument("Input directory contains a null character");
        const Summary summary = inspect(directory + "/" + benchmark);
        const std::string legacyDirectory = parserPath(directory);
        const std::string legacyBenchmark = parserPath(benchmark);
        ownsDatabase_ = true;
        db.Nodes.reserve(static_cast<std::size_t>(summary.modules - summary.fixed));
        db.Terminals.reserve(static_cast<std::size_t>(summary.fixed));
        db.Pins.reserve(static_cast<std::size_t>(summary.pins));
        db.Nets.reserve(static_cast<std::size_t>(summary.nets));
        db.SiteRows.reserve(static_cast<std::size_t>(summary.rows));
        db.readBookShelf(legacyDirectory, legacyBenchmark);
        arrange(db, region, summary);
    } catch (...) {
        clear();
        throw;
    }
}
}
