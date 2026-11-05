#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <shellapi.h>
#endif

#include "initial_place_io.h"
#include "placement_report.h"
#include <charconv>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {
struct Options {
    std::string inputDirectory;
    std::string benchmark;
    std::string method = "random";
    unsigned int seed = 42;
    fs::path output;
    fs::path report;
};

// 旧读取器会提前打印“解析完成”，且 Windows 路径使用系统编码。
// 暂存这些内部提示，正式成功信息等所有检查和保存完成后再输出。
struct ParserLogCapture {
    std::streambuf* previous;
    explicit ParserLogCapture(std::streambuf* buffer) : previous(std::cout.rdbuf(buffer)) {}
    ~ParserLogCapture() { std::cout.rdbuf(previous); }
};

// Windows 直接读取原始宽字符参数，避免系统代码页先把不能表示的字改成别的字。
// 内部文本统一为 UTF-8；这也让中文路径能够可靠地写入 JSON。
std::vector<std::string> utf8Arguments(int argc, char** argv) {
    std::vector<std::string> result;
#ifdef _WIN32
    (void)argc; (void)argv;
    int count = 0;
    wchar_t** wide = CommandLineToArgvW(GetCommandLineW(), &count);
    if (!wide) throw std::runtime_error("Cannot read Unicode command-line arguments");
    std::unique_ptr<void, decltype(&LocalFree)> storage(wide, &LocalFree);
    result.reserve(count);
    for (int i = 0; i < count; ++i) result.push_back(fs::path(wide[i]).u8string());
#else
    result.assign(argv, argv+argc);
#endif
    return result;
}

// writePl 保留原来的 string 接口；在 Windows 上传入系统能打开的文件名。
std::string outputFileName(const fs::path& path) {
#ifdef _WIN32
    const auto wide = path.wstring();
    BOOL replaced = FALSE;
    const UINT page = GetACP();
    const DWORD flags = page == CP_UTF8 ? 0 : WC_NO_BEST_FIT_CHARS;
    const int count = WideCharToMultiByte(page, flags, wide.data(), int(wide.size()),
                                         nullptr, 0, nullptr, page == CP_UTF8 ? nullptr : &replaced);
    if (!count || replaced) throw std::invalid_argument("Output path cannot be represented by the Windows file API");
    std::string result(count, '\0');
    WideCharToMultiByte(page, flags, wide.data(), int(wide.size()), &result[0], count,
                        nullptr, page == CP_UTF8 ? nullptr : &replaced);
    if (replaced) throw std::invalid_argument("Output path contains unsupported characters");
    return result;
#else
    return path.string();
#endif
}

void usage(std::ostream& out) {
    out << "Usage: myplace --input-dir DIR --benchmark NAME --output FILE\n"
           "               [--method random] [--seed UINT] [--report FILE]\n\n"
           "Reads NAME.nodes/.pl/.scl/.nets, keeps fixed modules unchanged,\n"
           "and writes a random initial placement. Overlaps are allowed.\n"
           "The default seed is 42. The default report is FILE.report.json.\n"
           "Only orientation N is supported in this stage.\n";
}

Options parseOptions(const std::vector<std::string>& arguments) {
    Options options;
    std::set<std::string> seen;
    for (std::size_t i = 1; i < arguments.size(); ++i) {
        const std::string& key = arguments[i];
        if (key != "--input-dir" && key != "--benchmark" && key != "--method"
                && key != "--seed" && key != "--output" && key != "--report")
            throw std::invalid_argument("Unknown option: " + key);
        if (!seen.insert(key).second) throw std::invalid_argument("Repeated option: " + key);
        if (++i == arguments.size()) throw std::invalid_argument("Missing value for " + key);
        const std::string& value = arguments[i];
        if (value.empty()) throw std::invalid_argument("Empty value for " + key);
        if (key == "--input-dir") options.inputDirectory = value;
        else if (key == "--benchmark") options.benchmark = value;
        else if (key == "--method") options.method = value;
        else if (key == "--output") options.output = fs::u8path(value);
        else if (key == "--report") options.report = fs::u8path(value);
        else {
            const auto converted = std::from_chars(value.data(), value.data()+value.size(), options.seed);
            if (converted.ec != std::errc{} || converted.ptr != value.data()+value.size())
                throw std::invalid_argument("Seed must be an unsigned 32-bit integer");
        }
    }
    if (options.inputDirectory.empty() || options.benchmark.empty() || options.output.empty())
        throw std::invalid_argument("--input-dir, --benchmark and --output are required");
    if (options.method != "random")
        throw std::invalid_argument("Only method random is implemented; cluster and quadratic are later work");
    if (options.report.empty()) options.report = fs::u8path(options.output.u8string()+".report.json");
    return options;
}

std::string pathKey(const fs::path& path) {
    std::string key = fs::weakly_canonical(fs::absolute(path)).generic_u8string();
#ifdef _WIN32
    // 本项目文件扩展名和参数为 ASCII；Windows 文件名不区分其大小写。
    for (char& c : key) if (c >= 'A' && c <= 'Z') c += 'a'-'A';
#endif
    return key;
}

bool sameFile(const fs::path& a, const fs::path& b) {
    if (fs::exists(a) && fs::exists(b) && fs::equivalent(a, b)) return true;
    return pathKey(a) == pathKey(b);
}

void prepareOutput(const Options& options) {
    if (sameFile(options.output, options.report))
        throw std::invalid_argument("Placement and report must use different output files");
    // 明确保护这份电路的原始文件，包含本阶段无需读取的 aux 和 wts。
    for (const char* extension : {".nodes", ".pl", ".scl", ".nets", ".aux", ".wts"}) {
        const fs::path source = fs::u8path(options.inputDirectory)
                              / fs::u8path(options.benchmark+extension);
        if (sameFile(source, options.output) || sameFile(source, options.report))
            throw std::invalid_argument("Output would overwrite an input file: " + source.u8string());
    }
    for (const auto& target : {options.output, options.report}) {
        if (fs::is_directory(target)) throw std::invalid_argument("Output is a directory: " + target.u8string());
        if (target.has_parent_path()) fs::create_directories(target.parent_path());
    }
    (void)outputFileName(options.output); // 布局开始前检查旧保存接口能否表示这个路径。
}

std::string jsonString(const std::string& value) {
    std::ostringstream escaped;
    escaped << '"';
    for (const unsigned char c : value) {
        if (c == '"' || c == '\\') escaped << '\\' << char(c);
        else if (c < 0x20) escaped << "\\u00" << std::hex << std::setw(2) << std::setfill('0') << int(c);
        else escaped << char(c);
    }
    escaped << '"';
    return escaped.str();
}

void writeReport(const Options& options, const initial_place::LoadedCircuit& circuit,
                 const initial_place::PlacementReport& result, double milliseconds) {
    std::ofstream out(options.report);
    if (!out) throw std::runtime_error("Cannot open report: " + options.report.u8string());
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "{\n  \"benchmark\": " << jsonString(options.benchmark)
        << ",\n  \"method\": \"random\",\n  \"seed\": " << options.seed
        << ",\n  \"parser_revision\": " << jsonString(initial_place::kParserRevision)
#ifdef _WIN32
        << ",\n  \"platform\": \"Windows\""
#else
        << ",\n  \"platform\": \"Linux/Unix\""
#endif
#ifdef __VERSION__
        << ",\n  \"compiler\": " << jsonString(__VERSION__)
#else
        << ",\n  \"compiler\": \"MSVC\""
#endif
        << ",\n  \"input_directory\": " << jsonString(fs::absolute(fs::u8path(options.inputDirectory)).u8string())
        << ",\n  \"output_file\": " << jsonString(fs::absolute(options.output).u8string())
        << ",\n  \"core_bounds\": [" << circuit.region.xmin << ", " << circuit.region.ymin
        << ", " << circuit.region.xmax << ", " << circuit.region.ymax << ']'
        << ",\n  \"total_modules\": " << result.totalModules
        << ",\n  \"movable_modules\": " << result.movableModules
        << ",\n  \"fixed_modules\": " << result.fixedModules
        << ",\n  \"nets\": " << circuit.db.Nets.size()
        << ",\n  \"pins\": " << circuit.db.Pins.size()
        << ",\n  \"rows\": " << circuit.db.SiteRows.size()
        << ",\n  \"out_of_bounds\": " << result.outOfBounds
        << ",\n  \"fixed_changed\": " << result.fixedChanged
        << ",\n  \"invalid_modules\": " << result.invalidModules
        << ",\n  \"hpwl\": " << result.hpwl
        << ",\n  \"placement_milliseconds\": " << milliseconds
        << ",\n  \"overlaps_allowed\": true\n}\n";
    out.close();
    if (!out) throw std::runtime_error("Failed to save report: " + options.report.u8string());
}
}

int main(int argc, char** argv) {
    try {
        const auto arguments = utf8Arguments(argc, argv);
        if (arguments.size() == 2 && arguments[1] == "--help") { usage(std::cout); return 0; }
        if (arguments.size() == 1) { usage(std::cerr); return 2; }
        const Options options = parseOptions(arguments);
        prepareOutput(options);
        initial_place::LoadedCircuit circuit;
        std::ostringstream parserMessages;
        {
            ParserLogCapture capture(parserMessages.rdbuf());
            circuit.load(options.inputDirectory, options.benchmark);
        }
        const auto before = initial_place::capturePlacement(circuit.db);

        // 只统计安排位置的时间，加载、检查和保存都在计时区间之外。
        const auto start = std::chrono::steady_clock::now();
        MyPlacer(&circuit.db).randomPlace(circuit.region, options.seed);
        const auto finish = std::chrono::steady_clock::now();
        const double milliseconds = std::chrono::duration<double, std::milli>(finish-start).count();
        const auto report = initial_place::inspectPlacement(circuit.db, circuit.region, before);
        if (report.outOfBounds || report.fixedChanged || report.invalidModules)
            throw std::runtime_error("Placement validation failed; result was not exported");
        MyPlacer(&circuit.db).writePl(outputFileName(options.output));
        writeReport(options, circuit, report, milliseconds);
        std::cout << "SUCCESS: method=random seed=" << options.seed
                  << " modules=" << report.totalModules << " movable=" << report.movableModules
                  << " fixed=" << report.fixedModules << " HPWL=" << std::setprecision(15) << report.hpwl
                  << " placement_ms=" << milliseconds << '\n'
                  << "Placement: " << options.output.u8string() << '\n'
                  << "Report: " << options.report.u8string() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
