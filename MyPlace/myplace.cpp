#include "myplace.h"
#include <Eigen/IterativeLinearSolvers>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <memory>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
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

std::vector<PlacementRange> collectPlacementRanges(const PlaceData& db, const CoreRegion& region)
{
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
    for (const std::vector<Module*>* group : {&db.Nodes, &db.Terminals}) {
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
    return ranges;
}

POS_2D randomCenter(const PlacementRange& range, std::mt19937& rng)
{
    std::uniform_real_distribution<double> randomX(range.x.first, range.x.second);
    std::uniform_real_distribution<double> randomY(range.y.first, range.y.second);
    // 固定先取 x 再取 y，保持相同输入和种子的随机布局结果。
    const float x = static_cast<float>(randomX(rng));
    const float y = static_cast<float>(randomY(rng));
    return POS_2D(x, y);
}

using ClusterId = std::size_t;
using ConnectionScore = std::int64_t;
constexpr std::size_t maxClusterSize = 8;
constexpr std::size_t maxScoredNetSize = 16;
constexpr ConnectionScore scorePerPoint = 360360; // 可被 1 到 15 整除。

struct Cluster {
    std::vector<std::size_t> members;
    double area;
    double x;
    double y;
    std::unordered_map<ClusterId, ConnectionScore> neighbors;
};

// 合并后原编号失效，给新组分配新编号。空指针表示旧组，便于识别过期候选。
using Clusters = std::vector<std::unique_ptr<Cluster>>;

struct MergeCandidate {
    ClusterId first;
    ClusterId second;
    ConnectionScore score;
};

struct MergeOrder {
    bool operator()(const MergeCandidate& a, const MergeCandidate& b) const
    {
        if (a.score != b.score) return a.score < b.score;
        if (a.first != b.first) return a.first > b.first;
        return a.second > b.second;
    }
};

using MergeQueue = std::priority_queue<MergeCandidate, std::vector<MergeCandidate>, MergeOrder>;

void scoreConnections(const PlaceData& db, const std::vector<PlacementRange>& ranges,
                      Clusters& clusters)
{
    std::unordered_map<const Module*, std::size_t> moduleIds;
    moduleIds.reserve(db.Nodes.size() + db.Terminals.size());
    for (std::size_t i = 0; i < ranges.size(); ++i) {
        moduleIds.emplace(ranges[i].module, i);
    }
    // 固定单元也各有编号，参与整条网的去重计数，但不加入可移动组。
    for (const std::vector<Module*>* group : {&db.Nodes, &db.Terminals}) {
        for (const Module* module : *group) {
            if (module->isFixed) moduleIds.emplace(module, moduleIds.size());
        }
    }

    std::unordered_set<const Net*> seenNets;
    seenNets.reserve(db.Nets.size());
    std::vector<std::size_t> connected;
    for (const Net* net : db.Nets) {
        if (!net || !seenNets.insert(net).second) {
            throw std::invalid_argument("Null or duplicate net in PlaceData::Nets");
        }
        connected.clear();
        for (const Pin* pin : net->netPins) {
            if (!pin || pin->net != net) {
                throw std::invalid_argument("Null pin or incorrect pin owner in Net::netPins");
            }
            const auto module = moduleIds.find(pin->module);
            if (module == moduleIds.end()) {
                throw std::invalid_argument("Pin references a module outside PlaceData");
            }
            connected.push_back(module->second);
        }
        std::sort(connected.begin(), connected.end());
        connected.erase(std::unique(connected.begin(), connected.end()), connected.end());
        if (connected.size() < 2 || connected.size() > maxScoredNetSize) continue;

        const auto movableEnd = std::lower_bound(connected.begin(), connected.end(), ranges.size());
        const auto points = scorePerPoint / static_cast<ConnectionScore>(connected.size() - 1);
        for (auto a = connected.begin(); a != movableEnd; ++a) {
            for (auto b = a + 1; b != movableEnd; ++b) {
                clusters[*a]->neighbors[*b] += points;
                clusters[*b]->neighbors[*a] += points;
            }
        }
    }
}

bool canMerge(const Cluster& a, const Cluster& b, double areaLimit)
{
    return a.members.size() + b.members.size() <= maxClusterSize
        && a.area + b.area <= areaLimit;
}

void mergeConnectedClusters(Clusters& clusters, double areaLimit)
{
    std::vector<MergeCandidate> initialCandidates;
    for (ClusterId id = 0; id < clusters.size(); ++id) {
        const Cluster& group = *clusters[id];
        for (const auto& [neighbor, score] : group.neighbors) {
            if (id < neighbor && canMerge(group, *clusters[neighbor], areaLimit)) {
                initialCandidates.push_back({id, neighbor, score});
            }
        }
    }
    MergeQueue candidates(MergeOrder{}, std::move(initialCandidates));
    while (!candidates.empty()) {
        const auto [first, second, score] = candidates.top();
        candidates.pop();
        if (!clusters[first] || !clusters[second]) continue; // 有一组已被合并。

        auto merged = std::move(clusters[first]);
        auto absorbed = std::move(clusters[second]);
        const ClusterId newId = clusters.size();
        const double combinedArea = merged->area + absorbed->area;
        merged->x = (merged->x * merged->area + absorbed->x * absorbed->area) / combinedArea;
        merged->y = (merged->y * merged->area + absorbed->y * absorbed->area) / combinedArea;
        merged->area = combinedArea;
        merged->members.insert(merged->members.end(), absorbed->members.begin(), absorbed->members.end());

        // 合并的两组之间不再计分。到同一邻组的分数相加，再更新邻组的反向记录。
        merged->neighbors.erase(second);
        for (const auto& [neighbor, connection] : absorbed->neighbors) {
            if (neighbor != first) merged->neighbors[neighbor] += connection;
        }
        for (const auto& [neighbor, connection] : merged->neighbors) {
            Cluster& other = *clusters[neighbor];
            other.neighbors.erase(first);
            other.neighbors.erase(second);
            other.neighbors.emplace(newId, connection);
            if (canMerge(*merged, other, areaLimit)) {
                candidates.push({neighbor, newId, connection});
            }
        }
        clusters.push_back(std::move(merged));
    }
}

void arrangeCluster(Cluster& group, const CoreRegion& region,
                    const std::vector<PlacementRange>& ranges, std::vector<POS_2D>& positions)
{
    const std::size_t count = group.members.size();
    if (count == 1) return; // 单个成员沿用临时中心。

    const double side = std::sqrt(2.0 * group.area);
    double width = side;
    double height = side;
    for (const std::size_t member : group.members) {
        width = std::max(width, static_cast<double>(ranges[member].module->width));
        height = std::max(height, static_cast<double>(ranges[member].module->height));
    }
    width = std::min(width, static_cast<double>(region.xmax) - region.xmin);
    height = std::min(height, static_cast<double>(region.ymax) - region.ymin);
    const double xmin = region.xmin;
    const double ymin = region.ymin;
    const double left = std::clamp(group.x - width / 2, xmin, std::max(xmin, region.xmax - width));
    const double bottom = std::clamp(group.y - height / 2, ymin, std::max(ymin, region.ymax - height));

    const auto columns = static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(count))));
    const auto rows = (count + columns - 1) / columns;
    const double cellWidth = width / columns;
    const double cellHeight = height / rows;
    std::sort(group.members.begin(), group.members.end());
    for (std::size_t i = 0; i < count; ++i) {
        const std::size_t member = group.members[i];
        const PlacementRange& range = ranges[member];
        const double x = left + (i % columns + 0.5) * cellWidth;
        const double y = bottom + (i / columns + 0.5) * cellHeight;
        positions[member] = POS_2D(
            static_cast<float>(std::clamp<double>(x, range.x.first, range.x.second)),
            static_cast<float>(std::clamp<double>(y, range.y.first, range.y.second)));
    }
}

constexpr int quadraticRounds = 3;
constexpr double quadraticMinDistance = 25.0;
constexpr double quadraticTolerance = 1e-8;
constexpr int quadraticIterationLimit = 1000;
constexpr int fixedPinRow = -1;

class MovableComponents {
public:
    explicit MovableComponents(std::size_t count) : parent(count, -1) {}

    int root(int row) {
        while (parent[row] >= 0) {
            if (parent[parent[row]] >= 0) parent[row] = parent[parent[row]];
            row = parent[row];
        }
        return row;
    }

    void join(int first, int second) {
        first = root(first);
        second = root(second);
        if (first == second) return;
        // 根记录组大小的负数，把较小的组并入较大的组。
        if (parent[first] > parent[second]) std::swap(first, second);
        parent[first] += parent[second];
        parent[second] = first;
    }

private:
    std::vector<int> parent;
};

struct QuadraticPin {
    int row;
    double x; // 可移动管脚存偏移，固定管脚存绝对坐标。
    double y;
};

struct QuadraticNet {
    std::size_t begin;
    std::size_t movableEnd;
    std::size_t end;
};

struct QuadraticTopology {
    std::vector<QuadraticPin> pins;
    std::vector<QuadraticNet> nets;
    std::vector<int> referenceRows;
    std::size_t tripletCount;
};

QuadraticTopology collectQuadraticTopology(const PlaceData& db,
                                           const std::vector<PlacementRange>& ranges)
{
    // Eigen 的默认稀疏存储用 int 保存下标和累计项数。
    const std::size_t indexLimit = std::numeric_limits<int>::max();
    if (ranges.size() > indexLimit) throw std::length_error("Too many movable modules for Eigen");
    QuadraticTopology topology{{}, {}, {}, ranges.size()};
    topology.pins.reserve(db.Pins.size());
    topology.nets.reserve(db.Nets.size());
    std::unordered_map<const Module*, int> rows;
    rows.reserve(db.Nodes.size() + db.Terminals.size());
    for (std::size_t row = 0; row < ranges.size(); ++row) {
        rows.emplace(ranges[row].module, static_cast<int>(row));
    }
    for (const auto* modules : {&db.Nodes, &db.Terminals}) {
        for (const Module* module : *modules) {
            if (module->isFixed) rows.emplace(module, fixedPinRow);
        }
    }

    MovableComponents components(ranges.size());
    std::vector<bool> touchesFixed(ranges.size(), false);
    std::vector<std::size_t> pinsPerRow(ranges.size(), 0);
    std::unordered_set<const Net*> seenNets;
    seenNets.reserve(db.Nets.size());
    std::unordered_set<const Pin*> seenPins;
    for (const Net* net : db.Nets) {
        if (!net || !seenNets.insert(net).second) {
            throw std::invalid_argument("Null or duplicate net in PlaceData::Nets");
        }
        const auto begin = topology.pins.size();
        seenPins.clear();
        for (const Pin* pin : net->netPins) {
            if (!pin || pin->net != net || !seenPins.insert(pin).second) {
                throw std::invalid_argument("Null, duplicate or misplaced pin in Net::netPins");
            }
            const auto entry = rows.find(pin->module);
            if (entry == rows.end()) {
                throw std::invalid_argument("Pin references a module outside PlaceData");
            }
            double x = pin->offset.x;
            double y = pin->offset.y;
            if (entry->second == fixedPinRow) {
                x += pin->module->center.x;
                y += pin->module->center.y;
            }
            if (!std::isfinite(x) || !std::isfinite(y)) {
                throw std::invalid_argument("Nonfinite pin coordinate for module: " + pin->module->name);
            }
            topology.pins.push_back({entry->second, x, y});
        }

        // 把固定管脚放在末尾，求解时直接略过所有固定–固定组合。
        const auto movableEnd = std::partition(topology.pins.begin() + begin, topology.pins.end(),
            [](const QuadraticPin& pin) { return pin.row != fixedPinRow; });
        const auto end = topology.pins.size();
        const auto split = static_cast<std::size_t>(movableEnd - topology.pins.begin());
        if (split == begin || end - begin < 2) {
            topology.pins.resize(begin);
            continue;
        }

        const int firstRow = topology.pins[begin].row;
        if (split < end) touchesFixed[firstRow] = true;
        for (auto index = begin; index < split; ++index) {
            const int row = topology.pins[index].row;
            components.join(firstRow, row);
            // 每个不同单元的移动管脚对只生成一条下三角记录。
            const auto newPairs = index - begin - pinsPerRow[row]++;
            if (newPairs > indexLimit - topology.tripletCount) {
                throw std::length_error("Too many pin pairs for Eigen sparse storage");
            }
            topology.tripletCount += newPairs;
        }
        for (auto index = begin; index < split; ++index) pinsPerRow[topology.pins[index].row] = 0;
        topology.nets.push_back({begin, split, end});
    }

    std::vector<bool> hasReference(ranges.size(), false);
    for (std::size_t row = 0; row < ranges.size(); ++row) {
        if (touchesFixed[row]) hasReference[components.root(static_cast<int>(row))] = true;
    }
    // 按行号从小到大选代表，给每个无固定参照的组补一个临时参照。
    for (std::size_t row = 0; row < ranges.size(); ++row) {
        const int root = components.root(static_cast<int>(row));
        if (hasReference[root]) continue;
        topology.referenceRows.push_back(static_cast<int>(row));
        hasReference[root] = true;
    }
    return topology;
}

using SparsePlacementMatrix = Eigen::SparseMatrix<double>;
using PlacementTriplet = Eigen::Triplet<double>;

struct AxisPin {
    int row;
    double offsetOrFixedCoordinate;
    double currentCoordinate;
};

void assembleQuadraticAxis(const QuadraticTopology& topology, const Eigen::VectorXd& current,
                           const Eigen::VectorXd& initial, char axis,
                           SparsePlacementMatrix& matrix, Eigen::VectorXd& rhs)
{
    std::vector<AxisPin> pins;
    pins.reserve(topology.pins.size());
    for (const QuadraticPin& pin : topology.pins) {
        const double coordinate = axis == 'x' ? pin.x : pin.y;
        const double base = pin.row == fixedPinRow ? 0.0 : current[pin.row];
        pins.push_back({pin.row, coordinate, base + coordinate});
    }
    Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(current.size());
    std::vector<PlacementTriplet> triplets;
    triplets.reserve(topology.tripletCount);
    for (const QuadraticNet& net : topology.nets) {
        const double pinCountMinusOne = static_cast<double>(net.end - net.begin - 1);
        for (auto i = net.begin; i < net.movableEnd; ++i) {
            const AxisPin& p = pins[i];
            for (auto j = i + 1; j < net.end; ++j) {
                const AxisPin& q = pins[j];
                if (p.row == q.row) continue;
                const double distance = std::abs(p.currentCoordinate - q.currentCoordinate);
                const double weight = 1.0 / (pinCountMinusOne * std::max(distance, quadraticMinDistance));
                const double right = weight * (q.offsetOrFixedCoordinate - p.offsetOrFixedCoordinate);
                diagonal[p.row] += weight;
                rhs[p.row] += right;
                if (q.row == fixedPinRow) continue;
                diagonal[q.row] += weight;
                rhs[q.row] -= right;
                triplets.emplace_back(std::max(p.row, q.row), std::min(p.row, q.row), -weight);
            }
        }
    }
    for (const int row : topology.referenceRows) {
        const double weight = diagonal[row] > 0.0 ? diagonal[row] : 1.0;
        diagonal[row] += weight;
        rhs[row] += weight * initial[row];
    }
    for (int row = 0; row < current.size(); ++row) {
        triplets.emplace_back(row, row, diagonal[row]);
    }
    // 同一矩阵位置可能来自不同管脚对，Eigen 压缩时会累加。
    matrix.setFromTriplets(triplets.begin(), triplets.end());
}

Eigen::VectorXd solveQuadraticAxis(const QuadraticTopology& topology, const Eigen::VectorXd& current,
                                   const Eigen::VectorXd& initial, QuadraticSolveStats& stats)
{
    using Clock = std::chrono::steady_clock;
    const auto started = Clock::now();
    SparsePlacementMatrix matrix(current.size(), current.size());
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(current.size());
    assembleQuadraticAxis(topology, current, initial, stats.axis, matrix, rhs);
    const auto assembled = Clock::now();

    Eigen::ConjugateGradient<SparsePlacementMatrix, Eigen::Lower> solver;
    solver.setTolerance(quadraticTolerance);
    solver.setMaxIterations(quadraticIterationLimit);
    solver.compute(matrix);
    const Eigen::VectorXd result = solver.solveWithGuess(rhs, current);
    const double residual = (matrix.selfadjointView<Eigen::Lower>() * result - rhs).norm()
                          / std::max(1.0, rhs.norm());
    if (solver.info() != Eigen::Success || !result.allFinite() || !std::isfinite(residual)
            || residual > quadraticTolerance) {
        std::ostringstream message;
        message << "quadratic solve failed at round " << stats.round << ", axis " << stats.axis
                << " (iterations=" << solver.iterations() << ", residual=" << residual << ')';
        throw std::runtime_error(message.str());
    }
    stats.triplets = topology.tripletCount;
    stats.storedCoefficients = static_cast<std::size_t>(matrix.nonZeros());
    stats.iterations = static_cast<int>(solver.iterations());
    stats.residual = residual;
    stats.assemblyMilliseconds = std::chrono::duration<double, std::milli>(assembled - started).count();
    stats.solveMilliseconds = std::chrono::duration<double, std::milli>(Clock::now() - assembled).count();
    return result;
}
}

void MyPlacer::randomPlace(const CoreRegion& region, unsigned int seed)
{
    if (!db) throw std::invalid_argument("randomPlace requires a PlaceData database");
    const auto ranges = collectPlacementRanges(*db, region);
    std::mt19937 rng(seed);
    for (const PlacementRange& range : ranges) {
        range.module->center = randomCenter(range, rng);
    }
}

void MyPlacer::clusterPlace(const CoreRegion& region, unsigned int seed,
                            std::vector<std::size_t>* groupIds)
{
    if (!db) throw std::invalid_argument("clusterPlace requires a PlaceData database");
    const auto ranges = collectPlacementRanges(*db, region);
    if (ranges.empty()) {
        if (groupIds) groupIds->clear();
        return;
    }

    Clusters clusters;
    clusters.reserve(ranges.size());
    std::vector<POS_2D> positions;
    positions.reserve(ranges.size());
    std::mt19937 rng(seed);
    double totalArea = 0;
    double largestArea = 0;
    for (std::size_t i = 0; i < ranges.size(); ++i) {
        const Module& module = *ranges[i].module;
        const double area = static_cast<double>(module.width) * module.height;
        const POS_2D start = randomCenter(ranges[i], rng);
        positions.push_back(start);
        clusters.push_back(std::make_unique<Cluster>(Cluster{{i}, area, start.x, start.y, {}}));
        totalArea += area;
        largestArea = std::max(largestArea, area);
    }
    const double areaLimit = std::max(largestArea, maxClusterSize * (totalArea / ranges.size()));
    scoreConnections(*db, ranges, clusters);
    mergeConnectedClusters(clusters, areaLimit);

    std::vector<std::size_t> ids(groupIds ? ranges.size() : 0);
    for (ClusterId id = 0; id < clusters.size(); ++id) {
        if (!clusters[id]) continue;
        Cluster& group = *clusters[id];
        arrangeCluster(group, region, ranges, positions);
        if (groupIds) {
            for (const std::size_t member : group.members) ids[member] = id;
        }
    }

    // 分组、边界检查及内存分配都完成后才改原坐标。
    if (groupIds) groupIds->swap(ids);
    for (std::size_t i = 0; i < ranges.size(); ++i) {
        ranges[i].module->center = positions[i];
    }
}

void MyPlacer::quadraticPlace(const CoreRegion& region, std::vector<QuadraticSolveStats>* solveStats)
{
    if (!db) throw std::invalid_argument("quadraticPlace requires a PlaceData database");
    const auto ranges = collectPlacementRanges(*db, region);
    if (ranges.empty()) {
        if (solveStats) solveStats->clear();
        return;
    }
    const auto topology = collectQuadraticTopology(*db, ranges);
    const double coreX = (static_cast<double>(region.xmin) + region.xmax) / 2.0;
    const double coreY = (static_cast<double>(region.ymin) + region.ymax) / 2.0;
    Eigen::VectorXd x(ranges.size()), y(ranges.size());
    for (std::size_t row = 0; row < ranges.size(); ++row) {
        const auto& range = ranges[row];
        x[row] = std::clamp<double>(coreX, range.x.first, range.x.second);
        y[row] = std::clamp<double>(coreY, range.y.first, range.y.second);
    }
    const Eigen::VectorXd initialX = x;
    const Eigen::VectorXd initialY = y;
    std::vector<QuadraticSolveStats> records;
    records.reserve(2 * quadraticRounds);
    for (int round = 1; round <= quadraticRounds; ++round) {
        QuadraticSolveStats xStats{};
        xStats.round = round;
        xStats.axis = 'x';
        x = solveQuadraticAxis(topology, x, initialX, xStats);
        QuadraticSolveStats yStats{};
        yStats.round = round;
        yStats.axis = 'y';
        y = solveQuadraticAxis(topology, y, initialY, yStats);
        records.push_back(xStats);
        records.push_back(yStats);
        for (std::size_t row = 0; row < ranges.size(); ++row) {
            const auto& range = ranges[row];
            x[row] = std::clamp<double>(x[row], range.x.first, range.x.second);
            y[row] = std::clamp<double>(y[row], range.y.first, range.y.second);
        }
    }
    if (solveStats) solveStats->swap(records);
    for (std::size_t row = 0; row < ranges.size(); ++row) {
        ranges[row].module->center = POS_2D(static_cast<float>(x[row]), static_cast<float>(y[row]));
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
