#pragma once
// ============================================================================
// 全局布局（对应任务5《全局布局设计》）
//
// 目标函数：  f(x,y) = min( Σ_{e∈E} WL(e;x,y) + λ·D(x,y) )
//   · WL(e;x,y) 为单个网络的线长（WA-LSE 平滑模型，见 ComputeWirelengthGradient）
//   · D(x,y)    为密度惩罚项（静电类比，见 ComputeDensityGradient）
//   · λ         为密度权重（初始由梯度量级比给出，迭代中按公式(37)自适应）
//
// 求解：带动量的梯度下降（步长按"单步位移上限"归一化）。
//
// 计算流程：
//   SetChipBoundary → BinInit → FillerInit → GradientVectorInitialization
//   → PenaltyFactorInit → { GetTotalGradient → MomentumUpdate → UpdateLambda
//                           → ComputeTauOverflow } 循环 → 清理 filler
// ============================================================================

#include "placedata.h"

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

struct GlobalPlacementStats
{
    bool converged{ false };
    int iterations{ 0 };
    int bestIteration{ -1 };
    double initialHPWL{ 0.0 };
    double finalHPWL{ 0.0 };
    double bestTau{ 0.0 };
};

// BinGrid：把布局区域均匀划分成 m×m 个密度网格，其中一个格子的数据。
// 各 *Density 字段保存的是"面积贡献"（长度单位²），不是 0~1 的占用比例；
// 占用比例需要再除以 bin 面积。
struct BinGrid
{
    double nodeDensity{ 0.0 };     // 可移动标准单元/宏单元的面积贡献
    double fillerDensity{ 0.0 };   // 虚拟填充单元(filler)的面积贡献
    double terminalDensity{ 0.0 }; // 固定终端(terminal)的面积贡献，初始化时统计
    double darkDensity{ 0.0 };     // 不可放置区域（SiteRow 覆盖范围之外）的面积贡献

    double x0{ 0.0 }, y0{ 0.0 };   // 网格左下角
    double x1{ 0.0 }, y1{ 0.0 };   // 网格右上角
    double area{ 0.0 };            // 网格面积 = (x1-x0)*(y1-y0)

    // 该网格中的"单元面积"总和（文档步骤 1-4-2-2 的 eDensity 分子）
    double CellArea() const
    {
        return nodeDensity + fillerDensity + terminalDensity + darkDensity;
    }
};

class GlobalPlacement
{
public:
    // 对 db 中已解析好的布局数据做全局布局，结果直接写回 Module::center。
    // 固定单元（terminal 或 .pl 中标记 /FIXED 的单元）不会被移动。
    GlobalPlacementStats run(PlaceData* db, const std::string& iterationLogPath,
                             unsigned int seed = 42);

private:
    // ------------------------------------------------------------------
    // 初始化阶段（文档 3.2 节 1. 初始化阶段）
    // ------------------------------------------------------------------
    void SetChipBoundary(PlaceData* db);           // 由 SiteRow 求芯片边界
    void FillerInit(PlaceData* db);                // 1-1 填充单元初始化
    void BinInit(PlaceData* db);                   // 1-2 密度 Bin 网格初始化
    void GradientVectorInitialization(PlaceData* db); // 1-3 梯度向量初始化
    void BuildNetPinIndex(PlaceData* db);          // 建立 pin→节点下标的快速索引
    double PenaltyFactorInit(PlaceData* db);       // 1-5 惩罚系数 lambda 初始化

    // ------------------------------------------------------------------
    // 梯度计算（文档步骤 1-4）
    // ------------------------------------------------------------------
    void UpdateBinDensities(PlaceData* db);        // 1-4-1 更新 bin 内单元/filler 密度
    void ComputeDensityGradient(PlaceData* db);    // 1-4-2 更新密度梯度 ∇D
    void ComputeWirelengthGradient(PlaceData* db); // 1-4-3 更新线长梯度 ∇WL
    void GetTotalGradient(PlaceData* db, double lambda); // 1-4-4 更新总梯度

    // ------------------------------------------------------------------
    // 位置更新与参数自适应（文档步骤 2-2 / 2-3）
    // ------------------------------------------------------------------
    void MomentumUpdate(PlaceData* db, double beta, double stepLimit);
    double UpdateLambda(double oldLambda, double deltaHPWL);
    void ClampModulesToChip(PlaceData* db);        // 把单元投影回芯片内部
    void ScatterInit(PlaceData* db);               // 退化初始解（全部堆在一点）的随机散布兜底
    bool NeedsScatterInit(const PlaceData* db) const;

    // ------------------------------------------------------------------
    // 评估（文档 3.3 节与停止判定）
    // ------------------------------------------------------------------
    double ComputeHPWL(PlaceData* db);
    double ComputeTauOverflow(PlaceData* db);

    // ------------------------------------------------------------------
    // 工具
    // ------------------------------------------------------------------
    void ComputeDensitySpectrum();                 // 步骤 1-4-2-3 DCT 系数 a_{j,k}
    void ComputeElectricField();                   // 步骤 1-4-2-5 电场 ξx, ξy
    double SampleField(const std::vector<double>& field, double bx, double by) const;
    double BinCoordX(double x) const;              // 物理坐标 → 连续 bin 坐标
    double BinCoordY(double y) const;
    double WirelengthSmoothingGamma() const;       // 公式(6)：γ 随密度溢出率自适应
    double BoundaryPenaltyGrad(double pos, double ll, double ur, double k) const;
    static bool IsValidFinite(double val);

    // ------------------------------------------------------------------
    // 成员变量
    // ------------------------------------------------------------------
    std::vector<Module*> fillerNodes;      // 本算法创建的虚拟填充单元（运行结束释放）
    std::vector<Module*> nodesAndFillers;  // 原始节点 + filler，密度向量按此下标索引
    std::vector<std::vector<BinGrid>> bins;

    int binM{ 64 };                        // 每维网格数（见 BinInit）
    double binWidth{ 0.0 };
    double binHeight{ 0.0 };
    double binArea{ 0.0 };
    POS_2D chipLL;                         // 芯片（core）左下角
    POS_2D chipUR;                         // 芯片（core）右上角

    // 梯度：线长梯度只对应原始单元，密度/总梯度对应单元+filler
    std::vector<double> gradWireX, gradWireY;
    std::vector<double> gradDenX, gradDenY;
    std::vector<double> gradTotalX, gradTotalY;
    std::vector<double> vx, vy;            // 梯度动量
    std::vector<double> momentumSample;    // 估计动量中位数时的抽样缓存

    // DCT 相关查表与频域/空域缓存（尺寸均为 binM*binM，按 [i*binM+j] 存放）
    std::vector<double> cosTable;          // cos(π·j·(x+0.5)/m)，正交 DCT-II 基
    std::vector<double> sinTable;          // sin(π·j·(x+0.5)/m)
    std::vector<double> densityGrid;       // ρ(x,y)
    std::vector<double> dctCoeff;          // a_{j,k}
    std::vector<double> fieldX, fieldY;    // 电场 ξx, ξy（在网格点上取值）
    std::vector<double> fieldWork;         // 计算电场的临时工作数组

    // 网络引脚 → 节点下标的 CSR 索引（-1 表示该 pin 属于固定单元）
    std::vector<std::size_t> netPinStart;
    std::vector<std::int32_t> netPinIndex;
    std::vector<double> pinPosX, pinPosY;  // 按网络临时收集的引脚坐标
    std::vector<std::int32_t> pinNodes;    // 与 pinPosX/Y 对应的节点下标
    std::vector<double> pinExp;            // 4*degree 的指数缓存（a±、b±），避免重复 exp

    // 面积统计与迭代状态
    double chipArea{ 0.0 };
    double stdCellArea{ 0.0 };
    double macroArea{ 0.0 };
    double terminalCoverage{ 0.0 };        // 固定终端在 core 内实际覆盖的面积
    double darkCoverage{ 0.0 };            // 不可放置区域（SiteRow 之外）的面积
    double lastTau{ 1.0 };
    double wireGradMax{ 0.0 };             // 线长梯度最大分量（步长归一化的量级参考）
    double initHPWLRef{ 0.0 };             // 初始 HPWL（lambda 更新时用来缩放 ΔHPWL_REF）
    std::vector<double> hpwlHistory;       // 最近若干次迭代的 HPWL，用于平滑 ΔHPWL
    bool binDensityDirty{ true };          // 单元位置未变时复用上一次的 bin 密度统计
    unsigned int randomSeed{ 42 };

    // 超参数
    double betaMomentum{ 0.90 };           // 动量系数 β
    double targetTau{ 0.10 };              // 目标密度溢出率 τ_min
    int minIter{ 30 };                     // 避免均匀随机起点造成“一轮虚假收敛”
    int maxIter{ 500 };                    // 最大迭代次数
    double targetBinDensity{ 0.75 };       // 目标密度 ρ_t
    double gammaBinScale{ 8.0 };           // 公式(6) 中的 8w_b
    double moveFraction{ 0.02 };           // 单次迭代允许的最大位移（相对芯片短边）
    bool useLengthScaling{ false };        // 文档 1-4-2-1 的 LengthScale（默认关闭，见 .cpp）
    double boundK{ 1200.0 };               // 软边界惩罚强度（安全保护）

    std::ofstream iterLogFile;             // 迭代日志，由主程序指定到输出目录
};
