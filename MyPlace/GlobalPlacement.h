#pragma once
#include "PlaceData.h"
#include <vector>

struct BinGrid
{
    double nodeDensity{ 0.0 };
    double fillerDensity{ 0.0 };
    double terminalDensity{ 0.0 };
    double darkDensity{ 0.0 };
    double x0{ 0.0 }, y0{ 0.0 };
    double x1{ 0.0 }, y1{ 0.0 };
};

class GlobalPlacement
{
public:
    void run(PlaceData* db);

private:
    // 初始化阶段
    void FillerInit(PlaceData* db);
    void BinInit(PlaceData* db);
    void GradientVectorInitialization(PlaceData* db);
    double PenaltyFactorInit(PlaceData* db);

    // 梯度计算
    void UpdateBinDensities(PlaceData* db);
    void ComputeDensityGradient(PlaceData* db);
    void ComputeWirelengthGradient(PlaceData* db);
    void GetTotalGradient(PlaceData* db, double lambda);

    // Nesterov NAG 更新
    void NesterovUpdate(PlaceData* db, double alpha, double beta);

    // 参数自适应
    double UpdateLambda(double oldLambda, double deltaHPWL);

    // 评估函数
    double ComputeHPWL(PlaceData* db);
    double ComputeTauOverflow(PlaceData* db);

    // 成员变量
    std::vector<Module*> fillerNodes;
    std::vector<Module*> nodesAndFillers;
    std::vector<std::vector<BinGrid>> bins;

    int binM{ 64 };
    POS_2D chipLL;
    POS_2D chipUR;

    // 梯度数组
    std::vector<double> gradWireX, gradWireY;
    std::vector<double> gradDenX, gradDenY;
    std::vector<double> gradTotalX, gradTotalY;
    std::vector<double> vx, vy; // Nesterov动量

    // 超参数
    double betaMomentum{ 0.80 };
    double targetTau{ 0.10 };
    int maxIter{ 1000 };
    double targetBinDensity{ 0.75 };
};
