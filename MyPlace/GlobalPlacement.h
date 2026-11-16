#pragma once
#include "placedata.h"
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
    // ��ʼ���׶�
    void FillerInit(PlaceData* db);
    void BinInit(PlaceData* db);
    void GradientVectorInitialization(PlaceData* db);
    double PenaltyFactorInit(PlaceData* db);

    // �ݶȼ���
    void UpdateBinDensities(PlaceData* db);
    void ComputeDensityGradient(PlaceData* db);
    void ComputeWirelengthGradient(PlaceData* db);
    void GetTotalGradient(PlaceData* db, double lambda);

    // Nesterov NAG ����
    void NesterovUpdate(PlaceData* db, double alpha, double beta);

    // ��������Ӧ
    double UpdateLambda(double oldLambda, double deltaHPWL);

    // ��������
    double ComputeHPWL(PlaceData* db);
    double ComputeTauOverflow(PlaceData* db);

    // ��Ա����
    std::vector<Module*> fillerNodes;
    std::vector<Module*> nodesAndFillers;
    std::vector<std::vector<BinGrid>> bins;

    int binM{ 64 };
    POS_2D chipLL;
    POS_2D chipUR;

    // �ݶ�����
    std::vector<double> gradWireX, gradWireY;
    std::vector<double> gradDenX, gradDenY;
    std::vector<double> gradTotalX, gradTotalY;
    std::vector<double> vx, vy; // Nesterov����

    // ������
    double betaMomentum{ 0.80 };
    double targetTau{ 0.10 };
    int maxIter{ 1000 };
    double targetBinDensity{ 0.75 };
};
