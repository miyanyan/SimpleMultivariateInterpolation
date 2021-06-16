#ifndef KRIGING_H
#define KRIGING_H

#include <map>
#include <utility>
#include "Eigen/Dense"
#include "interpolationfunction.h"

//reference:
//https://github.com/ByteShark/Kriging
//https://xg1990.com/blog/archives/222
//https://www.cnblogs.com/einyboy/p/3196679.html
class KrigingInterpolation : public InterpolationFunction
{
public:
    enum Model{
        Spherical   = 0, //球
        Exponential = 1, //指数
        Gaussian    = 2  //高斯
    };
    //TODO simpleKriging and any other kriging method
    enum KrigingType{
        OrdinaryKriging = 0,
    };

    KrigingInterpolation(std::vector<DataPoint> points, KrigingType krigingType = OrdinaryKriging, Model modelType = Spherical);

    double estimateValue(const DataPoint &point) override;

private:
    //计算已知点的半方差矩阵
    void calculateCovariogramMatrix();
    //计算未知点的半方差向量
    Eigen::VectorXd calculateCovariogramVector(DataPoint point);
    //计算半方差
    double calculateCovariogram(double distance);
    //普通 克里金
    double ordinaryKrigingForPoint(DataPoint point);

private:
    //kriging类型
    KrigingType m_krigingType;
    //变异函数类型
    Model m_ModelType;
    //已知点的半方差矩阵
    Eigen::MatrixXd m_covariogramMatrix;
    //块金
    double m_nugget;
    //基台
    double m_sill;
    //变成
    double m_range;
};

#endif // KRIGING_H
