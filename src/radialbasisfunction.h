#ifndef RADIALBASISFUNCTION_H
#define RADIALBASISFUNCTION_H

#include <algorithm>
#include <Eigen/Dense>
#include "interpolationfunction.h"

class RadialBasisFunction : public InterpolationFunction
{
public:
    //http://surferhelp.goldensoftware.com/griddata/idd_grid_data_radial_basis.htm?tocpath=Gridding%7CGridding%20Methods%7C_____13
    enum FunctionType{
        Multiquadric        = 0,
        InverseMultiquadric = 1,
        Multilog            = 2,
        NaturalCubicSpline  = 3,
        ThinPlateSpline     = 4
    };
    RadialBasisFunction(std::vector<DataPoint> points, FunctionType type = Multiquadric, double factor = 1);

    double estimateValue(const DataPoint &point) override;

private:
    double basisFunction(double r);
    //计算权重
    void calculateBasisFunctionMatrix();

private:
    FunctionType m_type;
    double m_factor;
    Eigen::MatrixXd m_basisFunctionMatrix;
    Eigen::VectorXd m_knowedPointVector;
    Eigen::VectorXd m_weight;
};

#endif // RADIALBASISFUNCTION_H
