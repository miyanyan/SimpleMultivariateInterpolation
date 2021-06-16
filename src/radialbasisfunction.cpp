#include "radialbasisfunction.h"

RadialBasisFunction::RadialBasisFunction(std::vector<DataPoint> points, FunctionType type, double factor)
    : InterpolationFunction(points),
      m_type(type),
      m_factor(factor)
{
    //注意到 m_basisFunctionMatrix 只是单纯的距离相关矩阵，与已知点的value无关，因此当已知点坐标不会变化时，可以只计算一次
    //AW = B --> W = A-1 * B
    //A is the basisfunction matrix
    calculateBasisFunctionMatrix();
    //B is the knowed point value vector
    m_knowedPointVector = Eigen::VectorXd(m_points.size());
    for(int i = 0; i < m_points.size(); ++i){
        m_knowedPointVector(i) = m_points[i].value;
    }
    //W is the weight vector
    m_weight = m_basisFunctionMatrix.inverse() * m_knowedPointVector;
}

/*!
 * \brief RadialBasisFunction::estimateValue 对未知点预测其值
 * \param point
 * \return
 * \time complexity O(n)
 */
double RadialBasisFunction::estimateValue(const DataPoint &point)
{
    //对于一个未知点，其径向基函数向量需要重新计算
    Eigen::RowVectorXd vector(m_points.size());
    for(int i = 0; i < m_points.size(); ++i){
        double r = calculateDistance(point, m_points[i]);
        vector(i) = basisFunction(r);
    }
    return vector * m_weight;
}

double RadialBasisFunction::basisFunction(double r)
{
    double ans = 0.0;
    double tmp = r * r + m_factor;
    switch (m_type) {
        case Multiquadric:
            ans = std::sqrt(tmp);
            break;
        case InverseMultiquadric:
            ans = 1 / std::sqrt(tmp);
        case Multilog:
            ans = log10(tmp);
        case NaturalCubicSpline:
            ans = std::pow(tmp, 1.5);
        case ThinPlateSpline:
            ans = tmp * log10(tmp);
        default:
            break;
    }
    return ans;
}

void RadialBasisFunction::calculateBasisFunctionMatrix()
{
    //AW = F -> W = A-1 F
    //基函数矩阵
    int n = m_points.size();
    m_basisFunctionMatrix = Eigen::MatrixXd(n, n);
    for(int i = 0; i < n; ++i){
        m_basisFunctionMatrix(i, i) = 0.0;
        for(int j = i+1; j < n; ++j){
            double r = calculateDistance(m_points[i], m_points[j]);
            m_basisFunctionMatrix(i, j) = basisFunction(r);
            m_basisFunctionMatrix(j, i) = m_basisFunctionMatrix(i, j);
        }
    }
}
