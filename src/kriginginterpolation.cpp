#include "kriginginterpolation.h"
#include <QTime>
#include <QDebug>

KrigingInterpolation::KrigingInterpolation(std::vector<DataPoint> points, KrigingType krigingType, KrigingInterpolation::Model modelType)
    : InterpolationFunction(points),
      m_krigingType(krigingType),
      m_ModelType(modelType),
      m_sill(0),
      m_nugget(0),
      m_range(0)
{
    //预计算
    int n = m_points.size();
    if(n == 0) return;
    double xMean = 0.0, yMean = 0.0, zMean;
    for(int i = 0; i < n; ++i){
        xMean += m_points[i].x;
        yMean += m_points[i].y;
        zMean += m_points[i].value;
    }
    xMean /= n;
    yMean /= n;
    zMean /= n;

    double numerator = 0.0, denominator = 0.0;
    for(int i = 0; i < n; ++i){
        //1. sill
        m_sill += (m_points[i].value - zMean) * (m_points[i].value - zMean);
        //2. nugget
        numerator   += (m_points[i].x - xMean) * (m_points[i].y - yMean);
        denominator += (m_points[i].x - xMean) * (m_points[i].x - xMean);
    }
    m_sill /= n;
    m_nugget = yMean - (numerator / denominator * xMean);
    //3. range是变成，即有变化的距离, 偷个懒，直接把它定义为点的最大辐射距离
    //不同的Model m_range有变化？
    m_range = 0.5;

    calculateCovariogramMatrix();
}

/*!
 * \brief Kriging::estimateValue 预测未知点的大小
 * \param point
 * \return
 */
double KrigingInterpolation::estimateValue(const DataPoint &point)
{
    double ans = 0.0;
    switch (m_krigingType) {
        case OrdinaryKriging:
            ans = ordinaryKrigingForPoint(point);
            break;
        default:
            break;
    }
    return ans;
}

double KrigingInterpolation::ordinaryKrigingForPoint(DataPoint point)
{
    //AW = b  ->  W = A-1 * b
    Eigen::VectorXd b = calculateCovariogramVector(point);
    Eigen::VectorXd weights =  m_covariogramMatrix * b;

    double estimatedZ = 0.0;
    for(int i = 0; i < m_points.size(); ++i){
        estimatedZ += weights(i) * m_points[i].value;
    }
    return estimatedZ;
}

/*!
 * \brief Kriging::calculateCovariogram 计算协方差
 * \param distance 给定距离
 * \return
 */
double KrigingInterpolation::calculateCovariogram(double distance)
{
    double covariogram = 0.0;
    switch (m_ModelType) {
        case Spherical:
            if(distance > m_range){
                covariogram = 0.0;
            }
            else if(distance > 0.0){
                double rate = distance / m_range;
                covariogram = m_sill * (1 - (1.5 * rate - 0.5 * std::pow(rate, 3.0)));
            }
            else{
                covariogram = m_sill;
            }
            break;
        case Exponential:
            if(distance > 0.0){
                covariogram = (m_sill - m_nugget) * (std::exp(-distance / m_range));
            }else{
                covariogram = m_sill;
            }
        case Gaussian:
            if(distance > 0){
                covariogram = (m_sill - m_nugget) * (std::exp(-std::pow(distance / m_range, 2.0)));
            }else{
                covariogram = m_sill;
            }
        default:
            break;
    }
    return covariogram;
}

/*!
 * \brief Kriging::calculateCovariogramMatrix 求解已知点的半方差矩阵，只需求一次即可
 * \return
 */
void KrigingInterpolation::calculateCovariogramMatrix()
{
    int n = m_points.size();
    m_covariogramMatrix = Eigen::MatrixXd(n+1, n+1);
    for(int i = 0; i < n; ++i){
        m_covariogramMatrix(i, i) = calculateCovariogram(0.0);
        for(int j = i + 1; j < n; ++j){
            double d = calculateDistance(m_points[i], m_points[j]);
            m_covariogramMatrix(i, j) = calculateCovariogram(d);
            m_covariogramMatrix(j, i) = m_covariogramMatrix(i, j);
        }
    }
    for(int i = 0; i < n; ++i){
        m_covariogramMatrix(i, n) = 1.0;
        m_covariogramMatrix(n, i) = 1.0;
    }
    m_covariogramMatrix(n , n) = 0.0;

    m_covariogramMatrix = m_covariogramMatrix.inverse();
}

/*!
 * \brief Kriging::calculateCovariogramVector 未知点的半方差矩阵
 * \param point 未知点
 * \return
 */
Eigen::VectorXd KrigingInterpolation::calculateCovariogramVector(DataPoint point)
{
    int n = m_points.size();
    Eigen::VectorXd distanceVector(n+1);
    for(int i = 0; i < n; ++i){
        double d = calculateDistance(point, m_points[i]);
        distanceVector(i) = calculateCovariogram(d);
    }
    distanceVector(n) = 1.0;

    return distanceVector;
}
