#ifndef INTERPOLATIONFUNCTION_H
#define INTERPOLATIONFUNCTION_H

#include <vector>

class DataPoint
{
public:
    DataPoint() : x(0.0), y(0.0), value(0.0){}

    DataPoint(double x_, double y_, double v_) : x(x_), y(y_), value(v_){}

    DataPoint(const DataPoint& sourceDataPoint)
    {
        *this = sourceDataPoint;
    }

    DataPoint& operator=(const DataPoint& sourceDataPoint)
    {
        if (this != &sourceDataPoint)
        {
            x = sourceDataPoint.x;
            y = sourceDataPoint.y;
            value = sourceDataPoint.value;
        }

        return *this;
    }

    ~DataPoint()
    {
    }

    double x;
    double y;
    double value;
};

class InterpolationFunction
{
public:
    InterpolationFunction(std::vector<DataPoint> points) : m_points(points){}
    //两点距离
    double calculateDistance(const DataPoint& point1, const DataPoint& point2){
        return std::sqrt(std::pow(point1.x - point2.x, 2) + std::pow(point1.y - point2.y, 2));
    }
    //预测未知点
    virtual double estimateValue(const DataPoint& point) = 0;

    std::vector<DataPoint> m_points;
};

#endif // INTERPOLATIONFUNCTION_H
