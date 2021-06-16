#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPointF>

#include "boost/polygon/voronoi.hpp"
#include "boost/polygon/segment_data.hpp"
#include "boost/polygon/voronoi_utils.hpp"

#include "../qcustomplot.h"


using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;

using Point = boost::polygon::point_data<double>;
using Segment = boost::polygon::segment_data<double>;
//using VoronoiUtils = boost::polygon::voronoi_utils<double>;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    void readChannelAxisFile(QString path);
    void visualization();

private:
    std::vector<Point> m_points;
    std::vector<Segment> m_segments;
    voronoi_diagram<double> m_voronoi;
    
    //VoronoiUtils m_voronoiUtils;

    QCustomPlot *m_customPlot;
};
#endif // MAINWINDOW_H
