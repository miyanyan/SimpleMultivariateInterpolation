#include "mainwindow.h"

#include <QFile>
#include <QTextStream>
#include <QPainter>
#include <QPen>
#include <QPixmap>

const double g_boxSize = 50;
const double g_boxX[5] = {-50, 50, 50, -50, -50};
const double g_boxY[5] = {50, 50, -50, -50, 50};

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    readChannelAxisFile(":/channelAxis.csv");
    //build voronoi
    boost::polygon::construct_voronoi(m_points.begin(), m_points.end(), &m_voronoi);
    
    visualization();
}

MainWindow::~MainWindow()
{
}

void MainWindow::readChannelAxisFile(QString path)
{
    //读csv文件
    QFile file(path);
    if(!file.exists()) return;
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
    //读取状态，1 通道数， 2 鼻尖坐标
    int flag = 0;
    QTextStream in(&file);
    while(!in.atEnd()){
        //按行读
        QString line = in.readLine();
        QStringList list = line.split(",");
        //通道数
        if(list[0].contains("POS")){
            flag = 1;
            continue;
        }
        //鼻尖
        else if(list[0].contains("NOSE")){
            flag = 2;
            continue;
        }
        //坐标
        if(flag == 1){
            m_points.push_back(Point(list[1].toFloat() * g_boxSize / 0.5, list[2].toFloat() * g_boxSize / 0.5));
        }
    }

//    m_segments.push_back(Segment(Point(-0.5, 0.5), Point(0.5, 0.5)));
//    m_segments.push_back(Segment(Point(-0.5, -0.5), Point(0.5, -0.5)));
    m_segments.push_back(Segment(Point(-g_boxSize, g_boxSize), Point(g_boxSize, g_boxSize)));
    m_segments.push_back(Segment(Point(-g_boxSize, -g_boxSize), Point(g_boxSize, -g_boxSize)));
}

void MainWindow::visualization()
{
    // configure axis rect:
    m_customPlot = new QCustomPlot(this);
    m_customPlot->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming
    m_customPlot->axisRect()->setupFullAxesBox(true);

    m_customPlot->xAxis->setRange(QCPRange(-g_boxSize, g_boxSize));
    m_customPlot->yAxis->setRange(QCPRange(-g_boxSize, g_boxSize));
    //
    m_customPlot->rescaleAxes();
    //
    this->setCentralWidget(m_customPlot);
    //box
    for(int i = 0; i < 4; ++i){
        QCPItemLine *lineLeft = new QCPItemLine(m_customPlot);
        lineLeft->start->setCoords(g_boxX[i], g_boxY[i]);
        lineLeft->end->setCoords(g_boxX[i+1], g_boxY[i+1]);
        lineLeft->setPen(QPen(Qt::black, 2));
    }
    //origin points
    for(auto& point : m_points){
        QCPItemEllipse* pointCircle = new QCPItemEllipse(m_customPlot);
        pointCircle->setAntialiased(true);
        double r = 0.01 * g_boxSize / 0.5;//每个通道的圆半径大小

        pointCircle->topLeft->setCoords(point.x() - r, point.y() + r);
        pointCircle->bottomRight->setCoords(point.x() + r, point.y() - r);
        pointCircle->setPen(QPen(Qt::green, 2));
    }

    for(voronoi_diagram<double>::const_edge_iterator it = m_voronoi.edges().begin();
        it != m_voronoi.edges().end(); ++it){
        //曲线边
        if(it->is_infinite()) continue;
        //直线边
        double x0 = it->vertex0()->x(), y0 = it->vertex0()->y();
        //if(x0 < -5 || x0 > 5 || y0 < -5 || y0 > 5) continue;
        double x1 = it->vertex1()->x(), y1 = it->vertex1()->y();
        //if(x1 < -5 || x1 > 5 || y1 < -5 || y1 > 5) continue;
        //画直线
        QCPItemLine *lineLeft = new QCPItemLine(m_customPlot);
        lineLeft->start->setCoords(x0, y0);
        lineLeft->end->setCoords(x1, y1);
        lineLeft->setPen(QPen(Qt::black, 2));
    }
    for (voronoi_diagram<double>::const_vertex_iterator it = m_voronoi.vertices().begin();
           it != m_voronoi.vertices().end(); ++it) {
        const voronoi_diagram<double>::vertex_type &vertex = *it;
        const voronoi_diagram<double>::edge_type *edge = vertex.incident_edge();
        // This is convenient way to iterate edges around Voronoi vertex.
        QCPItemEllipse *headMask = new QCPItemEllipse(m_customPlot);
        headMask->setAntialiased(true);
        double r = 0.01 * g_boxSize / 0.5;
        headMask->topLeft->setCoords(vertex.x() - r, vertex.y() + r);
        headMask->bottomRight->setCoords(vertex.x() + r, vertex.y() - r);
        headMask->setPen(QPen(Qt::red, 2));

      }
}


