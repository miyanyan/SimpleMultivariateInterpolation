#include "HotPlot.h"
#include <QTime>
#include <QMetaEnum>

#include "../../src/radialbasisfunction.h"
#include "../../src/kriginginterpolation.h"

//[nx, ny] 实际color映射
//[datax, datay] 预计算大小，然后双线性插值到 [nx, ny]
const int nx = 240;
const int ny = 240;
const int datax = 80;
const int datay = 80;
HotPlot::HotPlot(QWidget *parent, int contourCnt, InterpolationMethod method)
    : QWidget(parent),
      m_contourCnt(contourCnt),
      m_method(method)
{
    //初始化界面
    initWindow();
    //读取64通道电极点位置并决定是否显示通道、通道名
    readChannelAxisFile(":/channelAxis.csv");
    drawElements();
    showChannelCircle(true);
    showChannelName(true);
    //初始化要显示的通道数,默认显示所有的导联
    m_usedChannels.resize(m_channelAxis.size());
    for(int i = 0; i < m_usedChannels.size(); ++i){
        m_usedChannels[i] = i;
    }
    //预计算权重
    precompute();
    //开启测试
    test();
    //
    qDebug() << "interpolation method : " << QMetaEnum::fromType<HotPlot::InterpolationMethod>().valueToKey(m_method);
}

HotPlot::~HotPlot()
{
}

void HotPlot::setEnabledChannels(QStringList &channels)
{
    m_usedChannels.resize(channels.size());
    for(int i = 0; i < m_usedChannels.size(); ++i){
        m_usedChannels[i] = m_channelIndex[channels[i]];
    }
}

void HotPlot::setEnabledChannels(std::vector<int> &channels)
{
    if(channels.size() > m_channelAxis.size()) return;
    m_usedChannels = channels;
}

void HotPlot::setData(std::vector<double> &channels)
{
    calMatrix(channels);
    m_customPlot->replot();
}

void HotPlot::setColorRange(double lower, double upper)
{
    m_colorScale->axis()->setRange(lower, upper);
}

void HotPlot::readChannelAxisFile(QString path)
{
    //读csv文件
    QFile file(path);
    if(!file.exists()) return;
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
    //读取状态，1 通道数， 2 鼻尖坐标
    int flag = 0;
    //通道索引
    int index = 0;
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
        //hash坐标
        if(flag == 1){
            m_channelIndex[list[0]] = index++;
            m_channelAxis.push_back(QPointF(list[1].toFloat(), list[2].toFloat()));
        }
        else if(flag == 2){
            m_noseAxis[list[0]] = QPointF(list[1].toFloat(), list[2].toFloat());
        }
    }
}
/*!
 * \brief plot::drawChannelAxis 画出各个通道的位置，以点的形式
 */
void HotPlot::drawElements()
{
    m_channelCircle.resize(m_channelAxis.size());
    m_channelName.resize(m_channelAxis.size());
    //通道
    for(auto it = m_channelIndex.begin(); it != m_channelIndex.end(); ++it){
        int index = it.value();
        QString name = it.key();
        QPointF point = m_channelAxis[index];
        //每个通道的名字
        m_channelName[index] = new QCPItemText(m_customPlot);
        m_channelName[index]->setText(name);
        m_channelName[index]->position->setCoords(point.x(), point.y());
        m_channelName[index]->setColor(Qt::white);

        //每个通道的圆显示
        m_channelCircle[index] = new QCPItemEllipse(m_customPlot);
        m_channelCircle[index]->setAntialiased(true);
        double r = 0.02;//每个通道的圆半径大小

        m_channelCircle[index]->topLeft->setCoords(point.x() - r, point.y() + r);
        m_channelCircle[index]->bottomRight->setCoords(point.x() + r, point.y() - r);
        m_channelCircle[index]->setPen(QPen(Qt::green, 1));
    }
    //头,colormap 以矩阵形式存储，在轮廓处会有锯齿
    QCPItemEllipse *headMask = new QCPItemEllipse(m_customPlot);
    headMask->setAntialiased(true);
    headMask->topLeft->setCoords(-0.51, 0.51);
    headMask->bottomRight->setCoords(0.51, -0.51);
    headMask->setPen(QPen(Qt::white, 8));
    QCPItemEllipse *head = new QCPItemEllipse(m_customPlot);
    head->setAntialiased(true);
    head->topLeft->setCoords(-0.5, 0.5);
    head->bottomRight->setCoords(0.5, -0.5);
    head->setPen(QPen(Qt::black, 3));
    //鼻尖
    QCPItemLine *lineLeft = new QCPItemLine(m_customPlot);
    lineLeft->start->setCoords(m_noseAxis["left"].x(), m_noseAxis["left"].y());
    lineLeft->end->setCoords(m_noseAxis["top"].x(), m_noseAxis["top"].y());
    lineLeft->setPen(QPen(Qt::black, 2));
    QCPItemLine *lineRight = new QCPItemLine(m_customPlot);
    lineRight->start->setCoords(m_noseAxis["right"].x(), m_noseAxis["right"].y());
    lineRight->end->setCoords(m_noseAxis["top"].x(), m_noseAxis["top"].y());
    lineRight->setPen(QPen(Qt::black, 2.5));
}

void HotPlot::showChannelName(bool flag)
{
    for(int i = 0; i < m_channelName.size(); ++i){
        if(m_channelName[i] == nullptr) continue;
        m_channelName[i]->setVisible(flag);
    }
}

void HotPlot::showChannelCircle(bool flag)
{
    for(int i = 0; i < m_channelCircle.size(); ++i){
        if(m_channelCircle[i] == nullptr) continue;
        m_channelCircle[i]->setVisible(flag);
    }
}

void HotPlot::initWindow()
{
    // configure axis rect:
    m_customPlot = new QCustomPlot(this);
    m_customPlot->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming
    m_customPlot->axisRect()->setupFullAxesBox(true);
    m_customPlot->setOpenGl(true);
    //隐藏坐标，上下左右
    m_customPlot->xAxis->setVisible(true);
    m_customPlot->yAxis->setVisible(true);
    m_customPlot->xAxis2->setVisible(true);
    m_customPlot->yAxis2->setVisible(true);

    // set up the QCPColorMap:
    m_colorMap = new QCPColorMap(m_customPlot->xAxis, m_customPlot->yAxis);
    m_colorMap->setAntialiased(true);
    m_colorMap->data()->setSize(nx, ny); // we want the color map to have nx * ny data points
    m_colorMap->data()->setRange(QCPRange(-0.6, 0.6), QCPRange(-0.6, 0.6)); // and span the coordinate range -4..4 in both key (x) and value (y) dimensions
    // now we assign some data, by accessing the QCPColorMapData instance of the color map:
    //圆外透明度设为0
    m_inCircle = std::vector<std::vector<bool>>(nx, std::vector<bool>(ny, true));
    double x, y;
    for (int xIndex=0; xIndex<nx; ++xIndex)
    {
        for (int yIndex=0; yIndex<ny; ++yIndex)
        {
            m_colorMap->data()->cellToCoord(xIndex, yIndex, &x, &y);
            double r = qSqrt(x * x + y * y);
            if(r >= 0.5){
                m_colorMap->data()->setAlpha(xIndex, yIndex, 0);
                m_inCircle[xIndex][yIndex] = false;
            }
        }
    }

    // add a color scale:
    m_colorScale = new QCPColorScale(m_customPlot);
    m_customPlot->plotLayout()->addElement(0, 1, m_colorScale); // add it to the right of the main axis rect
    m_colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
    m_colorMap->setColorScale(m_colorScale); // associate the color map with the color scale
    m_colorScale->axis()->setLabel("Magnetic Field Strength");


    // set the color gradient of the color map to one of the presets:
    m_colorMap->setGradient(QCPColorGradient::gpJet);
    auto g = m_colorMap->gradient();
    g.setLevelCount(m_contourCnt);
    m_colorMap->setGradient(g);
    // we could have also created a QCPColorGradient instance and added own colors to
    // the gradient, see the documentation of QCPColorGradient for what's possible.

    // rescale the data dimension (color) such that all data points lie in the span visualized by the color gradient:
    m_colorMap->rescaleDataRange();

    // make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
    QCPMarginGroup *marginGroup = new QCPMarginGroup(m_customPlot);
    m_customPlot->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
    m_colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);

    m_colorScale->axis()->setRange(0.0, 4.0);

    // rescale the key (x) and value (y) axes so the whole color map is visible:
    m_customPlot->rescaleAxes();
    //
    QGridLayout *myLayout = new QGridLayout;
    myLayout->addWidget(m_customPlot);
    myLayout->setMargin(0);
    myLayout->setSpacing(0);
    this->setLayout(myLayout);
}
void HotPlot::test()
{
    m_timer.setInterval(100);
    connect(&m_timer, &QTimer::timeout, this, &HotPlot::handleTimeout);
    m_timer.start();
    qsrand(20);
}

void HotPlot::handleTimeout()
{
    std::vector<double> tmp(64, 0);
    for(int i = 0; i < tmp.size(); ++i){
        tmp[i] = qrand() * 1.0 / 10000;
    }

    setData(tmp);
}

/*!
 * \brief HotPlot::calMatrix 计算方法
 * 先缩小矩阵为50 * 50，此时由距离反比插值计算各坐标
 * 放大至200 * 200，采用双线性插值
 * \param channels
 */
void HotPlot::calMatrix(std::vector<double> &channels)
{
    QTime time;
    time.start();
    //缩放系数
    int ratex = nx/datax, ratey = ny/datay;
    //防止size过大
    int size = qMin(m_channelAxis.size(), channels.size());
    //
    if(m_method == IDW || m_method == ModifiedShepard){
        double x, y;
        for(int i = 0; i < datax; ++i){
            for(int j = 0; j < datay; ++j){
                //在圆外，透明度为0
                if(m_inCircle[i * ratex][j*ratey] == false) continue;
                //在圆内，64导数据叠加
                m_colorMap->data()->cellToCoord(i * ratex, j * ratey, &x, &y);
                double a = 0;
                for(int k = 0; k < size; ++k){
                    a += m_channelWeight[i][j][m_usedChannels[k]] * channels[k];
                }
                m_matrix1[i][j] = a;
            }
        }
    }
    else if(m_method == RBF){
        std::vector<DataPoint> points(m_channelAxis.size());
        for(int i = 0; i < points.size(); ++i){
            points[i] = DataPoint(m_channelAxis[i].x(), m_channelAxis[i].y(), 0);
        }
        for(int i = 0; i < size; ++i){
            points[m_usedChannels[i]] = DataPoint(m_channelAxis[m_usedChannels[i]].x(), m_channelAxis[m_usedChannels[i]].y(), channels[i]);
        }
        RadialBasisFunction RBF(points, RadialBasisFunction::Multiquadric, 0.001);
        double x, y;
        for(int i = 0; i < datax; ++i){
            for(int j = 0; j < datay; ++j){
                //在圆外，透明度为0
                if(m_inCircle[i * ratex][j*ratey] == false) continue;
                //在圆内，64导数据叠加
                m_colorMap->data()->cellToCoord(i * ratex, j * ratey, &x, &y);
                m_matrix1[i][j] = RBF.estimateValue(DataPoint(x, y, 0));
            }
        }
    }
    else if(m_method == Kriging){
        std::vector<DataPoint> points(size);
        for(int i = 0; i < size; ++i){
            points[i] = DataPoint(m_channelAxis[m_usedChannels[i]].x(), m_channelAxis[m_usedChannels[i]].y(), channels[i]);
        }
        KrigingInterpolation kriging(points, KrigingInterpolation::OrdinaryKriging, KrigingInterpolation::Spherical);
        double x, y;
        for(int i = 0; i < datax; ++i){
            for(int j = 0; j < datay; ++j){
                //在圆外，透明度为0
                if(m_inCircle[i * ratex][j*ratey] == false) continue;
                //在圆内，64导数据叠加
                m_colorMap->data()->cellToCoord(i * ratex, j * ratey, &x, &y);
                m_matrix1[i][j] = kriging.estimateValue(DataPoint(x, y, 0));
            }
        }
    }
    qDebug() << "interpolation cost : " << time.elapsed();
    time.start();
    //50*50 -> 200 * 200 双线性插值
    for (int i = 0; i < datax - 1; i++)
    {
        for (int j = 0; j < datay - 1; j++)
        {
            double V1 = m_matrix1[i][j];
            double V2 = m_matrix1[i + 1][j];
            double V3 = m_matrix1[i + 1][j + 1];
            double V4 = m_matrix1[i][j + 1];
            for (int m = 0; m < ratex; m++)
            {
                for (int n = 0; n < ratey; n++)
                {
                    int x = i * ratex + m, y = j * ratey + n;
                    if(m_inCircle[x][y] == false) continue;
                    m_matrix[x][y] = doubleLinear(m, n, ratex, ratey, V1, V2, V3, V4);
                    m_colorMap->data()->setCell(x, y, m_matrix[x][y]);
                }
            }
        }
    }
    qDebug() << "magnification cost : " << time.elapsed();
}
/*!
 * \brief HotPlot::precompute 计算每个电极点的权重, 提前计算，以减少复杂度
 */
void HotPlot::precompute()
{
    //初始化缓存矩阵
    m_matrix1 = std::vector<std::vector<double>>(datax, std::vector<double>(datay));
    m_matrix = std::vector<std::vector<double>>(nx, std::vector<double>(ny));
    //初始化权重
    m_channelWeight = std::vector<std::vector<std::vector<double>>>(datax, std::vector<std::vector<double>>(datay, std::vector<double>(m_channelAxis.size())));

    //开始计算
    double x, y;
    //缩放系数
    int ratex = nx/datax, ratey = ny/datay;
    for(int i = 0; i < datax; ++i){
        for(int j = 0; j < datay; ++j){
            m_colorMap->data()->cellToCoord(i * ratex, j * ratey, &x, &y);
            //计算权重值
            double sum_weight = 0.0;
            for(int k = 0; k < m_channelAxis.size(); ++k){
                //(x, y) 为 (i, j) 在 colormap坐标系下的映射
                auto& point = m_channelAxis[k];
                double d = sqrt( (point.x() - x) * (point.x() - x) + (point.y() - y) * (point.y() - y) );
                //插值方法不同，计算权重公式不同
                double weight = 0.0, R = 0.2;
                switch (m_method) {
                    case InterpolationMethod::IDW:
                        weight = 1 / (d * d);
                        break;
                    case InterpolationMethod::ModifiedShepard:
                        weight = qMax(0.0, R - d) / (R * d);
                        weight = weight * weight;
                        break;
                    default:
                        break;
                }
                sum_weight += weight;
                m_channelWeight[i][j][k] = weight;
            }

            for(int k = 0; k < m_channelAxis.size() && sum_weight > 0.0; ++k){
                m_channelWeight[i][j][k] /= sum_weight;
            }
        }
    }
}

double HotPlot::doubleLinear(int m, int n, int X, int Y, double V1, double V2, double V3, double V4)
{
    return (m * n * (V3 - V4 - V2 + V1) + X * n * (V4 - V1) + m * Y * (V2 - V1)) / (X * Y) + V1;
}
double HotPlot::singleLinear(int m, int X, double V1, double V2)
{
    return m * (V2 - V1) / X + V1;
}
