#include "HotPlot.h".h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    HotPlot w(nullptr, 100, HotPlot::InterpolationMethod::Kriging);
    w.resize(400, 300);
    w.show();

    return a.exec();
}
