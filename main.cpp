#include "gausseid.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    GausSeid w;
    w.show();

    return a.exec();
}
