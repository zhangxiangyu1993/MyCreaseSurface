#include <QtGui/QApplication>
#include "mycreasesurface.h"
#include "DataPretreating.h"
#include "CreaseSurface.h"
#include "glwidget.h"


int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	char *fileName = "single_fault.txt";
	DataPretreating DP(fileName);
	DP.InitialGrids(WITH_NOTHING, NON_RESAMPLE);

	ExtractSurface ES(DP.GetPoints());
	ES.ProcessAllCells();

	GLWidget GW;
	GW.show();

	return a.exec();
}
