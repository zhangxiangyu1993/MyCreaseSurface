#include "DefineData.h"
#include <fstream>


using namespace std;

const double threshold = 0.0;

class DataPretreating{
public:
	DataPretreating(char* filename);
	~DataPretreating(){}
	void ReadPoints_Resample(readMode _rm);
	void ReadPoints_WithoutResample(readMode _rm);
	void TransformInputPoints();
	void GetBoundary();
	void InitialGrids(readMode _rm, gridMode _gm);
	void GenerateGrid_Resample();
	void GenerateGrid_WithoutResample();
	double ComputeVol(vector<vector<vector<GridTemp>>>& t_Grid, 
		int idx, int idy, int idz, double max_dis);

	double CalcuFieldVolume(double x,double y,double z);
	void CalcuFieldGradient();
	void CalcuFieldHessian();

	double ComputeCubicBval(double x);
	double ComputeQuadraticBval(double x);
	double ComputeLinearBval(double x);

	vector<vector<vector<GridPoint>>>& GetPoints(){
		return m_PointsGrids;}


private:
	vector<NorPoint> m_Points;
	vector<vector<vector<GridPoint>>> m_PointsTemp;
	vector<vector<vector<GridPoint>>> m_PointsGrids;
	static char* filename;
	static double minX;
	static double maxX;
	static double minY;
	static double maxY;
	static double minZ;
	static double maxZ;
};