#ifndef _CREASESURFACE_H_INCLUDED_
#define _CREASESURFACE_H_INCLUDED_


#include "DefineData.h"
#include <teem/air.h>
#include <teem/seek.h>





class ExtractSurface{
public:
	ExtractSurface(vector<vector<vector<GridPoint>>>& Grids);
	~ExtractSurface(){}
	//初始化网格
	void InitialGrid(vector<vector<vector<GridPoint>>>& GridPoints);
	//核心函数，对所有的网格进行处理运算，保存结果
	void ProcessAllCells();
	//每个立方体棱上交点和退化点的求取
	void ProcessIntersectionPerCell(GridCell &cell);
	//判断棱上是否有交点
	bool IsIntersected(GridPoint& point1,GridPoint& point2);
	//计算棱的交点坐标
	void EvaluateEdgeInterCoord(GridCell& cell, int edgeid);
	//通过黑塞矩阵提取张量矩阵
	void ExtractT(Hessian &hessian, double T[9]);
	//计算棱上交点的法向量
	void ComputeEdgeNormal(GridCell &cell, int edgeId, int xi, int yi ,int zi);
	//法向量的具体运算
	void ComputeGradientLin(double *result, double *T, double *g, double *Txm, double *gxm, double *Txp, double *gxp,
		double *Tym, double *gym, double *Typ, double *gyp,
		double *Tzm, double *gzm, double *Tzp, double *gzp);

	//计算h的值
	void Computeh(GridPoint& point, double h[3]);
    //计算特征值以及对应的特征向量

	int DetermineCoordinate(GridCell &cell, int index, double& coord1, double& coord2);
	//计算面上退化点的坐标值

	void FindIntersectionPerFace(GridCell &cell, int faceId);
	//改进的牛顿迭代法
	bool NewtonMethod(GridCell &cell, int index, int sameCoord, double& coord1, double& coord2);
	//双线性插值法
	vector<double> BilinearInterpolation(GridCell &cell, int index, int sameCoord, double coord1, double coord2);
	//约束方程值
	vector<double> ExtractSurface::CF(Hessian &hessian);
	//面上退化点的法向量
	void ComputeFaceNormal(GridCell &cell, int faceId, int xi, int yi ,int zi);
	//找到每个面上的连接关系
	void FindConnectivity(char* pairs, double *bestval, int ct, char *idcs,
		int fixedct, double *coords, double *norms);
	//线线相交测试
	int lineIntersectionTest (double *l1p1, double *l1p2, double *l2p1, double *l2p2);
	//确定三角形的顶点方向
	int	checkTriOrientation (double *p1, double *p2, double *p3);
	//找到每个立方体所有退化点的连接关系
	void FindDegConnectivity(GridCell &cell, int x, int y, int z);
	//判断是否可以按此种方式连接，最后得到的向量积之和为多少
	double EvaluateDegConnection(double *p11, double *p12, double *p13, double *p14,
		double *p21, double *p22, double *p23, double *p24,
		double *p31, double *p32, double *p33, double *p34,
		double *n12, double *n13, double *n22, double *n23,
		double *n32, double *n33);
	//三角形相交测试
	int	triIntersectionTest (double *t1v1, double *t1v2, double *t1v3,
		double *t2v1, double *t2v2, double *t2v3);
	//将每个确定了连接关系的立方体进行三角化
	void TriangulatePerCell(GridCell &cell, int& pointNum, vector<int>& PointOrder, vector<double>& PointCoords,
		vector<double>& PointNormals);


private:
	vector<vector<vector<GridCell>>> m_Grids;
	vector<vector<vector<GridPoint>>> m_GridPoints;

};


#endif