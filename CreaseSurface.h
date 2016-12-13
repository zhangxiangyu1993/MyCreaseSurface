#ifndef _CREASESURFACE_H_INCLUDED_
#define _CREASESURFACE_H_INCLUDED_


#include "DefineData.h"
#include <teem/air.h>
#include <teem/seek.h>





class ExtractSurface{
public:
	ExtractSurface(vector<vector<vector<GridPoint>>>& Grids);
	~ExtractSurface(){}
	//��ʼ������
	void InitialGrid(vector<vector<vector<GridPoint>>>& GridPoints);
	//���ĺ����������е�������д������㣬������
	void ProcessAllCells();
	//ÿ�����������Ͻ�����˻������ȡ
	void ProcessIntersectionPerCell(GridCell &cell);
	//�ж������Ƿ��н���
	bool IsIntersected(GridPoint& point1,GridPoint& point2);
	//������Ľ�������
	void EvaluateEdgeInterCoord(GridCell& cell, int edgeid);
	//ͨ������������ȡ��������
	void ExtractT(Hessian &hessian, double T[9]);
	//�������Ͻ���ķ�����
	void ComputeEdgeNormal(GridCell &cell, int edgeId, int xi, int yi ,int zi);
	//�������ľ�������
	void ComputeGradientLin(double *result, double *T, double *g, double *Txm, double *gxm, double *Txp, double *gxp,
		double *Tym, double *gym, double *Typ, double *gyp,
		double *Tzm, double *gzm, double *Tzp, double *gzp);

	//����h��ֵ
	void Computeh(GridPoint& point, double h[3]);
    //��������ֵ�Լ���Ӧ����������

	int DetermineCoordinate(GridCell &cell, int index, double& coord1, double& coord2);
	//���������˻��������ֵ

	void FindIntersectionPerFace(GridCell &cell, int faceId);
	//�Ľ���ţ�ٵ�����
	bool NewtonMethod(GridCell &cell, int index, int sameCoord, double& coord1, double& coord2);
	//˫���Բ�ֵ��
	vector<double> BilinearInterpolation(GridCell &cell, int index, int sameCoord, double coord1, double coord2);
	//Լ������ֵ
	vector<double> ExtractSurface::CF(Hessian &hessian);
	//�����˻���ķ�����
	void ComputeFaceNormal(GridCell &cell, int faceId, int xi, int yi ,int zi);
	//�ҵ�ÿ�����ϵ����ӹ�ϵ
	void FindConnectivity(char* pairs, double *bestval, int ct, char *idcs,
		int fixedct, double *coords, double *norms);
	//�����ཻ����
	int lineIntersectionTest (double *l1p1, double *l1p2, double *l2p1, double *l2p2);
	//ȷ�������εĶ��㷽��
	int	checkTriOrientation (double *p1, double *p2, double *p3);
	//�ҵ�ÿ�������������˻�������ӹ�ϵ
	void FindDegConnectivity(GridCell &cell, int x, int y, int z);
	//�ж��Ƿ���԰����ַ�ʽ���ӣ����õ���������֮��Ϊ����
	double EvaluateDegConnection(double *p11, double *p12, double *p13, double *p14,
		double *p21, double *p22, double *p23, double *p24,
		double *p31, double *p32, double *p33, double *p34,
		double *n12, double *n13, double *n22, double *n23,
		double *n32, double *n33);
	//�������ཻ����
	int	triIntersectionTest (double *t1v1, double *t1v2, double *t1v3,
		double *t2v1, double *t2v2, double *t2v3);
	//��ÿ��ȷ�������ӹ�ϵ��������������ǻ�
	void TriangulatePerCell(GridCell &cell, int& pointNum, vector<int>& PointOrder, vector<double>& PointCoords,
		vector<double>& PointNormals);


private:
	vector<vector<vector<GridCell>>> m_Grids;
	vector<vector<vector<GridPoint>>> m_GridPoints;

};


#endif