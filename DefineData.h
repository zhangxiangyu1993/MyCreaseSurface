#ifndef _DEFINEDATA_H_INCLUDEED_
#define _DEFINEDATA_H_INCLUDEED_

#include <vector>
#include <map>

using namespace std;
const int WNUM = 128;
const int HNUM = 128;
const int DNUM = 128;

const int radius = 1;

const int INTERATION = 5; //ţ�ٵ�������������

const int WSampleNum = 50;
const int HSampleNum = 50;
const int DSampleNum = 50;

const int PointNum = 6710;

const double LITTLEVAL = 1e-3;//˫�������ж���ȵķ�Χ
const double DELTA = 0.1; //���޲�ַ��Ĳ���
const double THRESHOLD = 1;//�жϵ�����������ֵ
const double EVALDIFFTHRESH = 1;  //�����й�ʽ��4������ֵ��������������ֵ��ֵ


enum readMode{WITH_VOL_COORD = 0, WITH_ONLY_VOL, WITH_NOTHING};
enum gridMode{RESAMPLE = 0, NON_RESAMPLE};

struct Hessian{
	double fxx;
	double fxy;
	double fxz;
	double fyx;
	double fyy;
	double fyz;
	double fzx;
	double fzy;
	double fzz;
};

struct NorPoint{
	double vol;
	double x;
	double y;
	double z;
};

struct GridTemp{
	vector<int> m_relation;
};

struct GridPoint{
	double x;
	double y;
	double z;
	double gradient[3];
	Hessian hessian;
	double Volume;
};

const int EdgeConnection[12][2] = 
{
	{0,1}, {1,2}, {2,3}, {3,0},
	{4,5}, {5,6}, {6,7}, {7,4},
	{0,4}, {1,5}, {2,6}, {3,7},
};

const int FaceConnection[6][4] = {
	{0,3,2,1},{0,1,5,4},{0,3,7,4},
	{4,7,6,5},{3,2,6,7},{1,2,6,5},
};


const int FaceConnEdge[6][4] = {
	{3,2,1,0},{0,9,4,8},{3,11,7,8},
	{7,6,5,4},{2,10,6,11},{1,10,5,9},
};

struct Edge{
	double intersectionCoord[3];
	double intersectionNormal[3];
	int index;
	bool hasIntersection;
};

struct Face{
	double degenerationCoord[3];
	double degenerationNormal[3];
	int index;
	int interNum;
	char pairs[4]; // ���Ͻ����������Ӧ�ı߱�ţ����Ͻ������Ƕ�Ӧ�����ټ���12
	bool hasIntersection;
};

struct GridCell{
	vector<GridPoint> Points;
	vector<Edge> edges;
	vector<Face> faces;

	int degNum;
	char degpair[6];
};


#endif