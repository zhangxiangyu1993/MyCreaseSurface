#ifndef DATADEFINE_H
#define DATADEFINE_H


#include <QVector>
#include <QList>
#include <QSet>
#include <QPair>
#include <set>
#include <algorithm>


#define  INFINITE   -1;
struct  Vertex
{
	double x;
	double y;
	double z;
	double proper;
};
struct  VertexPro
{
	double x;
	double y;
	double z;
	double proper;
};

struct Cell
{
	int _v[4];              //��������ĸ�����������
	int now_position;       //�����嵱ǰλ��List�е�λ��
	Cell *cells[4];         //��������ĸ��ڽ�������
	bool infinite;          //�Ƿ�Ϊ��������Զ���������
};
struct Facet
{
	int vertex_index;    //���������������Ӧ�Ķ����
	Cell *cell;          //��������������
};
struct Outputcell
{
	int _v[ 4 ];         //�ĸ���������
	int _opp[ 4 ];       //�ĸ��ڽ���������������Ϊ-1����û���ڽ�������
};

//Ҫ����ĵ�������λ��
enum Location
{
	IN,           //�ڲ�
	BOUNDARY,     //�߽�����
	OUT           //�ⲿ
};

typedef QVector<Vertex> Vertices;
typedef QList<Cell> Cells;
typedef QVector<Cell> VCells;
typedef QVector<Facet> Facets;
typedef QVector<Outputcell> Outputcells;

#endif