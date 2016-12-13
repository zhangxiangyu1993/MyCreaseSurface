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
	int _v[4];              //四面体的四个顶点索引号
	int now_position;       //四面体当前位于List中的位置
	Cell *cells[4];         //四面体的四个邻接四面体
	bool infinite;          //是否为包含无穷远点的四面体
};
struct Facet
{
	int vertex_index;    //面所属的四面体对应的顶点号
	Cell *cell;          //面所属的四面体
};
struct Outputcell
{
	int _v[ 4 ];         //四个顶点索引
	int _opp[ 4 ];       //四个邻接四面体索引，若为-1，则没有邻接四面体
};

//要插入的点所属的位置
enum Location
{
	IN,           //内部
	BOUNDARY,     //边界面上
	OUT           //外部
};

typedef QVector<Vertex> Vertices;
typedef QList<Cell> Cells;
typedef QVector<Cell> VCells;
typedef QVector<Facet> Facets;
typedef QVector<Outputcell> Outputcells;

#endif