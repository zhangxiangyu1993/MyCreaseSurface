#include "CreaseSurface.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <time.h>

ExtractSurface::ExtractSurface(vector<vector<vector<GridPoint>>>& GridPoints){

	//InitialGrid(m_GridPoints);
	m_GridPoints = GridPoints;
	vector<vector<vector<GridPoint>>>().swap(GridPoints);

}

/* 初始化网格，最小单位为单元cell */
void ExtractSurface::InitialGrid(vector<vector<vector<GridPoint>>>& m_GridPoints){

	for (int i = 0;i < WSampleNum - 1;i++)
	{
		vector<vector<GridCell>> tempvec2d;
		for (int j = 0;j < HSampleNum - 1;j++)
		{
			vector<GridCell> tempvec1d;
			for (int k = 0;k < DSampleNum - 1;k++)
			{
				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[i][j][k]);
				tempCell.Points.push_back(m_GridPoints[i][j + 1][k]);
				tempCell.Points.push_back(m_GridPoints[i + 1][j + 1][k]);
				tempCell.Points.push_back(m_GridPoints[i + 1][j][k]);
				tempCell.Points.push_back(m_GridPoints[i][j][k + 1]);
				tempCell.Points.push_back(m_GridPoints[i][j + 1][k + 1]);
				tempCell.Points.push_back(m_GridPoints[i + 1][j + 1][k + 1]);
				tempCell.Points.push_back(m_GridPoints[i + 1][j][k + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				tempvec1d.push_back(tempCell);
			}

			tempvec2d.push_back(tempvec1d);
		}

		vector<vector<GridPoint>>().swap(m_GridPoints[i]);
		m_Grids.push_back(tempvec2d);
	}

	vector<vector<vector<GridPoint>>>().swap(m_GridPoints);

}




/* 处理所有立方体，找出每个立方体每条边的交点，每个面的退化点。
并确定每个面所有交点和退化点之间的连接性.函数包含了提取曲面算法的主要流程 */
void ExtractSurface::ProcessAllCells(){


	int PointNum = 0;
	int TriangleNum = 0;
	char * File1 = "points.txt";
	ofstream outfile1(File1);

	char * File2 = "norms.txt";
	ofstream outfile2(File2);

	char * File3 = "pointorders.txt";
	ofstream outfile3(File3);

	outfile1 << endl;
	outfile2 << endl;
	outfile3 << endl;

	vector<int> PointOrder;
	vector<double> PointCoords;
	vector<double> PointNormals;

	for (unsigned int i = 0; i < WSampleNum - 1; i++)
	{
		double dur = 0.0;
		clock_t start, end;
		start = clock();

		for (unsigned int j = 0; j < HSampleNum - 1; j++)
		{
			for (unsigned int k = 0; k < DSampleNum - 1; k++)
			{

				

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[i][j][k]);
				tempCell.Points.push_back(m_GridPoints[i][j + 1][k]);
				tempCell.Points.push_back(m_GridPoints[i + 1][j + 1][k]);
				tempCell.Points.push_back(m_GridPoints[i + 1][j][k]);
				tempCell.Points.push_back(m_GridPoints[i][j][k + 1]);
				tempCell.Points.push_back(m_GridPoints[i][j + 1][k + 1]);
				tempCell.Points.push_back(m_GridPoints[i + 1][j + 1][k + 1]);
				tempCell.Points.push_back(m_GridPoints[i + 1][j][k + 1]);

				tempCell.edges.resize(12);
				for (int ii = 0;ii < 12;ii++)
				{
					tempCell.edges[ii].index = ii;
					tempCell.edges[ii].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int ii = 0;ii < 6;ii++)
				{
					tempCell.faces[ii].index = ii;
					tempCell.faces[ii].hasIntersection = false;
					tempCell.faces[ii].interNum = 0;
				}

				
				ProcessIntersectionPerCell(tempCell);
				

				for (int faceid = 0; faceid < 6; faceid++)
				{
					int internum = 0;
					double coords[8];
					double normals[8];
					double temp1, temp2;
					int flag = DetermineCoordinate(tempCell, faceid, temp1, temp2);

					char temp[4];
					char pairs[4];
					int index[4];//通过index保存下点的1,2,3,4与实际边编号或者面编号的对应关系
					if (true == tempCell.faces[faceid].hasIntersection)
					{
						ComputeFaceNormal(tempCell, faceid, i, j, k);
						if (flag == 1)
						{
							coords[2 * internum] = tempCell.faces[faceid].degenerationCoord[1];
							coords[2 * internum + 1] = tempCell.faces[faceid].degenerationCoord[2];

							normals[2 * internum] = tempCell.faces[faceid].degenerationNormal[1];
							normals[2 * internum + 1] = tempCell.faces[faceid].degenerationNormal[2];


						}else if (flag == 2)
						{
							coords[2 * internum] = tempCell.faces[faceid].degenerationCoord[0];
							coords[2 * internum + 1] = tempCell.faces[faceid].degenerationCoord[2];

							normals[2 * internum] = tempCell.faces[faceid].degenerationNormal[0];
							normals[2 * internum + 1] = tempCell.faces[faceid].degenerationNormal[2];
						}else{
							coords[2 * internum] = tempCell.faces[faceid].degenerationCoord[0];
							coords[2 * internum + 1] = tempCell.faces[faceid].degenerationCoord[1];

							normals[2 * internum] = tempCell.faces[faceid].degenerationNormal[0];
							normals[2 * internum + 1] = tempCell.faces[faceid].degenerationNormal[1];
						}

						index[internum++] = faceid + 12;

					}

					for (int ii = 0; ii < 4; ii++)
					{
						int edgeid = FaceConnEdge[faceid][ii];

						if (true == tempCell.edges[edgeid].hasIntersection)
						{
							ComputeEdgeNormal(tempCell, edgeid, i, j, k);

							if (flag == 1)
							{
								coords[2 * internum] = tempCell.edges[edgeid].intersectionCoord[1];
								coords[2 * internum + 1] = tempCell.edges[edgeid].intersectionCoord[2] ;

								normals[2 * internum] = tempCell.edges[edgeid].intersectionNormal[1];
								normals[2 * internum + 1] = tempCell.faces[faceid].degenerationNormal[2];


							}else if (flag == 2)
							{
								coords[2 * internum] = tempCell.edges[edgeid].intersectionCoord[0];
								coords[2 * internum + 1] = tempCell.edges[edgeid].intersectionCoord[2];

								normals[2 * internum] = tempCell.edges[edgeid].intersectionNormal[0];
								normals[2 * internum + 1] = tempCell.faces[faceid].degenerationNormal[2];
							}else{
								coords[2 * internum] = tempCell.edges[edgeid].intersectionCoord[0];
								coords[2 * internum + 1] = tempCell.edges[edgeid].intersectionCoord[1];

								normals[2 * internum] = tempCell.edges[edgeid].intersectionNormal[0];
								normals[2 * internum + 1] = tempCell.faces[faceid].degenerationNormal[1];
							}


							index[internum++] = edgeid;

						}

					}

					tempCell.faces[faceid].interNum = internum;

					for (int jj = 0; jj < internum; jj++)
					{
						pairs[jj] = jj;
						temp[jj] = jj;
					}

					double bestscore = 1e20;

					if (internum > 0)
					{
						FindConnectivity(pairs, &bestscore, internum, temp, 0, coords, normals);

						for (int zz = 0; zz < internum; zz++)
						{
							pairs[zz] = index[pairs[zz]]; //pairs映射为实际的边编号或者面编号
						}

						memcpy(tempCell.faces[faceid].pairs, pairs, sizeof(pairs));

					}
				}

				FindDegConnectivity(tempCell, i, j, k);

				TriangulatePerCell(tempCell, PointNum, PointOrder, PointCoords, PointNormals);

				TriangleNum += PointOrder.size() / 3;

				int ii = 0;
				while (ii < PointCoords.size())
				{
					outfile1 << setiosflags(ios::fixed)<<setw(20)<<PointCoords[ii]<< " ";
					if (!((ii + 1) % 3))
					{
						outfile1 << endl;
					}

					ii++;
				}
		
				int jj = 0;
				while (jj < PointNormals.size())
				{
					outfile2 << setiosflags(ios::fixed)<<setw(20)<<PointNormals[jj]<< " ";
					if (!((jj + 1) % 3))
					{
						outfile2 << endl;
					}

					jj++;
				}

				int zz = 0;
				while (zz < PointOrder.size())
				{
					outfile3 << setiosflags(ios::fixed)<<setw(20)<<PointOrder[zz]<< " ";
					if (!((zz + 1) % 3))
					{
						outfile3 << endl;
					}

					zz++;
				}


				PointOrder.resize(0);
				PointCoords.resize(0);
				PointNormals.resize(0);
		

			}

		}

		end = clock();
		dur = (double)((end - start) / CLOCKS_PER_SEC);


	}

	outfile1.seekp(0,ios::beg);
	outfile2.seekp(0,ios::beg);
	outfile3.seekp(0,ios::beg);

	outfile1 << PointNum << endl;
	outfile2 << PointNum <<  endl;
	outfile3 << TriangleNum << endl;

	outfile1.close();
	outfile2.close();
	outfile3.close();


}


/*将每个立方体的交点连线进行三角化*/
void ExtractSurface::TriangulatePerCell(GridCell &cell, int& pointNum, vector<int>& PointOrder, vector<double>& PointCoords, vector<double>& PointNormals)
{
	int connections[36];
	memset(connections, -1, sizeof(connections));

	int orderTable[18];//将各个交点对应的在所有交点中的序号保存下来
	for (int faceid = 0; faceid < 6; faceid++)
	{
		for (int i = 0;i < cell.faces[faceid].interNum - 1; i += 2)
		{
			int index1 = cell.faces[faceid].pairs[i];
			int index2 = cell.faces[faceid].pairs[i + 1];

			if (connections[2 * index1] == -1)
				connections[2 * index1] = index2;
			else
				connections[2 * index1 + 1] = index2;

			if (connections[2 * index2] == -1)
				connections[2 * index2] = index1;
			else
				connections[2 * index2 + 1] = index1;

		}


	}

	for (int j = 0; j < cell.degNum; j += 2)
	{
		int index1 = cell.degpair[j] + 12;
		int index2 = cell.degpair[j + 1] + 12;

		if (connections[2 * index1] == -1)
			connections[2 * index1] = index2;
		else
			connections[2 * index1 + 1] = index2;

		if (connections[2 * index2] == -1)
			connections[2 * index2] = index1;
		else
			connections[2 * index2 + 1] = index1;
	}

	//从每条边开始搜索首尾相连的交点的形成的多边形
	for (int i=0; i < 12; i++) {
		if (connections[2*i] != -1) {
			signed char polygon[18];
			unsigned char polyct = 0;
			char thiz = i;
			char next = connections[2*i];
			polygon[polyct++] = i;
			connections[2*i] = -1;

			while (next != -1) {
				char helpnext;
				polygon[polyct++] = next;

				if (connections[2*next] == thiz) {
					helpnext = connections[2*next + 1];
				} else {
					helpnext = connections[2*next];
				}
				connections[2*next] = connections[2*next+1] = -1;
				thiz = next;
				next = helpnext;

				if (next == polygon[0])
					break; 
			}

			if (next != -1) { 
				int j;
				for (j=0; j < polyct; ++j) {
					if (polygon[j] < 12) { 
						orderTable[polygon[j]] = pointNum++;
						PointCoords.push_back(cell.edges[polygon[j]].intersectionCoord[0]);
						PointCoords.push_back(cell.edges[polygon[j]].intersectionCoord[1]);
						PointCoords.push_back(cell.edges[polygon[j]].intersectionCoord[2]);

						PointNormals.push_back(cell.edges[polygon[j]].intersectionNormal[0]);
						PointNormals.push_back(cell.edges[polygon[j]].intersectionNormal[1]);
						PointNormals.push_back(cell.edges[polygon[j]].intersectionNormal[2]);

					} else { 
						orderTable[polygon[j]] = pointNum++;
						PointCoords.push_back(cell.faces[polygon[j] - 12].degenerationCoord[0]);
						PointCoords.push_back(cell.faces[polygon[j] - 12].degenerationCoord[1]);
						PointCoords.push_back(cell.faces[polygon[j] - 12].degenerationCoord[2]);

						PointNormals.push_back(cell.faces[polygon[j] - 12].degenerationCoord[0]);
						PointNormals.push_back(cell.faces[polygon[j] - 12].degenerationCoord[1]);
						PointNormals.push_back(cell.faces[polygon[j] - 12].degenerationCoord[2]);
					}
				}


				if (polyct > 4) { 
					double tvertA[4], tvertAsum[4]={0,0,0,0},
						normsum[3]={0,0,0};
					int ovi;
					unsigned int vii[3];

					for (j = 0; j < polyct; j++) {
						if (polygon[j] < 12) {
							int eidx = polygon[j];
							ELL_3V_COPY(tvertA, cell.edges[eidx].intersectionCoord);
							tvertA[3] = 1.0;
							ELL_4V_INCR(tvertAsum,tvertA);

							if (ELL_3V_DOT(normsum, cell.edges[eidx].intersectionNormal) < 0)
								ELL_3V_SUB(normsum, normsum, cell.edges[eidx].intersectionNormal);
							else
								ELL_3V_INCR(normsum, cell.edges[eidx].intersectionNormal);

						} else {
							int fidx = polygon[j] - 12;
							ELL_3V_COPY(tvertA, cell.faces[fidx].degenerationCoord);
							tvertA[3] = 1.0;
							ELL_4V_INCR(tvertAsum,tvertA);
							if (ELL_3V_DOT(normsum,cell.faces[fidx].degenerationNormal) < 0)
								ELL_3V_SUB(normsum,normsum,cell.faces[fidx].degenerationNormal);
							else
								ELL_3V_INCR(normsum,cell.faces[fidx].degenerationNormal);
						}
					}
					//重心坐标的计算问题
					ELL_4V_HOMOG(tvertAsum, tvertAsum);

					ovi = pointNum++;

					PointCoords.push_back(tvertAsum[0]);
					PointCoords.push_back(tvertAsum[1]);
					PointCoords.push_back(tvertAsum[2]);


					double len = ELL_3V_LEN(normsum);//考虑double够不够用

					/* 归一化法向量 */
					ELL_3V_SCALE(normsum, 1.0/len, normsum);

					PointNormals.push_back(normsum[0]);
					PointNormals.push_back(normsum[1]);
					PointNormals.push_back(normsum[2]);

					vii[0] = ovi;

					double coord1[3], coord2[3], coord3[3];
					ELL_3V_COPY(coord1, tvertAsum);

					if (polygon[0] < 12){
						ELL_3V_COPY(coord2, cell.edges[polygon[0]].intersectionCoord);						
						vii[1] = orderTable[polygon[0]];
					}
					else {
						ELL_3V_COPY(coord2, cell.faces[polygon[0] - 12].degenerationCoord);
						vii[1] = orderTable[polygon[0]];
					}

					for (j = 0; j < polyct; ++j) {
						double edgeA[3], edgeB[3];							
						double norm[3];
						vii[2] = vii[1];
						ELL_3V_COPY(coord3, coord2);

						if (j == polyct - 1) {
							if (polygon[0] < 12){
								ELL_3V_COPY(coord2, cell.edges[polygon[0]].intersectionCoord);
								vii[1] = orderTable[polygon[0]];
							}
							else {
								ELL_3V_COPY(coord2, cell.faces[polygon[0] - 12].degenerationCoord);
								vii[1] = orderTable[polygon[0]];
							}
						} else {
							if (polygon[j + 1] < 12){
								ELL_3V_COPY(coord2, cell.edges[polygon[j + 1]].intersectionCoord);
								vii[1] = orderTable[polygon[j + 1]];
							}
							else {
								ELL_3V_COPY(coord2, cell.faces[polygon[j + 1] - 12].degenerationCoord);
								vii[1] = orderTable[polygon[j + 1]];
							}
						}


						ELL_3V_SUB(edgeA, coord2, coord1);
						ELL_3V_SUB(edgeB, coord3, coord1);
						ELL_3V_CROSS(norm, edgeA, edgeB);

						if (ELL_3V_DOT(norm, norm) != 0) {

							PointOrder.push_back(vii[0]);
							PointOrder.push_back(vii[1]);
							PointOrder.push_back(vii[2]);

						}

					}
				} else if (polyct > 2) {

					unsigned int vii[3];
					double coord1[3], coord2[3], coord3[3];
					if (polygon[0] < 12){
						ELL_3V_COPY(coord1, cell.edges[polygon[0]].intersectionCoord);
						vii[0] = orderTable[polygon[0]];

					}
					else {
						ELL_3V_COPY(coord1, cell.faces[polygon[0] - 12].degenerationCoord);
						vii[0] = orderTable[polygon[0]];
					}

					if (polygon[1] < 12){
						ELL_3V_COPY(coord2, cell.edges[polygon[1]].intersectionCoord);
						vii[1] = orderTable[polygon[1]];

					}
					else {
						ELL_3V_COPY(coord2, cell.faces[polygon[1] - 12].degenerationCoord);
						vii[1] = orderTable[polygon[1]];
					}

					if (polygon[2] < 12){
						ELL_3V_COPY(coord3, cell.edges[polygon[2]].intersectionCoord);
						vii[2] = orderTable[polygon[2]];

					}
					else {
						ELL_3V_COPY(coord3, cell.faces[polygon[2] - 12].degenerationCoord);
						vii[2] = orderTable[polygon[2]];
					}

					PointOrder.push_back(vii[0]);
					PointOrder.push_back(vii[1]);
					PointOrder.push_back(vii[2]);


					if (polyct == 4) {
						vii[1] = vii[2];
						ELL_3V_COPY(coord2, coord3);
						if (polygon[3] < 12){
							ELL_3V_COPY(coord3, cell.edges[polygon[3]].intersectionCoord);
							vii[2] = orderTable[polygon[3]];

						}
						else {
							ELL_3V_COPY(coord3, cell.faces[polygon[3] - 12].degenerationCoord);
							vii[2] = orderTable[polygon[3]];
						}

						PointOrder.push_back(vii[0]);
						PointOrder.push_back(vii[1]);
						PointOrder.push_back(vii[2]);
					}
				}

			}
		}
	}
}
/* 确定每个面上退化点，边上交点之间的连接性 */
void ExtractSurface::FindConnectivity(char* pairs, double *bestval, int ct, char* idcs, int fixedct, double *coords, double *norms){
	int i,j;
	if (fixedct==ct) {
		double weight=0;
		for (i=0; i<ct-1; i+=2) {
			double diff[2];
			ELL_2V_SUB(diff,coords+2*idcs[i],coords+2*idcs[i+1]);
			weight+=fabs(ELL_2V_DOT(diff,norms+2*idcs[i]))+
				fabs(ELL_2V_DOT(diff,norms+2*idcs[i+1]));
		}
		if (weight < *bestval) {
			*bestval = weight;
			memcpy(pairs,idcs,sizeof(char)*ct);;
		}
		return;
	}

	for (i=fixedct+1; i<ct; i++) {
		int intersect=0;
		char *idxnew;
		if (NULL == (idxnew = (char*) malloc (sizeof(char)*ct)))
			return;
		memcpy(idxnew,idcs,sizeof(char)*ct);

		idxnew[fixedct+1]=idcs[i];
		idxnew[i]=idcs[fixedct+1];

		for (j=0;j<fixedct;j+=2) {
			if (lineIntersectionTest(coords+2*idxnew[fixedct],
				coords+2*idxnew[fixedct+1],
				coords+2*idxnew[j],coords+2*idxnew[j+1])) {
					intersect=1;
					break;
			}
		}
		if (!intersect) {
			FindConnectivity(pairs, bestval, ct, idxnew, fixedct+2, coords, norms);
		}

		free(idxnew);

	}
}



/* 面上两条线段是否相交的测试 */
int ExtractSurface::lineIntersectionTest (double *l1p1, double *l1p2, double *l2p1, double *l2p2) {
	int or1 = checkTriOrientation(l1p1, l1p2, l2p1);
	int or2 = checkTriOrientation(l1p1, l1p2, l2p2);
	if (or1 != or2) {
		or1 = checkTriOrientation(l2p1, l2p2, l1p1);
		or2 = checkTriOrientation(l2p1, l2p2, l1p2);
		if (or1 != or2)
			return 1;
	}
	return 0;
}


/* 确定一个二维的三角形顶点顺序是顺时针方向还是逆时针方向，顺时针方向返回-1，逆时针方向返回1，共线则返回0.  */
int	ExtractSurface::checkTriOrientation (double *p1, double *p2, double *p3){
	double test = (((p2[0]-p1[0])*(p3[1]-p1[1])) - ((p3[0]-p1[0])*(p2[1]-p1[1])));
	if (test > 0) return 1;
	else if (test < 0) return -1;
	else return 0;
}

/* 确定每个立方体内退化点的连接关系 */
void ExtractSurface::FindDegConnectivity(GridCell &cell, int x, int y, int z){
	if (cell.degNum == 0)
	{
		return;
	}

	if (cell.degNum == 2)
	{
		int index = 0;
		for (int faceid = 0; faceid < 6; faceid++)
		{
			if (true == cell.faces[faceid].hasIntersection)
			{
				cell.degpair[index] = faceid;
				index++;
			}
		}

		memset(cell.degpair + 2, -1, sizeof(cell.degpair) - 2 * sizeof(char));
	}

	if (cell.degNum == 4)
	{
		char degfaceid[4];
		double* degcoord[4];
		double* degnormal[4];
		double* edgecoord[4];
		double sumval[3];
		int index = 0;

		for (int faceid = 0; faceid < 6; faceid++)
		{
			if (true == cell.faces[faceid].hasIntersection)
			{
				degfaceid[index] = faceid;
				degcoord[index] = cell.faces[faceid].degenerationCoord;
				degnormal[index] = cell.faces[faceid].degenerationNormal;
				edgecoord[index] = cell.edges[cell.faces[faceid].pairs[1]].intersectionCoord;//找到与面上的退化点相连接的边交点坐标

				index++;
			}
		}

		sumval[0] = EvaluateDegConnection(edgecoord[0], degcoord[0], degcoord[1], edgecoord[1],
			edgecoord[2], degcoord[2], degcoord[3], edgecoord[3], NULL, NULL, NULL, NULL,
			degnormal[0], degnormal[1], degnormal[2], degnormal[3], NULL, NULL);
		sumval[1] = EvaluateDegConnection(edgecoord[0], degcoord[0], degcoord[2], edgecoord[2],
			edgecoord[1], degcoord[1], degcoord[3], edgecoord[3], NULL, NULL, NULL, NULL,
			degnormal[0], degnormal[2], degnormal[1], degnormal[3], NULL, NULL);
		sumval[2] = EvaluateDegConnection(edgecoord[0], degcoord[0], degcoord[3], edgecoord[3],
			edgecoord[1], degcoord[1], degcoord[2], edgecoord[2], NULL, NULL, NULL, NULL,
			degnormal[0], degnormal[3], degnormal[1], degnormal[2], NULL, NULL);

		int maxidx = 0;
		for (int i = 0; i < 3; i++)
		{
			if (sumval[i] > sumval[maxidx])
			{
				maxidx = i;
			}
		}

		switch(maxidx){
			case 0:
				cell.degpair[0] = degfaceid[0];
				cell.degpair[1] = degfaceid[1];
				cell.degpair[2] = degfaceid[2];
				cell.degpair[3] = degfaceid[3];
				break;
			case 1:
				cell.degpair[0] = degfaceid[0];
				cell.degpair[1] = degfaceid[2];
				cell.degpair[2] = degfaceid[1];
				cell.degpair[3] = degfaceid[3];
				break;
			case 2:
				cell.degpair[0] = degfaceid[0];
				cell.degpair[1] = degfaceid[3];
				cell.degpair[2] = degfaceid[1];
				cell.degpair[3] = degfaceid[2];
				break;	
		}

		memset(cell.degpair + 4, -1, sizeof(cell.degpair) - 4 * sizeof(char));

	}

	if (cell.degNum == 6)
	{

		char pairings[15][6]={{0,1,2,3,4,5},{0,1,2,4,3,5},{0,1,2,5,3,4},
		{0,2,1,3,4,5},{0,2,1,4,3,5},{0,2,1,5,3,4},
		{0,3,1,2,4,5},{0,3,1,4,2,5},{0,3,1,5,2,4},
		{0,4,1,2,3,5},{0,4,1,3,2,5},{0,4,1,5,2,3},
		{0,5,1,2,3,4},{0,5,1,3,2,4},{0,5,1,4,2,3}};

		double *degcoord[6];
		double *degnormal[6];
		double *edgecoord[6];
		double sumval[15];


		for (int faceid = 0; faceid < 6; faceid++)
		{

			degcoord[faceid] = cell.faces[faceid].degenerationCoord;
			degnormal[faceid] = cell.faces[faceid].degenerationNormal;
			edgecoord[faceid] = cell.edges[cell.faces[faceid].pairs[1]].intersectionCoord;
		}

		for (int i = 0; i < 15; i++)
		{
			sumval[i] = EvaluateDegConnection(edgecoord[pairings[i][0]], degcoord[pairings[i][0]], degcoord[pairings[i][1]], edgecoord[pairings[i][1]],
				edgecoord[pairings[i][2]], degcoord[pairings[i][2]], degcoord[pairings[i][3]], edgecoord[pairings[i][3]],
				edgecoord[pairings[i][4]], degcoord[pairings[i][4]], degcoord[pairings[i][5]], edgecoord[pairings[i][5]],
				degnormal[pairings[i][0]], degnormal[pairings[i][1]], degnormal[pairings[i][2]], degnormal[pairings[i][3]],
				degnormal[pairings[i][4]], degnormal[pairings[i][5]]);
		}

		int maxidx = 0;
		for (int i = 0; i < 15; i++)
		{
			if (sumval[i] > sumval[maxidx])
			{
				maxidx = i;
			}
		}

		for (int j = 0; j < 6; j++)
		{
			cell.degpair[j] = pairings[maxidx][j];
		}

	}

}

//确定退化点的连接关系
double ExtractSurface::EvaluateDegConnection(double *p11, double *p12, double *p13, double *p14, double *p21, double *p22, double *p23, double *p24, double *p31, double *p32, double *p33, double *p34, double *n12, double *n13, double *n22, double *n23, double *n32, double *n33)
{
	double diff1[3], diff2[3], diff3[3], ret;
	/* first, perform intersection testing */
	if (triIntersectionTest(p11, p12, p13, p21, p22, p23) ||
		triIntersectionTest(p13, p14, p11, p21, p22, p23) ||
		triIntersectionTest(p11, p12, p13, p23, p24, p21) ||
		triIntersectionTest(p13, p14, p11, p23, p24, p21))
		return 1e20;
	if (p31 != NULL) { /* three pairs - some more to do */
		if (triIntersectionTest(p11, p12, p13, p31, p32, p33) ||
			triIntersectionTest(p11, p12, p13, p33, p34, p31) ||
			triIntersectionTest(p13, p14, p11, p31, p32, p33) ||
			triIntersectionTest(p13, p14, p11, p33, p34, p31) ||
			triIntersectionTest(p21, p22, p23, p31, p32, p33) ||
			triIntersectionTest(p21, p22, p23, p33, p34, p31) ||
			triIntersectionTest(p23, p24, p21, p31, p32, p33) ||
			triIntersectionTest(p23, p24, p21, p33, p34, p31))
			return 1e20;
	}
	ELL_3V_SUB(diff1,p13,p12);
	ELL_3V_SUB(diff2,p23,p22);
	ret=fabs(ELL_3V_DOT(diff1,n12))+fabs(ELL_3V_DOT(diff1,n13))+
		fabs(ELL_3V_DOT(diff2,n22))+fabs(ELL_3V_DOT(diff2,n23));
	if (p31 != NULL) {
		ELL_3V_SUB(diff3,p33,p32);
		ret+=fabs(ELL_3V_DOT(diff3,n32))+fabs(ELL_3V_DOT(diff3,n33));
	}
	return ret;
}

//三角形两两相交测试
int ExtractSurface::triIntersectionTest(double *t1v1, double *t1v2, double *t1v3, double *t2v1, double *t2v2, double *t2v3){
	double n1[3], n2[3], d1, d2;
	double diff1[3], diff2[3];
	double t2sd1, t2sd2, t2sd3;
	ELL_3V_SUB(diff1, t1v2, t1v1);
	ELL_3V_SUB(diff2, t1v3, t1v1);
	ELL_3V_CROSS(n1,diff1,diff2);
	d1=-ELL_3V_DOT(n1,t1v1);

	/* compute scaled signed distances of t2 to plane of t1 */
	t2sd1 = ELL_3V_DOT(n1, t2v1)+d1;
	t2sd2 = ELL_3V_DOT(n1, t2v2)+d1;
	t2sd3 = ELL_3V_DOT(n1, t2v3)+d1;

	if (t2sd1==0 && t2sd2==0 && t2sd3==0) {
		/* coplanar case: handle in 2D */
		double t1v12d[2], t1v22d[2], t1v32d[2], t2v12d[2], t2v22d[2], t2v32d[2];
		if (fabs(n1[0])>=fabs(n1[1]) && fabs(n1[0])>=fabs(n1[2])) {
			t1v12d[0]=t1v1[1]; t1v12d[1]=t1v1[2];
			t1v22d[0]=t1v2[1]; t1v22d[1]=t1v2[2];
			t1v32d[0]=t1v3[1]; t1v32d[1]=t1v3[2];
			t2v12d[0]=t2v1[1]; t2v12d[1]=t2v1[2];
			t2v22d[0]=t2v2[1]; t2v22d[1]=t2v2[2];
			t2v32d[0]=t2v3[1]; t2v32d[1]=t2v3[2];
		} else if (fabs(n1[1])>=fabs(n1[0]) && fabs(n1[1])>=fabs(n1[2])) {
			t1v12d[0]=t1v1[0]; t1v12d[1]=t1v1[2];
			t1v22d[0]=t1v2[0]; t1v22d[1]=t1v2[2];
			t1v32d[0]=t1v3[0]; t1v32d[1]=t1v3[2];
			t2v12d[0]=t2v1[0]; t2v12d[1]=t2v1[2];
			t2v22d[0]=t2v2[0]; t2v22d[1]=t2v2[2];
			t2v32d[0]=t2v3[0]; t2v32d[1]=t2v3[2];
		} else {
			t1v12d[0]=t1v1[0]; t1v12d[1]=t1v1[1];
			t1v22d[0]=t1v2[0]; t1v22d[1]=t1v2[1];
			t1v32d[0]=t1v3[0]; t1v32d[1]=t1v3[1];
			t2v12d[0]=t2v1[0]; t2v12d[1]=t2v1[1];
			t2v22d[0]=t2v2[0]; t2v22d[1]=t2v2[1];
			t2v32d[0]=t2v3[0]; t2v32d[1]=t2v3[1];
		}
		/* we may assume that none of the triangles is fully contained
		* within the other. Thus, it suffices to do a lot of 2D line-line
		* intersections */
		if (lineIntersectionTest(t1v12d, t1v22d, t2v12d, t2v22d) ||
			lineIntersectionTest(t1v22d, t1v32d, t2v12d, t2v22d) ||
			lineIntersectionTest(t1v32d, t1v12d, t2v12d, t2v22d) ||
			lineIntersectionTest(t1v12d, t1v22d, t2v22d, t2v32d) ||
			lineIntersectionTest(t1v22d, t1v32d, t2v22d, t2v32d) ||
			lineIntersectionTest(t1v32d, t1v12d, t2v22d, t2v32d) ||
			lineIntersectionTest(t1v12d, t1v22d, t2v32d, t2v12d) ||
			lineIntersectionTest(t1v22d, t1v32d, t2v32d, t2v12d) ||
			lineIntersectionTest(t1v32d, t1v12d, t2v32d, t2v12d))
			return 1;
		return 0;
	} else {
		/* pointers to the vertices on the same side / opposite side of plane */
		double *t2s11, *t2s12, *t2s2, t2s11sd, t2s12sd, t2s2sd;
		double t1sd1, t1sd2, t1sd3;
		double *t1s11, *t1s12, *t1s2, t1s11sd, t1s12sd, t1s2sd;
		double t1p11, t1p12, t1p2, t2p11, t2p12, t2p2;
		double D[3]; /* direction vector of line */
		double t1t1, t1t2, t2t1, t2t2;
		if (t2sd1*t2sd2>=0 && t2sd1*t2sd3<=0) {
			t2s11=t2v1; t2s12=t2v2; t2s2=t2v3; t2s11sd=t2sd1;
			t2s12sd=t2sd2; t2s2sd=t2sd3;
		} else if (t2sd1*t2sd3>=0 && t2sd1*t2sd2<=0) {
			t2s11=t2v1; t2s12=t2v3; t2s2=t2v2; t2s11sd=t2sd1;
			t2s12sd=t2sd3; t2s2sd=t2sd2;
		} else if (t2sd2*t2sd3>=0 && t2sd1*t2sd2<=0) {
			t2s11=t2v2; t2s12=t2v3; t2s2=t2v1; t2s11sd=t2sd2;
			t2s12sd=t2sd3; t2s2sd=t2sd1;
		} else
			return 0; /* all on the same side; no intersection */

		/* same game for triangle 2 */
		ELL_3V_SUB(diff1, t2v2, t2v1);
		ELL_3V_SUB(diff2, t2v3, t2v1);
		ELL_3V_CROSS(n2, diff1, diff2);
		d2=-ELL_3V_DOT(n2,t2v1);
		t1sd1 = ELL_3V_DOT(n2, t1v1)+d2;
		t1sd2 = ELL_3V_DOT(n2, t1v2)+d2;
		t1sd3 = ELL_3V_DOT(n2, t1v3)+d2;
		if (t1sd1*t1sd2>=0 && t1sd1*t1sd3<=0) {
			t1s11=t1v1; t1s12=t1v2; t1s2=t1v3; t1s11sd=t1sd1;
			t1s12sd=t1sd2; t1s2sd=t1sd3;
		} else if (t1sd1*t1sd3>=0 && t1sd1*t1sd2<=0) {
			t1s11=t1v1; t1s12=t1v3; t1s2=t1v2; t1s11sd=t1sd1;
			t1s12sd=t1sd3; t1s2sd=t1sd2;
		} else if (t1sd2*t1sd3>=0 && t1sd1*t1sd2<=0) {
			t1s11=t1v2; t1s12=t1v3; t1s2=t1v1; t1s11sd=t1sd2;
			t1s12sd=t1sd3; t1s2sd=t1sd1;
		} else
			return 0; /* all on the same side; no intersection */

		/* both planes intersect in a line; check if the intervals on that
		* line intersect */
		ELL_3V_CROSS(D,n1,n2);
		/* we are only interested in component magnitudes */
		D[0]=fabs(D[0]); D[1]=fabs(D[1]); D[2]=fabs(D[2]);
		if (D[0]>=D[1] && D[0]>=D[2]) {
			t1p11=t1s11[0]; t1p12=t1s12[0]; t1p2=t1s2[0];
			t2p11=t2s11[0]; t2p12=t2s12[0]; t2p2=t2s2[0];
		} else if (D[1]>=D[0] && D[1]>=D[2]) {
			t1p11=t1s11[1]; t1p12=t1s12[1]; t1p2=t1s2[1];
			t2p11=t2s11[1]; t2p12=t2s12[1]; t2p2=t2s2[1];
		} else {
			t1p11=t1s11[2]; t1p12=t1s12[2]; t1p2=t1s2[2];
			t2p11=t2s11[2]; t2p12=t2s12[2]; t2p2=t2s2[2];
		}
		/* compute interval boundaries */
		t1t1=t1p11+(t1p2-t1p11)*t1s11sd/(t1s11sd-t1s2sd);
		t1t2=t1p12+(t1p2-t1p12)*t1s12sd/(t1s12sd-t1s2sd);
		if (t1t1>t1t2) {
			double help=t1t1;
			t1t1=t1t2;
			t1t2=help;
		}
		t2t1=t2p11+(t2p2-t2p11)*t2s11sd/(t2s11sd-t2s2sd);
		t2t2=t2p12+(t2p2-t2p12)*t2s12sd/(t2s12sd-t2s2sd);
		if (t2t1>t2t2) {
			double help=t2t1;
			t2t1=t2t2;
			t2t2=help;
		}
		/* test for interval intersection */
		if (t2t1>t1t2 || t1t1>t2t2) return 0;
		return 1;
	}
}
/* 找出crease surface 与一个立方体的边上的交点以及面上的退化点 */
void ExtractSurface::ProcessIntersectionPerCell(GridCell &cell){
	for (int i = 0; i < 6; i++)
	{
		int countIntersection = 0;

		for (int j = 0; j < 4; j++)
		{

			int index1,index2;
			if (j == 3)
			{
				index1 = FaceConnection[i][j];
				index2 = FaceConnection[i][0];
			}else{
				index1 = FaceConnection[i][j];
				index2 = FaceConnection[i][j + 1];
			}

			GridPoint point1 = cell.Points[index1], point2 = cell.Points[index2];

			if (IsIntersected(point1,point2))
			{
				cell.edges[FaceConnEdge[i][j]].hasIntersection = true;
				EvaluateEdgeInterCoord(cell, FaceConnEdge[i][j]);
				countIntersection++;
			}

		}

		if ((countIntersection & 0x1) == 1)
		{
			FindIntersectionPerFace(cell, i);
		}
	}

}

void ExtractSurface::ComputeEdgeNormal(GridCell &cell, int edgeId, int xi, int yi ,int zi){
	double Txm[9], Txp[9], Tym[9], Typ[9], Tzm[9], Tzp[9], T[9], T0[9], T1[9],
		gxm[3], gxp[3], gym[3], gyp[3], gzm[3], gzp[3], g[3];
	double alpha = 0.5;
	GridPoint point1 = cell.Points[EdgeConnection[edgeId][0]], point2 = cell.Points[EdgeConnection[edgeId][1]];
	ExtractT(point1.hessian, T0);
	ExtractT(point2.hessian, T1);

	ELL_3M_LERP(T, alpha, T0, T1);
	ELL_3V_LERP(g, alpha, point1.gradient, point2.gradient);

	switch(edgeId){
		case 0:
			ExtractT(cell.Points[EdgeConnection[2][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[2][1]].hessian, T1);

			ELL_3M_LERP(Txp, alpha, T0, T1);
			ELL_3V_LERP(gxp, alpha, cell.Points[EdgeConnection[2][0]].gradient,
				cell.Points[EdgeConnection[2][1]].gradient);

			if (xi == 0)
			{
				ELL_3M_COPY(Txm, T); ELL_3V_COPY(gxm, g);
			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);

				ELL_3M_LERP(Txm, alpha, T0, T1);
				ELL_3V_LERP(gxm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[4][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[4][1]].hessian, T1);

			ELL_3M_LERP(Tzp, alpha, T0, T1);
			ELL_3V_LERP(gzp, alpha, cell.Points[EdgeConnection[4][0]].gradient,
				cell.Points[EdgeConnection[4][1]].gradient);

			if (zi == 0)
			{
				ELL_3M_COPY(Tzm, T); ELL_3V_COPY(gzm, g);

			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi- 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}


				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);

				ELL_3M_LERP(Txm, alpha, T0, T1);
				ELL_3V_LERP(gxm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(point1.hessian, T0);
			ExtractT(point1.hessian, T1);

			ELL_3M_COPY(Tym, T0); ELL_3V_COPY(gym, point1.gradient);
			ELL_3M_COPY(Typ, T1); ELL_3V_COPY(gyp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 1:

			ExtractT(cell.Points[EdgeConnection[3][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[3][1]].hessian, T1);

			ELL_3M_LERP(Typ, alpha, T0, T1);
			ELL_3V_LERP(gyp, alpha, cell.Points[EdgeConnection[3][0]].gradient,
				cell.Points[EdgeConnection[3][1]].gradient);

			if (yi == HSampleNum - 2)
			{
				ELL_3M_COPY(Tym, T); ELL_3V_COPY(gym, g);
			}else{
				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 2][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 2][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 2][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 2][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}


				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);

				ELL_3M_LERP(Tym, alpha, T0, T1);
				ELL_3V_LERP(gym, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[5][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[5][1]].hessian, T1);

			ELL_3M_LERP(Tzp, alpha, T0, T1);
			ELL_3V_LERP(gzp, alpha, cell.Points[EdgeConnection[5][0]].gradient,
				cell.Points[EdgeConnection[5][1]].gradient);

			if (zi == 0)
			{
				ELL_3M_COPY(Tzm, T); ELL_3V_COPY(gzm, g);

			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi- 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}


				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);

				ELL_3M_LERP(Tzm, alpha, T0, T1);
				ELL_3V_LERP(gzm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(point1.hessian, T0);
			ExtractT(point2.hessian, T1);

			ELL_3M_COPY(Txm, T0); ELL_3V_COPY(gxm, point1.gradient);
			ELL_3M_COPY(Txp, T1); ELL_3V_COPY(gxp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 2:

			ExtractT(cell.Points[EdgeConnection[0][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[0][1]].hessian, T1);

			ELL_3M_LERP(Txp, alpha, T0, T1);
			ELL_3V_LERP(gxp, alpha, cell.Points[EdgeConnection[0][0]].gradient,
				cell.Points[EdgeConnection[0][1]].gradient);

			if (xi == WSampleNum - 2)
			{
				ELL_3M_COPY(Txm, T); ELL_3V_COPY(gxm, g);
			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}


				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Txm, alpha, T0, T1);
				ELL_3V_LERP(gxm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[6][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[6][1]].hessian, T1);

			ELL_3M_LERP(Tzp, alpha, T0, T1);
			ELL_3V_LERP(gzp, alpha, cell.Points[EdgeConnection[6][0]].gradient,
				cell.Points[EdgeConnection[6][1]].gradient);

			if (zi == 0)
			{
				ELL_3M_COPY(Tzm, T); ELL_3V_COPY(gzm, g);

			}else{
				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi- 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}



				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tzm, alpha, T0, T1);
				ELL_3V_LERP(gzm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(point1.hessian, T0);
			ExtractT(point2.hessian, T1);

			ELL_3M_COPY(Tym, T0); ELL_3V_COPY(gym, point1.gradient);
			ELL_3M_COPY(Typ, T1); ELL_3V_COPY(gyp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 3:
			ExtractT(cell.Points[EdgeConnection[1][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[1][1]].hessian, T1);
			ELL_3M_LERP(Typ, alpha,T0,T1);
			ELL_3V_LERP(gyp, alpha, cell.Points[EdgeConnection[1][0]].gradient,
				cell.Points[EdgeConnection[1][1]].gradient);

			if (yi == 0)
			{
				ELL_3M_COPY(Tym, T); ELL_3V_COPY(gym, g);
			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi - 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi - 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi - 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi - 1][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian,T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian,T1);
				ELL_3M_LERP(Tym, alpha, T0, T1);
				ELL_3V_LERP(gym, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[7][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[7][1]].hessian, T1);
			ELL_3M_LERP(Tzp, alpha, T0 ,T1);
			ELL_3V_LERP(gzp, alpha, cell.Points[EdgeConnection[7][0]].gradient,
				cell.Points[EdgeConnection[7][1]].gradient);

			if (zi == 0)
			{
				ELL_3M_COPY(Tzm, T); ELL_3V_COPY(gzm, g);

			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi - 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi- 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tzm, alpha, T0, T1);
				ELL_3V_LERP(gzm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}
			ExtractT(point1.hessian, T0);
			ExtractT(point2.hessian, T1);
			ELL_3M_COPY(Txm, T0); ELL_3V_COPY(gxm, point1.gradient);
			ELL_3M_COPY(Txp, T1); ELL_3V_COPY(gxp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 4:

			ExtractT(cell.Points[EdgeConnection[6][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[6][1]].hessian, T1);
			ELL_3M_LERP(Txp, alpha, T0, T1);
			ELL_3V_LERP(gxp, alpha, cell.Points[EdgeConnection[6][0]].gradient,
				cell.Points[EdgeConnection[6][1]].gradient);

			if (xi == 0)
			{
				ELL_3M_COPY(Txm, T); ELL_3V_COPY(gxm, g);
			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Txm, alpha, T0, T1);
				ELL_3V_LERP(gxm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[0][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[0][1]].hessian, T1);
			ELL_3M_LERP(Tzp, alpha, T0, T1);
			ELL_3V_LERP(gzp, alpha, cell.Points[EdgeConnection[0][0]].gradient,
				cell.Points[EdgeConnection[0][1]].gradient);

			if (zi == DSampleNum - 2)
			{
				ELL_3M_COPY(Tzm, T); ELL_3V_COPY(gzm, g);

			}else{
				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 2]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tzm, alpha, T0, T1);
				ELL_3V_LERP(gzm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}
			ExtractT(point1.hessian, T0);
			ExtractT(point2.hessian, T1);
			ELL_3M_COPY(Tym, T0); ELL_3V_COPY(gym, point1.gradient);
			ELL_3M_COPY(Typ, T1); ELL_3V_COPY(gyp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 5:
			ExtractT(cell.Points[EdgeConnection[7][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[7][1]].hessian, T1);
			ELL_3M_LERP(Typ, alpha, T0, T1);
			ELL_3V_LERP(gyp, alpha, cell.Points[EdgeConnection[7][0]].gradient,
				cell.Points[EdgeConnection[7][1]].gradient);

			if (yi == HSampleNum - 2)
			{
				ELL_3M_COPY(Tym, T); ELL_3V_COPY(gym, g);
			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 2][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 2][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 2][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 2][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tym, alpha, T0, T1);
				ELL_3V_LERP(gym, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[1][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[1][1]].hessian, T1);
			ELL_3M_LERP(Tzp, alpha, T0, T1);
			ELL_3V_LERP(gzp, alpha, cell.Points[EdgeConnection[1][0]].gradient,
				cell.Points[EdgeConnection[1][1]].gradient);

			if (zi == DSampleNum - 2)
			{
				ELL_3M_COPY(Tzm, T); ELL_3V_COPY(gzm, g);

			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 2]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tzm, alpha, T0, T1);
				ELL_3V_LERP(gzm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(point1.hessian, T0);
			ExtractT(point2.hessian, T1);
			ELL_3M_COPY(Txm, T0); ELL_3V_COPY(gxm, point1.gradient);
			ELL_3M_COPY(Txp, T1); ELL_3V_COPY(gxp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 6:
			ExtractT(cell.Points[EdgeConnection[4][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[4][1]].hessian, T1);
			ELL_3M_LERP(Txp, alpha, T0, T1);
			ELL_3V_LERP(gxp, alpha, cell.Points[EdgeConnection[4][0]].gradient,
				cell.Points[EdgeConnection[4][1]].gradient);

			if (xi == WSampleNum - 2)
			{
				ELL_3M_COPY(Txm, T); ELL_3V_COPY(gxm, g);
			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Txm, alpha, T0, T1);
				ELL_3V_LERP(gxm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[2][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[2][1]].hessian, T1);
			ELL_3M_LERP(Tzp, alpha, T0, T1);
			ELL_3V_LERP(gzp, alpha, cell.Points[EdgeConnection[2][0]].gradient,
				cell.Points[EdgeConnection[2][1]].gradient);

			if (zi == DSampleNum - 2)
			{
				ELL_3M_COPY(Tzm, T); ELL_3V_COPY(gzm, g);

			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 2]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tzm, alpha, T0, T1);
				ELL_3V_LERP(gzm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}
			ExtractT(point1.hessian, T0);
			ExtractT(point2.hessian, T1);
			ELL_3M_COPY(Tym, T0); ELL_3V_COPY(gym, point1.gradient);
			ELL_3M_COPY(Typ, T1); ELL_3V_COPY(gyp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 7:
			ExtractT(cell.Points[EdgeConnection[5][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[5][1]].hessian, T1);
			ELL_3M_LERP(Typ, alpha, T0, T1);
			ELL_3V_LERP(gyp, alpha, cell.Points[EdgeConnection[5][0]].gradient,
				cell.Points[EdgeConnection[5][1]].gradient);

			if (yi == 0)
			{
				ELL_3M_COPY(Tym, T); ELL_3V_COPY(gym, g);
			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi - 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi - 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi - 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi - 1][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tym, alpha, T0, T1);
				ELL_3V_LERP(gym, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[3][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[3][1]].hessian, T1);
			ELL_3M_LERP(Tzp, alpha, T0, T1);
			ELL_3V_LERP(gzp, alpha, cell.Points[EdgeConnection[3][0]].gradient,
				cell.Points[EdgeConnection[3][1]].gradient);

			if (zi == DSampleNum - 2)
			{
				ELL_3M_COPY(Tzm, T); ELL_3V_COPY(gzm, g);

			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 2]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 2]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tzm, alpha, T0, T1);
				ELL_3V_LERP(gzm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}
			ExtractT(point1.hessian, T0);
			ExtractT(point2.hessian, T1);
			ELL_3M_COPY(Txm, T0); ELL_3V_COPY(gxm, point1.gradient);
			ELL_3M_COPY(Txp, T1); ELL_3V_COPY(gxp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 8:
			ExtractT(cell.Points[EdgeConnection[9][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[9][1]].hessian, T1);
			ELL_3M_LERP(Typ, alpha, T0, T1);
			ELL_3V_LERP(gyp, alpha, cell.Points[EdgeConnection[9][0]].gradient,
				cell.Points[EdgeConnection[9][1]].gradient);

			if (yi == 0)
			{
				ELL_3M_COPY(Tym, T); ELL_3V_COPY(gym, g);
			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi - 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi - 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi - 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi - 1][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tym, alpha, T0, T1);
				ELL_3V_LERP(gym, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[11][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[11][1]].hessian, T1);
			ELL_3M_LERP(Txp, alpha, T0, T1);
			ELL_3V_LERP(gxp, alpha, cell.Points[EdgeConnection[11][0]].gradient,
				cell.Points[EdgeConnection[11][1]].gradient);

			if (xi == 0)
			{
				ELL_3M_COPY(Txm, T); ELL_3V_COPY(gxm, g);

			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Txm, alpha, T0, T1);
				ELL_3V_LERP(gxm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}
			ExtractT(point1.hessian, T0);
			ExtractT(point2.hessian, T1);
			ELL_3M_COPY(Tzm, T0); ELL_3V_COPY(gzm, point1.gradient);
			ELL_3M_COPY(Tzp, T1); ELL_3V_COPY(gzp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 9:
			ExtractT(cell.Points[EdgeConnection[8][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[8][1]].hessian, T1);
			ELL_3M_LERP(Typ, alpha, T0, T1);
			ELL_3V_LERP(gyp, alpha, cell.Points[EdgeConnection[8][0]].gradient,
				cell.Points[EdgeConnection[8][1]].gradient);

			if (yi == HSampleNum - 2)
			{
				ELL_3M_COPY(Tym, T); ELL_3V_COPY(gym, g);
			}else{
				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 2][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 2][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 2][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 2][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tym, alpha, T0, T1);
				ELL_3V_LERP(gym, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[10][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[10][1]].hessian, T1);
			ELL_3M_LERP(Txp, alpha, T0, T1);
			ELL_3V_LERP(gxp, alpha, cell.Points[EdgeConnection[10][0]].gradient,
				cell.Points[EdgeConnection[10][1]].gradient);

			if (xi == 0)
			{
				ELL_3M_COPY(Txm, T); ELL_3V_COPY(gxm, g);

			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi - 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Txm, alpha, T0, T1);
				ELL_3V_LERP(gxm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}
			ExtractT(point1.hessian, T0);
			ExtractT(point2.hessian, T1);
			ELL_3M_COPY(Tzm, T0); ELL_3V_COPY(gzm, point1.gradient);
			ELL_3M_COPY(Tzp, T1); ELL_3V_COPY(gzp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 10:
			ExtractT(cell.Points[EdgeConnection[11][0]].hessian, T1);
			ExtractT(cell.Points[EdgeConnection[11][1]].hessian, T1);
			ELL_3M_LERP(Typ, alpha, T0, T1);
			ELL_3V_LERP(gyp, alpha, cell.Points[EdgeConnection[11][0]].gradient,
				cell.Points[EdgeConnection[11][1]].gradient);

			if (yi == HSampleNum - 2)
			{
				ELL_3M_COPY(Tym, T); ELL_3V_COPY(gym, g);
			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 2][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 2][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi + 2][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 2][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tym, alpha, T0, T1);
				ELL_3V_LERP(gym, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[9][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[9][1]].hessian, T1);
			ELL_3M_LERP(Txp, alpha, T0, T1);
			ELL_3V_LERP(gxp, alpha, cell.Points[EdgeConnection[9][0]].gradient,
				cell.Points[EdgeConnection[9][1]].gradient);

			if (xi == WSampleNum - 2)
			{
				ELL_3M_COPY(Txm, T); ELL_3V_COPY(gxm, g);

			}else{
				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Txm, alpha, T0, T1);
				ELL_3V_LERP(gxm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}
			ExtractT(point1.hessian, T0);
			ExtractT(point2.hessian, T1);
			ELL_3M_COPY(Tzm, T0); ELL_3V_COPY(gzm, point1.gradient);
			ELL_3M_COPY(Tzp, T1); ELL_3V_COPY(gzp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 11:
			ExtractT(cell.Points[EdgeConnection[10][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[10][1]].hessian, T1);
			ELL_3M_LERP(Typ, alpha, T0, T1);
			ELL_3V_LERP(gyp, alpha, cell.Points[EdgeConnection[10][0]].gradient,
				cell.Points[EdgeConnection[10][1]].gradient);

			if (yi == 0)
			{
				ELL_3M_COPY(Tym, T); ELL_3V_COPY(gym, g);
			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi][yi - 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi - 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi][yi - 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi - 1][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}
				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Tym, alpha, T0, T1);
				ELL_3V_LERP(gym, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}

			ExtractT(cell.Points[EdgeConnection[8][0]].hessian, T0);
			ExtractT(cell.Points[EdgeConnection[8][1]].hessian, T1);
			ELL_3M_LERP(Txp, alpha, T0, T1);
			ELL_3V_LERP(gxp, alpha, cell.Points[EdgeConnection[8][0]].gradient,
				cell.Points[EdgeConnection[8][1]].gradient);

			if (xi == WSampleNum - 2)
			{
				ELL_3M_COPY(Txm, T); ELL_3V_COPY(gxm, g);

			}else{

				GridCell tempCell;

				tempCell.degNum = 0;
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi + 1][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi][zi]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi + 1][zi + 1]);
				tempCell.Points.push_back(m_GridPoints[xi + 2][yi][zi + 1]);

				tempCell.edges.resize(12);
				for (int i = 0;i < 12;i++)
				{
					tempCell.edges[i].index = i;
					tempCell.edges[i].hasIntersection = false;
				}

				tempCell.faces.resize(6);
				for (int i = 0;i < 6;i++)
				{
					tempCell.faces[i].index = i;
					tempCell.faces[i].hasIntersection = false;
					tempCell.faces[i].interNum = 0;
				}	

				ExtractT(tempCell.Points[EdgeConnection[edgeId][0]].hessian, T0);
				ExtractT(tempCell.Points[EdgeConnection[edgeId][1]].hessian, T1);
				ELL_3M_LERP(Txm, alpha, T0, T1);
				ELL_3V_LERP(gxm, alpha, tempCell.Points[EdgeConnection[edgeId][0]].gradient,
					tempCell.Points[EdgeConnection[edgeId][1]].gradient);

			}
			ExtractT(point1.hessian, T0);
			ExtractT(point2.hessian, T1);
			ELL_3M_COPY(Tzm, T0); ELL_3V_COPY(gzm, point1.gradient);
			ELL_3M_COPY(Tzp, T1); ELL_3V_COPY(gzp, point2.gradient);

			ComputeGradientLin(cell.edges[edgeId].intersectionNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;	
	}
}

void ExtractSurface::ComputeGradientLin(double *result, double *T, double *g, double *Txm, double *gxm, double *Txp, double *gxp, double *Tym, double *gym, double *Typ, double *gyp, double *Tzm, double *gzm, double *Tzp, double *gzp)
{
	double Tder[9];
	double gder[3];
	double tmp[3], tmp1[3], tmp2[3];
	double derxv[3], deryv[3], derzv[3];
	ELL_3M_SUB(Tder,Txp,Txm);
	ELL_3V_SUB(gder,gxp,gxm);
	ell_3mv_mul_d(tmp,T,gder);
	ELL_3V_SUB(tmp,tmp,gder);
	ell_3mv_mul_d(derxv,Tder,g);
	ELL_3V_ADD2(derxv,derxv,tmp);

	ELL_3M_SUB(Tder,Typ,Tym);
	ELL_3V_SUB(gder,gyp,gym);
	ell_3mv_mul_d(tmp,T,gder);
	ELL_3V_SUB(tmp,tmp,gder);
	ell_3mv_mul_d(deryv,Tder,g);
	ELL_3V_ADD2(deryv,deryv,tmp);

	ELL_3M_SUB(Tder,Tzp,Tzm);
	ELL_3V_SUB(gder,gzp,gzm);
	ell_3mv_mul_d(tmp,T,gder);
	ELL_3V_SUB(tmp,tmp,gder);
	ell_3mv_mul_d(derzv,Tder,g);
	ELL_3V_ADD2(derzv,derzv,tmp);


	tmp1[0]=derxv[0]; tmp1[1]=deryv[0]; tmp1[2]=derzv[0];
	tmp2[0]=derxv[1]; tmp2[1]=deryv[1]; tmp2[2]=derzv[1];
	if (ELL_3V_DOT(tmp1,tmp2)<0)
		ELL_3V_SCALE(tmp2,-1.0,tmp2);
	ELL_3V_ADD2(tmp1,tmp1,tmp2);
	tmp2[0]=derxv[2]; tmp2[1]=deryv[2]; tmp2[2]=derzv[2];
	if (ELL_3V_DOT(tmp1,tmp2)<0)
		ELL_3V_SCALE(tmp2,-1.0,tmp2);
	ELL_3V_ADD2(result,tmp1,tmp2);

	double len = ELL_3V_LEN(result);

    if (len > 1.0)
    {
		ELL_3V_SCALE(result, 1.0/len, result);
    }

}



/* 判断 crease surface 与两点之间的棱是否有交点 */
bool ExtractSurface::IsIntersected(GridPoint& point1, GridPoint& point2){
	double h1[3], h2[3];
	Computeh(point1, h1);
	Computeh(point2, h2);


	return ELL_3V_DOT(h1,h2) < 0;
}

//计算边交点坐标
void ExtractSurface::EvaluateEdgeInterCoord(GridCell& cell, int edgeid){
	GridPoint point1 = cell.Points[EdgeConnection[edgeid][0]],
		point2 = cell.Points[EdgeConnection[edgeid][1]];

	cell.edges[edgeid].intersectionCoord[0]= (point1.x + point2.x)/2;
	cell.edges[edgeid].intersectionCoord[1] = (point1.y + point2.y)/2;
	cell.edges[edgeid].intersectionCoord[2] = (point1.z + point2.z)/2;

	return;
}


void ExtractSurface::Computeh(GridPoint& point, double d1[3]){
	double T[9];

	ExtractT(point.hessian, T);
	ell_3mv_mul_d(d1, T, point.gradient);
	ELL_3V_SUB(d1, d1, point.gradient);
}


void ExtractSurface::ExtractT(Hessian& hessian, double T[9]) {

	double lambdas[3]={0.0,0.0,0.0};
	double tmpMat[9], diag[9], evecsT[9];
	double hessnew[9], evecs[9], evals[3];

	hessnew[0] = hessian.fxx;
	hessnew[1] = hessian.fxy;
	hessnew[2] = hessian.fxz;
	hessnew[3] = hessian.fyx;
	hessnew[4] = hessian.fyy;
	hessnew[5] = hessian.fyz;
	hessnew[6] = hessian.fzx;
	hessnew[7] = hessian.fzy;
	hessnew[8] = hessian.fzz;

	ell_3m_eigensolve_d(evals, evecs, hessnew, AIR_TRUE);
	

	double diff = evals[1]-evals[2];
	lambdas[0] = lambdas[1] = 1.0;
	if (diff < EVALDIFFTHRESH)
		lambdas[2]=(1.0 - diff/EVALDIFFTHRESH)*(1.0 - diff/EVALDIFFTHRESH);
	else{
		lambdas[2] = 0.0;
	}

	ELL_3M_ZERO_SET(diag);
	ELL_3M_DIAG_SET(diag, lambdas[0], lambdas[1], lambdas[2]);
	ELL_3M_TRANSPOSE(evecsT, evecs);
	ELL_3M_MUL(tmpMat, diag, evecs);
	ELL_3M_MUL(T, evecsT, tmpMat);

}


/* 判断立方体的面是否有交点，若交点存在，则求取相应的交点坐标值 */
void ExtractSurface::FindIntersectionPerFace(GridCell &cell, int faceId){

	//cell.faces[i].hasIntersection = NewtonMethod(cell, i , sameCoord, centerPointCoord1, centerPointCoord2);

	cell.faces[faceId].hasIntersection = true;
	cell.degNum++;

	if (faceId % 3 == 0)
	{
		cell.faces[faceId].degenerationCoord[0] = (cell.Points[FaceConnection[faceId][0]].x 
			+ cell.Points[FaceConnection[faceId][1]].x) / 2;
		cell.faces[faceId].degenerationCoord[1] = (cell.Points[FaceConnection[faceId][0]].y 
			+ cell.Points[FaceConnection[faceId][3]].y) / 2;
		cell.faces[faceId].degenerationCoord[2] = cell.Points[FaceConnection[faceId][0]].z;
	}else if (faceId % 3 == 1)
	{
		cell.faces[faceId].degenerationCoord[0] = cell.Points[FaceConnection[faceId][0]].x;
		cell.faces[faceId].degenerationCoord[1] = (cell.Points[FaceConnection[faceId][0]].y 
			+ cell.Points[FaceConnection[faceId][1]].y) / 2;
		cell.faces[faceId].degenerationCoord[2] = (cell.Points[FaceConnection[faceId][0]].z 
			+ cell.Points[FaceConnection[faceId][3]].z) / 2;
	}else{
		cell.faces[faceId].degenerationCoord[0] = (cell.Points[FaceConnection[faceId][0]].x 
			+ cell.Points[FaceConnection[faceId][1]].x) / 2;
		cell.faces[faceId].degenerationCoord[1] = cell.Points[FaceConnection[faceId][0]].y;
		cell.faces[faceId].degenerationCoord[2] = (cell.Points[FaceConnection[faceId][0]].z
			+ cell.Points[FaceConnection[faceId][3]].z) / 2;
	}

}

/* 以每个面的中点为初值，通过改进的牛顿迭代法更新坐标值，若在规定迭代次数内收敛，
*则交点存在，否者视为无交点。Topological Lines in 3D Tensor Fields 文章中公式(10)、(11)中
提到的思路。
*/
// bool ExtractSurface::NewtonMethod(GridCell &cell, int index, int sameCoord, double& coord1, double& coord2){
// 	int iterNum = 0;
// 	double differ = DBL_MAX;
// 	while (differ > THRESHOLD && iterNum < INTERATION)
// 	{
// 		vector<double> CFvalVec1 = BilinearInterpolation(cell, index, sameCoord, coord1 - DELTA, coord2);
// 		vector<double> CFvalVec2 = BilinearInterpolation(cell, index, sameCoord, coord1, coord2 - DELTA);
// 		vector<double> CFvalVec3 = BilinearInterpolation(cell, index, sameCoord, coord1 + DELTA, coord2);
// 		vector<double> CFvalVec4 = BilinearInterpolation(cell, index, sameCoord, coord1, coord2 + DELTA);
// 
// 		MatrixXd m1(2,7);
// 		m1 << (CFvalVec3[0] - CFvalVec1[0]) / (2 * DELTA), (CFvalVec3[1] - CFvalVec1[1]) / (2 * DELTA), (CFvalVec3[2] - CFvalVec1[2]) / (2 * DELTA),
// 			(CFvalVec3[3] - CFvalVec1[3]) / (2 * DELTA), (CFvalVec3[4] - CFvalVec1[4]) / (2 * DELTA), (CFvalVec3[5] - CFvalVec1[5]) / (2 * DELTA),
// 			(CFvalVec3[6] - CFvalVec1[6]) / (2 * DELTA), (CFvalVec4[0] - CFvalVec2[0]) / (2 * DELTA), (CFvalVec4[1] - CFvalVec2[1]) / (2 * DELTA),
// 			(CFvalVec4[2] - CFvalVec2[2]) / (2 * DELTA), (CFvalVec4[3] - CFvalVec2[3]) / (2 * DELTA), (CFvalVec4[4] - CFvalVec2[4]) / (2 * DELTA),
// 			(CFvalVec4[5] - CFvalVec2[5]) / (2 * DELTA), (CFvalVec4[6] - CFvalVec2[6]) / (2 * DELTA);
// 
// 		MatrixXd m2(2,2);
// 		m2 = m1 * m1.transpose();
// 		m2 = m2.inverse();
// 
// 		VectorXd v(7);
// 		vector<double> CFvalVec = BilinearInterpolation(cell, index, sameCoord, coord1, coord2);
// 		v << CFvalVec[0], CFvalVec[1], CFvalVec[2], CFvalVec[3], CFvalVec[4], CFvalVec[5], CFvalVec[6];
// 
// 		VectorXd resultV(2);
// 		resultV = m2 * m1 * v;
// 
// 		coord1 -= resultV(0);
// 		coord2 -= resultV(1);
// 
// 		differ = abs(resultV(0) * resultV(0) + resultV(1) * resultV(1));
// 		iterNum++;
// 	}
// 
// 	if (differ < THRESHOLD)
// 	{
// 		return true;
// 	}
// 
// 	return false;
// }

/* 确定面上点坐标中相等的一维,返回值为1说明面与x轴垂直，返回值为2说明面与y轴垂直，返回值为3说明与z轴垂直 */
int ExtractSurface::DetermineCoordinate(GridCell &cell, int index, double& coord1, double& coord2){

	GridPoint point1 = cell.Points[FaceConnection[index][0]];
	GridPoint point2 = cell.Points[FaceConnection[index][1]];
	GridPoint point3 = cell.Points[FaceConnection[index][2]];
	GridPoint point4 = cell.Points[FaceConnection[index][3]];

	if (abs(point1.x - point3.x) < LITTLEVAL && abs(point2.x - point4.x) < LITTLEVAL)
	{
		coord1 = (point1.y + point2.y) / 2;
		coord2 = (point2.z + point3.z) / 2;
		return 1;
	}

	if (abs(point1.y - point3.y) < LITTLEVAL && abs(point2.y - point4.y) < LITTLEVAL)
	{
		coord1 = (point1.x + point2.x) / 2;
		coord2 = (point2.z + point3.z) / 2;
		return 2;
	}

	if (abs(point1.z - point3.z) < LITTLEVAL && abs(point2.z - point4.z) < LITTLEVAL)
	{
		coord1 = (point1.x + point2.x) / 2;
		coord2 = (point2.y + point3.y) / 2;
		return 3;
	}

	return 0;
} 

/* 通过双线性插值法，以及面的四个顶点值计算面内部的任意点的CF值。*/
vector<double> ExtractSurface::BilinearInterpolation(GridCell &cell, int index, int sameCoord, double coord1, double coord2){

	vector<double> valVec;
	GridPoint point1 = cell.Points[FaceConnection[index][0]];
	GridPoint point2 = cell.Points[FaceConnection[index][1]];
	GridPoint point3 = cell.Points[FaceConnection[index][2]];
	GridPoint point4 = cell.Points[FaceConnection[index][3]];
	vector<double> CFvec1 = CF(point1.hessian);
	vector<double> CFvec2 = CF(point2.hessian);
	vector<double> CFvec3 = CF(point3.hessian);
	vector<double> CFvec4 = CF(point4.hessian);

	for (int i = 0; i < 7; i++)
	{
		double a1 = 0.0, a2 = 0.0, b1 = 0.0, b2 = 0.0;
		if (sameCoord == 1)
		{
			a1 = point1.y;
			a2 = point2.y;
			b1 = point1.z;
			b2 = point3.z;
		}else if (sameCoord == 2)
		{
			a1 = point1.x;
			a2 = point2.x;
			b1 = point1.z;
			b2 = point3.z;
		}else{
			a1 = point1.x;
			a2 = point2.x;
			b1 = point1.y;
			b2 = point3.y;
		}

		double val = CFvec1[i] * (a2 - coord1) * (b2 - coord2) + CFvec2[i] * (coord1 - a1) * (b2 - coord2)
			+ CFvec3[i] * (coord1 - a1) * (coord2 - b1) + CFvec4[i] * (a2 - coord1) * (coord2 - b1);
		val = val / (a2 - a1) * (b2 - b1);
		valVec.push_back(val);

	}

	return valVec;		
}

/* 通过点的黑塞矩阵求取相应的约束方程向量值 */
vector<double> ExtractSurface::CF(Hessian &hessian){
	vector<double> valVec;
	double result = hessian.fxx * (hessian.fyy * hessian.fyy - hessian.fzz * hessian.fzz) +
		hessian.fxx * (hessian.fxy * hessian.fxy - hessian.fxz * hessian.fxz) +
		hessian.fyy * (hessian.fzz * hessian.fzz - hessian.fxx * hessian.fxx) +
		hessian.fyy * (hessian.fyz * hessian.fyz - hessian.fxy * hessian.fxy) +
		hessian.fzz * (hessian.fxx * hessian.fxx - hessian.fyy * hessian.fyy) +
		hessian.fzz * (hessian.fxz * hessian.fxz - hessian.fyz * hessian.fyz);
	valVec.push_back(result);

	result = hessian.fyz * (2 * (hessian.fyz * hessian.fyz - hessian.fxx * hessian.fxx) 
		- (hessian.fxz * hessian.fxz + hessian.fxy * hessian.fxy) 
		+ 2 * (hessian.fyy * hessian.fxx + hessian.fzz * hessian.fxx
		- hessian.fyy * hessian.fzz)) + hessian.fxy * hessian.fxz * (2 * hessian.fxx -
		hessian.fzz - hessian.fyy);
	valVec.push_back(result);

	result = hessian.fxz * (2 * (hessian.fxz * hessian.fxz - hessian.fyy * hessian.fyy) 
		- (hessian.fyz * hessian.fyz + hessian.fxy * hessian.fxy) 
		+ 2 * (hessian.fyy * hessian.fzz + hessian.fyy * hessian.fxx
		- hessian.fxx * hessian.fzz)) + hessian.fxy * hessian.fyz * (2 * hessian.fyy -
		hessian.fxx - hessian.fzz);
	valVec.push_back(result);

	result = hessian.fxy * (2 * (hessian.fxy * hessian.fxy - hessian.fzz * hessian.fzz) 
		- (hessian.fyz * hessian.fyz + hessian.fxz * hessian.fxz) 
		+ 2 * (hessian.fxx * hessian.fzz + hessian.fyy * hessian.fzz
		- hessian.fxx * hessian.fyy)) + hessian.fxz * hessian.fyz * (2 * hessian.fzz -
		hessian.fyy - hessian.fxx);
	valVec.push_back(result);

	result = hessian.fyz * (hessian.fxz * hessian.fxz + hessian.fxy * hessian.fxy)
		+ hessian.fxy * hessian.fxz * (hessian.fyy - hessian.fzz);
	valVec.push_back(result);

	result = hessian.fxz * (hessian.fxy * hessian.fxy + hessian.fyz * hessian.fyz)
		+ hessian.fyz * hessian.fxy * (hessian.fzz - hessian.fxx);
	valVec.push_back(result);

	result = hessian.fxy * (hessian.fyz * hessian.fyz + hessian.fxz * hessian.fxz)
		+ hessian.fxz * hessian.fyz * (hessian.fxx - hessian.fyy);
	valVec.push_back(result);


	return valVec;
}

void ExtractSurface::ComputeFaceNormal(GridCell &cell, int faceId, int xi, int yi ,int zi){
	double T[9], Txm[9], Txp[9], Tym[9], Typ[9], Tzm[9], Tzp[9],
		g[3], gxm[3], gxp[3], gym[3], gyp[3], gzm[3], gzp[3], T0[9], T1[9], T2[9], T3[9];
	double *coords = cell.faces[faceId].degenerationCoord;
	switch(faceId % 3){
		case 0:
			ExtractT(cell.Points[FaceConnection[faceId][0]].hessian, T0);
			ExtractT(cell.Points[FaceConnection[faceId][1]].hessian, T1);
			ExtractT(cell.Points[FaceConnection[faceId][2]].hessian, T2);
			ExtractT(cell.Points[FaceConnection[faceId][3]].hessian, T3);
			ELL_3M_LERP(Tym, coords[0], T0, T1);
			ELL_3V_LERP(gym, coords[0], cell.Points[FaceConnection[faceId][0]].gradient,
				cell.Points[FaceConnection[faceId][1]].gradient);
			ELL_3M_LERP(Typ, coords[0], T3, T2);
			ELL_3V_LERP(gyp, coords[0], cell.Points[FaceConnection[faceId][3]].gradient,
				cell.Points[FaceConnection[faceId][2]].gradient);

			ELL_3M_LERP(Txm, coords[1], T0, T3);
			ELL_3V_LERP(gxm, coords[1], cell.Points[FaceConnection[faceId][0]].gradient,
				cell.Points[FaceConnection[faceId][3]].gradient);
			ELL_3M_LERP(Txp, coords[1], T1, T2);
			ELL_3V_LERP(gxp, coords[1], cell.Points[FaceConnection[faceId][1]].gradient,
				cell.Points[FaceConnection[faceId][2]].gradient);

			ELL_3M_LERP(T, coords[0], Txm, Txp);
			ELL_3V_LERP(g, coords[0], gxm, gxp);

			ExtractT(cell.Points[FaceConnection[3][0]].hessian, T0);
			ExtractT(cell.Points[FaceConnection[3][1]].hessian, T1);
			ExtractT(cell.Points[FaceConnection[3][2]].hessian, T2);
			ExtractT(cell.Points[FaceConnection[3][3]].hessian, T3);
			if (faceId == 0)
			{
				double	tmTxm[9], tmTxp[9], tmgxm[3], tmgxp[3];

				ELL_3M_LERP(tmTxm, coords[1], T0, T3);
				ELL_3V_LERP(tmgxm, coords[1], cell.Points[FaceConnection[3][0]].gradient,
					cell.Points[FaceConnection[3][3]].gradient);
				ELL_3M_LERP(tmTxp, coords[1], T1, T2);
				ELL_3V_LERP(tmgxp, coords[1], cell.Points[FaceConnection[3][1]].gradient,
					cell.Points[FaceConnection[3][2]].gradient);

				ELL_3M_LERP(Tzp, coords[0], tmTxm, tmTxp);
				ELL_3V_LERP(gzp, coords[0], tmgxm, tmgxp);

				if (zi == 0)
				{
					ELL_3M_COPY(Tzm, T); ELL_3V_COPY(gzm, g);
				}else{

					double tmTxm[9], tmTxp[9], tmgxm[3], tmgxp[3];
					GridCell tempCell;

					tempCell.degNum = 0;
					tempCell.Points.push_back(m_GridPoints[xi][yi][zi - 1]);
					tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi - 1]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi - 1]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi- 1]);
					tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
					tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);

					tempCell.edges.resize(12);
					for (int i = 0;i < 12;i++)
					{
						tempCell.edges[i].index = i;
						tempCell.edges[i].hasIntersection = false;
					}

					tempCell.faces.resize(6);
					for (int i = 0;i < 6;i++)
					{
						tempCell.faces[i].index = i;
						tempCell.faces[i].hasIntersection = false;
						tempCell.faces[i].interNum = 0;
					}

					ExtractT(tempCell.Points[FaceConnection[faceId][0]].hessian, T0);
					ExtractT(tempCell.Points[FaceConnection[faceId][1]].hessian, T1);
					ExtractT(tempCell.Points[FaceConnection[faceId][2]].hessian, T2);
					ExtractT(tempCell.Points[FaceConnection[faceId][3]].hessian, T3);

					ELL_3M_LERP(tmTxm, coords[1], T0, T3);
					ELL_3V_LERP(tmgxm, coords[1], tempCell.Points[FaceConnection[faceId][0]].gradient,
						tempCell.Points[FaceConnection[faceId][3]].gradient);
					ELL_3M_LERP(tmTxp, coords[1], T1, T2);
					ELL_3V_LERP(tmgxp, coords[1], tempCell.Points[FaceConnection[faceId][1]].gradient,
						tempCell.Points[FaceConnection[faceId][2]].gradient);

					ELL_3M_LERP(Tzm, coords[0], tmTxm, tmTxp);
					ELL_3V_LERP(gzm, coords[0], tmgxm, tmgxp);

				}

			}else{
				double	tmTxm[9], tmTxp[9], tmgxm[3], tmgxp[3];

				ExtractT(cell.Points[FaceConnection[0][0]].hessian, T0);
				ExtractT(cell.Points[FaceConnection[0][1]].hessian, T1);
				ExtractT(cell.Points[FaceConnection[0][2]].hessian, T2);
				ExtractT(cell.Points[FaceConnection[0][3]].hessian, T3);

				ELL_3M_LERP(tmTxm, coords[1], T0, T3);
				ELL_3V_LERP(tmgxm, coords[1], cell.Points[FaceConnection[0][0]].gradient,
					cell.Points[FaceConnection[0][3]].gradient);
				ELL_3M_LERP(tmTxp, coords[1], T1, T2);
				ELL_3V_LERP(tmgxp, coords[1], cell.Points[FaceConnection[0][1]].gradient,
					cell.Points[FaceConnection[0][2]].gradient);

				ELL_3M_LERP(Tzp, coords[0], tmTxm, tmTxp);
				ELL_3V_LERP(gzp, coords[0], tmgxm, tmgxp);

				if (zi == DSampleNum - 2)
				{
					ELL_3M_COPY(Tzm, T); ELL_3V_COPY(gzm, g);
				}else{
					double	tmTxm[9], tmTxp[9], tmgxm[3], tmgxp[3];
					GridCell tempCell;

					tempCell.degNum = 0;
					tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 2]);
					tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 2]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 2]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 2]);

					tempCell.edges.resize(12);
					for (int i = 0;i < 12;i++)
					{
						tempCell.edges[i].index = i;
						tempCell.edges[i].hasIntersection = false;
					}

					tempCell.faces.resize(6);
					for (int i = 0;i < 6;i++)
					{
						tempCell.faces[i].index = i;
						tempCell.faces[i].hasIntersection = false;
						tempCell.faces[i].interNum = 0;
					}


					ExtractT(tempCell.Points[FaceConnection[faceId][0]].hessian, T0);
					ExtractT(tempCell.Points[FaceConnection[faceId][1]].hessian, T1);
					ExtractT(tempCell.Points[FaceConnection[faceId][2]].hessian, T2);
					ExtractT(tempCell.Points[FaceConnection[faceId][3]].hessian, T3);

					ELL_3M_LERP(tmTxm, coords[1], T0, T3);
					ELL_3V_LERP(tmgxm, coords[1], tempCell.Points[FaceConnection[faceId][0]].gradient,
						tempCell.Points[FaceConnection[faceId][3]].gradient);
					ELL_3M_LERP(tmTxp, coords[1], T1, T2);
					ELL_3V_LERP(tmgxp, coords[1], tempCell.Points[FaceConnection[faceId][1]].gradient,
						tempCell.Points[FaceConnection[faceId][2]].gradient);

					ELL_3M_LERP(Tzm, coords[0], tmTxm, tmTxp);
					ELL_3V_LERP(gzm, coords[0], tmgxm, tmgxp);				
				}
				ComputeGradientLin(cell.faces[faceId].degenerationNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
				break;
		case 1:
			ExtractT(cell.Points[FaceConnection[faceId][0]].hessian, T0);
			ExtractT(cell.Points[FaceConnection[faceId][1]].hessian, T1);
			ExtractT(cell.Points[FaceConnection[faceId][2]].hessian, T2);
			ExtractT(cell.Points[FaceConnection[faceId][3]].hessian, T3);
			ELL_3M_LERP(Tym, coords[2], T0, T3);
			ELL_3V_LERP(gym, coords[2], cell.Points[FaceConnection[faceId][0]].gradient,
				cell.Points[FaceConnection[faceId][3]].gradient);
			ELL_3M_LERP(Typ, coords[2], T1 ,T2);
			ELL_3V_LERP(gyp, coords[2], cell.Points[FaceConnection[faceId][1]].gradient,
				cell.Points[FaceConnection[faceId][2]].gradient);

			ELL_3M_LERP(Tzm, coords[1], T0, T1);
			ELL_3V_LERP(gzm, coords[1], cell.Points[FaceConnection[faceId][0]].gradient,
				cell.Points[FaceConnection[faceId][1]].gradient);
			ELL_3M_LERP(Tzp, coords[1], T3, T2);
			ELL_3V_LERP(gzp, coords[1], cell.Points[FaceConnection[faceId][3]].gradient,
				cell.Points[FaceConnection[faceId][2]].gradient);

			ELL_3M_LERP(T, coords[1], Tym, Typ);
			ELL_3V_LERP(g, coords[1], gym, gyp);

			if (faceId == 1)
			{
				double	tmTym[9], tmTyp[9], tmgym[3], tmgyp[3];

				ExtractT(cell.Points[FaceConnection[4][0]].hessian, T0);
				ExtractT(cell.Points[FaceConnection[4][1]].hessian, T1);
				ExtractT(cell.Points[FaceConnection[4][2]].hessian, T2);
				ExtractT(cell.Points[FaceConnection[4][3]].hessian, T3);
				ELL_3M_LERP(tmTym, coords[2], T0, T3);
				ELL_3V_LERP(tmgym, coords[2], cell.Points[FaceConnection[4][0]].gradient,
					cell.Points[FaceConnection[4][3]].gradient);
				ELL_3M_LERP(tmTyp, coords[2], T1, T2);
				ELL_3V_LERP(tmgyp, coords[2], cell.Points[FaceConnection[4][1]].gradient,
					cell.Points[FaceConnection[4][2]].gradient);

				ELL_3M_LERP(Txp, coords[1], tmTym, tmTyp);
				ELL_3V_LERP(gxp, coords[1], tmgym, tmgyp);



				if (xi == 0)
				{
					ELL_3M_COPY(Txm, T); ELL_3V_COPY(gxm, g);
				}else{
					double	tmTxm[9], tmTxp[9], tmgxm[3], tmgxp[3];

					GridCell tempCell;

					tempCell.degNum = 0;
					tempCell.Points.push_back(m_GridPoints[xi - 1][yi][zi]);
					tempCell.Points.push_back(m_GridPoints[xi - 1][yi + 1][zi]);
					tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
					tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
					tempCell.Points.push_back(m_GridPoints[xi - 1][yi][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi - 1][yi + 1][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);

					tempCell.edges.resize(12);
					for (int i = 0;i < 12;i++)
					{
						tempCell.edges[i].index = i;
						tempCell.edges[i].hasIntersection = false;
					}

					tempCell.faces.resize(6);
					for (int i = 0;i < 6;i++)
					{
						tempCell.faces[i].index = i;
						tempCell.faces[i].hasIntersection = false;
						tempCell.faces[i].interNum = 0;
					}

					ExtractT(tempCell.Points[FaceConnection[faceId][0]].hessian, T0);
					ExtractT(tempCell.Points[FaceConnection[faceId][1]].hessian, T1);
					ExtractT(tempCell.Points[FaceConnection[faceId][2]].hessian, T2);
					ExtractT(tempCell.Points[FaceConnection[faceId][3]].hessian, T3);

					ELL_3M_LERP(tmTym, coords[2], T0, T3);
					ELL_3V_LERP(tmgym, coords[2], tempCell.Points[FaceConnection[faceId][0]].gradient,
						tempCell.Points[FaceConnection[faceId][3]].gradient);
					ELL_3M_LERP(tmTyp, coords[2], T1, T2);
					ELL_3V_LERP(tmgyp, coords[2], tempCell.Points[FaceConnection[faceId][1]].gradient,
						tempCell.Points[FaceConnection[faceId][2]].gradient);

					ELL_3M_LERP(Txm, coords[1], tmTym, tmTyp);
					ELL_3V_LERP(gxm, coords[1], tmgym, tmgyp);
				}

			}else{
				double	tmTym[9], tmTyp[9], tmgym[3], tmgyp[3];

				ExtractT(cell.Points[FaceConnection[1][0]].hessian, T0);
				ExtractT(cell.Points[FaceConnection[1][1]].hessian, T1);
				ExtractT(cell.Points[FaceConnection[1][2]].hessian, T2);
				ExtractT(cell.Points[FaceConnection[1][3]].hessian, T3);
				ELL_3M_LERP(tmTym, coords[2], T0, T3);
				ELL_3V_LERP(tmgym, coords[2], cell.Points[FaceConnection[1][0]].gradient,
					cell.Points[FaceConnection[1][3]].gradient);
				ELL_3M_LERP(tmTyp, coords[2], T1, T2);
				ELL_3V_LERP(tmgyp, coords[2], cell.Points[FaceConnection[1][1]].gradient,
					cell.Points[FaceConnection[1][2]].gradient);

				ELL_3M_LERP(Txp, coords[1], tmTym, tmTyp);
				ELL_3V_LERP(gxp, coords[1], tmgym, tmgyp);



				if (xi == WSampleNum - 2)
				{
					ELL_3M_COPY(Txm, T); ELL_3V_COPY(gxm, g);
				}else{
					double	tmTxm[9], tmTxp[9], tmgxm[3], tmgxp[3];
					GridCell tempCell;

					tempCell.degNum = 0;
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
					tempCell.Points.push_back(m_GridPoints[xi + 2][yi + 1][zi]);
					tempCell.Points.push_back(m_GridPoints[xi + 2][yi][zi]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi + 2][yi + 1][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi + 2][yi][zi + 1]);

					tempCell.edges.resize(12);
					for (int i = 0;i < 12;i++)
					{
						tempCell.edges[i].index = i;
						tempCell.edges[i].hasIntersection = false;
					}

					tempCell.faces.resize(6);
					for (int i = 0;i < 6;i++)
					{
						tempCell.faces[i].index = i;
						tempCell.faces[i].hasIntersection = false;
						tempCell.faces[i].interNum = 0;
					}

					ExtractT(tempCell.Points[FaceConnection[faceId][0]].hessian, T0);
					ExtractT(tempCell.Points[FaceConnection[faceId][1]].hessian, T1);
					ExtractT(tempCell.Points[FaceConnection[faceId][2]].hessian, T2);
					ExtractT(tempCell.Points[FaceConnection[faceId][3]].hessian, T3);

					ELL_3M_LERP(tmTym, coords[2], T0, T3);
					ELL_3V_LERP(tmgym, coords[2], tempCell.Points[FaceConnection[faceId][0]].gradient,
						tempCell.Points[FaceConnection[faceId][3]].gradient);
					ELL_3M_LERP(tmTyp, coords[2], T1, T2);
					ELL_3V_LERP(tmgyp, coords[2], tempCell.Points[FaceConnection[faceId][1]].gradient,
						tempCell.Points[FaceConnection[faceId][2]].gradient);

					ELL_3M_LERP(Txm, coords[1], tmTym, tmTyp);
					ELL_3V_LERP(gxm, coords[1], tmgym, tmgyp);
				}			
			}
			ComputeGradientLin(cell.faces[faceId].degenerationNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;
		case 2:

			ExtractT(cell.Points[FaceConnection[faceId][0]].hessian, T0);
			ExtractT(cell.Points[FaceConnection[faceId][1]].hessian, T1);
			ExtractT(cell.Points[FaceConnection[faceId][2]].hessian, T2);
			ExtractT(cell.Points[FaceConnection[faceId][3]].hessian, T3);

			ELL_3M_LERP(Txm, coords[2], T0, T3);
			ELL_3V_LERP(gxm, coords[2], cell.Points[FaceConnection[faceId][0]].gradient,
				cell.Points[FaceConnection[faceId][3]].gradient);
			ELL_3M_LERP(Txp, coords[2], T1, T2);
			ELL_3V_LERP(gxp, coords[2], cell.Points[FaceConnection[faceId][1]].gradient,
				cell.Points[FaceConnection[faceId][2]].gradient);
			ExtractT(cell.Points[FaceConnection[faceId][0]].hessian, T0);
			ExtractT(cell.Points[FaceConnection[faceId][1]].hessian, T1);
			ExtractT(cell.Points[FaceConnection[faceId][2]].hessian, T2);
			ExtractT(cell.Points[FaceConnection[faceId][3]].hessian, T3);

			ELL_3M_LERP(Tzm, coords[0], T0, T1);
			ELL_3V_LERP(gzm, coords[0], cell.Points[FaceConnection[faceId][0]].gradient,
				cell.Points[FaceConnection[faceId][1]].gradient);
			ELL_3M_LERP(Tzp, coords[0], T3, T2);
			ELL_3V_LERP(gzp, coords[0], cell.Points[FaceConnection[faceId][3]].gradient,
				cell.Points[FaceConnection[faceId][2]].gradient);

			ELL_3M_LERP(T, coords[0], Txm, Txp);
			ELL_3V_LERP(g, coords[0], gxm, gxp);

			if (faceId == 2)
			{
				double	tmTxm[9], tmTxp[9], tmgxm[3], tmgxp[3];
				ExtractT(cell.Points[FaceConnection[5][0]].hessian, T0);
				ExtractT(cell.Points[FaceConnection[5][1]].hessian, T1);
				ExtractT(cell.Points[FaceConnection[5][2]].hessian, T2);
				ExtractT(cell.Points[FaceConnection[5][3]].hessian, T3);

				ELL_3M_LERP(tmTxm, coords[2], T0, T3);
				ELL_3V_LERP(tmgxm, coords[2], cell.Points[FaceConnection[5][0]].gradient,
					cell.Points[FaceConnection[5][3]].gradient);
				ELL_3M_LERP(tmTxp, coords[2], T1, T2);
				ELL_3V_LERP(tmgxp, coords[2], cell.Points[FaceConnection[5][1]].gradient,
					cell.Points[FaceConnection[5][2]].gradient);

				ELL_3M_LERP(Typ, coords[0], tmTxm, tmTxp);
				ELL_3V_LERP(gyp, coords[0], tmgxm, tmgxp);



				if (yi == 0)
				{
					ELL_3M_COPY(Tym, T); ELL_3V_COPY(gym, g);
				}else{
					double	tmTxm[9], tmTxp[9], tmgxm[3], tmgxp[3];

					GridCell tempCell;

					tempCell.degNum = 0;
					tempCell.Points.push_back(m_GridPoints[xi][yi - 1][zi]);
					tempCell.Points.push_back(m_GridPoints[xi][yi][zi]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi - 1][zi]);
					tempCell.Points.push_back(m_GridPoints[xi][yi - 1][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi][yi][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi - 1][zi + 1]);

					tempCell.edges.resize(12);
					for (int i = 0;i < 12;i++)
					{
						tempCell.edges[i].index = i;
						tempCell.edges[i].hasIntersection = false;
					}

					tempCell.faces.resize(6);
					for (int i = 0;i < 6;i++)
					{
						tempCell.faces[i].index = i;
						tempCell.faces[i].hasIntersection = false;
						tempCell.faces[i].interNum = 0;
					}

					ExtractT(tempCell.Points[FaceConnection[faceId][0]].hessian, T0);
					ExtractT(tempCell.Points[FaceConnection[faceId][1]].hessian, T1);
					ExtractT(tempCell.Points[FaceConnection[faceId][2]].hessian, T2);
					ExtractT(tempCell.Points[FaceConnection[faceId][3]].hessian, T3);

					ELL_3M_LERP(tmTxm, coords[2], T0, T3);
					ELL_3V_LERP(tmgxm, coords[2], tempCell.Points[FaceConnection[faceId][0]].gradient,
						tempCell.Points[FaceConnection[faceId][3]].gradient);
					ELL_3M_LERP(tmTxp, coords[2], T1, T2);
					ELL_3V_LERP(tmgxp, coords[2], tempCell.Points[FaceConnection[faceId][1]].gradient,
						tempCell.Points[FaceConnection[faceId][2]].gradient);

					ELL_3M_LERP(Tym, coords[0], tmTxm, tmTxp);
					ELL_3V_LERP(gym, coords[0], tmgxm, tmgxp);
				}

			}else{
				double	tmTxm[9], tmTxp[9], tmgxm[3], tmgxp[3];
				ExtractT(cell.Points[FaceConnection[2][0]].hessian, T0);
				ExtractT(cell.Points[FaceConnection[2][1]].hessian, T1);
				ExtractT(cell.Points[FaceConnection[2][2]].hessian, T2);
				ExtractT(cell.Points[FaceConnection[2][3]].hessian, T3);

				ELL_3M_LERP(tmTxm, coords[2], T0, T3);
				ELL_3V_LERP(tmgxm, coords[2], cell.Points[FaceConnection[2][0]].gradient,
					cell.Points[FaceConnection[2][3]].gradient);
				ELL_3M_LERP(tmTxp, coords[2], T1, T2);
				ELL_3V_LERP(tmgxp, coords[2], cell.Points[FaceConnection[2][1]].gradient,
					cell.Points[FaceConnection[2][2]].gradient);

				ELL_3M_LERP(Typ, coords[0], tmTxm, tmTxp);
				ELL_3V_LERP(gyp, coords[0], tmgxm, tmgxp);


				if (yi == HSampleNum - 2)
				{
					ELL_3M_COPY(Tym, T); ELL_3V_COPY(gym, g);
				}else{
					double	tmTxm[9], tmTxp[9], tmgxm[3], tmgxp[3];
					GridCell tempCell;

					tempCell.degNum = 0;
					tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi]);
					tempCell.Points.push_back(m_GridPoints[xi][yi + 2][zi]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 2][zi]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi]);
					tempCell.Points.push_back(m_GridPoints[xi][yi + 1][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi][yi + 2][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 2][zi + 1]);
					tempCell.Points.push_back(m_GridPoints[xi + 1][yi + 1][zi + 1]);

					tempCell.edges.resize(12);
					for (int i = 0;i < 12;i++)
					{
						tempCell.edges[i].index = i;
						tempCell.edges[i].hasIntersection = false;
					}

					tempCell.faces.resize(6);
					for (int i = 0;i < 6;i++)
					{
						tempCell.faces[i].index = i;
						tempCell.faces[i].hasIntersection = false;
						tempCell.faces[i].interNum = 0;
					}



					ExtractT(tempCell.Points[FaceConnection[faceId][0]].hessian, T0);
					ExtractT(tempCell.Points[FaceConnection[faceId][1]].hessian, T1);
					ExtractT(tempCell.Points[FaceConnection[faceId][2]].hessian, T2);
					ExtractT(tempCell.Points[FaceConnection[faceId][3]].hessian, T3);


					ELL_3M_LERP(tmTxm, coords[2], T0, T3);
					ELL_3V_LERP(tmgxm, coords[2], tempCell.Points[FaceConnection[faceId][0]].gradient,
						tempCell.Points[FaceConnection[faceId][3]].gradient);
					ELL_3M_LERP(tmTxp, coords[2], T1, T2);
					ELL_3V_LERP(tmgxp, coords[2], tempCell.Points[FaceConnection[faceId][1]].gradient,
						tempCell.Points[FaceConnection[faceId][2]].gradient);

					ELL_3M_LERP(Tym, coords[0], tmTxm, tmTxp);
					ELL_3V_LERP(gym, coords[0], tmgxm, tmgxp);
				}			
			}
			ComputeGradientLin(cell.faces[faceId].degenerationNormal, T, g,Txm, gxm, Txp, gxp,Tym, gym, Typ, gyp,Tzm, gzm, Tzp, gzp);
			break;

			}

	}
}
