#include "DataPretreating.h"
#include <teem/air.h>
#include <teem/seek.h>

double DataPretreating::maxX = 1;
double DataPretreating::minX = -1;
double DataPretreating::maxY = 1;
double DataPretreating::minY = -1;
double DataPretreating::maxZ = 1;
double DataPretreating::minZ = -1;

char* DataPretreating::filename = NULL;

DataPretreating::DataPretreating(char* _filename){
	filename = _filename;
}



void DataPretreating::ReadPoints_Resample(readMode _rm){
	ifstream infile(filename,ios_base::in);

	if (!infile)
	{
		return ;
	}

	if (_rm == WITH_VOL_COORD)
	{
		m_Points.resize(PointNum);

		for (int i = 0; i < PointNum; i++)
		{

			double t_x, t_y, t_z, t_vol;		
			infile >> t_vol;
			infile >> t_x >> t_y >> t_z;

			if (i == 0)
			{
				maxX = t_x;
				minX = t_x;
				maxY = t_y;
				minY = t_y;
				maxZ = t_z;
				minZ = t_z;
			}else{
				if (maxX < t_x)
				{
					maxX = t_x;
				}

				if (minX > t_x)
				{
					minX = t_x;
				}

				if (maxY < t_y)
				{
					maxY = t_y;
				}

				if (minY > t_y)
				{
					minY = t_y;
				}

				if (maxZ < t_z)
				{
					maxZ = t_z;
				}

				if (minZ > t_z)
				{
					minZ = t_z;
				}
			}

			NorPoint t_point;
			t_point.vol = t_vol;
			t_point.x = t_x;
			t_point.y = t_y;
			t_point.z = t_z;

			m_Points[i] = t_point;
		}
	}else{

		m_Points.resize(PointNum);

		for (int i = 0; i < PointNum; i++)
		{
			int t_x, t_y, t_z;	
			t_x = i / (HNUM * DNUM);
			t_y = (i - t_x * HNUM * DNUM) / DNUM;
			t_z = i - t_x * HNUM * DNUM - t_y * DNUM;

			NorPoint t_point;

			infile >> t_point.vol;
			t_point.x = AIR_AFFINE(0, t_x, WNUM - 1, minX, maxX);
			t_point.y = AIR_AFFINE(0, t_y, HNUM - 1, minY, maxY);
			t_point.z = AIR_AFFINE(0, t_z, DNUM - 1, minZ, maxZ);

			m_Points[i] = t_point;
		}
	}

}

void DataPretreating::ReadPoints_WithoutResample(readMode _rm){
	ifstream infile(filename,ios_base::in);

	if (!infile)
	{
		return ;
	}

	m_PointsGrids.resize(WSampleNum);
	m_PointsTemp.resize(WSampleNum);

	for (int i = 0; i < m_PointsGrids.size(); i++)
	{
		m_PointsGrids[i].resize(HSampleNum);
		m_PointsTemp[i].resize(WSampleNum);

	}

	for (int i = 0; i < m_PointsGrids.size(); i++)
	{
		for (int j = 0; j < m_PointsGrids[0].size(); j++)
		{
			m_PointsGrids[i][j].resize(DSampleNum);
			m_PointsTemp[i][j].resize(DSampleNum);

		}
	}

	if (_rm == WITH_VOL_COORD)
	{
		for (int i = 0; i < WSampleNum; i++)
		{
			for (int j = 0; j < HSampleNum; j++)
			{
				for (int k = 0; k < DSampleNum; k++)
				{
					GridPoint t_point;
					infile >> t_point.Volume >> t_point.x >> t_point.y >> t_point.z;
					m_PointsTemp[i][j][k] = t_point;

					if (i == 0 && j == 0 && k == 0)
					{
						maxX = t_point.x;
						minX = t_point.x;
						maxY = t_point.y;
						minY = t_point.y;
						maxZ = t_point.z;
						minZ = t_point.z;
					}else{
						if (maxX < t_point.x)
						{
							maxX = t_point.x;
						}

						if (minX > t_point.x)
						{
							minX = t_point.x;
						}

						if (maxY < t_point.y)
						{
							maxY = t_point.y;
						}

						if (minY > t_point.y)
						{
							minY = t_point.y;
						}

						if (maxZ < t_point.z)
						{
							maxZ = t_point.z;
						}

						if (minZ > t_point.z)
						{
							minZ = t_point.z;
						}
					}
	
					
				}
			}
		}
	}else if(_rm == WITH_ONLY_VOL){

		for (int i = 0; i < WSampleNum; i++)
		{
			for (int j = 0; j < HSampleNum; j++)
			{
				for (int k = 0; k < DSampleNum; k++)
				{
					GridPoint t_point;
					infile >> t_point.Volume;
						
					t_point.x = AIR_AFFINE(0, i, WSampleNum - 1, minX, maxX);
					t_point.y = AIR_AFFINE(0, j, HSampleNum - 1, minY, maxY);
					t_point.z = AIR_AFFINE(0, k, DSampleNum - 1, minZ, maxZ);
					m_PointsTemp[i][j][k] = t_point;
				}
			}
		}

	}else{

		ofstream outfile("genvol.txt");

		for (int i = 0; i < WSampleNum; i++)
		{
			for (int j = 0; j < HSampleNum; j++)
			{
				for (int k = 0; k < DSampleNum; k++)
				{
					GridPoint t_point;

					if (k == -j + DSampleNum)
					{
						t_point.Volume = 1.0;
					}else{
						t_point.Volume = 0.0;
					}

					t_point.x = AIR_AFFINE(0, i, WSampleNum - 1, minX, maxX);
					t_point.y = AIR_AFFINE(0, j, HSampleNum - 1, minY, maxY);
					t_point.z = AIR_AFFINE(0, k, DSampleNum - 1, minZ, maxZ);
					outfile << t_point.Volume << " " << t_point.x 
						<< " " << t_point.y << " " << t_point.z << endl;
					m_PointsTemp[i][j][k] = t_point;
				}
			}
		}

		outfile.close();
	}
}
// void DataPretreating::GetBoundary()
// {
// 	minX = m_PointsGrids[0][0][0].x, maxX = m_PointsGrids[0][0][0].x;
// 	minY = m_PointsGrids[0][0][0].y, maxY = m_PointsGrids[0][0][0].y;
// 	minZ = m_PointsGrids[0][0][0].z, maxZ = m_PointsGrids[0][0][0].z;
// 	for (int i = 0;i < WNUM;++i)
// 	{
// 		for (int j = 0; j < HNUM; j++)
// 		{
// 			for (int k = 0; k < DNUM; k++)
// 			{
// 				if (minX > m_PointsGrids[i][j][k].x)
// 				{
// 					minX = m_PointsGrids[i][j][k].x;
// 				}
// 				if (minY > m_PointsGrids[i][j][k].y)
// 				{
// 					minY = m_PointsGrids[i][j][k].y;
// 				}
// 				if (minZ > m_PointsGrids[i][j][k].z)
// 				{
// 					minZ = m_PointsGrids[i][j][k].z;
// 				}
// 				if (maxX < m_PointsGrids[i][j][k].x)
// 				{
// 					maxX = m_PointsGrids[i][j][k].x;
// 				}
// 				if (maxY < m_PointsGrids[i][j][k].y)
// 				{
// 					maxY = m_PointsGrids[i][j][k].y;
// 				}
// 				if (maxZ < m_PointsGrids[i][j][k].z)
// 				{
// 					maxZ = m_PointsGrids[i][j][k].z;
// 				}
// 			}
// 		}
// 
// 	}
// 
// }
// 
// 
// void DataPretreating::TransformInputPoints(){
// 	if (m_PointsGrids.size() <= 0)
// 	{
// 		return ;
// 	}
// 
// 	//初始化设置输入点不需要进行坐标转换
// 	bool pointTransformOrNot = false;
// 	for (int i = 0;i < WNUM;++i)
// 	{
// 		for (int j = 0; j < HNUM; j++)
// 		{
// 			for (int k = 0; k < DNUM; k++)
// 			{
// 				if ((abs(m_PointsGrids[i][j][k].x)>1) || (abs(m_PointsGrids[i][j][k].y)>1) || (abs(m_PointsGrids[i][j][k].z)>1))
// 				{
// 					pointTransformOrNot = true;
// 					break;
// 				}
// 			}
// 		}
// 		
// 	}
// 
// 	if (pointTransformOrNot)
// 	{
// 		//转换
// 		for (int i = 0;i < WNUM;++i)
// 		{
// 			for (int j = 0; j < HNUM; j++)
// 			{
// 				for (int k = 0; k < DNUM; k++)
// 				{
// 					m_PointsGrids[i][j][k].x = (m_PointsGrids[i][j][k].x - (maxX+minX)/2.0) / ((maxX-minX)/2.0);
// 					m_PointsGrids[i][j][k].y = (m_PointsGrids[i][j][k].y - (maxY+minY)/2.0) / ((maxY-minY)/2.0);
// 					m_PointsGrids[i][j][k].z = (m_PointsGrids[i][j][k].z - (maxZ+minZ)/2.0) / ((maxZ-minZ)/2.0);
// 				}
// 			}
// 
// 		}
// 		
// 	}
// }


/* 初始化网格点，计算其相应的标量值、梯度值和黑塞矩阵值,整个预处理类的核心函数 */
void DataPretreating::InitialGrids(readMode _rm, gridMode _gm){

	if (_gm == RESAMPLE)
	{
		ReadPoints_Resample(_rm);
		GenerateGrid_Resample();
	}else{
		ReadPoints_WithoutResample(_rm);
		GenerateGrid_WithoutResample();
	}
	
	CalcuFieldGradient();
	CalcuFieldHessian();

}




double DataPretreating::ComputeVol(vector<vector<vector<GridTemp>>>& t_Grid, 
								   int idx, int idy, int idz, double max_dis){

	double result = 0.0;
	
	for (int i = 0; i < t_Grid[idx][idy][idz].m_relation.size(); i++)
	{
		double distance = 0.0;
		double v1[3], v2[3], v3[3];

		v1[0] = m_PointsGrids[idx][idy][idz].x;
		v1[1] = m_PointsGrids[idx][idy][idz].y;
		v1[2] = m_PointsGrids[idx][idy][idz].z;

		v2[0] = m_Points[t_Grid[idx][idy][idz].m_relation[i]].x;
		v2[1] = m_Points[t_Grid[idx][idy][idz].m_relation[i]].y;
		v2[2] = m_Points[t_Grid[idx][idy][idz].m_relation[i]].z;

		ELL_3V_SUB(v3, v2, v1);
		distance = ELL_3V_LEN(v3);

		result += ComputeCubicBval(2 * distance / max_dis) * m_Points[t_Grid[idx][idy][idz].m_relation[i]].vol;

	}

	vector<int>().swap(t_Grid[idx][idy][idz].m_relation);

	return result;
}


void DataPretreating::GenerateGrid_Resample(){

	vector<vector<vector<GridTemp>>> t_Grid;

	t_Grid.resize(WSampleNum);
	m_PointsGrids.resize(WSampleNum);

	for (int i = 0; i < WSampleNum; i++)
	{
		t_Grid[i].resize(HSampleNum);
		m_PointsGrids[i].resize(HSampleNum);
	}

	for (int i = 0; i < WSampleNum; i++)
	{
		for (int j = 0; j < HSampleNum; j++)
		{
			t_Grid[i][j].resize(DSampleNum);
			m_PointsGrids[i][j].resize(DSampleNum);
		}
	}

	double t_gap_x = (maxX - minX) / (WSampleNum - 1);
	double t_gap_y = (maxY - minY) / (HSampleNum - 1);
	double t_gap_z = (maxZ - minZ) / (DSampleNum - 1);

	for (int k = 0; k < m_Points.size(); k++)
	{
		double t_x, t_y, t_z;
		t_x = m_Points[k].x;
		t_y = m_Points[k].y;
		t_z = m_Points[k].z;

		double indx = (t_x - minX) / t_gap_x;
		double indy = (t_y - minY) / t_gap_y;
		double indz = (t_z - minZ) / t_gap_z;

		for (int i = (int)(indx - radius) + 1; i <= (int)(indx + radius); i++) 
		{
			for (int j = (int)(indy - radius) + 1; j <= (int)(indy + radius); j++)
			{
				for(int p = (int)(indz - radius) + 1; p <= (int)(indz + radius); p++)
				{
					
					if (i >= 0 && i < WSampleNum && j >= 0 && j < HSampleNum 
						&& p >= 0 && p < DSampleNum)
					{
						t_Grid[i][j][p].m_relation.push_back(k);
						
					}

				}
			}
		}

	}

	double len[3] = {radius*t_gap_x, radius*t_gap_y, radius*t_gap_z};
	double max_dis = ELL_3V_LEN(len);

	for (int i = 0; i < WSampleNum; i++)
	{
		for (int j = 0; j < HSampleNum; j++)
		{
			for (int k = 0; k < DSampleNum; k++)
			{
				
				m_PointsGrids[i][j][k].x = minX + i * t_gap_x;
				m_PointsGrids[i][j][k].y = minY + j * t_gap_y;
				m_PointsGrids[i][j][k].z = minZ + k * t_gap_z;

				m_PointsGrids[i][j][k].Volume = ComputeVol(t_Grid, i, j, k, max_dis);

			}
		}
	}

	vector<NorPoint>().swap(m_Points);
}

void DataPretreating::GenerateGrid_WithoutResample(){

	double delta_x = (maxX - minX) / (WSampleNum - 1);
	double delta_y = (maxY - minY) / (HSampleNum - 1);
	double delta_z = (maxZ - minZ) / (DSampleNum - 1);
	double max_dis = 0.0, len[3] = {radius * delta_x, radius * delta_y, radius * delta_z};
	max_dis = ELL_3V_LEN(len);

	for (int i = 0; i < WSampleNum; i++)
	{
		for (int j = 0; j < HSampleNum; j++)
		{
			for (int k = 0; k < DSampleNum; k++)
			{

				double result = 0.0;
				double v1[3];

				v1[0] = m_PointsTemp[i][j][k].x;
				v1[1] = m_PointsTemp[i][j][k].y;
				v1[2] = m_PointsTemp[i][j][k].z;

				for (int ii = i - radius; ii <= i + radius; ii++)
				{
					for (int jj = j - radius; jj <= j +radius; jj++ )
					{
						for(int kk = k - radius; kk <= k + radius; kk++)
						{
							if (ii >= 0 && ii < WSampleNum && jj >= 0 && jj < HSampleNum 
								&& kk >= 0 && kk < DSampleNum)
							{
								double v2[3], v3[3], distance;
								v2[0] = m_PointsTemp[ii][jj][kk].x;
								v2[1] = m_PointsTemp[ii][jj][kk].y;
								v2[2] = m_PointsTemp[ii][jj][kk].z;

								ELL_3V_SUB(v3, v2, v1);
								distance = ELL_3V_LEN(v3);

								result += ComputeCubicBval(2 * distance / max_dis) * m_PointsTemp[ii][jj][kk].Volume;
							}

						}
					}
				
				}

				m_PointsGrids[i][j][k].Volume = result;
				m_PointsGrids[i][j][k].x = m_PointsTemp[i][j][k].x;
				m_PointsGrids[i][j][k].y = m_PointsTemp[i][j][k].y;
				m_PointsGrids[i][j][k].z = m_PointsTemp[i][j][k].z;


			}
		}
	}


	vector<vector<vector<GridPoint>>>().swap(m_PointsTemp);

}

/* 计算点在标量场的梯度值 */
void DataPretreating::CalcuFieldGradient( ){

	double gap_x = (double)(maxX - minX)/(WSampleNum - 1);
	double gap_y = (double)(maxY - minY)/(HSampleNum - 1);
	double gap_z = (double)(maxZ - minZ)/(DSampleNum - 1);

	
	for (int i = 0; i < WSampleNum; i++)
	{
		for (int j = 0; j < HSampleNum; j++)
		{
			for (int k = 0; k < DSampleNum; k++)
			{
				//计算x方向的导数
				if (i == 0)
				{
					m_PointsGrids[0][j][k].gradient[0] = (m_PointsGrids[1][j][k].Volume - m_PointsGrids[0][j][k].Volume) / gap_x;

				}else if (i == WSampleNum - 1){

					m_PointsGrids[WSampleNum - 1][j][k].gradient[0] = (m_PointsGrids[WSampleNum - 1][j][k].Volume - m_PointsGrids[WSampleNum - 2][j][k].Volume) / gap_x;

				}else{
					m_PointsGrids[i][j][k].gradient[0] = (m_PointsGrids[i + 1][j][k].Volume - m_PointsGrids[i - 1][j][k].Volume) / (2 * gap_x);
				}


				// 计算y方向导数
				if (j == 0)
				{
					m_PointsGrids[i][0][k].gradient[1] = (m_PointsGrids[i][1][k].Volume - m_PointsGrids[i][0][k].Volume) / gap_y;
				}else if (j == HSampleNum - 1){
					m_PointsGrids[i][HSampleNum - 1][k].gradient[1] = (m_PointsGrids[i][HSampleNum - 1][k].Volume - m_PointsGrids[i][HSampleNum - 2][k].Volume) / gap_y;
				}else{
					m_PointsGrids[i][j][k].gradient[1] = (m_PointsGrids[i][j + 1][k].Volume - m_PointsGrids[i][j - 1][k].Volume) / (2 * gap_y);
				}


				//计算z方向导数
				if(k == 0)
				{
					m_PointsGrids[i][j][0].gradient[2] = (m_PointsGrids[i][j][1].Volume - m_PointsGrids[i][j][0].Volume) / gap_z;
					
				}else if (k == DSampleNum - 1)
				{
					m_PointsGrids[i][j][DSampleNum - 1].gradient[2] = (m_PointsGrids[i][j][DSampleNum - 1].Volume - m_PointsGrids[i][j][DSampleNum - 2].Volume) / gap_z;
				}else{
					m_PointsGrids[i][j][k].gradient[2] = (m_PointsGrids[i][j][k + 1].Volume - m_PointsGrids[i][j][k - 1].Volume) / (2 * gap_z);
				}

			}
		}
	}
}


/* 计算点在标量场中的黑塞矩阵值 */
void DataPretreating::CalcuFieldHessian(){

	double gap_x = (double)(maxX - minX)/(WSampleNum - 1);
	double gap_y = (double)(maxY - minY)/(HSampleNum - 1);
	double gap_z = (double)(maxZ - minZ)/(DSampleNum - 1);

	for (int i = 0; i < WSampleNum; i++)
	{
		for (int j = 0; j< HSampleNum; j++)
		{
			for (int k = 0; k < DSampleNum; k++)
			{
				//计算二阶偏导数fxx, fyx, fzx
				if (i == 0)
				{
					m_PointsGrids[0][j][k].hessian.fxx = (m_PointsGrids[1][j][k].gradient[0] -m_PointsGrids[0][j][k].gradient[0]) / gap_x;
					m_PointsGrids[0][j][k].hessian.fyx = (m_PointsGrids[1][j][k].gradient[1] -m_PointsGrids[0][j][k].gradient[1]) / gap_x;
					m_PointsGrids[0][j][k].hessian.fzx = (m_PointsGrids[1][j][k].gradient[2] -m_PointsGrids[0][j][k].gradient[2]) / gap_x;
		
				}else if (i == WSampleNum - 1){
					m_PointsGrids[WSampleNum - 1][j][k].hessian.fxx = (m_PointsGrids[WSampleNum - 1][j][k].gradient[0]
						- m_PointsGrids[WSampleNum - 2][j][k].gradient[0]) / gap_x;
					m_PointsGrids[WSampleNum - 1][j][k].hessian.fyx = (m_PointsGrids[WSampleNum - 1][j][k].gradient[1]
						- m_PointsGrids[WSampleNum - 2][j][k].gradient[1]) / gap_x;
					m_PointsGrids[WSampleNum - 1][j][k].hessian.fzx = (m_PointsGrids[WSampleNum - 1][j][k].gradient[2]
						- m_PointsGrids[WSampleNum - 2][j][k].gradient[2]) / gap_x;
				}else{
					m_PointsGrids[i][j][k].hessian.fxx = (m_PointsGrids[i + 1][j][k].gradient[0] - m_PointsGrids[i - 1][j][k].gradient[0]) / (2 * gap_x);
					m_PointsGrids[i][j][k].hessian.fyx = (m_PointsGrids[i + 1][j][k].gradient[1] - m_PointsGrids[i - 1][j][k].gradient[1]) / (2 * gap_x);
					m_PointsGrids[i][j][k].hessian.fzx = (m_PointsGrids[i + 1][j][k].gradient[2] - m_PointsGrids[i - 1][j][k].gradient[2]) / (2 * gap_x);

				}

				//计算二阶偏导数fxy, fyy, fzy
				if (j == 0)
				{
					m_PointsGrids[i][0][k].hessian.fxy = (m_PointsGrids[i][1][k].gradient[0] - m_PointsGrids[i][0][k].gradient[0]) / gap_y;
					m_PointsGrids[i][0][k].hessian.fyy = (m_PointsGrids[i][1][k].gradient[1] - m_PointsGrids[i][0][k].gradient[1]) / gap_y;
					m_PointsGrids[i][0][k].hessian.fzy = (m_PointsGrids[i][1][k].gradient[2] - m_PointsGrids[i][0][k].gradient[2]) / gap_y;

				}else if (j == HSampleNum - 1){
					m_PointsGrids[i][HSampleNum - 1][k].hessian.fxy = (m_PointsGrids[i][HSampleNum - 1][k].gradient[0] 
						- m_PointsGrids[i][HSampleNum - 2][k].gradient[0]) / gap_y;
					m_PointsGrids[i][HSampleNum - 1][k].hessian.fyy = (m_PointsGrids[i][HSampleNum - 1][k].gradient[1] 
						- m_PointsGrids[i][HSampleNum - 2][k].gradient[1]) / gap_y;
					m_PointsGrids[i][HSampleNum - 1][k].hessian.fzy = (m_PointsGrids[i][HSampleNum - 1][k].gradient[2] 
						- m_PointsGrids[i][HSampleNum - 2][k].gradient[2]) / gap_y;
				}else{
					m_PointsGrids[i][j][k].hessian.fxy = (m_PointsGrids[i][j + 1][k].gradient[0] - m_PointsGrids[i][j - 1][k].gradient[0]) / (2 * gap_y);
					m_PointsGrids[i][j][k].hessian.fyy = (m_PointsGrids[i][j + 1][k].gradient[1] - m_PointsGrids[i][j - 1][k].gradient[1]) / (2 * gap_y);
					m_PointsGrids[i][j][k].hessian.fzy = (m_PointsGrids[i][j + 1][k].gradient[2] - m_PointsGrids[i][j - 1][k].gradient[2]) / (2 * gap_y);
				}

				//计算二阶偏导数fxz, fyz, fzz
				if (k == 0)
				{
					m_PointsGrids[i][j][0].hessian.fxz = (m_PointsGrids[i][j][1].gradient[0] - m_PointsGrids[i][j][0].gradient[0]) / gap_z;
					m_PointsGrids[i][j][0].hessian.fyz = (m_PointsGrids[i][j][1].gradient[1] - m_PointsGrids[i][j][0].gradient[1]) / gap_z;
					m_PointsGrids[i][j][0].hessian.fzz = (m_PointsGrids[i][j][1].gradient[2] - m_PointsGrids[i][j][0].gradient[2]) / gap_z;

				}else if (k == DSampleNum - 1)
				{
					m_PointsGrids[i][j][DSampleNum - 1].hessian.fxz = (m_PointsGrids[i][j][DSampleNum - 1].gradient[0] 
						- m_PointsGrids[i][j][DSampleNum - 2].gradient[0]) / gap_z;
					m_PointsGrids[i][j][DSampleNum - 1].hessian.fyz = (m_PointsGrids[i][j][DSampleNum - 1].gradient[1]
						- m_PointsGrids[i][j][DSampleNum - 2].gradient[1]) / gap_z;
					m_PointsGrids[i][j][DSampleNum - 1].hessian.fzz = (m_PointsGrids[i][j][DSampleNum - 1].gradient[2]
						- m_PointsGrids[i][j][DSampleNum - 2].gradient[2]) / gap_z;
				}else{
					m_PointsGrids[i][j][k].hessian.fxz = (m_PointsGrids[i][j][k + 1].gradient[0] + m_PointsGrids[i][j][k - 1].gradient[0]) / (2 * gap_z);
					m_PointsGrids[i][j][k].hessian.fyz = (m_PointsGrids[i][j][k + 1].gradient[1] + m_PointsGrids[i][j][k - 1].gradient[1]) / (2 * gap_z);
					m_PointsGrids[i][j][k].hessian.fzz = (m_PointsGrids[i][j][k + 1].gradient[2] + m_PointsGrids[i][j][k - 1].gradient[2]) / (2 * gap_z);

				}
			}
		}
	}
}

/* 计算三次B样条函数值 */
double DataPretreating::ComputeCubicBval(double x){
	x = abs(x);
	double result = 0;

	if (x < 1)
	{
		result = (x - 2) * x * x / 2 + 2.0 / 3.0;
	}else if (x < 2)
	{
		result = pow((2 - x),3) / 6;
	}else{
		result = 0;
	}

	return result;
}

/* 计算二次B样条函数值 */
double DataPretreating::ComputeQuadraticBval(double x){
	x = abs(x);
	double result  = 0;
	if (x < 1 / 2)
	{
		result = 3 / 4 - x * x;
	}else if (x < -3 / 2)
	{
		result = x * x / 2 + 3 * x / 2 + 9.0/8.0;
	}else{
			result = 0;
		}

	return result;
}

/* 计算一次B样条函数值 */
double DataPretreating::ComputeLinearBval(double x){
	x = abs(x);
	double result = 0;
	if (x < 1)
	{
		result = 1 - x;
	}else{
		result = 0;
	}
	return result;
}