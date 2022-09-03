#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <time.h>
#include <vector>
#include "counttime.cpp"
using namespace std;
const int timelimit = 5; //最优解
const int level = 2;
const int MaxN = 3000;
const double alpha = 0.8;
int shuffle[level] = {0, 0};
// int size[level] = {377,94};//每一层节点个数
// int movable[level] = {263,65};//每一层不可动节点个数（相应节点编号为0~n-1）
int size[level] = {0, 0};			   //每一层节点个数
int movable[level] = {0, 0};		   //每一层不可动节点个数（相应节点编号为0~n-1）
int pos[level][MaxN] = {0};			   //每个位置上的节点编号
int loc[level][MaxN] = {0};			   //每个编号的节点位置
int evmatrix[level][MaxN][MaxN] = {0}; //评估矩阵，【i】【j】表示i在j前的交叉数
int delta[level][MaxN][MaxN] = {0};	   //移动的评估矩阵，【i】【j】表示j放到i前（后）的交叉变化数，优化则为正

int cpoint[level][MaxN][MaxN] = {0}; //向右的连接节点编号
int nc[level][MaxN] = {0};			 //向右的连接节点个数

int bestnow = 0x7fffffff;
int besthistory = bestnow;
int bestpos[level][MaxN] = {0};
int bestloc[level][MaxN] = {0};
int maxevmatrix[level][MaxN][MaxN] = {0}; //最优解记录

int Data[MaxN][MaxN] = {-1};
int Right[MaxN][MaxN] = {-1};
int num_connect[MaxN] = {0};
int right_connect[MaxN] = {0}; //以上四个为二层时方便处理的单元集合
int max(int a, int b, int c)
{
	if (a >= b && a >= c)
		return a;
	else if (b >= a && b >= c)
		return b;
	else
		return c;
}
void orig()
{
	//目的是更新两个pos、loc
	//再赋值evmatrix
	//	for(int i = movable[0];i<size[0];i++) //取出左点
	//	{
	//		pos[0][i] = -1;
	//	}
	int minnum[size[0]] = {0};
	int minplace[size[0]] = {0};
	for (int key = movable[0]; key < size[0]; key++)
	{
		int count = 0;
		for (int i = 0; i < movable[0]; i++) //初始化最初位置交叉数
		{
			count += evmatrix[0][key][i];
		}
		int min = count;
		int min_place = 0;
		//		printf("\nkey = %d\n",key);
		for (int i = 0; i < movable[0]; i++) //计算每个位置的交叉数
		{
			count -= evmatrix[0][key][i];
			count += evmatrix[0][i][key];
			if (count < min)
			{
				//				printf("min = %d, place = %d\n",count,i+1);
				min = count;
				min_place = i + 1;
			}
		}
		//		printf("\n");
		minnum[key] = min;
		minplace[key] = min_place;
	}
	for (int key = movable[0]; key < size[0]; key++)
	{
		printf("%d %d\n", minnum[key], minplace[key]);
	}

	getchar();
}

int active(int id, int t)
{
	if (pos[t][id] >= movable[t])
		return 1;
	return 0;
}
int findsuit(int t, int j, int i)
{
	if (i == j)
		return -1;

	if (active(i, t) && active(j, t))
	{
		int tmp = max(delta[t][j][i], delta[t][i][j], delta[t][j][i] + delta[t][i - 1][j]);
		if (tmp < 0)
			return -1;
		if (tmp == delta[t][j][i])
			return 0;
		if (tmp == delta[t][i][j])
			return 1;
		if (tmp == delta[t][j][i] + delta[t][i - 1][j])
			return 2;
	}
	if (active(i, t))
	{
		if (delta[t][j][i] >= 0)
		{
			return 0;
		}
	}
	if (active(j, t))
	{
		if (delta[t][i][j] >= 0)
			return 1;
	}
	return -1;
}
void renewEvArray(int t, int j, int i, int dir)
{
	if (t == 0)
	{
		// insert i before j;
		if (dir == 0)
		{
			//			printf("\n\n\n以下开始更新矩阵 j = %d,i = %d",j,i);
			//更新矩阵
			int counter[size[1 - t]] = {0};
			for (int k = j; k < i; k++) // k是位置
			{
				for (int q = 0; q < num_connect[pos[t][k]]; q++) //对每一个夹点的邻接点
				{
					counter[Data[pos[t][k]][q + 2]]++; //被夹点连接的点记录数自增
				}
			}
			for (int k = 0; k < num_connect[pos[t][i]]; k++) //对移动点的每一个节点
			{
				for (int q = 0; q < size[1 - t]; q++) //对右侧每一个被夹点连接的点
				{
					if (counter[q] != 0)
					{
						evmatrix[1 - t][Data[pos[t][i]][k + 2]][q] -= counter[q];
						evmatrix[1 - t][q][Data[pos[t][i]][k + 2]] += counter[q];
					}
				}
			}
		}
		else if (dir == 1)
		{
			//更新矩阵
			int counter[size[1 - t]] = {0};
			for (int k = j + 1; k <= i; k++) // k是位置
			{
				for (int q = 0; q < num_connect[pos[t][k]]; q++) //对每一个夹点的邻接点
				{
					counter[Data[pos[t][k]][q + 2]]++; //被夹点连接的点记录数自增
				}
			}
			for (int k = 0; k < num_connect[pos[t][j]]; k++) //对移动点的每一个节点
			{
				for (int q = 0; q < size[1 - t]; q++) //对右侧每一个被夹点连接的点
				{
					if (counter[q] != 0)
					{
						evmatrix[1 - t][Data[pos[t][j]][k + 2]][q] += counter[q];
						evmatrix[1 - t][q][Data[pos[t][j]][k + 2]] -= counter[q];
					}
				}
			}
		}
	}
	else if (t == 1)
	{
		//		printf("**********************\nj = %d ,i = %d ********************\n");
		if (dir == 0)
		{
			int counter[size[1 - t]] = {0};
			for (int k = j; k < i; k++) // k是位置
			{
				for (int q = 0; q < right_connect[pos[t][k]]; q++) //对每一个夹点的邻接点
				{

					counter[Right[pos[t][k]][q + 2]]++; //被夹点连接的点记录数自增
				}
			}
			for (int k = 0; k < right_connect[pos[t][i]]; k++) //对移动点的每一个节点
			{
				for (int q = 0; q < size[1 - t]; q++) //对右侧每一个被夹点连接的点
				{
					if (counter[q] != 0)
					{
						evmatrix[1 - t][Right[pos[t][i]][k + 2]][q] -= counter[q];
						evmatrix[1 - t][q][Right[pos[t][i]][k + 2]] += counter[q];
					}
				}
			}
		}
		else if (dir == 1)
		{
			//更新矩阵
			int counter[size[1 - t]] = {0};
			for (int k = j + 1; k <= i; k++) // k是位置
			{
				for (int q = 0; q < right_connect[pos[t][k]]; q++) //对每一个夹点的邻接点
				{
					counter[Right[pos[t][k]][q + 2]]++; //被夹点连接的点记录数自增
				}
			}
			for (int k = 0; k < right_connect[pos[t][j]]; k++) //对移动点的每一个节点
			{
				for (int q = 0; q < size[1 - t]; q++) //对右侧每一个被夹点连接的点
				{
					if (counter[q] != 0)
					{
						evmatrix[1 - t][Right[pos[t][j]][k + 2]][q] += counter[q];
						evmatrix[1 - t][q][Right[pos[t][j]][k + 2]] -= counter[q];
					}
				}
			}
		}
	}
}
void initialize(string filepath)
{
	srand((int)time(NULL));
	ifstream fq;
	fq.open(filepath, ios::in);
	if (!fq.is_open())
		cout << "open file failure" << endl; // f:\\incgraph_2_0.06_5_30_1.20_1.txt
	int g1 = 0;
	int a;

	// line start from 0
	int tmp1 = 0;
	//读取
	int left_num = 10;
	int right_num = 10;

	while (!fq.eof())
	{
		char buffer[256] = {'\0'};
		fq.getline(buffer, 500); //读入每行

		int tmp2 = 0;
		for (int i = 0; i < 255; i++)
		{
			if (buffer[i] >= '0' && buffer[i] <= '9')
			{
				if (i == 0 || (buffer[i - 1] < '0' || buffer[i - 1] > '9'))
				{
					if (buffer[i + 1] >= '0' && buffer[i + 1] <= '9' && buffer[i + 2] >= '0' && buffer[i + 2] <= '9' &&
						buffer[i + 3] >= '0' && buffer[i + 3] <= '9')
						a = (int)(buffer[i] - '0') * 1000 + (int)(buffer[i + 1] - '0') * 100 + (int)(buffer[i + 2] - '0') * 10 + (int)(buffer[i + 3] - '0');
					else if (buffer[i + 1] >= '0' && buffer[i + 1] <= '9' && buffer[i + 2] >= '0' && buffer[i + 2] <= '9')
						a = (int)(buffer[i] - '0') * 100 + (int)(buffer[i + 1] - '0') * 10 + (int)(buffer[i + 2] - '0');
					else if (buffer[i + 1] >= '0' && buffer[i + 1] <= '9')
						a = (int)(buffer[i] - '0') * 10 + (int)(buffer[i + 1] - '0');
					else
						a = (int)(buffer[i] - '0');
					// read line 2
					if (tmp1 == 1)
					{
						size[tmp2] = a;
						tmp2++;
					}
					// when read line 3 set left and right num
					if (tmp1 == 2)
					{
						left_num = size[0];
						right_num = size[1];
					}
					// when you read the 3rd to the left_num+2(included) lines.
					if (tmp1 > 1 && tmp1 < left_num + 2)
					{
						Data[tmp1 - 2][tmp2] = a; // directly put all data into the ``Data[]`` array
						tmp2++;
					}
					// when you read the left_num+3(included) line and so on.
					if (tmp1 >= left_num + 2)
					{
						Right[tmp1 - left_num - 2][tmp2] = a; // directly put all data into the ``Right[]`` array
						tmp2++;
					}
				}
			}
		}
		// after read each line, set all useless element in Data to -1
		if (tmp2 < MaxN && tmp1 < (left_num + 2) && tmp1 > 1)
		{
			for (int i = tmp2; i < MaxN; i++)
				Data[tmp1 - 2][i] = -1;
		}
		tmp1++;
	}

	fq.close();
	// cout<<"leftnum: "<<left_num<<endl;
	// cout<<"rightnum: "<<right_num<<endl;
	// for each i, Data[i][0]=0 means that vertex can move, =1 means that vertex can not move

	int ct = 0; // ct is the number of 1s in Data[:][0]
	for (int i = 0; i < left_num; i++)
	{
		if (Data[i][0] == 1)
			ct++;
	}
	// number of 1s in Data[:][0]
	movable[0] = ct;
	ct = 0;
	for (int i = 0; i < right_num; i++)
	{
		if (Right[i][0] == 1)
			ct++;
		for (int j = 2; j < right_num; j++)
		{
			Right[i][j] = -1;
		}
	}

	movable[1] = ct;
	// the number of vertices in the left layer * some ratio
	shuffle[0] = int((left_num - movable[0]) * alpha);
	shuffle[1] = int((right_num - movable[1]) * alpha);

	// int Data1[left_num][left_num] = {-1};
	vector<vector<int>> Data1(left_num,vector<int>(left_num,-1));
	for (int i = 0; i < left_num; i++)
	{
		for (int j = 0; j < left_num; j++)
			Data1[Data[i][1]][j] = Data[i][j];
	}
	for (int i = 0; i < left_num; i++)
	{
		for (int j = 0; j < left_num; j++)
			Data[i][j] = Data1[i][j];
	}
	for (int i = 0; i < left_num; i++)
	{
		int j = MaxN - 1;
		for (; j >= 0; j--)
		{
			if (Data[i][j] != -1)
				break;
		}
		nc[0][i] = j - 1;
		num_connect[i] = j - 1;
	}
	//左侧完成
	//左侧完成更新
	//以上实现读取数据
	//以下进行数据规定化,包括归零化，计算另一边的邻接点
	int set[left_num + right_num] = {0};

	for (int i = 0; i < right_num; i++)
	{
		if (Right[i][0] == 0 && Right[i][1] > left_num + movable[1] - 1)
			Right[i][1] -= left_num;
		//		set[Right[i][1]] = i+left_num;
	}

	for (int i = 0; i < left_num; i++)
	{
		for (int j = 2; j < 2 + num_connect[i]; j++)
		{
			Data[i][j] -= left_num;
			Right[Data[i][j]][2 + right_connect[Data[i][j]]] = i;
			right_connect[Data[i][j]]++;
		}
	}

	int s = 0;

	s = 0;
	vector<vector<int>> Right1(right_num,vector<int>(right_num,0));
	// int Right1[right_num][right_num] = {0};
	for (int i = 0; i < right_num; i++)
	{
		for (int j = 0; j < right_connect[i] + 2; j++)
		{

			Right1[Right[i][1]][j] = Right[i][j];
		}
	}
	for (int i = 0; i < left_num; i++)
	{
		for (int j = 2; j < num_connect[i] + 2; j++)
		{
			Data[i][j] = Right[Data[i][j]][1];
		}
	}
	int rightem[right_num] = {0};
	for (int i = 0; i < right_num; i++)
	{
		rightem[Right[i][1]] = i;
	}
	for (int i = 0; i < right_num; i++)
	{
		int j = 0;

		for (; j < right_connect[rightem[i]] + 2; j++)
		{

			Right[i][j] = Right1[i][j];
		}
		for (; j < right_num; j++)
			Right[i][j] = -1;
	}
	for (int i = 0; i < right_num; i++)
	{
		int j = 0;
		while (Right[i][j + 2] != -1)
			j++;
		right_connect[i] = j;
	}

	for (int i = 0; i < left_num - movable[0]; i++)
	{
		int a = rand() % movable[0];
		int tmp = pos[0][i + movable[0]];
		for (int k = a; k < i + movable[0]; k++)
		{
			pos[0][k + 1] = pos[0][k];
		}
		pos[0][a] = tmp;
	}
	for (int i = 0; i < right_num - movable[1]; i++)
	{
		int a = rand() % movable[1];
		int tmp = pos[1][i + movable[1]];
		for (int k = a; k < i + movable[1]; k++)
		{
			pos[1][k + 1] = pos[1][k];
		}
		pos[1][a] = tmp;
	}
	for (int i = 0; i < left_num; i++)
	{
		pos[0][i] = i;
		loc[0][i] = i;
	}

	for (int i = 0; i < right_num; i++)
	{
		pos[1][i] = i;
		loc[1][i] = i;
	}
}
void datainit()
{
	int left_num = size[0];
	int right_num = size[1];
	//以下计算位置评估矩阵，表示编号i在编号j节点前两点边交叉个数
	for (int i = 0; i < right_num; i++)
	{
		for (int j = 0; j < right_num; j++)
		{
			if (i == j)
				evmatrix[1][i][j] = 0;
			else
			{
				for (int m = 0; m < right_connect[i]; m++)
				{
					for (int n = 0; n < right_connect[j]; n++)
					{
						if ((loc[0][Right[i][m + 2]] - loc[0][Right[j][n + 2]]) > 0) //若满足位置关系
						{
							evmatrix[1][i][j]++; //
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < left_num; i++)
	{
		for (int j = 0; j < left_num; j++)
		{
			if (i == j)
				evmatrix[0][i][j] = 0;
			else
			{
				for (int m = 0; m < num_connect[i]; m++)
				{
					for (int n = 0; n < num_connect[j]; n++)
					{
						if ((loc[1][Data[i][m + 2]] - loc[1][Data[j][n + 2]]) > 0) //若满足位置关系
						{
							evmatrix[0][i][j]++; //
						}
					}
				}
			}
		}
	}
}
void dp(int t) // t即type b表示层数
{
	int left_num = size[0];
	int right_num = size[1];
	//计算交换delta矩阵

	for (int j = 0; j < size[t]; j++)
	{
		delta[t][0][j] = 0;
		for (int i = 0; i <= j; i++)
		{
			delta[t][0][j] += (evmatrix[t][pos[t][i]][pos[t][j]] - evmatrix[t][pos[t][j]][pos[t][i]]);
		}
	}
	//		for(int j = 0;j < left_num;j++)
	//		printf("%d  ",delta_matrix[0][j]);
	for (int j = 1; j < size[t]; j++)
	{
		for (int i = 1; i <= j; i++)
		{
			delta[t][i][j] = delta[t][i - 1][j] + evmatrix[t][pos[t][j]][pos[t][i - 1]] - evmatrix[t][pos[t][i - 1]][pos[t][j]];
		}
	}
	for (int j = size[t] - 1; j >= 0; j--)
	{
		delta[t][size[t] - 1][j] = 0;
		for (int i = j; i <= size[t] - 1; i++)
		{
			delta[t][size[t] - 1][j] += (evmatrix[t][pos[t][j]][pos[t][i]] - evmatrix[t][pos[t][i]][pos[t][j]]);
		}
	}
	for (int j = 0; j < size[t]; j++)
	{
		for (int i = size[t] - 1 - 1; i >= j; i--)
		{
			delta[t][i][j] = delta[t][i + 1][j] + evmatrix[t][pos[t][i + 1]][pos[t][j]] - evmatrix[t][pos[t][j]][pos[t][i + 1]];
		}
	}
	// ofstream fout("Ddpls.out");
    // for(int i=0;i<left_num;i++)
    //     for(int j=0;j<left_num;j++)
    //         fout<<-delta[0][i][j]<<endl;
    // cout<<"finish\n";int aaa;cin>>aaa;
	//下面进行动态规划，获得最优运动
	int dir[size[t]] = {0};
	int behind[size[t]] = {0};
	int cross_delta[size[t]] = {0};
	int count = 0;

	for (int i = 1; i < size[t]; i++) //内部动态规划 !!!!0-30，from 1 effective
	{

		int rem = -1;
		int in_tmp_delta = -10000;
		int in_tmpj = 0;
		int max = 0;
		int tmp_dir = 0;
		int in_dir = 0;
		for (int j = 0; j <= i; j++)
		{
			tmp_dir = findsuit(t, j, i);
			switch (tmp_dir)
			{
			case -1:
				max = 0;
				break;
			case 0:
				max = delta[t][j][i];
				break;
			case 1:
				max = delta[t][i][j];
				break;
			case 2:
				max = delta[t][j][i] + delta[t][i - 1][j];
				break;
			}
			if (j == 0)
				count = max;
			else
				count = cross_delta[j - 1] + max;

			// transact(j+1,i)
			if (count >= in_tmp_delta)
			{
				in_tmp_delta = count; //记录一重值
				in_tmpj = j;
				in_dir = tmp_dir;
			}
		}
		cross_delta[i] = in_tmp_delta;
		behind[i] = in_tmpj;
		dir[i] = in_dir;
	}

	int i = size[t] - 1;
	while (i > 0) //实施序列变换，并更新矩阵
	{
		int j = behind[i];
		int tmp = pos[t][i];
		// insert i before j
		if (dir[i] == 0)
		{
			tmp = pos[t][i];
			renewEvArray(t, j, i, dir[i]);
			for (int k = i; k > j; k--)
			{
				pos[t][k] = pos[t][k - 1];
			}
			pos[t][j] = tmp;
		}

		else if (dir[i] == 1)
		{
			tmp = pos[t][j];
			renewEvArray(t, j, i, dir[i]);
			for (int k = j; k < i; k++)
			{
				pos[t][k] = pos[t][k + 1];
			}
			pos[t][i] = tmp;
		}
		else if (dir[i] == 2)
		{
			//			printf("\n\n***********************\nnotice!!!!\ncross_delta=%d\n****************************\n",cross_delta[size[t]-1]);
			tmp = pos[t][i];
			renewEvArray(t, j, i, 0);
			for (int k = i; k > j; k--)
			{
				pos[t][k] = pos[t][k - 1];
			}
			pos[t][j] = tmp;

			tmp = pos[t][j + 1];
			renewEvArray(t, j + 1, i, 1);
			for (int k = j + 1; k < i; k++)
			{
				pos[t][k] = pos[t][k + 1];
			}
			pos[t][i] = tmp;
		}
		i = j - 1;
	}
	for (int i = 0; i < size[t]; i++) //计算loc
	{
		loc[t][pos[t][i]] = i;
	}
	for (int i = 0; i < size[1 - t]; i++) //计算loc
	{
		loc[1 - t][pos[1 - t][i]] = i;
	}

	int crosssum = 0;
	for (int i = 0; i < size[t]; i++)
	{
		for (int j = i; j < size[t]; j++)
		{
			crosssum += evmatrix[t][pos[t][i]][pos[t][j]];
		}
	}
	//	printf("\n%d\n",crosssum);
	int crosssum1 = 0;
	for (int i = 0; i < size[1 - t]; i++)
	{
		for (int j = i; j < size[1 - t]; j++)
		{
			crosssum1 += evmatrix[1 - t][pos[1 - t][i]][pos[1 - t][j]];
		}
	}
	int s = 0;

	bestnow = crosssum;
	//	getchar();
}
void extract()
{
	int left_num = size[0];
	int right_num = size[1];
	int tmppoint[MaxN] = {0}; //记录位置
	int count = 0;
	for (int i = 0; i < left_num; i++)
	{
		if (pos[0][i] >= movable[0])
		{
			tmppoint[count] = i;
			count++;
		}
	}
	int ta = 0;
	int tb = 0;
	int tmp = 0;
	if (count >= 2)
	{
		for (int mm = 0; mm < shuffle[0]; mm++)
		{
			int t = 0;
			do
			{
				ta = rand() % count;
				tb = rand() % count;
			} while (ta == tb);

			if (ta > tb) //让ta较小
			{
				int q = ta;
				ta = tb;
				tb = q;
			}
			int i = tmppoint[tb];
			int j = tmppoint[ta];
			//进行交换和更新
			int tmp = pos[t][i];
			renewEvArray(t, j, i, 0);
			for (int k = i; k > j; k--)
			{
				pos[t][k] = pos[t][k - 1];
			}
			pos[t][j] = tmp;

			tmp = pos[t][j + 1];
			renewEvArray(t, j + 1, i, 1);
			for (int k = j + 1; k < i; k++)
			{
				pos[t][k] = pos[t][k + 1];
			}
			pos[t][i] = tmp;
			//
		}
		for (int cpp = 0; cpp < 5; cpp++)
		{
			int miku = rand() % left_num;

			do
			{
				miku = rand() % left_num;
			} while (pos[0][miku] < movable[0]);
			int mov = miku;
			int bi = rand() % left_num;
			int hatsune = 0;
			if (bi >= mov)
				hatsune = 1;

			if (hatsune == 0)
			{
				renewEvArray(0, bi, mov, hatsune);
				int tmp = pos[0][mov];
				for (int k = mov; k > bi; k--)
				{
					pos[0][k] = pos[0][k - 1];
				}
				pos[0][bi] = tmp;
			}
		}
		for (int i = 0; i < left_num; i++)
		{
			loc[0][pos[0][i]] = i;
		}
	}

	//右侧
	count = 0;

	for (int i = 0; i < right_num; i++)
	{
		if (pos[1][i] >= movable[1])
		{
			tmppoint[count] = i;
			count++;
		}
	}
	if (count >= 2)
	{
		for (int mm = 0; mm < shuffle[1]; mm++)
		{
			int t = 1;
			do
			{
				ta = rand() % count;
				tb = rand() % count;
			} while (ta == tb);

			if (ta > tb) //让ta较小
			{
				int q = ta;
				ta = tb;
				tb = q;
			}
			int i = tmppoint[tb];
			int j = tmppoint[ta];
			//进行交换和更新
			int tmp = pos[t][i];
			renewEvArray(t, j, i, 0);
			for (int k = i; k > j; k--)
			{
				pos[t][k] = pos[t][k - 1];
			}
			pos[t][j] = tmp;

			tmp = pos[t][j + 1];
			renewEvArray(t, j + 1, i, 1);
			for (int k = j + 1; k < i; k++)
			{
				pos[t][k] = pos[t][k + 1];
			}
			pos[t][i] = tmp;
			//
		}
		for (int cpp = 0; cpp < 5; cpp++)
		{
			int miku = rand() % right_num;

			do
			{
				miku = rand() % right_num;
			} while (pos[1][miku] < movable[1]);
			int mov = miku;
			int bi = rand() % right_num;
			int hatsune = 0;
			if (bi >= mov)
				hatsune = 1;

			if (hatsune == 0)
			{
				renewEvArray(1, bi, mov, hatsune);
				int tmp = pos[1][mov];
				for (int k = mov; k > bi; k--)
				{
					pos[1][k] = pos[1][k - 1];
				}
				pos[1][bi] = tmp;
			}
		}
		for (int i = 0; i < left_num; i++)
		{
			loc[0][pos[0][i]] = i;
		}
	}

	int crosssum = 0;

	crosssum = 0;
}
void storebest()
{
	besthistory = bestnow;
	for (int i = 0; i <= 1; i++)
	{
		for (int j = 0; j < MaxN; j++)
		{
			bestpos[i][j] = pos[i][j];
			for (int k = 0; k < MaxN; k++)
			{
				maxevmatrix[i][j][k] = evmatrix[i][j][k];
			}
		}
	}
}
void renew()
{

	bestnow = 0x7fffffff;
	besthistory = 0x7fffffff;

	for (int i = 0; i < 2; i++)
	{
		shuffle[i] = 0;
		size[i] = 0;
		movable[i] = 0;
		for (int j = 0; j < MaxN; j++)
		{
			num_connect[j] = 0;
			right_connect[j] = 0;
			pos[i][j] = 0;
			loc[i][j] = 0;
			nc[i][j] = 0;
			bestpos[i][j] = 0;
			bestloc[i][j] = 0;
			for (int k = 0; k < MaxN; k++)
			{
				Data[j][k] = -1;
				Right[j][k] = -1;
				evmatrix[i][j][k] = 0;
				delta[i][j][k] = 0;
				cpoint[i][j][k] = 0;
				maxevmatrix[i][j][k] = 0;
				;
			}
		}
	}
}
void restorebest()
{
	for (int i = 0; i <= 1; i++)
	{
		for (int j = 0; j < MaxN; j++)
		{
			pos[i][j] = bestpos[i][j];
			for (int k = 0; k < MaxN; k++)
			{
				evmatrix[i][j][k] = maxevmatrix[i][j][k];
			}
		}
	}
}
int main(int argc,char ** argv)
{
	// ifstream fq;
	// //	string source;
	// //	cin>>source;
	// fq.open("f://batalpha.txt", ios::in);
	// if (!fq.is_open())
	// 	cout << "open file failure" << endl; // f:\\incgraph_2_0.06_5_30_1.20_1.txt
	string filepath(argv[1]);
	st = system_clock::now();
	int t = 1;
	while (t--)
	{

		// char a[256] = {'\0'};
		// fq.getline(a, 500);
		// const int len = sizeof(a) + 1;
		// char *b = (char *)malloc(len);
		// memcpy(b, a, sizeof(a));
		// b[len - 1] = '\0';
		// //		printf("%s\n",b);
		srand((unsigned)time(NULL));
		renew();
		initialize(filepath);
		clock_t start, finish, bestfinish;
		start = clock();
		finish = clock();
		bestfinish = clock();
		datainit();
		int cac = bestnow;
		int runtimes=0;
		while (counttime(st, ed) < timelimit)
		{
			int a = 10;
			int cnt=0;
			while (a != bestnow)
			{
				a = bestnow;
				dp(0);
				dp(1);
				cnt+=2;
				ed = system_clock::now();
				if(counttime(st, ed) >= timelimit)	break;
			}
			if (a <= besthistory)
			{

				if (a < besthistory)
				{

					bestfinish = clock();
					storebest();
					// cout<<a<<endl;
				}
				else if (rand() % 3 == 0)
					storebest();
			}
			else
			{
				restorebest();
			}

			// ed=system_clock::now();
			// cout<<counttime()<<endl;
			// cout<<"cnt: "<<cnt<<endl;
			// break;

			extract();
			runtimes++;
			// finish = clock();
			ed = system_clock::now();
		}
		// cout<<"run times: "<<runtimes<<endl;
		printf("%d+%f\n", besthistory, counttime(st, ed));
	}
	// ed = system_clock::now();
	// cout << counttime(st, ed) << endl;
	ofstream fout("result2.out");
	// fout << "left: \n";
	for (int i = 0; i < size[0]; i++)
	{
		fout << bestpos[0][i] << endl;
	}
	// fout << "\nright: \n";
	for (int i = 0; i < size[1]; i++)
	{
		fout << bestpos[1][i] << endl;
	}
}
