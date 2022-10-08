#include<iostream>
#include<cstring>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>
#include<string>
#include<cmath>
#include<math.h>
#include<time.h>
using namespace std;
const int timelimit = 5000;//最优解
const int level = 2;
const int MaxN = 400;
const double alpha = 0.8;
int shuffle[level] = { 0,0 };
//int size[level] = {377,94};//每一层节点个数
//int movable[level] = {263,65};//每一层不可动节点个数（相应节点编号为0~n-1）
int size1[level] = { 0,0 };//每一层节点个数
int immovable[level] = { 0,0 };//每一层不可动节点个数（相应节点编号为0~n-1）
int pos[level][MaxN] = { 0 };//每个位置上的节点编号
int loc[level][MaxN] = { 0 };//每个编号的节点位置
int evmatrix[level][MaxN][MaxN] = { 0 };//评估矩阵，【i】【j】表示i在j前的交叉数
int delta[level][MaxN][MaxN] = { 0 }; //移动的评估矩阵，【i】【j】表示j放到i前（后）的交叉变化数，优化则为正

// 这俩个数组实际没有发挥作用
int cpoint[level][MaxN][MaxN] = { 0 };//向右的连接节点编号
int nc[level][MaxN] = { 0 };//向右的连接节点个数


int bestnow = 10000000;
int besthistory = 10000000;
int bestpos[level][MaxN] = { 0 };
int bestloc[level][MaxN] = { 0 };
int maxevmatrix[level][MaxN][MaxN] = { 0 };//最优解记录

int Data[MaxN][MaxN] = { -1 };
int Right[MaxN][MaxN] = { -1 };
int num_connect[MaxN] = { 0 };
int right_connect[MaxN] = { 0 };//以上四个为二层时方便处理的单元集合

int tenure = 1;   // 禁忌步长
const int length = 1000000;  // 禁忌列表长度(该值很大影响程度也有限，因为受步长约束）
int List_len[level] = { 0 }; //左、右两侧禁忌列表中的元素个数
int Tabulist0[length][3] = { -1 }; //设置禁忌列表
int Tabulist1[length][3] = { -1 };
int times = 0; //单层图的搜索次数

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
    int minnum[MaxN] = { 0 };
    int minplace[MaxN] = { 0 };
    for (int key = immovable[0]; key < size1[0]; key++)
    {
        int count = 0;
        for (int i = 0; i < immovable[0]; i++)//初始化最初位置交叉数
        {
            count += evmatrix[0][key][i];
        }
        int min = count;
        int min_place = 0;
        for (int i = 0; i < immovable[0]; i++)//计算每个位置的交叉数
        {
            count -= evmatrix[0][key][i];
            count += evmatrix[0][i][key];
            if (count < min)
            {
                min = count;
                min_place = i + 1;
            }
        }
        minnum[key] = min;
        minplace[key] = min_place;
    }
    for (int key = immovable[0]; key < size1[0]; key++)
    {
        printf("%d %d\n", minnum[key], minplace[key]);
    }



    getchar();
}

// 判断编号为id,在第t层的节点是不是增量节点
int active(int id, int t)
{
    if (pos[t][id] >= immovable[t])
        return 1;
    return 0;
}
// 返回邻域动作增效最大的一个
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
// 根据insert或是swap动作，更新评估矩阵的值
// dir: 确定由j到i还是由i到j
void renewEvArray(int t, int j, int i, int dir)
{
    if (t == 0)
    {
        if (dir == 0)
        {
            //更新矩阵
            int counter[MaxN] = { 0 };
            for (int k = j; k < i; k++)//k是位置
            {
                for (int q = 0; q < num_connect[pos[t][k]]; q++)//对每一个夹点的邻接点
                {
                    counter[Data[pos[t][k]][q + 2]]++;//被夹点连接的点记录数自增
                }
            }
            for (int k = 0; k < num_connect[pos[t][i]]; k++)//对移动点的每一个节点
            {
                for (int q = 0; q < size1[1 - t]; q++)//对右侧每一个被夹点连接的点
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
            int counter[MaxN] = { 0 };
            for (int k = j + 1; k <= i; k++)//k是位置
            {
                for (int q = 0; q < num_connect[pos[t][k]]; q++)//对每一个夹点的邻接点
                {
                    counter[Data[pos[t][k]][q + 2]]++;//被夹点连接的点记录数自增
                }
            }
            for (int k = 0; k < num_connect[pos[t][j]]; k++)//对移动点的每一个节点
            {
                for (int q = 0; q < size1[1 - t]; q++)//对右侧每一个被夹点连接的点
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
        if (dir == 0)
        {
            int counter[MaxN] = { 0 };
            for (int k = j; k < i; k++)//k是位置
            {
                for (int q = 0; q < right_connect[pos[t][k]]; q++)//对每一个夹点的邻接点
                {

                    counter[Right[pos[t][k]][q + 2]]++;//被夹点连接的点记录数自增
                }
            }
            for (int k = 0; k < right_connect[pos[t][i]]; k++)//对移动点的每一个节点
            {
                for (int q = 0; q < size1[1 - t]; q++)//对右侧每一个被夹点连接的点
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
            int counter[MaxN] = { 0 };
            for (int k = j + 1; k <= i; k++)//k是位置
            {
                for (int q = 0; q < right_connect[pos[t][k]]; q++)//对每一个夹点的邻接点
                {
                    counter[Right[pos[t][k]][q + 2]]++;//被夹点连接的点记录数自增
                }
            }
            for (int k = 0; k < right_connect[pos[t][j]]; k++)//对移动点的每一个节点
            {
                for (int q = 0; q < size1[1 - t]; q++)//对右侧每一个被夹点连接的点
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
void  initialize(char* string)
{
    srand((int)time(NULL));
    ifstream fq;
    cout << string << endl;
    fq.open(string, ios::in);
    if (!fq.is_open())
        cout << "open file failure 2" << endl; //f:\\incgraph_2_0.06_5_30_1.20_1.txt
    int g1 = 0;
    int a;

    int tmp1 = 0;
    //读取
    int left_num = 10;
    int right_num = 10;

    while (!fq.eof())
    {
        char buffer[256] = { '\0' };
        fq.getline(buffer, 500);  //读入每行

        int tmp2 = 0;
        for (int i = 0; i < 255; i++)
        {
            if (buffer[i] >= '0' && buffer[i] <= '9')
            {
                if (i == 0 || (buffer[i - 1] < '0' || buffer[i - 1] > '9'))
                {
                    if (buffer[i + 1] >= '0' && buffer[i + 1] <= '9' && buffer[i + 2] >= '0' && buffer[i + 2] <= '9')
                        a = (int)(buffer[i] - '0') * 100 + (int)(buffer[i + 1] - '0') * 10 + (int)(buffer[i + 2] - '0');
                    else if (buffer[i + 1] >= '0' && buffer[i + 1] <= '9')
                        a = (int)(buffer[i] - '0') * 10 + (int)(buffer[i + 1] - '0');
                    else
                        a = (int)(buffer[i] - '0');

                    if (tmp1 == 1)
                    {
                        size1[tmp2] = a;
                        tmp2++;
                    }
                    if (tmp1 == 2)
                    {
                        left_num = size1[0];
                        right_num = size1[1];
                    }
                    if (tmp1 > 1 && tmp1 < left_num + 2)
                    {
                        Data[tmp1 - 2][tmp2] = a;//
                        tmp2++;
                    }
                    if (tmp1 >= left_num + 2)
                    {
                        Right[tmp1 - left_num - 2][tmp2] = a;
                        tmp2++;
                    }
                }
            }
        }
        if (tmp2 < MaxN && tmp1 < (left_num + 2) && tmp1 > 1)
        {
            // 将范围外的设为-1
            for (int i = tmp2; i < MaxN; i++)
                Data[tmp1 - 2][i] = -1;
        }
        tmp1++;
    }

    fq.close();

    int ct = 0;
    for (int i = 0; i < left_num; i++)
    {
        if (Data[i][0] == 1)
            ct++;
    }

    immovable[0] = ct;
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

    immovable[1] = ct;
    shuffle[0] = int((left_num - immovable[0]) * alpha);
    shuffle[1] = int((right_num - immovable[1]) * alpha);

    int Data1[MaxN][MaxN] = { -1 };
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
    int set[MaxN * 2] = { 0 };

    for (int i = 0; i < right_num; i++)
    {
        if (Right[i][0] == 0 && Right[i][1] > left_num + immovable[1] - 1)
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

    int Right1[MaxN][MaxN] = { 0 };
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
    int rightem[MaxN] = { 0 };
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

    for (int i = 0; i < left_num - immovable[0]; i++)
    {
        int a = rand() % immovable[0];
        int tmp = pos[0][i + immovable[0]];
        for (int k = a; k < i + immovable[0]; k++)
        {
            pos[0][k + 1] = pos[0][k];
        }
        pos[0][a] = tmp;
    }
    for (int i = 0; i < right_num - immovable[1]; i++)
    {
        int a = rand() % immovable[1];
        int tmp = pos[1][i + immovable[1]];
        for (int k = a; k < i + immovable[1]; k++)
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
    int left_num = size1[0];
    int right_num = size1[1];
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
                        if ((loc[0][Right[i][m + 2]] - loc[0][Right[j][n + 2]]) > 0)//若满足位置关系
                        {
                            evmatrix[1][i][j]++;//
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
                        if ((loc[1][Data[i][m + 2]] - loc[1][Data[j][n + 2]]) > 0)//若满足位置关系
                        {
                            evmatrix[0][i][j]++;//
                        }
                    }
                }
            }
        }
    }

}
void dp(int t)//t即type b表示层数
{
    int left_num = size1[0];
    int right_num = size1[1];
    //计算交换delta矩阵

    // 代码开头定义delta和evmatrix的位置有这两个数组的明确含义
    // 根据含义理解下面通过evmatrix计算delta数组的方法
    for (int j = 0; j < size1[t]; j++)
    {
        // 这个变量记录的是j移动到0前面的交叉变化数，相当于j和[0,j-1]的点的前后位置关系都发生了变化，后面几个循环同理
        delta[t][0][j] = 0;
        for (int i = 0; i <= j; i++)
        {
            delta[t][0][j] += (evmatrix[t][pos[t][i]][pos[t][j]] - evmatrix[t][pos[t][j]][pos[t][i]]);
        }
    }
    for (int j = 1; j < size1[t]; j++)
    {
        for (int i = 1; i <= j; i++)
        {
            delta[t][i][j] = delta[t][i - 1][j] + evmatrix[t][pos[t][j]][pos[t][i - 1]] - evmatrix[t][pos[t][i - 1]][pos[t][j]];
        }
    }
    for (int j = size1[t] - 1; j >= 0; j--)
    {
        delta[t][size1[t] - 1][j] = 0;
        for (int i = j; i <= size1[t] - 1; i++)
        {
            delta[t][size1[t] - 1][j] += (evmatrix[t][pos[t][j]][pos[t][i]] - evmatrix[t][pos[t][i]][pos[t][j]]);
        }
    }
    for (int j = 0; j < size1[t]; j++)
    {
        for (int i = size1[t] - 1 - 1; i >= j; i--)
        {
            delta[t][i][j] = delta[t][i + 1][j] + evmatrix[t][pos[t][i + 1]][pos[t][j]] - evmatrix[t][pos[t][j]][pos[t][i + 1]];
        }
    }

    //下面进行动态规划，获得最优运动
    int dir[MaxN] = { 0 };
    int behind[MaxN] = { 0 };
    int cross_delta[MaxN] = { 0 };
    int count = 0;

    for (int i = 1; i < size1[t]; i++)//内部动态规划 !!!!0-30，from 1 effective
    {
        int rem = -1;
        int in_tmp_delta = -10000;
        int in_tmpj = 0;
        int max = 0;
        int tmp_dir = 0;
        int in_dir = 0;
        for (int j = 0; j <= i; j++)
        {
            tmp_dir = findsuit(t, j, i); // 确定三个邻域动作中那个增益最大
            switch (tmp_dir)
            {
            case -1:max = 0; break;
            case 0:max = delta[t][j][i]; break;
            case 1:max = delta[t][i][j]; break;
            case 2:max = delta[t][j][i] + delta[t][i - 1][j]; break;
            }
            if (j == 0)
                count = max;
            else
                count = cross_delta[j - 1] + max; // 即式（12）

            //transact(j+1,i)
            // 优于此前最优解，记录下各变量的值
            if (count >= in_tmp_delta)
            {
                in_tmp_delta = count;//记录一重值
                in_tmpj = j;
                in_dir = tmp_dir;
            }
        }
        cross_delta[i] = in_tmp_delta;  //最优动作改进量
        behind[i] = in_tmpj;            //确定节点j(和相应的i对应)
        dir[i] = in_dir;                //确定最优动作
    }

    int i = size1[t] - 1;
    int cd = cross_delta[i];
    // 上面已经通过动态规划确定了这一步最优的邻域动作
    // 下面就是实施这一动作，对位置数组已经评估矩阵的值进行更新
    // 每个条件语句里的代码都是完成了某个插入操作，也就是将数组中的某个值移动了位置
    if (cd >= 0) {
        while (i > 0)//实施序列变换，并更新矩阵
        {
            int j = behind[i];
            int tmp = pos[t][i];

            // i移动到j前
            if (dir[i] == 0)
            {
                // 将i位置的值将保存
                tmp = pos[t][i];
                renewEvArray(t, j, i, dir[i]);
                // 后面的值挨个前移
                for (int k = i; k > j; k--)
                {
                    pos[t][k] = pos[t][k - 1];
                }
                // 此前i位置的值放到j位置
                pos[t][j] = tmp;

            }
            // j移动到i前
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
            // 交换, 用两次移动表示
            else if (dir[i] == 2)
            {
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
    }
    // 根据pos更新loc
    for (int i = 0; i < size1[t]; i++)//计算loc
    {
        loc[t][pos[t][i]] = i;
    }
    for (int i = 0; i < size1[1 - t]; i++)//计算loc
    {
        loc[1 - t][pos[1 - t][i]] = i;
    }

    // 计算交叉数，算了两遍，程序正确的情况下crossnum和crossnum1的结果应当是一样的
    int crosssum = 0;
    for (int i = 0; i < size1[t]; i++)
    {
        for (int j = i; j < size1[t]; j++)
        {
            crosssum += evmatrix[t][pos[t][i]][pos[t][j]];
        }
    }
    int crosssum1 = 0;
    for (int i = 0; i < size1[1 - t]; i++)
    {
        for (int j = i; j < size1[1 - t]; j++)
        {
            crosssum1 += evmatrix[1 - t][pos[1 - t][i]][pos[1 - t][j]];
        }
    }
    int s = 0;

    //printf("\nt = %d,cross1 = %d; cross2 = %d,delta = %d", t, crosssum, crosssum1, cd);
    bestnow = crosssum;
    if (bestnow < besthistory) {
        //printf("\n当前besthistory = %d", bestnow);
    }
}

void dpts(int t)//t即type b表示层数
{
    int left_num = size1[0];
    int right_num = size1[1];
    //计算交换delta矩阵

    // 代码开头定义delta和evmatrix的位置有这两个数组的明确含义
    // 根据含义理解下面通过evmatrix计算delta数组的方法
    for (int j = 0; j < size1[t]; j++)
    {
        // 这个变量记录的是j移动到0前面的交叉变化数，相当于j和[0,j-1]的点的前后位置关系都发生了变化，后面几个循环同理
        delta[t][0][j] = 0;
        for (int i = 0; i <= j; i++)
        {
            delta[t][0][j] += (evmatrix[t][pos[t][i]][pos[t][j]] - evmatrix[t][pos[t][j]][pos[t][i]]);
        }
    }
    for (int j = 1; j < size1[t]; j++)
    {
        for (int i = 1; i <= j; i++)
        {
            delta[t][i][j] = delta[t][i - 1][j] + evmatrix[t][pos[t][j]][pos[t][i - 1]] - evmatrix[t][pos[t][i - 1]][pos[t][j]];
        }
    }
    for (int j = size1[t] - 1; j >= 0; j--)
    {
        delta[t][size1[t] - 1][j] = 0;
        for (int i = j; i <= size1[t] - 1; i++)
        {
            delta[t][size1[t] - 1][j] += (evmatrix[t][pos[t][j]][pos[t][i]] - evmatrix[t][pos[t][i]][pos[t][j]]);
        }
    }
    for (int j = 0; j < size1[t]; j++)
    {
        for (int i = size1[t] - 1 - 1; i >= j; i--)
        {
            delta[t][i][j] = delta[t][i + 1][j] + evmatrix[t][pos[t][i + 1]][pos[t][j]] - evmatrix[t][pos[t][j]][pos[t][i + 1]];
        }
    }

    //下面进行动态规划，获得最优运动
    int dir[MaxN] = { 0 };
    int behind[MaxN] = { 0 };
    int cross_delta[MaxN] = { 0 };
    int count = 0;

    //做两次动态规划来判断是否突破禁忌
    for (int i = 1; i < size1[t]; i++)//内部动态规划 !!!!0-30，from 1 effective
    {
        int rem = -1;
        int in_tmp_delta = -10000;
        int in_tmpj = 0;
        int max = 0;
        int tmp_dir = 0;
        int in_dir = 0;
        for (int j = 0; j <= i; j++)
        {
            tmp_dir = findsuit(t, j, i); // 确定三个邻域动作中那个增益最大
            switch (tmp_dir)
            {
            case -1:max = 0; break;
            case 0:max = delta[t][j][i]; break;
            case 1:max = delta[t][i][j]; break;
            case 2:max = delta[t][j][i] + delta[t][i - 1][j]; break;
            }
            if (j == 0)
                count = max;
            else
                count = cross_delta[j - 1] + max; // 即式（12）

            //transact(j+1,i)
            // 优于此前最优解，记录下各变量的值
            if (count >= in_tmp_delta)
            {
                in_tmp_delta = count;//记录一重值
                in_tmpj = j;
                in_dir = tmp_dir;
            }
        }
        cross_delta[i] = in_tmp_delta;  //最优动作改进量
        behind[i] = in_tmpj;            //确定节点j(和相应的i对应)
        dir[i] = in_dir;                //确定最优动作
    }

    int i = size1[t] - 1;
    int cd0 = cross_delta[i];
    //printf("\n非禁忌的改进量为 = %d\n", cd0);
    if (cd0 >= 0) {  //如果能改进最优解，无视禁忌
    // 上面已经通过动态规划确定了这一步最优的邻域动作
    // 下面就是实施这一动作，对位置数组已经评估矩阵的值进行更新
    // 每个条件语句里的代码都是完成了某个插入操作，也就是将数组中的某个值移动了位置
        while (i > 0)//实施序列变换，并更新矩阵
        {
            int j = behind[i];
            int tmp = pos[t][i];

            if (t == 0) {              //加入禁忌点
                int K = List_len[t];
                if (K < length) {
                    Tabulist0[K][0] = j;
                    Tabulist0[K][1] = i;
                    Tabulist0[K][2] = times;
                }
                else {
                    Tabulist0[length - 1][0] = j;
                    Tabulist0[length - 1][1] = i;
                    Tabulist0[length - 1][2] = times;
                    for (i = 1; i < length; i++) {
                        Tabulist0[i - 1][0] = Tabulist0[i][0];
                        Tabulist0[i - 1][1] = Tabulist0[i][1];
                        Tabulist0[i - 1][2] = Tabulist0[i][2];
                    }
                }
                List_len[t] += 1;
            }
            if (t == 1) {              //加入禁忌点
                int K = List_len[t];
                if (K < length) {
                    Tabulist1[K][0] = j;
                    Tabulist1[K][1] = i;
                    Tabulist1[K][2] = times;
                }
                else {
                    Tabulist1[length - 1][0] = j;
                    Tabulist1[length - 1][1] = i;
                    Tabulist1[length - 1][2] = times;
                    for (i = 1; i < length; i++) {
                        Tabulist1[i - 1][0] = Tabulist1[i][0];
                        Tabulist1[i - 1][1] = Tabulist1[i][1];
                        Tabulist1[i - 1][2] = Tabulist1[i][2];
                    }
                }
                List_len[t] += 1;
            }

            // i移动到j前
            if (dir[i] == 0)
            {
                // 将i位置的值将保存
                tmp = pos[t][i];
                renewEvArray(t, j, i, dir[i]);
                // 后面的值挨个前移
                for (int k = i; k > j; k--)
                {
                    pos[t][k] = pos[t][k - 1];
                }
                // 此前i位置的值放到j位置
                pos[t][j] = tmp;

            }
            // j移动到i前
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
            // 交换, 用两次移动表示
            else if (dir[i] == 2)
            {
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
    }
    else { //如果不能改进最优解，动态规划时考虑禁忌，挑选非禁忌的最优解
        for (int i = 1; i < size1[t]; i++)//内部动态规划 !!!!0-30，from 1 effective
        {
            int rem = -1;
            int in_tmp_delta = -10000;
            int in_tmpj = 0;
            int max = 0;
            int tmp_dir = 0;
            int in_dir = 0;
            for (int j = 0; j <= i; j++)
            {
                tmp_dir = findsuit(t, j, i); // 确定三个邻域动作中那个增益最大
                switch (tmp_dir)
                {
                case -1:max = 0; break;
                case 0:max = delta[t][j][i]; break;
                case 1:max = delta[t][i][j]; break;
                case 2:max = delta[t][j][i] + delta[t][i - 1][j]; break;
                }

                //*如果是禁忌节点对，相应动作赋值为0
                int judge = 0;  //判断是否被禁忌的条件
                judge = 0;

                if (t == 0) {
                    for (int n = 0; n < min(List_len[t], length); n++) {
                        if (j == Tabulist0[n][0] && i == Tabulist0[n][1] && times - Tabulist0[n][2] <= tenure) {
                            judge = judge + 1;
                        }
                    }
                }
                if (t == 1) {
                    for (int n = 0; n < min(List_len[t], length); n++) {
                        if (j == Tabulist1[n][0] && i == Tabulist1[n][1] && times - Tabulist1[n][2] <= tenure) {
                            judge = judge + 1;
                        }
                    }
                }
                if (judge != 0) { //如果是禁忌的，该节点对之间的动作改变量赋值为0
                    max = 0;
                }

                if (j == 0)
                    count = max;
                else
                    count = cross_delta[j - 1] + max; // 即式（12）

                //transact(j+1,i)
                // 优于此前最优解，记录下各变量的值
                if (count >= in_tmp_delta)
                {
                    in_tmp_delta = count;//记录一重值
                    in_tmpj = j;
                    in_dir = tmp_dir;
                }
            }
            cross_delta[i] = in_tmp_delta;  //最优动作改进量
            behind[i] = in_tmpj;            //确定节点j(和相应的i对应)
            dir[i] = in_dir;                //确定最优动作
        }
        int i = size1[t] - 1;
        while (i > 0)//实施序列变换，并更新矩阵
        {
            int j = behind[i];
            int tmp = pos[t][i];

            if (t == 0) {              //加入禁忌点
                int K = List_len[t];
                if (K < length) {
                    Tabulist0[K][0] = j;
                    Tabulist0[K][1] = i;
                    Tabulist0[K][2] = times;
                }
                else {
                    Tabulist0[length - 1][0] = j;
                    Tabulist0[length - 1][1] = i;
                    Tabulist0[length - 1][2] = times;
                    for (i = 1; i < length; i++) {
                        Tabulist0[i - 1][0] = Tabulist0[i][0];
                        Tabulist0[i - 1][1] = Tabulist0[i][1];
                        Tabulist0[i - 1][2] = Tabulist0[i][2];
                    }
                }
                List_len[t] += 1;
            }
            if (t == 1) {              //加入禁忌点
                int K = List_len[t];
                if (K < length) {
                    Tabulist1[K][0] = j;
                    Tabulist1[K][1] = i;
                    Tabulist1[K][2] = times;
                }
                else {
                    Tabulist1[length - 1][0] = j;
                    Tabulist1[length - 1][1] = i;
                    Tabulist1[length - 1][2] = times;
                    for (i = 1; i < length; i++) {
                        Tabulist1[i - 1][0] = Tabulist1[i][0];
                        Tabulist1[i - 1][1] = Tabulist1[i][1];
                        Tabulist1[i - 1][2] = Tabulist1[i][2];
                    }
                }
                List_len[t] += 1;
            }

            // i移动到j前
            if (dir[i] == 0)
            {
                // 将i位置的值将保存
                tmp = pos[t][i];
                renewEvArray(t, j, i, dir[i]);
                // 后面的值挨个前移
                for (int k = i; k > j; k--)
                {
                    pos[t][k] = pos[t][k - 1];
                }
                // 此前i位置的值放到j位置
                pos[t][j] = tmp;

            }
            // j移动到i前
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
            // 交换, 用两次移动表示
            else if (dir[i] == 2)
            {
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
    }

    i = size1[t] - 1;
    int cd = cross_delta[i];
    // 上面已经通过动态规划确定了这一步最优的邻域动作
    // 下面就是实施这一动作，对位置数组已经评估矩阵的值进行更新
    // 每个条件语句里的代码都是完成了某个插入操作，也就是将数组中的某个值移动了位置

    // 根据pos更新loc
    for (int i = 0; i < size1[t]; i++)//计算loc
    {
        loc[t][pos[t][i]] = i;
    }
    for (int i = 0; i < size1[1 - t]; i++)//计算loc
    {
        loc[1 - t][pos[1 - t][i]] = i;
    }

    // 计算交叉数，算了两遍，程序正确的情况下crossnum和crossnum1的结果应当是一样的
    int crosssum = 0;
    for (int i = 0; i < size1[t]; i++)
    {
        for (int j = i; j < size1[t]; j++)
        {
            crosssum += evmatrix[t][pos[t][i]][pos[t][j]];
        }
    }
    int crosssum1 = 0;
    for (int i = 0; i < size1[1 - t]; i++)
    {
        for (int j = i; j < size1[1 - t]; j++)
        {
            crosssum1 += evmatrix[1 - t][pos[1 - t][i]][pos[1 - t][j]];
        }
    }
    int s = 0;

    //printf("\nt = %d,cross1 = %d; cross2 = %d,delta = %d", t, crosssum, crosssum1, cd);
    bestnow = crosssum;
    if (bestnow < besthistory) {
        //printf("\n当前besthistory = %d", bestnow);
    }
}

// 随机扰乱
void extract()
{
    int left_num = size1[0];
    int right_num = size1[1];
    int tmppoint[MaxN] = { 0 };//记录位置
    int count = 0;
    for (int i = 0; i < left_num; i++)
    {
        if (pos[0][i] >= immovable[0])
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

            if (ta > tb)//让ta较小
            {
                int q = ta;
                ta = tb;
                tb = q;
            }
            int i = tmppoint[tb];
            int j = tmppoint[ta];
            //进行交换和更新
            int	tmp = pos[t][i];
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

            do {
                miku = rand() % left_num;
            } while (pos[0][miku] < immovable[0]);
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
        if (pos[1][i] >= immovable[1])
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

            if (ta > tb)//让ta较小
            {
                int q = ta;
                ta = tb;
                tb = q;
            }
            int i = tmppoint[tb];
            int j = tmppoint[ta];
            //进行交换和更新
            int	tmp = pos[t][i];
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

            do {
                miku = rand() % right_num;
            } while (pos[1][miku] < immovable[1]);
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
// 记录下最优解的位置和评估矩阵信息
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
// 重置数组, 把用到的数组都初始化，根据需要清零或是设为某个值
void renew()
{
    bestnow = 10000000;
    besthistory = 10000000;

    for (int i = 0; i < 2; i++)
    {
        shuffle[i] = 0;
        size1[i] = 0;
        immovable[i] = 0;
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
                maxevmatrix[i][j][k] = 0;;
            }
        }
    }
}
// 恢复最优解，就是把记录的最优解的位置和评估矩阵信息复原到dp使用的数组中
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

int main()
{
    // ifstream fq;
    // //fstream iofile;
    // fq.open("E:\\bat alpha4.txt", ios::in);
    // if (!fq.is_open())
    //     cout << "open file failure 1" << endl; //f:\\incgraph_2_0.06_5_30_1.20_1.txt
    // // 循环读取运行不同的测试集
    //iofile.open("result1.txt", ios::out);//创建txt文件，并以写入的模式打开
    // while (1)
    {
        // char a[256] = { '\0' };
        // fq.getline(a, 500); // 读取测试集的路径
        // if (strlen(a) == 0) break;
        // const int len = sizeof(a) + 1;
        // char* b = (char*)malloc(len);
        // memcpy(b, a, sizeof(a));
        // b[len - 1] = '\0';
        srand((unsigned)time(NULL)); // 改变随机数种子，使每次运行结果不同
        renew(); // 重置数组
        initialize("/home/congyu/dpls/instances/GB_1_rnd1_01/GB_1_rnd1_01_0001_30.txt"); // 读入数据及一些简单的初始化
        datainit(); // 计算评估矩阵

        times = 0;
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < 3; j++) {
                Tabulist0[i][j] = { -1 };
                Tabulist1[i][j] = { -1 };
            }
        }
        List_len[0] = 0;
        List_len[1] = 0;

        // 初始化了三个计时的变量
        clock_t start, finish, bestfinish;
        start = clock();
        finish = clock();
        bestfinish = clock();

        // 下面是论文算法的主体, 对应于论文中的伪代码算法1
        // 终止条件(与论文中用次数做条件稍有不同)
        // 当间隔了timelimit时间都没有搜索到更优解的时候停止搜索

        int n = 0;
        while (difftime(finish, bestfinish) < timelimit)
        {
            if (n > 1000 or difftime(finish, bestfinish) > timelimit) {
                break;
            }
            int a = 10;
            // 使用while循环因为在dp的时候最优值bestnow可能会变化
            // 持续dp直到左右两层的邻域内不存在更优邻域解, 对应论文算法2
            while (a != bestnow)
            {
                a = bestnow;
                // 搜索左边一侧的邻域更优解
                dpts(0);
                if (a < besthistory)
                {
                    bestfinish = clock();
                    storebest();
                    n = 0;
                }
                // 搜索右边一侧的邻域更优解
                dpts(1);
                if (a < besthistory)
                {
                    bestfinish = clock();
                    storebest();
                    n = 0;
                }
                times += 1;
            }
            //printf("\n当前a的值为%d\n", a);
            if (a <= besthistory)
            {
                // 优于此前最优解，记录下最优解
                if (a < besthistory)
                {
                    bestfinish = clock();
                    storebest();
                    //n = 0;
                }
                // 等于此前最优解，随机确定采用原值还是现值
                else if (rand() % 3 == 0)
                    storebest();
            }
            else
            {
                // 不如此前最优解，回到最优解的状态
                restorebest();
            }

            // 随机扰乱
            extract();
            for (int i = 0; i < length; i++) {
                for (int j = 0; j < 3; j++) {
                    Tabulist0[i][j] = { -1 };
                    Tabulist1[i][j] = { -1 };
                }
            }
            List_len[0] = 0;
            List_len[1] = 0;

            n = n + 1;
            //printf("\n我扰动%d次啦！", n);

            // 一轮结束时间记录
            finish = clock();
        }
        printf("%d %f\n", besthistory, 0.001 * (difftime(bestfinish, start) + 1));
        for(int i=0;i<size1[0];i++) cout<<bestpos[0][i]<<endl;
        for(int i=0;i<size1[1];i++) cout<<bestpos[1][i]<<endl;

        //double time1 = 0.001 * (difftime(bestfinish, start) + 1);
        //iofile << besthistory <<"  " << time1 << endl;
    }
    //iofile.close();
    return 0;
}
