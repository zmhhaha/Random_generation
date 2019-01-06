// QSGS_3D_optimization.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream> //输入输出流头文件
#include <stdlib.h> //标准库头文件
#include <string> //字符串头文件
#include <ctime> //时间函数头文件
#include <fstream> //文件流函数头文件

using namespace std; //为标准程序库提供全局命名空间

const int N = 300; //网格数（不一定非要是方形的）

					 /*四参数参数值修改点*/
const float phi = 0.7; //孔隙率
const float p_cd = 0.01; //生成核概率
const float z = 0.2, f = 0.05; //核生长概率

const float epi = 0.1; //粗放生长误差
const float epi_exact = 0.01; //二次生长误差

/*基于随机抽样方法的孔隙率修正，在孔隙率较大的时候易造成计算量过大*/
const float p_r = 0.8; //随机抽样概率

/*基于生长概率等比例减小的孔隙率修正,在孔隙率较小时容易造成误差过大*/
const float g_r = 0.05;

/*
//生长方向参考D3Q15
Oe ,
nX , pX , nY , pY , nZ , pZ , //主轴方向
nXnYnZ , nXnYpZ , nXpYnZ , nXpYpZ , pXnYnZ , pXnYpZ , pXpYnZ , pXpYpZ //格子正方体的8最长个对角线
*/
const float p_d1 = z, p_d2 = z, p_d3 = z, p_d4 = z, p_d5 = z, p_d6 = z;
const float p_d7 = f, p_d8 = f, p_d9 = f, p_d10 = f, p_d11 = f, p_d12 = f, p_d13 = f, p_d14 = f;

//全局变量
int Solid[N][N][N], Solid_p[N][N][N];
int Grow_Times = 0;

//初始化_开始没有固体
void init()
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				Solid[i][j][k] = 0;
			}
		}
	}
}

void grow_part(int i, int j, int k, float r)
{
	float p_r_d1 = r*p_d1;
	float p_r_d2 = r*p_d2;
	float p_r_d3 = r*p_d3;
	float p_r_d4 = r*p_d4;
	float p_r_d5 = r*p_d5;
	float p_r_d6 = r*p_d6;
	float p_r_d7 = r*p_d7;
	float p_r_d8 = r*p_d8;
	float p_r_d9 = r*p_d9;
	float p_r_d10 = r*p_d10;
	float p_r_d11 = r*p_d11;
	float p_r_d12 = r*p_d12;
	float p_r_d13 = r*p_d13;
	float p_r_d14 = r*p_d14;

	if (i == 0)
	{
		if (j == 0)
		{
			if (k == 0)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d7) Solid_p[i + 1][j + 1][k + 1] = 1;
				}
			}
			else if (k == N - 1)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d8) Solid_p[i + 1][j + 1][k - 1] = 1;
				}
			}
			else
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d7) Solid_p[i + 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d8) Solid_p[i + 1][j + 1][k - 1] = 1;
				}
			}
		}
		else if (j == N - 1)
		{
			if (k == 0)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d9) Solid_p[i + 1][j - 1][k + 1] = 1;
				}
			}
			else if (k == N - 1)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d10) Solid_p[i + 1][j - 1][k - 1] = 1;
				}
			}
			else
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d9) Solid_p[i + 1][j - 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d10) Solid_p[i + 1][j - 1][k - 1] = 1;
				}
			}
		}
		else
		{
			if (k == 0)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d7) Solid_p[i + 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d9) Solid_p[i + 1][j - 1][k + 1] = 1;
				}
			}
			else if (k == N - 1)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d8) Solid_p[i + 1][j + 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d10) Solid_p[i + 1][j - 1][k - 1] = 1;
				}
			}
			else
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d7) Solid_p[i + 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d8) Solid_p[i + 1][j + 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d9) Solid_p[i + 1][j - 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d10) Solid_p[i + 1][j - 1][k - 1] = 1;
				}
			}
		}
	} //i=0的四个点四条边和一个面
	else if (i == N - 1)
	{
		if (j == 0)
		{
			if (k == 0)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d11) Solid_p[i - 1][j + 1][k + 1] = 1;
				}
			}
			else if (k == N - 1)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d12) Solid_p[i - 1][j + 1][k - 1] = 1;
				}
			}
			else
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d11) Solid_p[i - 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d12) Solid_p[i - 1][j + 1][k - 1] = 1;
				}
			}
		}
		else if (j == N - 1)
		{
			if (k == 0)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d13) Solid_p[i - 1][j - 1][k + 1] = 1;
				}
			}
			else if (k == N - 1)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d14) Solid_p[i - 1][j - 1][k - 1] = 1;
				}
			}
			else
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d13) Solid_p[i - 1][j - 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d14) Solid_p[i - 1][j - 1][k - 1] = 1;
				}
			}
		}
		else
		{
			if (k == 0)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d11) Solid_p[i - 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d13) Solid_p[i - 1][j - 1][k + 1] = 1;
				}
			}
			else if (k == N - 1)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d12) Solid_p[i - 1][j + 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d14) Solid_p[i - 1][j - 1][k - 1] = 1;
				}
			}
			else
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d11) Solid_p[i - 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d12) Solid_p[i - 1][j + 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d13) Solid_p[i - 1][j - 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d14) Solid_p[i - 1][j - 1][k - 1] = 1;
				}
			}
		}
	} //i=N-1的四个点四条边和一个面
	else
	{
		if (j == 0)
		{
			if (k == 0)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d7) Solid_p[i + 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d11) Solid_p[i - 1][j + 1][k + 1] = 1;
				}
			}
			else if (k == N - 1)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d8) Solid_p[i + 1][j + 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d12) Solid_p[i - 1][j + 1][k - 1] = 1;
				}
			}
			else
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d7) Solid_p[i + 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d8) Solid_p[i + 1][j + 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d11) Solid_p[i - 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d12) Solid_p[i - 1][j + 1][k - 1] = 1;
				}
			}
		}
		else if (j == N - 1)
		{
			if (k == 0)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d9) Solid_p[i + 1][j - 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d13) Solid_p[i - 1][j - 1][k + 1] = 1;
				}
			}
			else if (k == N - 1)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d10) Solid_p[i + 1][j - 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d14) Solid_p[i - 1][j - 1][k - 1] = 1;
				}
			}
			else
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d9) Solid_p[i + 1][j - 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d10) Solid_p[i + 1][j - 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d13) Solid_p[i - 1][j - 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d14) Solid_p[i - 1][j - 1][k - 1] = 1;
				}
			}
		}
		else
		{
			if (k == 0)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d7) Solid_p[i + 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d9) Solid_p[i + 1][j - 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d11) Solid_p[i - 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d13) Solid_p[i - 1][j - 1][k + 1] = 1;
				}
			}
			else if (k == N - 1)
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d8) Solid_p[i + 1][j + 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d10) Solid_p[i + 1][j - 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d12) Solid_p[i - 1][j + 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d14) Solid_p[i - 1][j - 1][k - 1] = 1;
				}
			}
			else
			{
				if (Solid[i][j][k] == 1)
				{
					if ((rand() / double(RAND_MAX)) < p_r_d1) Solid_p[i + 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d2) Solid_p[i - 1][j][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d3) Solid_p[i][j + 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d4) Solid_p[i][j - 1][k] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d5) Solid_p[i][j][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d6) Solid_p[i][j][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d7) Solid_p[i + 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d8) Solid_p[i + 1][j + 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d9) Solid_p[i + 1][j - 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d10) Solid_p[i + 1][j - 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d11) Solid_p[i - 1][j + 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d12) Solid_p[i - 1][j + 1][k - 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d13) Solid_p[i - 1][j - 1][k + 1] = 1;
					if ((rand() / double(RAND_MAX)) < p_r_d14) Solid_p[i - 1][j - 1][k - 1] = 1;
				}
			}
		}
	} //剩下的四个边和四个面以及内部节点
}

//粗放生长
void grow()
{
	Grow_Times += 1;
	float r = 1;

	for (int i = 0; i < N; i++) //复制上一时间步的状况
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				Solid_p[i][j][k] = Solid[i][j][k];
			}
		}
	}

	//向邻近方向生长
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				grow_part(i, j, k, r);
			}
		}
	}

	//将本时间步的变化赋值给固体节点
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				Solid[i][j][k] = Solid_p[i][j][k];
			}
		}
	}
}

//精细修正

/*基于生长概率等比例减小的孔隙率修正,在孔隙率较小时容易造成误差过大*/
void grow_exact()
{
	float r = g_r;

	for (int i = 0; i < N; i++) //复制上一时间步的状况
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				Solid_p[i][j][k] = Solid[i][j][k];
			}
		}
	}

	//向邻近方向生长
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				grow_part(i, j, k, r);
			}
		}
	}

	//将本时间步的变化赋值给固体节点
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				Solid[i][j][k] = Solid_p[i][j][k];
			}
		}
	}
}


/*基于随机抽样方法的孔隙率修正，在孔隙率较大的时候易造成计算量过大*/
void grow_exact_exact()
{
	float r = 1;
	//复制固体现有状况
	for (int i = 0; i < N; i++) 
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				Solid_p[i][j][k] = Solid[i][j][k];
			}
		}
	}
	for ( int t = 0; t < N*p_r ; t++)
	{
		int i = rand() % N ;
		int j = rand() % N ;
		int k = rand() % N ;
		grow_part(i, j, k, r);
	}

	//将变化赋值给固体节点
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				Solid[i][j][k] = Solid_p[i][j][k];
			}
		}
	}
}


int main()
{
	float phi_p;
	
	//打开孔隙率记录文件
	ofstream fout;
	fout.open("phi.dat", ios::trunc); //ios::trunc 如果文件已存在则先删除文件
	fout.close();
	fout.open("phi.dat", ios::app); //ios::app 所有输出附加在文件末尾

	//生成随机种子
	srand((unsigned)time(NULL));

	//初始化
	init();

	//生成生长核
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				if ((rand() / double(RAND_MAX)) < p_cd)
				{
					Solid[i][j][k] = 1;
				}
			}
		}
	}

	//生长过程
	do
	{
		grow();

		//计算孔隙率
		float t_temp = 0.0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					t_temp += Solid[i][j][k];
				}
			}
		}
		phi_p = 1.0 - t_temp / (N * N * N);
	} while (phi_p - epi > phi);
	cout << "实际孔隙率=" << phi_p << endl;
	fout << "实际孔隙率=" << phi_p << endl;

	//孔隙率临界值修正
	do
	{
		grow_exact();

		//计算孔隙率
		float t_temp = 0.0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					t_temp += Solid[i][j][k];
				}
			}
		}
		phi_p = 1.0 - t_temp / (N * N * N);
	} while (phi_p - epi_exact > phi);

	cout << "实际孔隙率=" << phi_p << endl;
	fout << "实际孔隙率=" << phi_p << endl;

	//随机抽样细微生长
	do
	{
		grow_exact_exact();

		//计算孔隙率
		float t_temp = 0.0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					t_temp += Solid[i][j][k];
				}
			}
		}
		phi_p = 1.0 - t_temp / (N * N * N);
	} while (phi_p > phi);
	
	cout << "实际孔隙率=" << phi_p << endl;
	fout << "实际孔隙率=" << phi_p << endl;
	fout.close();

	//输出数据
	fout.open("data.dat", ios::trunc); //ios::trunc 如果文件已存在则先删除文件
	fout.close();
	fout.open("data.dat", ios::app); //ios::app 所有输出附加在文件末尾
	fout << "TITLE = \"contour\"\nvariables = \"x\", \"y\", \"z\", \"solid\"\nZone I = 300, J = 300, K = 300 F = POINT" << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				fout << i << "\t" << j << "\t" << k << "\t" << Solid[i][j][k] << endl;
			}
		}
	}
	fout.close();

	system("pause");

	return 0;
}

