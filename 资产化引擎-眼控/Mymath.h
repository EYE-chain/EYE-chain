//"Mymath.h"//常用到的一些数学函数
#ifndef _Mymath_h_
#define  _Mymath_h_

#define _SIMULATION_FUN_
#include <math.h>
bool MtxInv( double *A, double *AI,  int n); // A(n*n) -> AI(n*n)// 矩阵求逆

void MtxMlt( double *A,  double *B, double *AB, int m,  int n,  int s); // A(m*n) * B(n*s) = AB(m*s)// 矩阵相乘(均为按行存储)

//连乘
void TreeMtxMlt(double *A, double *B, double *C, double *RE, int m, int n, int l, int s );//A(m*n) *B(n*l) * C(l*s)

void MtxTrs(double *A, double *AT, int m, int n); // A(m*n) -> AT(n*m)// 矩阵转置


// 奇异值分解
int bmuav(double *a,int m,int n,double *u,double *v,double eps,int ka);//矩阵奇异值分解 M(m*n),ka=max(m,n)+1
//计算后a存储奇异值，u为左正交矩阵，v为右正交矩阵，如果u*a*v则得到原矩阵，v每一行对应矩阵at*a的特征向量，
//对应的特征值亦即奇异值的平方从大到小排列


void ppp(double *a,double *e,double *s,double *v,int m,int n);
void sss(double *fg,double *cs);

// 最小二乘解线性方程组 // 成功返回true，失败返回false // M(m*n) * X(n*1) = B(m*1) -> X
bool LS(double *M, double *X, double *B, int m, int n);

int agmiv(double *a,int m,int n,double *b,double *x,double *aa,double eps,double *u,double *v,int ka);
int bginv(double *a,int m,int n,double *aa,double eps,double *u,double *v,int ka);


 void cstrq(double *a,//输入：储存要转化的矩阵
            	int n,//输入：阶数
        	double *q,//输出：乘积矩阵
	        double *b,//输出：长度为n,返回对角元素
	       double *c//输出：返回次对角元素
		   );//矩阵对称三对角化



 int csstq(int n,//输入：阶数
	       double *b,//存放对角元素，返回特征值
		   double *c,//存放次对角元素
		   double *q,//存放乘积矩阵，返回特征向量，顺序与特征值对应
		   double eps,//精度
		   int l//最大迭代次数
		   );//求解对称三对角阵的特征值和特征向量

 int aldle(double *a,int n,int m,double *c);//求解对称方程组的分解法

 int mrab1(int a,int b,int *r);//生成a\b区间的一个均匀分布的随机整数要求种子r>=1的奇数
  void mrabs(int a,int b,int *r,int *p,int n);//生成a\b区间的一个均匀分布的随机整数序列,长度为n,要求种子r>=0

 int brank(double *a,int m,int n);//求m*n矩阵秩，返回结果为秩大小

 double sqar(int n,double *a);//求向量的模平方

 void MatrixAdd(double *A, double *B, double *SUM, int m, int n);//两矩阵相加
 void MatrixMatNum(double *A, double a,int m, int n);//矩阵数乘

 void VectorMat(double *a, double *A);//三维向量对应的斜对称矩阵，a^b=Ab;^表示向量叉乘
 void VectorMatilply(double *a, double *b, double *c);//两三维向量叉乘
 double  VectorPoinMut(double *a, double *b,int n);//n维向量内积

 void hpir1(double *x,double *y,int n,double *a,int m,double *dt);//最小二乘曲线拟合 x,y为点坐标，拟合为m-1次多项式，dt[3]分别为：误差平方和、误差绝对值和、最大误差绝对值
 
 double mgrn1(double u, double g,double *r);//产生正态分布的随机数r，均值为u,标准差为g

 void get_weights(double* samples, int nsize);
 void get_weights_huber(double* samples, int nsize, double scale = 5.0);

#endif