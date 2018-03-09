#ifndef ROTATION_H
#define ROTATION_H
#include <cv.h>
using namespace cv;

/********************************************************************
created:	2013/09/16
created:	2013/09/15   20:13
author:	    张跃强
purpose:	旋转矩阵同欧拉角相互转换函数库
*********************************************************************/

struct myAngle 
{
	double ax;
	double ay;
	double az;
	myAngle(){ax = ay = az = 0;};
	myAngle(double x, double y, double z){ax = x; ay = y; az = z;};
	myAngle operator - ()
	{
		return myAngle( -ax, -ay , -az);
	};
	friend myAngle operator + (const myAngle &a1, const myAngle &a2)
	{
		return myAngle( a1.ax + a2.ax, a1.ay + a2.ay,  a1.az + a2.az);
	};
	friend myAngle operator / (const myAngle &a1, const double &b)
	{
		return myAngle( a1.ax/(b+1e-5), a1.ay/(b+1e-5),  a1.az/(b+1e-5));
	};
	friend myAngle operator - (const myAngle &a1, const myAngle &a2)
	{
		return myAngle( a1.ax - a2.ax, a1.ay - a2.ay,  a1.az - a2.az);
	};
};
struct myPosi 
{
	double tx;
	double ty;
	double tz;
	myPosi(){tx = ty = tz = 0;};
	myPosi(double x, double y, double z){tx = x; ty = y; tz = z;};
	myPosi operator - ()
	{
		return myPosi( -tx, -ty , -tz);
	};
	friend myPosi operator + (const myPosi &a1, const myPosi &a2)
	{
		return myPosi( a1.tx + a2.tx, a1.ty + a2.ty,  a1.tz + a2.tz);
	};
	friend myPosi operator - (const myPosi &a1, const myPosi &a2)
	{
		return myPosi( a1.tx - a2.tx, a1.ty - a2.ty,  a1.tz - a2.tz);
	};
	friend myPosi operator / (const myPosi &a1, const double &b)
	{
		return myPosi( a1.tx/(b+1e-5), a1.ty/(b+1e-5),  a1.tz/(b+1e-5));
	};
};

cv::Mat Get_RFromQuaternion(const cv::Mat& q);

cv::Mat Get_QuaternionFromR(const cv::Mat& R);

//旋转顺序x->y->z, 沿坐标轴正方向看过去顺时针为正
void EulerAngle2Rm(const double* angle , cv::Mat& R);

//相同顶定义下的欧拉角求解, 旋转顺序x->y->z, 沿坐标轴正方向看过去顺时针为正
void Rm2EulerAngle(const cv::Mat& rm , double* angle);

//欧拉角到旋转矩阵
/***********************************************    
       1:      xyz
       2:      yzx
	   3:      zxy
	   4:      xzy
       5:      yxz
       6:      zyx                
*************************************************/
cv::Mat Get_Rx(double ax);

cv::Mat Get_Ry(double ay);

cv::Mat Get_Rz(double az);

cv::Mat Get_R(const double* ang = NULL,  int order = 1);

cv::Mat Get_R(myAngle ang , int order = 1);

myAngle Get_EularAngle(const cv::Mat& R, int order = 1);


//右手，正向
myAngle  Get_Angle_xyz(const cv::Mat& R);

myAngle  Get_Angle_yzx(const cv::Mat& R);

myAngle  Get_Angle_zxy(const cv::Mat& R);
//右手要反向的 三个,全都以此函数为例
myAngle  Get_Angle_xzy(const cv::Mat& R);

myAngle  Get_Angle_yxz(const cv::Mat& R);

myAngle Get_Angle_zyx(const cv::Mat& R);
//cv::mat -> myPosi
myPosi Get_Posi(const cv::Mat& t);
cv::Mat Get_T(const myPosi posi);
cv::Mat Get_T(double* posi);
myAngle  Rm2EulerAngle2(const cv::Mat& rm);

double uniform_random(void);
double gaussian_random(void);
double uniform_random2(void);
double evaluate_gaussian(double val, double sigma);
#endif