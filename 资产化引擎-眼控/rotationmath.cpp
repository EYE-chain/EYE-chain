#include "stdafx.h"
#include "rotationmath.h"
//沿x轴正方向，顺时针为正
/********************************************************************
created:	2013/09/16
created:	18:5:2010   20:13
author:		张跃强	
purpose:	标准化角度使其在[head+0　head+360]之间
*********************************************************************/
double standardizeAngle(double a,double head = 0)												// [head,head + 360)
{ 
	if( a >= head )	
		a -= int( (a-head)/360 ) * 360;
	else			
		a += int( -(a-head)/360 + 1 ) * 360;
	return a;
}

/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强	
purpose:	标准化角度使其在[head+0　head+360]之间
*********************************************************************/
 myAngle Get_StandardizeAngle(myAngle angle,double head /*= 0*/)
 {
 	angle.ax = standardizeAngle(angle.ax,head);
 	angle.ay = standardizeAngle(angle.ay,head);
 	angle.az = standardizeAngle(angle.az,head);
 	return angle;
 }

/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强		
purpose:	前后两角∈[-180,180)中间第二个角∈[-90,90]
*********************************************************************/
myAngle Get_StandardizeAngle(myAngle angle)
{
	// 	if( !VM_ISDEGREE)
	// 		angle = Mrad2deg(angle);

	angle.ay = standardizeAngle(angle.ay,-90);
	if(angle.ay > 90)
	{
		angle.ax += 180;
		angle.ay  = 180 - angle.ay;
		angle.az += 180;
	}
	angle.ax = standardizeAngle(angle.ax,-180);
	angle.az = standardizeAngle(angle.az,-180);

	return angle;
}

/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强	
purpose:	x角∈[-180,180), 顺时针为正，得到旋转矩阵
*********************************************************************/
cv::Mat Get_Rx(double ax)
{
	//角度转换为弧度
	double aax = ax*CV_PI/180.0;

	return (Mat_<float>(3,3) <<  	
		1,	       0,	        0,  
		0,	cos(aax),	-sin(aax),
		0,	sin(aax),	 cos(aax) );
}

/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强		
purpose:	y角∈[-90,90), 顺时针为正，得到旋转矩阵
*********************************************************************/
cv::Mat Get_Ry(double ay)
{
	//角度转换为弧度
	double aay = ay*CV_PI/180.0;

	return (Mat_<float>(3,3) <<  	
		cos(aay),	0,	sin(aay),
		0,    1,	      0,
		-sin(aay),	0,	cos(aay));
}

/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强			
purpose:	z角∈[-180,180), 顺时针为正，得到旋转矩阵
*********************************************************************/
cv::Mat Get_Rz(double az)
{
	//角度转换为弧度
	double aaz = az*CV_PI/180.0;

	return (Mat_<float>(3,3) <<  
		cos(aaz),	-sin(aaz),	0,
		sin(aaz),	 cos(aaz),	0,
		0,		   0,   1);
}
/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强	
purpose:    从欧拉角到旋转矩阵，ang - 欧拉角；order-旋转顺序
*********************************************************************/
cv::Mat Get_R(const double* ang , int order)
{
	double ax = ang[0];
	double ay = ang[1];
	double az = ang[2];
	cv::Mat R;
	switch (order)
	{
		case 1:      //xyz
			R = Get_Rz(az)*Get_Ry(ay)*Get_Rx(ax);
			break;
		case 2:      //yzx
			R = Get_Rx(ax)*Get_Rz(az)*Get_Ry(ay);
			break;
		case 3:      //zxy
			R = Get_Ry(ay)*Get_Rx(ax)*Get_Rz(az);
			break;
		case 4:      //xzy
			R = Get_Ry(ay)*Get_Rz(az)*Get_Rx(ax);
			break;
		case 5:      //yxz
			R = Get_Rz(az)*Get_Rx(ax)*Get_Ry(ay);
			break;
		case 6:      //zyx
			R = Get_Rx(ax)*Get_Ry(ay)*Get_Rz(az);
			break;
		default:
			break;
	}
	return R;
}

cv::Mat Get_R(myAngle ang , int order)
{
	double angle[3] = {ang.ax, ang.ay, ang.az};
	return Get_R(angle, order);
}
/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强	
purpose:    从旋转矩阵到欧拉角，myAngle - 欧拉角；order-旋转顺序
*********************************************************************/
myAngle Get_EularAngle(const cv::Mat& R, int order)
{
	myAngle ang;
	switch (order)
	{
	case 1:      //xyz
		ang=Get_Angle_xyz(R);
		break;
	case 2:      //yzx
		ang=Get_Angle_yzx(R);
		if (fabs(ang.az)>90)
			ang=myAngle(ang.ax, -ang.ay, -ang.az);
		break;
	case 3:      //zxy
		ang=Get_Angle_zxy(R);
		if (fabs(ang.ax)>90)
			ang=myAngle(-ang.ax, -ang.ay, ang.az);
		break;
	case 4:      //xzy
		ang=Get_Angle_xzy(R);
		if (fabs(ang.az)>90)
			ang=myAngle(ang.ax, -ang.ay, -ang.az);
		break;
	case 5:      //yxz
		ang=Get_Angle_yxz(R);
		if (fabs(ang.ax)>90)
			ang=myAngle(-ang.ax, -ang.ay, ang.az);
		break;
	case 6:      //zyx
		ang=Get_Angle_zyx(R);
		break;
	default:
		break;
	}
	return ang;

}
/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强	
purpose:    四元数<->旋转矩阵
*********************************************************************/
cv::Mat Get_RFromQuaternion(const cv::Mat& q)
{
	double q1 = q.at<float>(0, 0), q2 = q.at<float>(1, 0), q3 = q.at<float>(2, 0), q4 = q.at<float>(3, 0);
	double pR[9] = {
		q1*q1+q2*q2-q3*q3-q4*q4,	2*(q2*q3-q1*q4),					2*(q2*q4+q1*q3),
		2*(q2*q3+q1*q4),			q1*q1+q3*q3-q2*q2-q4*q4,	        2*(q3*q4-q1*q2),
		2*(q2*q4-q1*q3),			2*(q3*q4+q1*q2),			        q1*q1+q4*q4-q2*q2-q3*q3};
	cv::Mat AA(3, 3,CV_32F);
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			AA.at<float>(i, j) = pR[i*3+j];
		}
	}
		//cout << "AA: " << AA << endl;
	return AA;
}
/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强	
purpose:    四元数<->旋转矩阵
*********************************************************************/
cv::Mat Get_QuaternionFromR(const cv::Mat& R)
{
	//四元数与转轴转角的关系 顺时针旋转 q = [ cos(θ/2), sin(θ/2) * N]
	// q1=cos(θ/2) → cosθ= 2*q1^2 - 1 = (tr+1)/2 - 1 = (tr-1) / 2  ， θ用迹表示
	cv::Scalar rt = cv::trace(R);
	double q1     = sqrt(rt[0]+1) / 2.0; // 第一个规定大于0
	double pQ[4]  = {q1,                                         (R.at<float>(2,1)-R.at<float>(1,2))/(4*q1), 
		            (R.at<float>(0,2)-R.at<float>(2,0))/(4*q1),  (R.at<float>(1,0)-R.at<float>(0,1))/(4*q1)};
	cv::Mat AA(3, 3,CV_32F);
	for (int i=0; i<4; i++)
	{
		AA.at<float>(i, 0) = pQ[i];
		
	}
	//cout << "AA: " << AA << endl;
	return AA;
}
/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强	
purpose:    旋转顺序x->y->z, 沿坐标轴正方向看过去顺时针为正
*********************************************************************/
void EulerAngle2Rm(const double* angle , cv::Mat& R)
{
	double alpha = angle[0], beta = angle[1], gamma = angle[2];
	alpha =  alpha*CV_PI/180.0;
	beta  =  beta*CV_PI/180.0;
	gamma =  gamma*CV_PI/180.0;
	float sx = sin(alpha);
	float cx = cos(alpha);
	float sy = sin(beta);
	float cy = cos(beta);
	float sz = sin(gamma);
	float cz = cos(gamma);
	R.at<float>(0,0) = cy*cz; R.at<float>(0,1) = sx*sy*cz-cx*sz; R.at<float>(0,2) = cx*sy*cz+sx*sz;
	R.at<float>(1,0) = cy*sz; R.at<float>(1,1) = sx*sy*sz+cx*cz; R.at<float>(1,2) = cx*sy*sz-sx*cz;
	R.at<float>(2,0) = -sy;   R.at<float>(2,1) = sx*cy;          R.at<float>(2,2) = cx*cy;
}
/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强	
purpose:    旋转顺序x->y->z, 沿坐标轴正方向看过去顺时针为正
*********************************************************************/
double*  Rm2EulerAngle(const cv::Mat& rm)
{
	double R[3][3];
	double angle[3];
	for (int i=0; i<3;i++)
	{
		for (int j=0; j<3; j++)
		{
			R[i][j] = rm.at<float>(i, j);
		}
	}
	double Rz = atan2(R[1][0],R[0][0]); 
	double Sz = sin(Rz);
	double Cz = cos(Rz);
	double Ry = atan2(-R[2][0], (R[0][0]*Cz+R[1][0]*Sz));
	double Rx = atan2((R[0][2]*Sz-R[1][2]*Cz), (R[1][1]*Cz-R[0][1]*Sz));
	angle[0] = Rx/CV_PI*180;
	angle[1] = Ry/CV_PI*180;
	angle[2] = Rz/CV_PI*180;
	return angle;
}

myAngle  Rm2EulerAngle2(const cv::Mat& rm)
{
	double R[3][3];
	double angle[3];
	/*R必须判断其正负*/
	if(cv::determinant(rm) < 0)
	{
		rm *= -1;
	}
	for (int i=0; i<3;i++)
	{
		for (int j=0; j<3; j++)
		{
			R[i][j] = rm.at<float>(i, j);
		}
	}
	double Rz = atan2(R[1][0],R[0][0]); 
	double Sz = sin(Rz);
	double Cz = cos(Rz);
	double Ry = atan2(-R[2][0], (R[0][0]*Cz+R[1][0]*Sz));
	double Rx = atan2((R[0][2]*Sz-R[1][2]*Cz), (R[1][1]*Cz-R[0][1]*Sz));
	angle[0] = Rx/CV_PI*180;
	angle[1] = Ry/CV_PI*180;
	angle[2] = Rz/CV_PI*180;
	
	return myAngle(angle[0], angle[1], angle[2]);
}
/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强	
purpose:    旋转顺序x->y->z, 沿坐标轴正方向看过去顺时针为正 
            输出:y∈[-90,90],x,z∈[-180,180),输出角度
*********************************************************************/
myAngle Get_Angle(const cv::Mat& R)
{
	double angle[3] = {0};	
	double x, y, z;

	/*R必须判断其正负*/
	if(cv::determinant(R) < 0)
	{
		R *= -1;
	}
	y = asin( -R.at<float>(2,0) );	//[-90,90], cy >=0
	if( fabs(y*180/CV_PI-90 )<1e-6 )
	{
		z = 0;					// 不妨设z=0
		x = atan( R.at<float>(0,1)/R.at<float>(0,2) );// [-90,90]
	}
	else if( fabs(y*180/CV_PI+90 )<1e-6 )
	{
		z = 0;					// 不妨设z=0
		x = atan( R.at<float>(0,1) / R.at<float>(0,2) );// [-90,90]
	}
	else
	{
		x = atan2( R.at<float>(2,1) , R.at<float>(2,2) );   //[-180,180]

		z = atan2( R.at<float>(1,0) , R.at<float>(0,0) );   //[-180,180]
	}	
	angle[0] = x*180/CV_PI;
	angle[1] = y*180/CV_PI;
	angle[2] = z*180/CV_PI;

	/*double rm[3][3];
	double angle[3];
	for (int i=0; i<3;i++)
	{
		for (int j=0; j<3; j++)
		{
			rm[i][j] = R.at<float>(i, j);
		}
	}
	double Rz = atan(rm[1][0]/rm[0][0]); 
	double Sz = sin(Rz);
	double Cz = cos(Rz);
	double Ry = atan(-rm[2][0]/(rm[0][0]*Cz+rm[1][0]*Sz));
	double Rx = atan((rm[0][2]*Sz-rm[1][2]*Cz)/(rm[1][1]*Cz-rm[0][1]*Sz));
	angle[0] = Rx/CV_PI*180;
	angle[1] = Ry/CV_PI*180;
	angle[2] = Rz/CV_PI*180;*/
	return myAngle(angle[0], angle[1], angle[2]);	
}
/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强	
purpose:    旋转顺序x->y->z, 沿坐标轴正方向看过去顺时针为正 
            输出:y∈[-90,90],x,z∈[-180,180),输出角度
*********************************************************************/
myAngle Get_Angle_xyz(const cv::Mat& R)
{
	return Get_Angle(R);	
}

myAngle  Get_Angle_yzx(const cv::Mat& R)
{
	//       ↑zw          ↑yc
	//  ①xw↙→yw  ②zc←↙xc 

	cv::Mat Rpq = (Mat_<float>(3,3) <<   
		0, 1, 0,
		0, 0, 1,
		1, 0, 0);
	myAngle ang = Get_Angle( Rpq * R * Rpq.t());
	return Get_StandardizeAngle(myAngle(ang.az, ang.ax, ang.ay));
}
/********************************************************************
created:	2013/09/16
created:	2013/09/16   10:10
author:		张跃强	
purpose:    旋转顺序z->x->y, 沿坐标轴正方向看过去顺时针为正 
输出:y∈[-90,90],x,z∈[-180,180),输出角度
*********************************************************************/
myAngle Get_Angle_zxy(const cv::Mat& R)
{
	//       ↑zw           ↑xc
	//  ①xw↙→yw    ② yc↙→zc
	cv::Mat Rpq = (Mat_<float>(3,3) <<    
		0, 0, 1,
		1, 0, 0,
		0, 1, 0);
	myAngle ang = Get_Angle( Rpq * R * Rpq.t());
	return Get_StandardizeAngle(myAngle(ang.ay, ang.az, ang.ax));
}

/********************************************************************
右手要反向的 三个,全都以此函数为例
Xw' = Rww' Xw;
现在需要求xw zw yw的顺序转动的角度，都要将轴变成xc yc zc,先不管正负号
xc        xw             1 0  0
yc = Rwc  yw ，显然Rwc = 0 0 -1 
zc        zw             0 1  0
c坐标系中旋转cx cy cz→c'→w' 
Xw' = Rc'w' Rcc' Rwc Xw
Rww' = ~Rwc Rcc' Rwc
也就是待求的wx wz wy对应cx cy cz
xw   1 0  0 xc                      wx   1      cx
yw = 0 0  1 xy交换一下yw和zw顺序得  wz =   -1   cy
zw   0 -1 0 xz                      wy        1 cz
*********************************************************************/
myAngle Get_Angle_xzy(const cv::Mat& R)
{
	//       ↑zw         ↑2        ↑yc
	//  ①xw↙→yw  ②3←↙1  即zc←↙xc   ③ c→c' ④ c'→w'
	cv::Mat Rpq = (Mat_<float>(3,3) <<
		1, 0, 0,
		0, 0, -1,
		0, 1, 0);
	/*Rww' = ~Rwc Rcc' Rwc得Rcc' = Rwc * Rww' * ~Rwc*/
	myAngle aq = Get_Angle( Rpq * R * Rpq.t());
	aq.ay = -aq.ay; 	//角度位置对应，只差一个正负
	return Get_StandardizeAngle(myAngle(aq.ax, aq.az, aq.ay));
}
// 从上面开始都是最后一行跳上去
myAngle Get_Angle_yxz(const cv::Mat&  R)
{
	//       ↑zw          →xc
	//  ①xw↙→yw  ②yc↙↓zc 
	cv::Mat Rpq = (Mat_<float>(3,3) <<
		0, 1, 0,
		1, 0, 0,
		0, 0, -1);
	myAngle ap = Get_Angle(Rpq * R * Rpq.t());
	ap.az = -ap.az;
	return Get_StandardizeAngle(myAngle(ap.ay, ap.ax, ap.az));
}
myAngle Get_Angle_zyx(const cv::Mat& R)
{
	//       ↑zw     xc↑↗zc
	//  ①xw↙→yw  ②   →yc
	cv::Mat Rpq = (Mat_<float>(3,3) <<
		0, 0, 1,
		0, 1, 0,
		-1, 0, 0 );
	myAngle ap = Get_Angle(Rpq * R * Rpq.t());
	ap.az = -ap.az;
	return Get_StandardizeAngle(myAngle(ap.az, ap.ay, ap.ax));
}

myPosi Get_Posi(const cv::Mat& t)
{
	return myPosi(t.at<float>(0, 0), t.at<float>(1, 0), t.at<float>(2, 0));
}
cv::Mat Get_T(const myPosi posi)
{
	return (Mat_<float>(3, 1)<<posi.tx, posi.ty, posi.tz);
}

cv::Mat Get_T(double* posi)
{
	return Get_T(myPosi(posi[0], posi[1], posi[2]));
}

/* The following are some utility routines. */

/* This is the worst random number routine in the known universe, but
   I use it for portability. Feel free to replace it. */
double uniform_random(void)
{
  return (double) rand() / (double) RAND_MAX;
}

double uniform_random2(void)
{
	return 2.0*uniform_random()-1.0;
}


/* This Gaussian routine is stolen from Numerical Recipes and is their
   copyright. */

double gaussian_random(void)
{
  static int next_gaussian = 0;
  static double saved_gaussian_value;

  double fac, rsq, v1, v2;

  if (next_gaussian == 0) {
    do {
      v1 = 2.0*uniform_random()-1.0;
      v2 = 2.0*uniform_random()-1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    saved_gaussian_value=v1*fac;
    next_gaussian=1;
    return v2*fac;
  } else {
    next_gaussian=0;
    return saved_gaussian_value;
  }
}

double evaluate_gaussian(double val, double sigma)
{
  /* This private definition is used for portability */
  static const double Condense_PI = 3.14159265358979323846;

  return 1.0/(sqrt(2.0*Condense_PI) * sigma) *
    exp(-0.5 * (val*val / (sigma*sigma)));
}
