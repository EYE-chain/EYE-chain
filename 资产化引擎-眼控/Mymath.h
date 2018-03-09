//"Mymath.h"//���õ���һЩ��ѧ����
#ifndef _Mymath_h_
#define  _Mymath_h_

#define _SIMULATION_FUN_
#include <math.h>
bool MtxInv( double *A, double *AI,  int n); // A(n*n) -> AI(n*n)// ��������

void MtxMlt( double *A,  double *B, double *AB, int m,  int n,  int s); // A(m*n) * B(n*s) = AB(m*s)// �������(��Ϊ���д洢)

//����
void TreeMtxMlt(double *A, double *B, double *C, double *RE, int m, int n, int l, int s );//A(m*n) *B(n*l) * C(l*s)

void MtxTrs(double *A, double *AT, int m, int n); // A(m*n) -> AT(n*m)// ����ת��


// ����ֵ�ֽ�
int bmuav(double *a,int m,int n,double *u,double *v,double eps,int ka);//��������ֵ�ֽ� M(m*n),ka=max(m,n)+1
//�����a�洢����ֵ��uΪ����������vΪ�������������u*a*v��õ�ԭ����vÿһ�ж�Ӧ����at*a������������
//��Ӧ������ֵ�༴����ֵ��ƽ���Ӵ�С����


void ppp(double *a,double *e,double *s,double *v,int m,int n);
void sss(double *fg,double *cs);

// ��С���˽����Է����� // �ɹ�����true��ʧ�ܷ���false // M(m*n) * X(n*1) = B(m*1) -> X
bool LS(double *M, double *X, double *B, int m, int n);

int agmiv(double *a,int m,int n,double *b,double *x,double *aa,double eps,double *u,double *v,int ka);
int bginv(double *a,int m,int n,double *aa,double eps,double *u,double *v,int ka);


 void cstrq(double *a,//���룺����Ҫת���ľ���
            	int n,//���룺����
        	double *q,//������˻�����
	        double *b,//���������Ϊn,���ضԽ�Ԫ��
	       double *c//��������شζԽ�Ԫ��
		   );//����Գ����Խǻ�



 int csstq(int n,//���룺����
	       double *b,//��ŶԽ�Ԫ�أ���������ֵ
		   double *c,//��ŴζԽ�Ԫ��
		   double *q,//��ų˻����󣬷�������������˳��������ֵ��Ӧ
		   double eps,//����
		   int l//����������
		   );//���Գ����Խ��������ֵ����������

 int aldle(double *a,int n,int m,double *c);//���ԳƷ�����ķֽⷨ

 int mrab1(int a,int b,int *r);//����a\b�����һ�����ȷֲ����������Ҫ������r>=1������
  void mrabs(int a,int b,int *r,int *p,int n);//����a\b�����һ�����ȷֲ��������������,����Ϊn,Ҫ������r>=0

 int brank(double *a,int m,int n);//��m*n�����ȣ����ؽ��Ϊ�ȴ�С

 double sqar(int n,double *a);//��������ģƽ��

 void MatrixAdd(double *A, double *B, double *SUM, int m, int n);//���������
 void MatrixMatNum(double *A, double a,int m, int n);//��������

 void VectorMat(double *a, double *A);//��ά������Ӧ��б�Գƾ���a^b=Ab;^��ʾ�������
 void VectorMatilply(double *a, double *b, double *c);//����ά�������
 double  VectorPoinMut(double *a, double *b,int n);//nά�����ڻ�

 void hpir1(double *x,double *y,int n,double *a,int m,double *dt);//��С����������� x,yΪ�����꣬���Ϊm-1�ζ���ʽ��dt[3]�ֱ�Ϊ�����ƽ���͡�������ֵ�͡����������ֵ
 
 double mgrn1(double u, double g,double *r);//������̬�ֲ��������r����ֵΪu,��׼��Ϊg

 void get_weights(double* samples, int nsize);
 void get_weights_huber(double* samples, int nsize, double scale = 5.0);

#endif