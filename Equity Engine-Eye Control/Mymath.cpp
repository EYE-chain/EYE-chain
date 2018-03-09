#include "stdafx.h"
#include "Mymath.h"
#include "stdio.h"
#include "stdlib.h"
//*************************************************************************//
//矩阵奇异值分解
int bmuav(double *a,int m,int n,double *u,double *v,double eps,int ka)
{
	int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
    double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh,fg[2],cs[2];
    double *s,*e,*w;
    s=(double*)malloc(ka*sizeof(double));
    e=(double*)malloc(ka*sizeof(double));
    w=(double*)malloc(ka*sizeof(double));
    it=60; k=n;
    if (m-1<n) k=m-1;
    l=m;
    if (n-2<m) l=n-2;
    if (l<0) l=0;
    ll=k;
    if (l>k) ll=l;
    if (ll>=1)
      { for (kk=1; kk<=ll; kk++)
          { if (kk<=k)
              { d=0.0;
                for (i=kk; i<=m; i++)
                  { ix=(i-1)*n+kk-1; d=d+a[ix]*a[ix];}
                s[kk-1]=sqrt(d);
                if (s[kk-1]!=0.0)
                  { ix=(kk-1)*n+kk-1;
                    if (a[ix]!=0.0)
                      { s[kk-1]=fabs(s[kk-1]);
                        if (a[ix]<0.0) s[kk-1]=-s[kk-1];
                      }
                    for (i=kk; i<=m; i++)
                      { iy=(i-1)*n+kk-1;
                        a[iy]=a[iy]/s[kk-1];
                      }
                    a[ix]=1.0+a[ix];
                  }
                s[kk-1]=-s[kk-1];
              }
            if (n>=kk+1)
              { for (j=kk+1; j<=n; j++)
                  { if ((kk<=k)&&(s[kk-1]!=0.0))
                      { d=0.0;
                        for (i=kk; i<=m; i++)
                          { ix=(i-1)*n+kk-1;
                            iy=(i-1)*n+j-1;
                            d=d+a[ix]*a[iy];
                          }
                        d=-d/a[(kk-1)*n+kk-1];
                        for (i=kk; i<=m; i++)
                          { ix=(i-1)*n+j-1;
                            iy=(i-1)*n+kk-1;
                            a[ix]=a[ix]+d*a[iy];
                          }
                      }
                    e[j-1]=a[(kk-1)*n+j-1];
                  }
              }
            if (kk<=k)
              { for (i=kk; i<=m; i++)
                  { ix=(i-1)*m+kk-1; iy=(i-1)*n+kk-1;
                    u[ix]=a[iy];
                  }
              }
            if (kk<=l)
              { d=0.0;
                for (i=kk+1; i<=n; i++)
                  d=d+e[i-1]*e[i-1];
                e[kk-1]=sqrt(d);
                if (e[kk-1]!=0.0)
                  { if (e[kk]!=0.0)
                      { e[kk-1]=fabs(e[kk-1]);
                        if (e[kk]<0.0) e[kk-1]=-e[kk-1];
                      }
                    for (i=kk+1; i<=n; i++)
                      e[i-1]=e[i-1]/e[kk-1];
                    e[kk]=1.0+e[kk];
                  }
                e[kk-1]=-e[kk-1];
                if ((kk+1<=m)&&(e[kk-1]!=0.0))
                  { for (i=kk+1; i<=m; i++) w[i-1]=0.0;
                    for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                        w[i-1]=w[i-1]+e[j-1]*a[(i-1)*n+j-1];
                    for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                        { ix=(i-1)*n+j-1;
                          a[ix]=a[ix]-w[i-1]*e[j-1]/e[kk];
                        }
                  }
                for (i=kk+1; i<=n; i++)
                  v[(i-1)*n+kk-1]=e[i-1];
              }
          }
      }
    mm=n;
    if (m+1<n) mm=m+1;
    if (k<n) s[k]=a[k*n+k];
    if (m<mm) s[mm-1]=0.0;
    if (l+1<mm) e[l]=a[l*n+mm-1];
    e[mm-1]=0.0;
    nn=m;
    if (m>n) nn=n;
    if (nn>=k+1)
      { for (j=k+1; j<=nn; j++)
          { for (i=1; i<=m; i++)
              u[(i-1)*m+j-1]=0.0;
            u[(j-1)*m+j-1]=1.0;
          }
      }
    if (k>=1)
      { for (ll=1; ll<=k; ll++)
          { kk=k-ll+1; iz=(kk-1)*m+kk-1;
            if (s[kk-1]!=0.0)
              { if (nn>=kk+1)
                  for (j=kk+1; j<=nn; j++)
                    { d=0.0;
                      for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+kk-1;
                          iy=(i-1)*m+j-1;
                          d=d+u[ix]*u[iy]/u[iz];
                        }
                      d=-d;
                      for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+j-1;
                          iy=(i-1)*m+kk-1;
                          u[ix]=u[ix]+d*u[iy];
                        }
                    }
                  for (i=kk; i<=m; i++)
                    { ix=(i-1)*m+kk-1; u[ix]=-u[ix];}
                  u[iz]=1.0+u[iz];
                  if (kk-1>=1)
                    for (i=1; i<=kk-1; i++)
                      u[(i-1)*m+kk-1]=0.0;
              }
            else
              { for (i=1; i<=m; i++)
                  u[(i-1)*m+kk-1]=0.0;
                u[(kk-1)*m+kk-1]=1.0;
              }
          }
      }
    for (ll=1; ll<=n; ll++)
      { kk=n-ll+1; iz=kk*n+kk-1;
        if ((kk<=l)&&(e[kk-1]!=0.0))
          { for (j=kk+1; j<=n; j++)
              { d=0.0;
                for (i=kk+1; i<=n; i++)
                  { ix=(i-1)*n+kk-1; iy=(i-1)*n+j-1;
                    d=d+v[ix]*v[iy]/v[iz];
                  }
                d=-d;
                for (i=kk+1; i<=n; i++)
                  { ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
                    v[ix]=v[ix]+d*v[iy];
                  }
              }
          }
        for (i=1; i<=n; i++)
          v[(i-1)*n+kk-1]=0.0;
        v[iz-n]=1.0;
      }
    for (i=1; i<=m; i++)
    for (j=1; j<=n; j++)
      a[(i-1)*n+j-1]=0.0;
    m1=mm; it=60;
    while (1==1)
      { if (mm==0)
          { ppp(a,e,s,v,m,n);
            free(s); free(e); free(w); return(1);
          }
        if (it==0)
          { ppp(a,e,s,v,m,n);
            free(s); free(e); free(w); return(-1);
          }
        kk=mm-1;
	while ((kk!=0)&&(fabs(e[kk-1])!=0.0))
          { d=fabs(s[kk-1])+fabs(s[kk]);
            dd=fabs(e[kk-1]);
            if (dd>eps*d) kk=kk-1;
            else e[kk-1]=0.0;
          }
        if (kk==mm-1)
          { kk=kk+1;
            if (s[kk-1]<0.0)
              { s[kk-1]=-s[kk-1];
                for (i=1; i<=n; i++)
                  { ix=(i-1)*n+kk-1; v[ix]=-v[ix];}
              }
            while ((kk!=m1)&&(s[kk-1]<s[kk]))
              { d=s[kk-1]; s[kk-1]=s[kk]; s[kk]=d;
                if (kk<n)
                  for (i=1; i<=n; i++)
                    { ix=(i-1)*n+kk-1; iy=(i-1)*n+kk;
                      d=v[ix]; v[ix]=v[iy]; v[iy]=d;
                    }
                if (kk<m)
                  for (i=1; i<=m; i++)
                    { ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
                      d=u[ix]; u[ix]=u[iy]; u[iy]=d;
                    }
                kk=kk+1;
              }
            it=60;
            mm=mm-1;
          }
        else
          { ks=mm;
            while ((ks>kk)&&(fabs(s[ks-1])!=0.0))
              { d=0.0;
                if (ks!=mm) d=d+fabs(e[ks-1]);
                if (ks!=kk+1) d=d+fabs(e[ks-2]);
                dd=fabs(s[ks-1]);
                if (dd>eps*d) ks=ks-1;
                else s[ks-1]=0.0;
              }
            if (ks==kk)
              { kk=kk+1;
                d=fabs(s[mm-1]);
                t=fabs(s[mm-2]);
                if (t>d) d=t;
                t=fabs(e[mm-2]);
                if (t>d) d=t;
                t=fabs(s[kk-1]);
                if (t>d) d=t;
                t=fabs(e[kk-1]);
                if (t>d) d=t;
                sm=s[mm-1]/d; sm1=s[mm-2]/d;
                em1=e[mm-2]/d;
                sk=s[kk-1]/d; ek=e[kk-1]/d;
                b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
                c=sm*em1; c=c*c; shh=0.0;
                if ((b!=0.0)||(c!=0.0))
                  { shh=sqrt(b*b+c);
                    if (b<0.0) shh=-shh;
                    shh=c/(b+shh);
                  }
                fg[0]=(sk+sm)*(sk-sm)-shh;
                fg[1]=sk*ek;
                for (i=kk; i<=mm-1; i++)
                  { sss(fg,cs);
                    if (i!=kk) e[i-2]=fg[0];
                    fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
                    e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
                    fg[1]=cs[1]*s[i];
                    s[i]=cs[0]*s[i];
                    if ((cs[0]!=1.0)||(cs[1]!=0.0))
                      for (j=1; j<=n; j++)
                        { ix=(j-1)*n+i-1;
                          iy=(j-1)*n+i;
                          d=cs[0]*v[ix]+cs[1]*v[iy];
                          v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                          v[ix]=d;
                        }
                    sss(fg,cs);
                    s[i-1]=fg[0];
                    fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
                    s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
                    fg[1]=cs[1]*e[i];
                    e[i]=cs[0]*e[i];
                    if (i<m)
                      if ((cs[0]!=1.0)||(cs[1]!=0.0))
                        for (j=1; j<=m; j++)
                          { ix=(j-1)*m+i-1;
                            iy=(j-1)*m+i;
                            d=cs[0]*u[ix]+cs[1]*u[iy];
                            u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                            u[ix]=d;
                          }
                  }
                e[mm-2]=fg[0];
                it=it-1;
              }
            else
              { if (ks==mm)
                  { kk=kk+1;
                    fg[1]=e[mm-2]; e[mm-2]=0.0;
                    for (ll=kk; ll<=mm-1; ll++)
                      { i=mm+kk-ll-1;
                        fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        if (i!=kk)
                          { fg[1]=-cs[1]*e[i-2];
                            e[i-2]=cs[0]*e[i-2];
                          }
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=n; j++)
                            { ix=(j-1)*n+i-1;
                              iy=(j-1)*n+mm-1;
                              d=cs[0]*v[ix]+cs[1]*v[iy];
                              v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                              v[ix]=d;
                            }
                      }
                  }
                else
                  { kk=ks+1;
                    fg[1]=e[kk-2];
                    e[kk-2]=0.0;
                    for (i=kk; i<=mm; i++)
                      { fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        fg[1]=-cs[1]*e[i-1];
                        e[i-1]=cs[0]*e[i-1];
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=m; j++)
                            { ix=(j-1)*m+i-1;
                              iy=(j-1)*m+kk-2;
                              d=cs[0]*u[ix]+cs[1]*u[iy];
                              u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                              u[ix]=d;
                            }
                      }
                  }
              }
          }
      }
    return(1);
  }

void ppp(double *a,double *e,double *s,double *v,int m,int n)
  {
	int i,j,p,q;
    double d;
    if (m>=n) i=n;
    else i=m;
    for (j=1; j<=i-1; j++)
      { a[(j-1)*n+j-1]=s[j-1];
        a[(j-1)*n+j]=e[j-1];
      }
    a[(i-1)*n+i-1]=s[i-1];
    if (m<n) a[(i-1)*n+i]=e[i-1];
    for (i=1; i<=n-1; i++)
    for (j=i+1; j<=n; j++)
      { p=(i-1)*n+j-1; q=(j-1)*n+i-1;
        d=v[p]; v[p]=v[q]; v[q]=d;
      }
    return;
  }

void sss(double *fg,double *cs)
  { double r,d;
    if ((fabs(fg[0])+fabs(fg[1]))==0.0)
      { cs[0]=1.0; cs[1]=0.0; d=0.0;}
    else 
      { d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
        if (fabs(fg[0])>fabs(fg[1]))
          { d=fabs(d);
            if (fg[0]<0.0) d=-d;
          }
        if (fabs(fg[1])>=fabs(fg[0]))
          { d=fabs(d);
            if (fg[1]<0.0) d=-d;
          }
        cs[0]=fg[0]/d; cs[1]=fg[1]/d;
      }
    r=1.0;
    if (fabs(fg[0])>fabs(fg[1])) r=cs[1];
    else
      if (cs[0]!=0.0) r=1.0/cs[0];
    fg[0]=d; fg[1]=r;
    return;
  }


// 矩阵求逆
bool MtxInv( double *A, double *AI,  int n) // A(n*n) -> AI(n*n)
{
	int t, i, j, k, u, v, n2 = n * n;
	int *is = new int[n];
	int *js = new int[n];
	double d, p, **MatIO = new double*[n];
	for(t = 0; t < n2; t++)	AI[t] = A[t];
	for (MatIO[0] = AI, t = 1; t < n; t++)	MatIO[t] = MatIO[t - 1] + n;
	for (k = 0; k < n; k++)
	{
		d = 0;
        for (i = k; i < n; i++)
			for (j = k; j < n; j++)
			{
				p = fabs(MatIO[i][j]);
				if (p > d)
				{
					d = p;	is[k] = i;	js[k] = j;
				}
			}
        if (fabs(d) < 1E-15)	break;
        if (is[k] != k)
			for (j = 0; j < n; j++)
			{
				v = is[k];	p = MatIO[k][j];
				if( (v >= n) || (v < 0) )	
				{
					delete[] is;
					delete[] js;
					delete[] MatIO;	
					return false;
				}
				MatIO[k][j]	= MatIO[v][j];	MatIO[v][j]	= p;
			}
        if (js[k] != k)
			for (i = 0; i < n; i++)
			{
				u = js[k];	p = MatIO[i][k];
				if( (u >= n) || (u < 0) )
				{
					delete[] is;
					delete[] js;
					delete[] MatIO;	
					return false;
				}
				MatIO[i][k]	= MatIO[i][u];	MatIO[i][u] = p;
			}
		if(fabs(MatIO[k][k]) > 0)
			MatIO[k][k] = 1 / MatIO[k][k];
		else
			MatIO[k][k] = 999999.0;
        for (j = 0; j < n; j++)	
			if (j != k) 
				MatIO[k][j] *= MatIO[k][k];
		for (i = 0; i < n; i++)
			if (i != k) 
			{
				for (j = 0; j < n; j++) 
					if (j != k) 
						MatIO[i][j] -= MatIO[i][k] * MatIO[k][j];
			}
        for (i = 0; i < n; i++)	
			if (i != k) 
				MatIO[i][k] *= -MatIO[k][k];
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j < n; j++)
			{
				v = js[k];	p = MatIO[k][j];
				if( (v >= n) || (v < 0) )
				{
					delete[] is;
					delete[] js;
					delete[] MatIO;	
					return false;
				}
				MatIO[k][j]	= MatIO[v][j];	MatIO[v][j]	= p;
			}
		if (is[k] != k)
			for (i = 0; i < n; i++)
			{
				u = is[k];	p = MatIO[i][k];
				if( (u >= n) || (u < 0) )
				{
					delete[] is;
					delete[] js;
					delete[] MatIO;	
					return false;
				}
				MatIO[i][k]	= MatIO[i][u];	MatIO[i][u]	= p;
			}
		if(k == 0)
			break;
	}
	delete[] is;
	delete[] js;
	delete[] MatIO;	
	return true;
}

// 矩阵相乘(按行储存)
void MtxMlt( double *A,  double *B, double *AB, int m,  int n,  int s) // A(m*n) * B(n*s) = AB(m*s)
{
	int i, j, l, k, is, in;
    for (i = 0; i < m; i++)
	{
		is = i * s;	in = i * n;
		for (j = 0; j < s; j++)
		{ 
			k = is + j;	AB[k] = 0.0;
			for (l = 0; l < n; l++)
                AB[k] += A[in + l] * B[l * s + j];
		}
	}

}

// 矩阵转置
void MtxTrs(double *A, double *AT, int m, int n) // A(m*n) -> AT(n*m)
{
	int i, j, jn;
	for (j = 0; j < n; j++)
	{
		jn = j * m;
		for (i = 0; i < m; i++)	AT[jn + i] = A[i * n + j];
	}
}

// 最小二乘解线性方程组 // 成功返回true，失败返回false // M(m*n) * X(n*1) = B(m*1) -> X
bool LS( double *M, double *X, double *B, int m, int n)
{
	int nm = n * m, nn = n * n;
	double *MT = new double[nm];
	double *MT_M = new double[nn];
	double *MT_M_A = new double[nn];
	double *MT_M_A_MT = new double[nm];
	MtxTrs(M, MT, n, m);
	MtxMlt(MT, M, MT_M, n, m, n);
	if(!MtxInv(MT_M, MT_M_A, n) )
	{
		delete[] MT;
		delete[] MT_M;
		delete[] MT_M_A;
		delete[] MT_M_A_MT;	
		return false;
	}
	MtxMlt(MT_M_A, MT, MT_M_A_MT, n, n, m);
	MtxMlt(MT_M_A_MT, B, X, n, m, 1);
	delete[] MT;
	delete[] MT_M;
	delete[] MT_M_A;
	delete[] MT_M_A_MT;
	return true;
}

int agmiv(double *a,int m,int n,double *b,double *x,double *aa,double eps,double *u,double *v,int ka)
{ 
	int i,j;
    i=bginv(a,m,n,aa,eps,u,v,ka);
    if (i<0) return(-1);
    for (i=0; i<=n-1; i++)
	{ x[i]=0.0;
	for (j=0; j<=m-1; j++)
		x[i]=x[i]+aa[i*m+j]*b[j];
	}
    return(1);
}

int bginv(double *a,int m,int n,double *aa,double eps,double *u,double *v,int ka)
{ 
	int i,j,k,l,t,p,q,f;
    i=bmuav(a,m,n,u,v,eps,ka);
    if (i<0) return(-1);
    j=n;
    if (m<n) j=m;
    j=j-1;
    k=0;
    while ((k<=j)&&(a[k*n+k]!=0.0)) k=k+1;
    k=k-1;
    for (i=0; i<=n-1; i++)
		for (j=0; j<=m-1; j++)
		{ t=i*m+j; aa[t]=0.0;
        for (l=0; l<=k; l++)
		{ f=l*n+i; p=j*m+l; q=l*n+l;
		aa[t]=aa[t]+v[f]*u[p]/a[q];
		}
		}
	return(1);
}
void TreeMtxMlt(double *A, double *B, double *C, double *RE, int m, int n, int l, int s )//A(m*n) *B(n*l) * C(l*s)
{
	double *AB=new double[m*l];
	MtxMlt(A,B,AB,m,n,l);
	MtxMlt(AB,C,RE,m,l,s);
	delete[]AB;
	return;

}

//矩阵化为对称三对角矩阵
void cstrq(double *a,int n,double *q,double *b,double *c)
{
	int i,j,k,u;
    double h,f,g,h2;
    for (i=0; i<=n-1; i++)
    for (j=0; j<=n-1; j++)
      { u=i*n+j; q[u]=a[u];}
    for (i=n-1; i>=1; i--)
      { h=0.0;
        if (i>1)
          for (k=0; k<=i-1; k++)
            { u=i*n+k; h=h+q[u]*q[u];}
        if (h+1.0==1.0)
          { c[i]=0.0;
            if (i==1) c[i]=q[i*n+i-1];
            b[i]=0.0;
          }
        else
          { c[i]=sqrt(h);
            u=i*n+i-1;
            if (q[u]>0.0) c[i]=-c[i];
            h=h-q[u]*c[i];
            q[u]=q[u]-c[i];
            f=0.0;
            for (j=0; j<=i-1; j++)
              { q[j*n+i]=q[i*n+j]/h;
                g=0.0;
                for (k=0; k<=j; k++)
                  g=g+q[j*n+k]*q[i*n+k];
                if (j+1<=i-1)
                  for (k=j+1; k<=i-1; k++)
                    g=g+q[k*n+j]*q[i*n+k];
                c[j]=g/h;
                f=f+g*q[j*n+i];
              }
            h2=f/(h+h);
            for (j=0; j<=i-1; j++)
              { f=q[i*n+j];
                g=c[j]-h2*f;
                c[j]=g;
                for (k=0; k<=j; k++)
                  { u=j*n+k;
                    q[u]=q[u]-f*c[k]-g*q[i*n+k];
                  }
              }
            b[i]=h;
          }
      }
    for (i=0; i<=n-2; i++) c[i]=c[i+1];
    c[n-1]=0.0;
    b[0]=0.0;
    for (i=0; i<=n-1; i++)
      { if ((b[i]!=0.0)&&(i-1>=0))
          for (j=0; j<=i-1; j++)
            { g=0.0;
              for (k=0; k<=i-1; k++)
                g=g+q[i*n+k]*q[k*n+j];
              for (k=0; k<=i-1; k++)
                { u=k*n+j;
                  q[u]=q[u]-g*q[k*n+i];
                }
            }
        u=i*n+i;
        b[i]=q[u]; q[u]=1.0;
        if (i-1>=0)
          for (j=0; j<=i-1; j++)
            { q[i*n+j]=0.0; q[j*n+i]=0.0;}
      }
    return;
}

//求对称三对角阵的特征值和特征向量
int csstq(int n,double *b,double *c,double *q,double eps,int l)
{
	int i,j,k,m,it,u,v;
    double d,f,h,g,p,r,e,s;
    c[n-1]=0.0; d=0.0; f=0.0;
    for (j=0; j<=n-1; j++)
      { it=0;
        h=eps*(fabs(b[j])+fabs(c[j]));
        if (h>d) d=h;
        m=j;
        while ((m<=n-1)&&(fabs(c[m])>d)) m=m+1;
        if (m!=j)
          { do
              { if (it==l)
                  { printf("fail\n");
                    return(-1);
                  }
                it=it+1;
                g=b[j];
                p=(b[j+1]-g)/(2.0*c[j]);
                r=sqrt(p*p+1.0);
                if (p>=0.0) b[j]=c[j]/(p+r);
                else b[j]=c[j]/(p-r);
                h=g-b[j];
                for (i=j+1; i<=n-1; i++)
                  b[i]=b[i]-h;
                f=f+h; p=b[m]; e=1.0; s=0.0;
                for (i=m-1; i>=j; i--)
                  { g=e*c[i]; h=e*p;
                    if (fabs(p)>=fabs(c[i]))
                      { e=c[i]/p; r=sqrt(e*e+1.0);
                        c[i+1]=s*p*r; s=e/r; e=1.0/r;
                      }
                    else
		      { e=p/c[i]; r=sqrt(e*e+1.0);
                        c[i+1]=s*c[i]*r;
                        s=1.0/r; e=e/r;
                      }
                    p=e*b[i]-s*g;
                    b[i+1]=h+s*(e*g+s*b[i]);
                    for (k=0; k<=n-1; k++)
                      { u=k*n+i+1; v=u-1;
                        h=q[u]; q[u]=s*q[v]+e*h;
                        q[v]=e*q[v]-s*h;
                      }
                  }
                c[j]=s*p; b[j]=e*p;
              }
            while (fabs(c[j])>d);
          }
        b[j]=b[j]+f;
      }
    for (i=0; i<=n-1; i++)
      { k=i; p=b[i];
        if (i+1<=n-1)
          { j=i+1;
            while ((j<=n-1)&&(b[j]<=p))
              { k=j; p=b[j]; j=j+1;}
          }
        if (k!=i)
          { b[k]=b[i]; b[i]=p;
            for (j=0; j<=n-1; j++)
              { u=j*n+i; v=j*n+k;
                p=q[u]; q[u]=q[v]; q[v]=p;
              }
          }
      }
    return(1);
}

int aldle(double *a,int n,int m,double *c)
{
	int i,j,l,k,u,v,w,k1,k2,k3;
    double p;
    if (fabs(a[0])+1.0==1.0)
      { printf("fail\n"); return(-2);}
    for (i=1; i<=n-1; i++)
      { u=i*n; a[u]=a[u]/a[0];}
    for (i=1; i<=n-2; i++)
      { u=i*n+i;
        for (j=1; j<=i; j++)
          { v=i*n+j-1; l=(j-1)*n+j-1;
            a[u]=a[u]-a[v]*a[v]*a[l];
          }
        p=a[u];
        if (fabs(p)+1.0==1.0)
          { printf("fail\n"); return(-2);}
        for (k=i+1; k<=n-1; k++)
          { u=k*n+i;
            for (j=1; j<=i; j++)
              { v=k*n+j-1; l=i*n+j-1; w=(j-1)*n+j-1;
		a[u]=a[u]-a[v]*a[l]*a[w];
              }
            a[u]=a[u]/p;
          }
      }
    u=n*n-1;
    for (j=1; j<=n-1; j++)
      { v=(n-1)*n+j-1; w=(j-1)*n+j-1;
        a[u]=a[u]-a[v]*a[v]*a[w];
      }
    p=a[u];
    if (fabs(p)+1.0==1.0)
      { printf("fail\n"); return(-2);}
    for (j=0; j<=m-1; j++)
    for (i=1; i<=n-1; i++)
      { u=i*m+j;
        for (k=1; k<=i; k++)
          { v=i*n+k-1; w=(k-1)*m+j;
            c[u]=c[u]-a[v]*c[w];
          }
      }
    for (i=1; i<=n-1; i++)
      { u=(i-1)*n+i-1;
        for (j=i; j<=n-1; j++)
          { v=(i-1)*n+j; w=j*n+i-1;
            a[v]=a[u]*a[w];
          }
      }
    for (j=0; j<=m-1; j++)
      { u=(n-1)*m+j;
        c[u]=c[u]/p;
        for (k=1; k<=n-1; k++)
          { k1=n-k; k3=k1-1; u=k3*m+j;
            for (k2=k1; k2<=n-1; k2++)
              { v=k3*n+k2; w=k2*m+j;
                c[u]=c[u]-a[v]*c[w];
              }
            c[u]=c[u]/a[k3*n+k3];
          }
      }
    return(2);
}

  int mrab1(int a,int b,int *r)
  { int k,l,m,i,p;
    k=b-a+1; l=2;
    while (l<k) l=l+l;
    m=4*l; k=*r; i=1;
    while (i<=1)
      { k=k+k+k+k+k;
	k=k%m; l=k/4+a;
        if (l<=b) { p=l; i=i+1;}
      }
    *r=k;
    return(p);
  }

  void mrabs(int a,int b,int *r,int *p,int n)
  { int k,l,m,i;
    k=b-a+1; l=2;
    while (l<k) l=l+l;
    m=4*l; k=*r; i=0;
    while (i<=n-1)
      { k=k+k+k+k+k;
	k=k%m; l=k/4+a;
        if (l<=b) { p[i]=l; i=i+1;}
      }
    *r=k;
    return;
  }


 int brank(double *a,int m,int n)
 { 
	  int i,j,k,nn,is,js,l,ll,u,v;
    double q,d;
    nn=m;
    if (m>=n) nn=n;
    k=0;
    for (l=0; l<=nn-1; l++)
      { q=0.0;
        for (i=l; i<=m-1; i++)
        for (j=l; j<=n-1; j++)
          { ll=i*n+j; d=fabs(a[ll]);
	    if (d>q) { q=d; is=i; js=j;}
          }
        if (q+1.0==1.0) return(k);
        k=k+1;
        if (is!=l)
          { for (j=l; j<=n-1; j++)
              { u=l*n+j; v=is*n+j;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        if (js!=l)
          { for (i=l; i<=m-1; i++)
              { u=i*n+js; v=i*n+l;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        ll=l*n+l;
        for (i=l+1; i<=n-1; i++)
          { d=a[i*n+l]/a[ll];
            for (j=l+1; j<=n-1; j++)
              { u=i*n+j;
                a[u]=a[u]-d*a[l*n+j];
              }
          }
      }
    return(k);
 }

  double sqar(int n,double *a)//求向量的模平方
  {
	  int i;
	  double sqr=0;
	  for(i=0;i<n;i++)
	  {
		  sqr+=(a[i]*a[i]);

	  }
	  return sqr;
  }

   void MatrixAdd(double *A, double *B, double *SUM, int m, int n)//两矩阵相加
   {
	   int i,j;
	   for(i=0;i<m;i++)
	   {
		   for(j=0;j<n;j++)
		   {
			   SUM[i*m+j]=A[i*m+j]+B[i*m+j];
		   }
	   }
   }

    void MatrixMatNum(double *A, double a,int m,int n)//矩阵数乘
	{
		int i,j;
		for(i=0;i<m;i++)
		{
			for(j=0;j<n;j++)
			{
				A[i*m+j]*=a;
			}
		}
	}

	void VectorMat(double *a, double *A)
	{
		A[0]=A[4]=A[8]=0;
		A[1]=-a[2];
		A[2]=a[1];
		A[3]=a[2];
		A[5]=-a[0];
		A[6]=-a[1];
		A[7]=a[0];
	}
	 void VectorMatilply(double *a, double *b, double *c)//两三维向量叉乘
	 {
		 double M[9];
		 VectorMat(a,M);
		 MtxMlt(M,b,c,3,3,1);

	 }


 double  VectorPoinMut(double *a, double *b,int n)//n维向量内积
 {
	 int i;
	 double mut=0;
	 for(i=0;i<n;i++)
	 {
		 mut+=(a[i]*b[i]);
	 }
	 return mut;
 }

 void hpir1(double *x,double *y,int n,double *a,int m,double *dt)
  { int i,j,k;
    double z,p,c,g,q,d1,d2,s[20],t[20],b[20];
    for (i=0; i<=m-1; i++) a[i]=0.0;
    if (m>n) m=n;
    if (m>20) m=20;
    z=0.0;
    for (i=0; i<=n-1; i++) z=z+x[i]/(1.0*n);
    b[0]=1.0; d1=1.0*n; p=0.0; c=0.0;
    for (i=0; i<=n-1; i++)
      { p=p+(x[i]-z); c=c+y[i];}
    c=c/d1; p=p/d1;
    a[0]=c*b[0];
    if (m>1)
      { t[1]=1.0; t[0]=-p;
        d2=0.0; c=0.0; g=0.0;
        for (i=0; i<=n-1; i++)
          { q=x[i]-z-p; d2=d2+q*q;
            c=c+y[i]*q;
            g=g+(x[i]-z)*q*q;
          }
        c=c/d2; p=g/d2; q=d2/d1;
        d1=d2;
        a[1]=c*t[1]; a[0]=c*t[0]+a[0];
      }
    for (j=2; j<=m-1; j++)
      { s[j]=t[j-1];
        s[j-1]=-p*t[j-1]+t[j-2];
        if (j>=3)
          for (k=j-2; k>=1; k--)
            s[k]=-p*t[k]+t[k-1]-q*b[k];
        s[0]=-p*t[0]-q*b[0];
        d2=0.0; c=0.0; g=0.0;
        for (i=0; i<=n-1; i++)
          { q=s[j];
            for (k=j-1; k>=0; k--)
              q=q*(x[i]-z)+s[k];
            d2=d2+q*q; c=c+y[i]*q;
            g=g+(x[i]-z)*q*q;
          }
        c=c/d2; p=g/d2; q=d2/d1;
        d1=d2;
        a[j]=c*s[j]; t[j]=s[j];
        for (k=j-1; k>=0; k--)
          { a[k]=c*s[k]+a[k];
            b[k]=t[k]; t[k]=s[k];
          }
      }
    dt[0]=0.0; dt[1]=0.0; dt[2]=0.0;
    for (i=0; i<=n-1; i++)
      { q=a[m-1];
        for (k=m-2; k>=0; k--)
          q=a[k]+q*(x[i]-z);
        p=q-y[i];
        if (fabs(p)>dt[2]) dt[2]=fabs(p);
        dt[0]=dt[0]+p*p;
        dt[1]=dt[1]+fabs(p);
      }
    return;
  }

 double mgrn1(double u, double g,double *r)
 {
	int i,m;
    double s,w,v,t;
    s=65536.0; w=2053.0; v=13849.0;
    t=0.0;
    for (i=1; i<=12; i++)
      { *r=(*r)*w+v; m=(int)(*r/s);
        *r=*r-m*s; t=t+(*r)/s;
      }
    t=u+g*(t-6.0);
    return(t);
  }
 #define TUKEY_C 4.6851
 double tukey_weight(double x)
 {
	 if (fabs(x)>=TUKEY_C)
	 {
		 return 0;
	 }
	 x = x/TUKEY_C;
	 return (1-x*x)*(1-x*x);
 }

 static int weight_cmp( const void* p1, const void* p2 )
 {
	 double* _p1 = (double*)p1;
	 double* _p2 = (double*)p2;

	 if( (*_p1) > (*_p2) )
		 return -1;
	 if( (*_p1) < (*_p2) )
		 return 1;
	 return 0;
 }
 //获得中值
#include <string>
 double get_median(double* samples, int nsize)
 {
	 double med_weight = 10;
	 double* weights = new double[nsize];
	 memset(weights, 0, sizeof(double)*nsize);
	 memcpy(weights, samples, sizeof(double)*nsize);
	 qsort(weights, nsize, sizeof(double),&weight_cmp);
	 //获得中值
	 med_weight = (weights[nsize/2]+weights[nsize/2-1])/2.0;
	 delete[] weights;
	 return med_weight;
 }

 void get_weights(double* samples, int nsize)
 {
	 double med_dis     = 10;                              //残差中值
	 double med_detar   = 10;                             //detar中值
	 double detar       = 10;
	 double* detars    = new double[nsize];
	 memset(detars,    0, sizeof(double)*nsize);
	 med_dis = get_median(samples, nsize);
	 //距离残差使用中值进行归一化
	 for (int k2=0; k2<nsize; k2++)
	 {
		 detars[k2] = samples[k2]-med_dis;
	 }
	 //获取归一化距离中值
	 med_detar = get_median(detars, nsize);
	 //计算偏差绝对值中值
	 for (int k2=0; k2<nsize; k2++)
	 {
		 detars[k2] = fabs(detars[k2]-med_detar);
	 }
	 //计算数据sigmar值
	 med_detar = get_median(detars, nsize);
	 detar     = med_detar/1.48;
	 for (int k2=0; k2<nsize; k2++)
	 {
		 //samples[k2] = tukey_weight((samples[k2]-med_dis)/(detar+1e-5));    //小bug，2014-8-26
		 //不进行加权
		 samples[k2] = 1;
	 }
	 delete[] detars;
 }

 void get_weights_huber(double* samples, int nsize, double scale)
 {
	 for (int k2=0; k2<nsize; k2++)
	 {
		 samples[k2] = 1/(1+fabs(samples[k2])/scale);
	 }
 }