#include "stdafx.h"
#include "SaliantObjectDetectRef.h"
#include <queue>

#define DIS_CV
#include <opencv2/opencv.hpp>
#include <opencv2/legacy/compat.hpp>
#include "ransac_ellipse.h"
using namespace cv;
using namespace std;



bool check_in_image_range(int x, int y, int w, int h, int r = 0)
{
	if (x-r<0||x>w-1-r)
	{
		return false;
	}
	if (y-r<0||y>h-1-r)
	{
		return false;
	}
	return true;
}

//��Ե�������Բ���

box2d getBox2d(double* ep)
{
	box2d box;
	box.center.x = ep[0];
	box.center.y = ep[1];
	box.size.x   = ep[2]*2;
	box.size.y   = ep[3]*2;
	box.angle    = ep[4];
	return box;
}

//�����������
bool getNextPoint8p(unsigned char* imgcanny, int w, int h, int* rnow,int* cnow)
{
	int roff[9]={1 , 0 , -1 , 0 , 1 , 1 , -1 , -1, 0};   //x����ƫ��,����ͨ��
	int coff[9]={0 , 1 , 0 , -1 , 1 , -1 , -1 , 1, 0};   //y����ƫ��
	int r[9];
	int c[9];
	//������������
	for (int i=0;i<9;i++)
	{
		r[i]=*rnow+roff[i];
		c[i]=*cnow+coff[i];
	}
	int ff[8] = {0};
	int n = 0;
	for (int i=0;i<8;i++)
	{
		*rnow=r[i];
		*cnow=c[i];
		if ((*rnow>=0)&&(*rnow<w)&&(*cnow>=0)&&(*cnow<h))
		{
			if (imgcanny[(*cnow)*w+*rnow]!=0)
			{
				n++;
				ff[i] = 1;
			}
		}
	}
	//if ((n==0)||(n>2) )
	{
		////*rnow=r[8];
		//*cnow=c[8];
		//return false;
	}
	//else
	{
		for (int i=0; i<8; i++)
		{
			if (ff[i])
			{
				*rnow=r[i];
				*cnow=c[i];
				return true;
			}
		}
	}
	return false;
}

bool getNextPoint8p2(unsigned char* imgcanny, int w, int h, int* rnow,int* cnow, unsigned char* img1b)
{
	int roff[9]={1 , 0 , -1 , 0 , 1 , 1 , -1 , -1, 0};   //x����ƫ��,����ͨ��
	int coff[9]={0 , 1 , 0 , -1 , 1 , -1 , -1 , 1, 0};   //y����ƫ��
	int r[9];
	int c[9];
	//������������
	for (int i=0;i<9;i++)
	{
		r[i]=*rnow+roff[i];
		c[i]=*cnow+coff[i];
	}
	int ff[8] = {0};
	int n = 0;
	for (int i=0;i<8;i++)
	{
		*rnow=r[i];
		*cnow=c[i];
		if ((*rnow>=0)&&(*rnow<w)&&(*cnow>=0)&&(*cnow<h))
		{
			if ((imgcanny[(*cnow)*w+*rnow]!=0)&&(img1b[(*cnow)*w+*rnow]!=0))
			{
				n++;
				ff[i] = 1;
			}
		}
	}
	//if ((n==0)||(n>2) )
	{
		////*rnow=r[8];
		//*cnow=c[8];
		//return false;
	}
	//else
	{
		for (int i=0; i<8; i++)
		{
			if (ff[i])
			{
				*rnow=r[i];
				*cnow=c[i];
				return true;
			}
		}
	}
	return false;
}

//����ͨ��
bool getNextPoint4p(unsigned char* imgcanny, int w, int h, int* rnow,int* cnow)
{
	int roff[5]={1 , 0 , -1 , 0 ,  0};   //x����ƫ��,����ͨ��
	int coff[5]={0 , 1 , 0 , -1 ,  0};   //y����ƫ��
	int r[5];
	int c[5];
	//������������
	for (int i=0;i<5;i++)
	{
		r[i]=*rnow+roff[i];
		c[i]=*cnow+coff[i];
	}
	for (int i=0;i<4;i++)
	{
		*rnow=r[i];
		*cnow=c[i];
		if ((*rnow>=0)&&(*rnow<w)&&(*cnow>=0)&&(*cnow<h))
		{
			if (imgcanny[(*cnow)*w+*rnow]!=0)
				return true;
		}
	}
	*rnow=r[4];
	*cnow=c[4];
	return false;
}

void CalcCirclePara2(point* pts,  int nsize, CvBox2D& box)
{
	CvPoint2D32f* edge= new CvPoint2D32f[nsize];
	//ͼ������ϵתΪ�ѿ�������ϵ
	for (int i = 0; i<nsize; i++)
	{
		edge[i]=cvPoint2D32f(pts[i].x, pts[i].y);
	}
	cvFitEllipse(edge, nsize, &box);
	delete[] edge;
}

void CalcCirclePara2(vector<CvPoint2D32f> pts,  CvBox2D& box)
{
	CvPoint2D32f* edge= new CvPoint2D32f[pts.size()];
	//ͼ������ϵתΪ�ѿ�������ϵ
	for (int i = 0; i<pts.size(); i++)
	{
		edge[i]=pts[i];
	}
	cvFitEllipse(edge, pts.size(), &box);
	delete[] edge;
}


float dis(point p0, point p1)
{
	point p;
	p.x = p0.x-p1.x;
	p.y = p1.y-p1.y;
	return p.x*p.x+p.y*p.y;
}

void trackedge(unsigned char* imgCanny, int w, int h,int rstart,int cstart,int minedgelength, point** edgeSegment, int &nsize)
{
	int r=rstart;
	int c=cstart;
	point pt;
	pt.x = rstart;
	pt.y = cstart;
	std::vector<point> pts;
	pts.push_back(pt);
	bool thereIsAPoint=false;
	int ang = 0;
	thereIsAPoint=getNextPoint8p(imgCanny,w, h, &rstart,&cstart);
	//������,��һ������

	while (thereIsAPoint)
	{
		pt.x = rstart;
		pt.y = cstart;
		pts.push_back(pt);
		imgCanny[cstart*w+rstart]=0;
		thereIsAPoint=getNextPoint8p(imgCanny,w, h, &rstart,&cstart);
	}
	//������
	rstart=r;
	cstart=c;
	thereIsAPoint=getNextPoint8p(imgCanny, w, h, &rstart,&cstart);
	while (thereIsAPoint)
	{
		pt.x = rstart;
		pt.y = cstart;
		pts.insert(pts.begin(),pt);
		imgCanny[cstart*w+rstart]=0;
		thereIsAPoint=getNextPoint8p(imgCanny,w, h, &rstart,&cstart);
	}
	if (pts.size()<minedgelength)
	{
		for (int i=0;i<pts.size();i++)
		{ 
			pt = pts[i];
			r  = pt.x;
			c  = pt.y;
			imgCanny[c*w+r] = 255;
		}
		pts.clear();
	}
	else
	{
		nsize = pts.size();
		*edgeSegment = new point[nsize];
		for (int i=0;i<pts.size();i++)
		{ 
			(*edgeSegment)[i] = pts[i];
		}
		pts.clear();
	}
}

void trackedge2(unsigned char* imgCanny, int w, int h,int rstart,int cstart,int minedgelength, point** edgeSegment, int &nsize, unsigned char* img1b)
{
	int r=rstart;
	int c=cstart;
	point pt;
	pt.x = rstart;
	pt.y = cstart;
	std::vector<point> pts;
	pts.push_back(pt);
	bool thereIsAPoint=false;
	int ang = 0;
	thereIsAPoint=getNextPoint8p2(imgCanny,w, h, &rstart,&cstart, img1b);
	//������,��һ������

	while (thereIsAPoint)
	{
		pt.x = rstart;
		pt.y = cstart;
		pts.push_back(pt);
		imgCanny[cstart*w+rstart]=0;
		thereIsAPoint=getNextPoint8p2(imgCanny,w, h, &rstart,&cstart, img1b);
		if (pts.size()>40)
		{
			break;
		}
	}
	//������
	rstart=r;
	cstart=c;
	thereIsAPoint=getNextPoint8p2(imgCanny, w, h, &rstart,&cstart, img1b);
	while (thereIsAPoint)
	{
		pt.x = rstart;
		pt.y = cstart;
		pts.insert(pts.begin(),pt);
		imgCanny[cstart*w+rstart]=0;
		thereIsAPoint=getNextPoint8p2(imgCanny,w, h, &rstart,&cstart, img1b);
		if (pts.size()>80)
		{
			break;
		}
	}
	if (pts.size()<minedgelength)
	{
		for (int i=0;i<pts.size();i++)
		{ 
			pt = pts[i];
			r  = pt.x;
			c  = pt.y;
			imgCanny[c*w+r] = 255;
		}
		pts.clear();
	}
	else
	{
		nsize = pts.size();
		*edgeSegment = new point[nsize];
		for (int i=0;i<pts.size();i++)
		{ 
			(*edgeSegment)[i] = pts[i];
		}
		pts.clear();
	}
}

//����������ֵ��Եͼ�õ���Ե������
point** getedge(IplImage* img, int minedgelength,int low_threshole,int high_threshild, int** ess, int *ne)
{
	int w = img->width;
	int h = img->height;
	unsigned char* imgCanny = 0;
	unsigned char* imgMask = 0;
	unsigned int grayNum=0;
	std::vector<point*> edges;
	std::vector<int> ns;
	IplImage* imgcv = cvCreateImage(cvSize(w, h), 8, 1);
	IplImage* img1b = cvCreateImage(cvSize(w, h), 8, 1);
	cvCanny(img, imgcv, low_threshole, high_threshild);
	cvThreshold(img, img1b, 200, 255, CV_THRESH_BINARY_INV);
	cvErode(img1b, img1b);
#ifdef DIS_CV
	//cvShowImage("img1b", img1b);
#endif
	imgCanny = (unsigned char*)imgcv->imageData;
    imgMask  = (unsigned char*)img1b->imageData;
	for (int j=0;j<h;j++)
	{
		for (int i=0;i<w;i++)
		{
			grayNum=imgCanny[j*w+i];
			if (grayNum!=0)
			{
				point* edgeSegment = 0;
				int   nsize = 0;
				trackedge2(imgCanny, w, h, i,j, minedgelength, &edgeSegment, nsize, imgMask);
				if (nsize > minedgelength)
				{
					if(dis(edgeSegment[0], edgeSegment[nsize/2])<7)
						continue;
					edges.push_back(edgeSegment);
					ns.push_back(nsize);
				}
			}
		}
	}
	cvReleaseImage(&imgcv);
	cvReleaseImage(&img1b);
	if (ns.size())
	{
		point** img_edges = new point*[ns.size()];
		*ess  = new int[ns.size()];
		*ne   = ns.size();
		for (int i=0; i<ns.size(); i++)
		{
			img_edges[i] = edges[i];
			(*ess)[i]      = ns[i];
		}
		return img_edges;
	}
	return 0;
}  
//�ֲ���Բ���
double* detect_ellipse_local_cv(IplImage* img, CvRect roi)
{
    double*  elines = NULL;
	//ȫͼ������Բ���
	CvBox2D box0;
	int w = img->width;
	int h = img->height;
	//roi
	int ww = roi.width;
	int hh = roi.height;
	IplImage* imgroi = cvCreateImage(cvSize(ww, hh), 8, 1);
	cvSetImageROI(img, roi);
	cvCopy(img, imgroi);
	cvResetImageROI(img);
	int s1 = 2.0;
	int s2 = s1*s1;
	int www = ww/s1;
	int hhh = hh/s1;
	IplImage* img1u0   = cvCreateImage(cvSize(www, hhh), 8, 1);;
	IplImage* img1u1   = cvCreateImage(cvSize(www, hhh), 8, 1);;
	cvResize(imgroi, img1u0);    
	//opencv,��ֵ�˲�
	cvSmooth(img1u0, img1u1, CV_MEDIAN, 15/s1, 15/s1);
	int npara1 = 30,  npara2 = 50;                             //canny��Ե������
	int maxr   = 15000/s2, minr = 2000/s2;                    //���������ֵ������ͫ�׳��������������            
	int* nes = 0;
	int  ns  = 0;
	point** edges = getedge(img1u1, 100/s1, npara1, npara2, &nes, &ns);   //�б�Եͼ��õ���Ե����
	if (edges)
	{
		int* nboxs = new int[ns];
		int ee = 0;
		int ntotal = 0;
		for (int e=0; e<ns; e++)
		{
			//opencv ��Բ���
			CalcCirclePara2(edges[e], nes[e], box0);                                       //��Ե�����Բ��
			double area = box0.size.width*box0.size.height;
			if ((area<maxr)&&(area>minr))
			{
				if (fabs(max(box0.size.height, box0.size.width)/min(box0.size.height, box0.size.width)-1.0)<2.0)   //��Բ��������Բ�����
				{
					nboxs[ee++] = e;                     //�洢����Ҫ������еı�Ե����
					ntotal += nes[e];
					//break;
				}
			}	
		}
		int eee = ee;
		if (eee)
		{
			ee = 0;
			point* epts = new point[ntotal];
			for (int e=0; e<eee; e++)
			{
				//opencv ��Բ���
				for (int i=0 ;i<nes[nboxs[e]]; i++)
				{
					epts[ee++] = edges[nboxs[e]][i];  
				}
			}
			if (ntotal)
			{
				CalcCirclePara2(epts, ntotal, box0);                                        //ͫ�ױ�Ե���������Բ
				box0.center.x = box0.center.x*s1;
				box0.center.y = box0.center.y*s1;
				box0.size.width   = box0.size.width*s1;
				box0.size.height  = box0.size.height*s1;
				elines = new double[5];
				elines[0] = box0.center.x + roi.x;
				elines[1] = box0.center.y + roi.y;
				elines[2] = box0.size.width/2.0;
				elines[3] = box0.size.height/2.0;
				elines[4] = box0.angle*3.14159/180.0;
			}
			delete[] epts;
		}
	    //�ͷ��ڴ�
		for (int i=0;i<ns;i++)
		{
			delete[] edges[i];
		}
		delete[] edges;
		delete[] nes;
	}
	cvReleaseImage(&imgroi);
	cvReleaseImage(&img1u0);
	cvReleaseImage(&img1u1);
	return elines;
}
float getRegionMeanValue(IplImage* img, CvPoint pt, int nR)
{
	// Mean values
	float meanvalue= 0;
	int   count=0;
	for (int i=-nR; i<=nR; i++){
		for(int j=-nR; j<=nR; j++)
		{
			meanvalue += ((uchar*)img->imageData+(j+pt.y)*img->widthStep)[i+pt.x]; //��
			count++;
		}
	}
	meanvalue/=count;
	return meanvalue;
}

void RegionGrow(IplImage* org, IplImage* dst, point seedPt, int nThreshold)
{
	int nDx[]={-1,0,1,0};
	int nDy[]={0,1,0,-1};
	unsigned char* pimg = (unsigned char*)org->imageData;
	unsigned char* pout = (unsigned char*)dst->imageData;
	int w = org->width;
	int h = org->height;
	int nWidth             = w;
	int nHeight            = h;
	// ͼ�����ڴ���ÿһ������ռ�õ�ʵ�ʿռ�
	int nSaveWidth = w;
	memset(pout, 0, w*h*sizeof(unsigned char));
	// ���ӵ�
	int nSeedX, nSeedY;
	// �������ӵ�Ϊͼ�������
	nSeedX = seedPt.x;
	nSeedY = seedPt.y;
	// �����ջ���洢����
	int * pnGrowQueX ;
	int * pnGrowQueY ;
	// ����ռ�
	pnGrowQueX = new int [nWidth*nHeight];
	pnGrowQueY = new int [nWidth*nHeight];
	// ͼ�����ݵ�ָ��
	unsigned char *  pUnchInput  = pimg;
	unsigned char *  pUnchOutput = pout;
	// �����ջ�������յ�
	// ��nStart=nEnd, ��ʾ��ջ��ֻ��һ����
	int nStart ;
	int nEnd   ;
	//��ʼ��
	nStart = 0 ;
	nEnd   = 0 ;
	// �����ӵ������ѹ��ջ
	pnGrowQueX[nEnd] = nSeedX;
	pnGrowQueY[nEnd] = nSeedY;
	// ��ǰ���ڴ��������
	int nCurrX ;
	int nCurrY ;
	// ѭ�����Ʊ���
	int k ;
	// ͼ��ĺ�������,�����Ե�ǰ���ص�������б���
	int xx;
	int yy;
	int v0 = pUnchInput[nSeedY*nWidth+nSeedX];
	while (nStart<=nEnd)

	{
		// ��ǰ���ӵ������
		nCurrX = pnGrowQueX[nStart];
		nCurrY = pnGrowQueY[nStart];
		// �Ե�ǰ���������б���
		for (k=0; k<4; k++)   
		{   
			// 4�������ص�����
			xx = nCurrX + nDx[k];
			yy = nCurrY + nDy[k];
			// �ж�����(xx��yy) �Ƿ���ͼ���ڲ�
			// �ж�����(xx��yy) �Ƿ��Ѿ������
			// pUnRegion[yy*nWidth+xx]==0 ��ʾ��û�д���
			// �����������ж�����(xx��yy)�͵�ǰ����(nCurrX,nCurrY) ����ֵ��ľ���ֵ
			int nf = pUnchOutput[yy*nWidth+xx];
			//�ж����������Ƿ�����
			int v1 = pUnchInput[yy*nWidth+xx];
			int v2 = pUnchInput[nCurrY*nWidth+nCurrX];
			if ( (xx < nWidth) && (xx>=0) && (yy>=0)&&(yy<nHeight)
				&& (nf ==0)&&(fabs(float(v1-v2))<nThreshold)
				&& (fabs(float(v0-v1))<1.05*nThreshold))
			{
				// ��ջ��β��ָ�����һλ
				nEnd++;
				// ����(xx��yy) ѹ��ջ
				pnGrowQueX[nEnd] = xx;
				pnGrowQueY[nEnd] = yy;
				// ������(xx��yy)���ó��߼�����
				// ͬʱҲ���������ش����
				pUnchOutput[yy*nWidth+xx] = 255 ;
			}
		}
		nStart++;
	}
	// �ͷ��ڴ�
	delete []pnGrowQueX;
	delete []pnGrowQueY;
	pnGrowQueX = NULL ;
	pnGrowQueY = NULL ;
}

double* detect_ellipse_global_region(IplImage* pimg1u)
{
	double*  elines = NULL;
	//ȫͼ������Բ���
	CvBox2D box0, box1;
	//�������
	CvSeq* c = NULL;
	int w = pimg1u->width;
	int h = pimg1u->height;
	int s1 = 3.0;
	int s2 = s1*s1;
	int ww = w/s1;
	int hh = h/s1;
	IplImage* pimg1u0   = cvCreateImage(cvSize(ww, hh), 8, 1);;
	cvResize(pimg1u, pimg1u0);
	CvMemStorage* mem = cvCreateMemStorage(0);
	CvSeq* regions    = NULL;
	vector<CvPoint2D32f> pts;
	cvSmooth( pimg1u0, pimg1u0, CV_MEDIAN, MAX(9/s1, 3), MAX(9/s1, 3));
	IplImage* pimg1b = cvCreateImage(cvGetSize(pimg1u0), 8, 1);
	double area_max = 22000/s2;
	double area_min = 2000/s2;
	for (int i=0 ;i<2; i++)
	{
		pts.clear();
		int nthreshold = 20-i*5;
		//nthreshold = otsu2(pimg1u);
		cvThreshold(pimg1u0, pimg1b, nthreshold, 255, CV_THRESH_BINARY_INV);
		//������ͨ
		//cvDilate(pimg1b, pimg1b);
		//cvErode(pimg1b, pimg1b);
#ifdef DIS_CV
		cvShowImage("pimg1b", pimg1b);
		cvWaitKey(1);
#endif
		cvFindContours(pimg1b, mem, &regions, sizeof(CvContour), CV_RETR_EXTERNAL, CV_CHAIN_APPROX_TC89_KCOS);
		
		for(c=regions; c!=NULL;c=c->h_next)
		{
			double  region=cvContourArea(c);
			if ((fabs(region)>area_max)||(fabs(region)<area_min))
				continue;
			//�ж������Բ�� l^2/4piS
			double  len   = cvContourPerimeter(c);
			double dcircle = len*len/(4*CV_PI*region);
			if (dcircle> 2.0)
				continue;
			CvBox2D box=cvMinAreaRect2(c,mem);
			double area = 3.141516*box.size.height*box.size.width/4.0;
			//if (fabs(area-region)/area>0.5)
				//continue;
			if (box.center.x<50/s1||box.center.x>pimg1u0->width-50/s1)
				continue;
			if (fabs(MAX(box.size.height, box.size.width)/MIN(box.size.height, box.size.width)-1.0)<1.0)
			{
				int onetourlength = c->total;   
				CvSeqReader reader;       //-- ������һ����������
				CvPoint pt = cvPoint(0,0);   
				cvStartReadSeq(c, &reader);       //��ʼ��ȡ   
				for(int i = 0 ;i < onetourlength; i++){   
					CV_READ_SEQ_ELEM(pt,reader);     //--������һ�������е�һ��Ԫ�ص�
					pts.push_back(cvPointTo32f(pt));
				} 
				//opencv ��Բ���
				CalcCirclePara2(pts, box0);
				elines = new double[5];     // by Max������ò�ƻ�����ڴ�й©��
				elines[0] = box0.center.x*s1;
				elines[1] = box0.center.y*s1;
				elines[2] = box0.size.width*s1/2.0;
				elines[3] = box0.size.height*s1/2.0;
				elines[4] = box0.angle*3.14159/180.0;
				cout<<i<<endl;
				goto stop;
				
			}
		}
	}
stop:
	cvReleaseImage(&pimg1b);
	cvReleaseMemStorage(&mem);
	cvReleaseImage(&pimg1u0);
	return elines;
}

double* detect_ellipse_local_cv(IplImage* img, CvRect roi, CvBox2D ep)
{
	double*  elines = NULL;
	//ȫͼ������Բ���
	CvBox2D box0, box;
	int w = img->width;
	int h = img->height;
	//roi
	int ww = roi.width;
	int hh = roi.height;
	IplImage* imgroi = cvCreateImage(cvSize(ww, hh), 8, 1);
	cvSetImageROI(img, roi);
	cvCopy(img, imgroi);
	cvResetImageROI(img);
	int s1 = 1.0;
	int s2 = s1*s1;
	int www = ww/s1;
	int hhh = hh/s1;
	IplImage* img1u0   = cvCreateImage(cvSize(www, hhh), 8, 1);;
	IplImage* img1u1   = cvCreateImage(cvSize(www, hhh), 8, 1);;
	cvResize(imgroi, img1u0);    
	//opencv,��ֵ�˲�
	cvSmooth(img1u0, img1u1, CV_MEDIAN, 5/s1, 5/s1);
	int npara1 = 30,  npara2 = 60;                             //canny��Ե������    
	int* nes = 0;
	int  ns  = 0;
	point** edges = getedge(img1u1, 20/s1, npara1, npara2, &nes, &ns);   //�б�Եͼ��õ���Ե����
	//��ʼ����ȷ��������
	CvSeq* contour;
	int npts = 36;
	CvPoint* epts = new CvPoint[npts];
	ep.center = cvPoint2D32f(roi.width/s1/2.0, roi.height/s1/2.0);
	ep.size.width  /= s1; 
	ep.size.height /= s1;
	cvEllipse2Poly(cvPointFrom32f(ep.center), cvSize(ep.size.width/2.0, ep.size.height/2.0), ep.angle, 0, 350, epts, 360/npts);
	int nR = 8/s1;
	CvMemStorage* storage = cvCreateMemStorage(0);
	contour = cvCreateSeq( CV_SEQ_POLYGON, /* sequence of integer elements */
		sizeof(CvContour), /* header size - no extra fields */
		sizeof(CvPoint), /* element size */
		storage /* the container storage */ );
	
	for (int i=0; i<npts; i++)
	{
		cvSeqPush(contour,&epts[i]);
	}
	if (edges)
	{
		int* nboxs = new int[ns];
		int ee = 0;
		int ntotal = 0;
		for (int e=0; e<ns; e++)
		{
			//�жϱ�Ե�˵㲿�ֻҶȾ�ֵ
			CvPoint p1 = cvPoint(edges[e][0].x, edges[e][0].y);
			CvPoint p2 = cvPoint(edges[e][nes[e]-1].x, edges[e][nes[e]-1].y);
			CvPoint p0 = cvPoint(edges[e][nes[e]/2].x, edges[e][nes[e]/2].y);
			//float v1 = getRegionMeanValue(img1u1, p0, 1);
			//if(v1>100) continue;
			double ret1 = cvPointPolygonTest( contour, cvPointTo32f(p1), 1); //���ֶ�ָ������ʱ������1 ��cvFindContoursʱ0,1�Կ�
			double ret2 = cvPointPolygonTest( contour, cvPointTo32f(p2), 1); //���ֶ�ָ������ʱ������1 ��cvFindContoursʱ0,1�Կ�
			if (fabs(ret1)<nR&&(fabs(ret2)<nR))
			{
				nboxs[ee++] = e;                     //�洢����Ҫ������еı�Ե����
				ntotal += nes[e];
			}	
		}

		int eee = ee;
		if (eee)
		{
			ee = 0;
			point* epts = new point[ntotal];
			for (int e=0; e<eee; e++)
			{
				//opencv ��Բ���
				for (int i=0 ;i<nes[nboxs[e]]; i++)
				{
					epts[ee++] = edges[nboxs[e]][i];  
				}
			}
			if (ntotal)
			{
				CalcCirclePara2(epts, ntotal, box);                                        //ͫ�ױ�Ե���������Բ
				box0 = box;
				box0.center.x = box0.center.x*s1;
				box0.center.y = box0.center.y*s1;
				box0.size.width   = box0.size.width*s1;
				box0.size.height  = box0.size.height*s1;
				elines = new double[5];
				elines[0] = box0.center.x + roi.x;
				elines[1] = box0.center.y + roi.y;
				elines[2] = box0.size.width/2.0;
				elines[3] = box0.size.height/2.0;
				elines[4] = box0.angle*3.14159/180.0;
			}
			delete[] epts;
		}
		//for test
#ifdef DIS_CV
		cvZero(img1u1);
		for (int e=0; e<ns; e++)
		{
			CvPoint p0 = cvPoint(edges[e][nes[e]/2].x, edges[e][nes[e]/2].y);
			float v1 = getRegionMeanValue(img1u1, p0, 1);
			if(v1>100)
				continue;
			for (int j=0; j<nes[e]; j++)
			{
				CvPoint pt = cvPoint(edges[e][j].x, edges[e][j].y);
				((uchar*)img1u1->imageData+pt.y*img1u1->widthStep)[pt.x] = 255;
			}
		}
#endif
		//�ͷ��ڴ�
		for (int i=0;i<ns;i++)
		{
			delete[] edges[i];
		}
		delete[] edges;
		delete[] nes;
	}

	//for test
#ifdef DIS_CV
	CvPoint* epts1 = new CvPoint[npts];
	CvPoint* epts2 = new CvPoint[npts];
	CvPoint* epts3 = new CvPoint[npts];
	cvEllipse2Poly(cvPointFrom32f(ep.center), cvSize(ep.size.width/2.0+nR, ep.size.height/2.0+nR), ep.angle, 0, 350, epts1, 360/npts);
	cvEllipse2Poly(cvPointFrom32f(ep.center), cvSize(ep.size.width/2.0-nR, ep.size.height/2.0-nR), ep.angle, 0, 350, epts2, 360/npts);
	cvEllipse2Poly(cvPointFrom32f(box.center), cvSize(box.size.width/2.0, box.size.height/2.0), box.angle, 0, 350, epts3, 360/npts);
	IplImage* imgroi3u = cvCreateImage(cvGetSize(img1u1), 8, 3);
	//cvCanny(img1u1, img1u1, npara1, npara2);
	cvCvtColor(img1u1, imgroi3u, CV_GRAY2RGB);
	cvPolyLine(imgroi3u, &epts, &npts, 1, 1, CV_RGB(255, 0, 0), 1, 16);
	cvPolyLine(imgroi3u, &epts1, &npts, 1, 1, CV_RGB(255, 255, 0), 1, 16);
	cvPolyLine(imgroi3u, &epts2, &npts, 1, 1, CV_RGB(255, 255, 0), 1, 16);
	cvPolyLine(imgroi3u, &epts3, &npts, 1, 1, CV_RGB(0, 255, 0), 1, 16);
	cvShowImage("imgroi3u", imgroi3u);
	cvWaitKey(1);
	cvReleaseImage(&imgroi3u);
	delete[] epts1;
	delete[] epts2;
	delete[] epts3;
#endif

	delete[] epts;
	cvReleaseImage(&imgroi);
	cvReleaseImage(&img1u0);
	cvReleaseImage(&img1u1);
	cvReleaseMemStorage(&storage);
	return elines;
}

double* detect_ellipse_local_Starburst(IplImage* img, CvRect roi, CvBox2D ep)
{
	double*  elines = NULL;
	//ȫͼ������Բ���
	CvBox2D box0;
	int w = img->width;
	int h = img->height;
	//roi
	int ww = roi.width;
	int hh = roi.height;
	IplImage* imgroi = cvCreateImage(cvSize(ww, hh), 8, 1);
	cvSetImageROI(img, roi);
	cvCopy(img, imgroi);
	cvResetImageROI(img);
	int s1 = 1.0;
	int s2 = s1*s1;
	int www = ww/s1;
	int hhh = hh/s1;
	IplImage* img1u0   = cvCreateImage(cvSize(www, hhh), 8, 1);;
	IplImage* img1u1   = cvCreateImage(cvSize(www, hhh), 8, 1);;
	cvResize(imgroi, img1u0);    
	//opencv,��ֵ�˲�
	cvSmooth(img1u0, img1u1, CV_MEDIAN, 5/s1, 5/s1);
	//parameters for the algorithm
	int edge_threshold = 20;		//threshold of pupil edge points detection
	int rays = 18;				//number of rays to use to detect feature points
	int min_feature_candidates = 10;	//minimum number of pupil feature candidates
	//starburst pupil contour detection�� ��ȡͫ�����ĵ㣬˼·�ǲ��ϴӳ�ʼͫ�����ķ������ߣ��������߷��������ݶ����ֵ�㣬�������Բ����������ֱ��ͫ�����Ĳ���
	starburst_pupil_contour_detection((UINT8*)img1u1->imageData, img1u1->width, img1u1->height,
		edge_threshold, rays, min_feature_candidates);
	int *inliers_index;
	inliers_num = 0;
	inliers_index = pupil_fitting_inliers((UINT8*)img1u1->imageData, img1u1->width, img1u1->height, inliers_num);
#ifdef DIS_CV
	IplImage* img3u = cvCreateImage(cvGetSize(img1u0), 8, 3);
	cvCvtColor(img1u1, img3u, CV_GRAY2RGB);
	bool is_inliers = 0;
	for (int i = 0; i < edge_point.size(); i++) {
		is_inliers = 0;
		for (int j = 0; j < inliers_num; j++) {
			if (i == inliers_index[j])
				is_inliers = 1;
		}
		stuDPoint *edge = edge_point.at(i);
		if (is_inliers)
			cvCircle(img3u, cvPoint((int)edge->x, (int)edge->y), 1, CV_RGB(0, 0, 255));
		else
			cvCircle(img3u, cvPoint((int)edge->x, (int)edge->y), 1, CV_RGB(0, 255, 0));
	}
	cvShowImage("img3u", img3u);
	cvWaitKey(1);
	cvReleaseImage(&img3u);
#endif 
	free(inliers_index);

	box0.center.x = pupil_param[2]*s1;
	box0.center.y = pupil_param[3]*s1;
	box0.size.width   = pupil_param[0]*2*s1;
	box0.size.height  = pupil_param[1]*2*s1;
	box0.angle        = pupil_param[4]*180/3.14159;
	elines = new double[5];
	elines[0] = box0.center.x + roi.x;
	elines[1] = box0.center.y + roi.y;
	elines[2] = box0.size.width/2.0;
	elines[3] = box0.size.height/2.0;
	elines[4] = box0.angle*3.14159/180.0;

	cvReleaseImage(&imgroi);
	cvReleaseImage(&img1u0);
	cvReleaseImage(&img1u1);

	return elines;
}

void detect_pupil_obj_fast(IplImage* pimg1u, box2d &pupil, bool bfit)
{
	int rp = 80;
	double* ept = NULL;
	double* ept2 = NULL;
	bool bfirst = true;
	//����һ֡�����������Բ��⣬����ȫͼ����
	int w = pimg1u->width;
	int h = pimg1u->height;
	CvRect roi = cvRect(0, 0, w, h);
	//�ӱ߽��ж�
	if (pupil.center.x != 0)
	{
		if (check_in_image_range(pupil.center.x, pupil.center.y, w, h, rp))
		{
			roi = cvRect(pupil.center.x-rp, pupil.center.y-rp, 2*rp, 2*rp);
		}
		bfirst = false;
	}
	ept = detect_ellipse_global_region(pimg1u);
	if (ept)
	{
		pupil.center.x = ept[0];
		pupil.center.y = ept[1];
		pupil.size.x   = ept[2]*2;
		pupil.size.y   = ept[3]*2;
		pupil.angle    = ept[4]*180/3.14159;     //���ȵ���
		if (bfit)
		{
			if (check_in_image_range(pupil.center.x, pupil.center.y, w, h, rp))
			{
				CvBox2D roi2d;
				roi2d.center.x = ept[0];
				roi2d.center.y = ept[1];
				roi2d.size.width    = ept[2]*2;
				roi2d.size.height   = ept[3]*2;
				roi2d.angle         = ept[4]*180/3.14159;     //���ȵ���
				roi = cvRect(pupil.center.x-rp, pupil.center.y-rp, 2*rp, 2*rp);
				ept2 = detect_ellipse_local_Starburst(pimg1u, roi, roi2d);
				if (ept2)
				{
					if (((ept[0]-ept2[0])*(ept[0]-ept2[0])+(ept[1]-ept2[1])*(ept[1]-ept2[1]))<20*20)
					{
						pupil.center.x = ept2[0];
				        pupil.center.y = ept2[1];
						pupil.size.x   = ept2[2]*2;
						pupil.size.y   = ept2[3]*2;
						pupil.angle    = ept2[4]*180/3.14159;     //���ȵ���
					}
				}
			}
		}		
	}



}



