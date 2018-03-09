#include "stdafx.h"
#include "IRLightDetectionRef.h"

#include <opencv2/legacy/compat.hpp>
//#define  DIS_CV
void AbsAngle(const Mat& cmplx32FC2, Mat& mag32FC1, Mat& ang32FC1)
{
	CV_Assert(cmplx32FC2.type() == CV_32FC2);
	mag32FC1.create(cmplx32FC2.size(), CV_32FC1);
	ang32FC1.create(cmplx32FC2.size(), CV_32FC1);

	for (int y = 0; y < cmplx32FC2.rows; y++)
	{
		const float* cmpD = cmplx32FC2.ptr<float>(y);
		float* dataA = ang32FC1.ptr<float>(y);
		float* dataM = mag32FC1.ptr<float>(y);
		for (int x = 0; x < cmplx32FC2.cols; x++, cmpD += 2)
		{
			dataA[x] = atan2(cmpD[1], cmpD[0]);
			dataM[x] = sqrt(cmpD[0] * cmpD[0] + cmpD[1] * cmpD[1]);
		}
	}
}

void GetCmplx(const Mat& mag32F, const Mat& ang32F, Mat& cmplx32FC2)
{
	CV_Assert(mag32F.type() == CV_32FC1 && ang32F.type() == CV_32FC1 && mag32F.size() == ang32F.size());
	cmplx32FC2.create(mag32F.size(), CV_32FC2);
	for (int y = 0; y < mag32F.rows; y++)
	{
		float* cmpD = cmplx32FC2.ptr<float>(y);
		const float* dataA = ang32F.ptr<float>(y);
		const float* dataM = mag32F.ptr<float>(y);
		for (int x = 0; x < mag32F.cols; x++, cmpD += 2)
		{
			cmpD[0] = dataM[x] * cos(dataA[x]);
			cmpD[1] = dataM[x] * sin(dataA[x]);
		}
	}
}

Mat GetSR(const Mat &img)
{
	Size sz(64, 64);
	Mat img1f[2], sr1f, cmplxSrc2f, cmplxDst2f;
	img1f[1] = img.clone();
	resize(img1f[1], img1f[0], sz, 0, 0, CV_INTER_AREA); 

	img1f[1] = Mat::zeros(sz, CV_32F);
	merge(img1f, 2, cmplxSrc2f);
	dft(cmplxSrc2f, cmplxDst2f);
	AbsAngle(cmplxDst2f, img1f[0], img1f[1]);

	log(img1f[0], img1f[0]);
	blur(img1f[0], sr1f, Size(3, 3));
	sr1f = img1f[0] - sr1f;

	exp(sr1f, sr1f);
	GetCmplx(sr1f, img1f[1], cmplxDst2f);
	dft(cmplxDst2f, cmplxSrc2f, DFT_INVERSE | DFT_SCALE);
	split(cmplxSrc2f, img1f);

	pow(img1f[0], 2, img1f[0]);
	pow(img1f[1], 2, img1f[1]);
	img1f[0] += img1f[1];

	GaussianBlur(img1f[0], img1f[0], Size(3, 3), 0);
	normalize(img1f[0], img1f[0], 0, 1, NORM_MINMAX);
	resize(img1f[0], img1f[1], img.size(), 0, 0, INTER_CUBIC);

	return img1f[1];
}

// --- ---以上为显著度检测




//至少五点
void CalcEllipsePara(vector<CvPoint2D32f> pts,  CvBox2D& box)
{
	//x^2+by^2+ax+by+cxy=d
	CvPoint2D32f* edge= new CvPoint2D32f[pts.size()];
	//图像坐标系转为笛卡尔坐标系
	for (int i = 0; i<pts.size(); i++)
	{
		edge[i]=pts[i];
	}
	cvFitEllipse(edge, pts.size(), &box);
	delete[] edge;
}

//至少三个点
void CalcCirclePara(vector<CvPoint2D32f> pts,  CvBox2D& box)
{
	//x^2+y^2-ax-by=c
	//图像坐标系转为笛卡尔坐标系
	int nUnknown = 3;
	int EquNum = pts.size();
	float  *fCoeff = new float[EquNum * nUnknown];
	float  *fConst = new float[EquNum];

	for (int i=0; i<pts.size(); i++)
	{
		fCoeff[0+i*3] = pts[i].x;	fCoeff[1+i*3] = pts[i].y;	    fCoeff[2+i*3] = 1; fConst[i]= -pts[i].x*pts[i].x-pts[i].y*pts[i].y;   //plane1
	}

	cv::Mat A(EquNum, 3, CV_32F, fCoeff);
	cv::Mat B(EquNum, 1, CV_32F, fConst);
	cv::Mat x(3, 1, CV_32F);
	cv::solve(A, B,  x, DECOMP_SVD);
	float a = x.at<float>(0, 0);
	float b = x.at<float>(1, 0);
	float c = x.at<float>(2, 0);
	float r = a*a+b*b-4*c; 
	r = sqrt(r)/2.0;
	box.angle = 0;
	box.center = cvPoint2D32f(-a/2.0, -b/2.0);
	box.size   = cvSize2D32f(2*r, 2*r);
}

#define DIS_CV


//  获得灰度平均值
float get_mean_value(IplImage* img, CvPoint pt, int nR)
{
	// Mean values
	float meanvalue= 0;
	int   count=0;
	for (int i=-nR; i<=nR; i++){
		int inx = i+pt.x;
		if (inx<0)
		{
			inx = 0;
		}
		if (inx>img->width-1)
		{
			inx = img->width-1;
		}
		for(int j=-nR; j<=nR; j++)
		{
			int iny = j+pt.y;
			if (iny<0)
			{
				iny = 0;
			}
			if (iny>img->height-1)
			{
				iny=img->height-1;
			}
			meanvalue += ((uchar*)img->imageData+(iny)*img->widthStep)[inx]; //右
			count++;
		}
	}
	meanvalue/=count;
	return meanvalue;
}


// 求取灰度加权质心，近似亚像素拟合
CvPoint2D32f get_center(IplImage* img, CvPoint pt, int nR)
{
	// Mean values
	double xx = 0;
	double yy = 0;
	int   count=0;
	int   gg  = 0;
	int ggsum = 0;
	for (int i=-nR; i<=nR; i++){
		int inx = i+pt.x;
		if (inx<0)
		{
			continue;
		}
		if (inx>img->width-1)
		{
			continue;
		}
		for(int j=-nR; j<=nR; j++)
		{
			int iny = j+pt.y;
			if (iny<0)
			{
				continue;
			}
			if (iny>img->height-1)
			{
				continue;
			}
			gg = ((uchar*)img->imageData+(iny)*img->widthStep)[inx]; //右
			//if (gg>ggmax)
			{
				ggsum += gg;
			}
			xx = xx + inx*gg;
			yy = yy + iny*gg;
			count++;
		}
	}
	xx /= ggsum;
	yy /= ggsum;

	return cvPoint2D32f(xx, yy);
}

// 在原始图范围像上根据瞳孔大体范围，寻找亮点圆心
vector<CvPoint2D32f> detect_IRlignts2(IplImage* pimg1u, box2d pupil, CvBox2D &cirle)
{
	//double nstart, nend1, nend2, nend3;
    //nstart = cv::getTickCount();
	double rad  = (pupil.size.x/2.0 + pupil.size.y/2.0)/2.0;
	double x0   = pupil.center.x;
	double y0   = pupil.center.y;
	int rr   = 2*int(rad);
	rad = 1.8*rad;
	
	IplImage* pimg = cvCreateImage(cvSize(rr, rr), 8, 1);
	IplImage* pimg0 = cvCreateImage(cvSize(rr, rr), 8, 1);
	int roi_x = int(x0)-rr/2;
	int roi_y = int(y0)-rr/2;
	cvSetImageROI(pimg1u, cvRect(roi_x, roi_y, rr, rr));
	cvCopy(pimg1u, pimg0);
	cvCopy(pimg0, pimg);
#ifdef DIS_CV
	//cvShowImage("pim1u", pimg);
	//cvWaitKey(1);
#endif
	cvResetImageROI(pimg1u);
	cvSmooth(pimg, pimg, CV_GAUSSIAN, 3, 3, 0.8, 0.8);
	cvThreshold(pimg, pimg, 160, 255, CV_THRESH_BINARY);
	cvDilate(pimg, pimg);
#ifdef DIS_CV
	cvShowImage("pim1b", pimg);
	cvWaitKey(1);
#endif
	//连通域判断
	//nend1 = cv::getTickCount();
	CvBox2D box;
	CvSeq* c = NULL;
	CvMemStorage* mem = cvCreateMemStorage(0);
	CvSeq* regions    = NULL;
	vector<CvPoint2D32f> pts0;
	vector<CvPoint2D32f> pts1;
	cvFindContours(pimg, mem, &regions, sizeof(CvContour), CV_RETR_EXTERNAL, CV_CHAIN_APPROX_TC89_KCOS);  //二值区域轮廓检测
	
//	//显著性检测
//	cv::Mat img1u(pimg);
//	cv::Mat img1f;
//	img1u.convertTo(img1f, CV_32F);
//	cv::Mat sal_img1f = GetSR(img1f);
//	sal_img1f = sal_img1f*255.0;
//	cv::Mat sal_img1u;
//	sal_img1f.convertTo(sal_img1u, CV_8U);
//	IplImage* pim1u = &sal_img1u.operator IplImage();
//#ifdef DIS_CV
//	cvShowImage("pim1s", pim1u);
//	cvWaitKey(0);
//#endif

	double rad2 = rad*rad;
	//nend2 = cv::getTickCount();
	for(c=regions; c!=NULL;c=c->h_next)
	{
		double  region=cvContourArea(c);
		if ((fabs(region)>200)||(fabs(region)<5))
			continue;
		CvBox2D box=cvMinAreaRect2(c,mem);
		//亮点应该在虹膜范围内
		if (((box.center.x-rr/2)*(box.center.x-rr/2)+(box.center.y-rr/2)*(box.center.y-rr/2))>rad2)
			continue;
		float mV = get_mean_value(pimg0, cvPointFrom32f(box.center), 20);
		if (mV>150)
		{
			continue;
		}
		if (fabs(MAX(box.size.height, box.size.width)/MIN(box.size.height, box.size.width)-1.0)<2.5)
		{
			CvPoint2D32f pt = get_center(pimg0, cvPointFrom32f(box.center), 5);
			pts0.push_back(cvPoint2D32f(pt.x+roi_x, pt.y+roi_y));
		}
	}
   if (pts0.size()>2)
   {
	    CalcCirclePara(pts0, box);
		cirle = box;
#ifdef DIS_CV
		//IplImage* pimg3u = cvCreateImage(cvGetSize(pimg1u), 8, 3);
		//cvCvtColor(pimg1u, pimg3u, CV_GRAY2RGB);
		////cvCircle(pimg3u, cvPointFrom32f(box.center), box.size.width/2.0, CV_RGB(255, 0, 0), 2, 16);
		//cvShowImage("pim3u", pimg3u);
		//cvWaitKey(1);
		//cvReleaseImage(&pimg3u);
#endif
   }

	cvReleaseMemStorage(&mem);
	//nend3 = cv::getTickCount();
	//cout<<"pre-img:"<<1000*(nend1-nstart)/cv::getTickFrequency()<<"  ms"<<endl;
	//cout<<"pre-con:"<<1000*(nend2-nend1)/cv::getTickFrequency()<<"  ms"<<endl;
	//cout<<"pro-con:"<<1000*(nend3-nend2)/cv::getTickFrequency()<<"  ms"<<endl;
	//int h = pimg->height;
	//int w = pimg->width;
	//BYTE** img = new BYTE*[h];
	//BYTE* imgdata  = (BYTE*)pimg1u->imageData;
	//for (int i = 0; i<h; i++)
	//{
	//	img[i] = new BYTE[w];
	//}
	//for (int i = 0; i<h; i++)
	//{
	//	for (int j=0; j<w; j++)
	//	{
	//		img[i][j] = imgdata[i*w+j];
	//	}
	//}
	////亚像素位置
	//for (int n=0 ;n<pts0.size(); n++)
	//{
	//	double x = pts0[n].x;
	//	double y = pts0[n].y;
	//	//if (FitSubPix_QuadSurfFit(img, w, h, x, y, 0, 5))
	//	pts1.push_back(cvPoint2D32f(x, y));
	//}
	//for (int i = 0; i<h; i++)
	//{
	//	delete [] img[i];
	//}
	//delete[] img;
	cvReleaseImage(&pimg);
	cvReleaseImage(&pimg0);
	return pts0;
}