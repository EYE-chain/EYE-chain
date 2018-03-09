#include "stdafx.h"
#include "PupilandBlobDetectDef.h"
#include "IRLightDetectionRef.h"
#include <queue>
extern bool  bfirst;
vector<double> itpts_angs;
extern vector<int>  inds;
CvBox2D circle0;
#define DIS_CV
void detect_pupil_irpts_fast(IplImage* img1u, box2d &pupil, vector<CvPoint2D32f> &irpts, CvPoint2D32f &center, bool bfit)
{
	//初始化
	inds.clear();
	irpts.clear();
	box2d pt;
	CvBox2D circle;
	vector<CvPoint2D32f> pts;
	detect_pupil_obj_fast(img1u, pupil, bfit);
	pts = detect_IRlignts2(img1u, pupil, circle);
	center = circle.center;
	//对检测的亮点进行排序，按照亮点在圆周上的分布
	vector<double> angs;
	if (bfirst)     //此处first指的是离线标定采集的含有八个亮点的图像
	{
		/*pts.clear();
		pts.push_back(cvPoint2D32f(265, 275));          
		pts.push_back(cvPoint2D32f(282, 238));
		pts.push_back(cvPoint2D32f(325, 226));
		pts.push_back(cvPoint2D32f(357, 253));
		pts.push_back(cvPoint2D32f(360, 294));
		pts.push_back(cvPoint2D32f(340, 324));
		pts.push_back(cvPoint2D32f(307, 331));
		pts.push_back(cvPoint2D32f(276, 309));*/
		circle0 = circle;
	}
	for (int i=0; i<pts.size(); i++)
	{
		double x = pts[i].x-circle.center.x;
		double y = pts[i].y-circle.center.y;
		double ang = atan2(y, x); ang = ang*180/CV_PI;
		angs.push_back(ang);
	}
	if (bfirst)
	{
		itpts_angs.clear();
		typedef pair< float, unsigned int> PDI;
		priority_queue< PDI > Q;
		for (int i=0; i<angs.size(); i++)
		{
			Q.push(PDI(-angs[i], i));
		}
		for (int i=0; i<angs.size(); i++)
		{
			itpts_angs.push_back(angs[Q.top().second]);
			irpts.push_back(pts[Q.top().second]);
			inds.push_back(i);
			Q.pop();
		}
	}else
	{
		for (int i=0; i<angs.size(); i++)
		{
			irpts.push_back(pts[i]);

			double ang_min = 100;
			int ind = 0;
			for (int j=0 ;j<itpts_angs.size(); j++)
			{
				double ang = fabs(angs[i]-itpts_angs[j]);
				
				if (ang<ang_min)
				{
					ang_min = ang;
					ind = j;
				}
			}
			inds.push_back(ind);
		}
	}
	//如果拟合的圆不稳定，则亮点位置不稳定，因为不需要每帧都进行亮点提取，
	//可以使用最近的帧的亮点质心位置作为当前帧的亮点质心
	//同样在线计算单应只有当亮点质心位置较为稳定时才进行
	if (MAX((circle.size.width*circle.size.height),(circle0.size.width*circle0.size.height))/
		MIN((circle.size.width*circle.size.height),(circle0.size.width*circle0.size.height))>1.4)
	{
		center = cvPoint2D32f(0, 0);
		irpts.clear();
		
	}
	if ((fabs(circle.center.x-circle0.center.x)+fabs(circle.center.y-circle0.center.y))>50)
	{
		center = cvPoint2D32f(0, 0);
		irpts.clear();
	}
	#ifdef DIS_CV
	if (center.x>0)
	{
		IplImage* pimg3u = cvCreateImage(cvGetSize(img1u), 8, 3);
		cvCvtColor(img1u, pimg3u, CV_GRAY2RGB);
		cvCircle(pimg3u, cvPointFrom32f(circle.center), circle.size.width/2.0, CV_RGB(255, 0, 0), 2, 16);
		cvShowImage("pim3u", pimg3u);
		cvWaitKey(1);
		cvReleaseImage(&pimg3u);

	}
	#endif
}