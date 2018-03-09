//  瞳孔中亮点检测，拟合圆



#ifndef PROCESS_IRLP_H
#define PROCESS_IRLP_H
#include <opencv2/opencv.hpp>
#include "SaliantObjectDetectRef.h"
using namespace cv;
using namespace std;

// 三个以上的点拟合圆
// box为输出（圆心，半径）
void CalcCirclePara(vector<CvPoint2D32f> pts,  CvBox2D& box);

//在原始图像中找亮点（亮点要求在虹膜范围内（pupil *1.8））拟合圆
//circle为输出（圆心，半径）
vector<CvPoint2D32f> detect_IRlignts2(IplImage* pimg1u, box2d pupil, CvBox2D &cirle);

vector<CvPoint2D32f> detect_IRlignts(IplImage* pimg1u);
#endif
