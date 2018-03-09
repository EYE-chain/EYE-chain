#ifndef PROCESS_PUPIL_IRLP_H
#define PROCESS_PUPIL_IRLP_H
#include <string>
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include "SaliantObjectDetectRef.h"
using namespace cv;
using namespace std;


// 瞳孔检测和亮点检测；
// imglu 输入原始图像
// pupil 输出瞳孔位置（圆点，半径）
// irpts 所有亮点坐标（pixel），亮点质心；
// center 所有亮点拟合的圆心
// bfit 是否拟合亮点的圆心；

void detect_pupil_irpts_fast(IplImage* img1u, box2d &pupil, vector<CvPoint2D32f> &irpts, CvPoint2D32f &center, bool bfit = true);
#endif