#ifndef PROCESS_PUPIL_IRLP_H
#define PROCESS_PUPIL_IRLP_H
#include <string>
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include "SaliantObjectDetectRef.h"
using namespace cv;
using namespace std;


// ͫ�׼��������⣻
// imglu ����ԭʼͼ��
// pupil ���ͫ��λ�ã�Բ�㣬�뾶��
// irpts �����������꣨pixel�����������ģ�
// center ����������ϵ�Բ��
// bfit �Ƿ���������Բ�ģ�

void detect_pupil_irpts_fast(IplImage* img1u, box2d &pupil, vector<CvPoint2D32f> &irpts, CvPoint2D32f &center, bool bfit = true);
#endif