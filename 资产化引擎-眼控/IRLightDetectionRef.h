//  ͫ���������⣬���Բ



#ifndef PROCESS_IRLP_H
#define PROCESS_IRLP_H
#include <opencv2/opencv.hpp>
#include "SaliantObjectDetectRef.h"
using namespace cv;
using namespace std;

// �������ϵĵ����Բ
// boxΪ�����Բ�ģ��뾶��
void CalcCirclePara(vector<CvPoint2D32f> pts,  CvBox2D& box);

//��ԭʼͼ���������㣨����Ҫ���ں�Ĥ��Χ�ڣ�pupil *1.8�������Բ
//circleΪ�����Բ�ģ��뾶��
vector<CvPoint2D32f> detect_IRlignts2(IplImage* pimg1u, box2d pupil, CvBox2D &cirle);

vector<CvPoint2D32f> detect_IRlignts(IplImage* pimg1u);
#endif
