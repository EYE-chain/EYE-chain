#ifndef PROCESS_PUPIL_H
#define PROCESS_PUPIL_H
#include <opencv2/opencv.hpp>
struct point {int x,y;};
struct pointf {double x,y;};
struct box2d
{
	pointf center;  /* Center of the box.                          */
	point  size;    /* Box width and length.                       */
	float angle;          /* Angle between the horizontal axis           */
	/* and the first side (i.e. length) in degrees */
};

// ¼ì²âÍ«¿×³ÌÐò
void detect_pupil_obj_fast(IplImage* pimg1u, box2d &pupil, bool bfit = true);
#endif