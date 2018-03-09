// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>



// TODO: 在此处引用程序需要的其他头文件
#include "opencv2/core/version.hpp"
#define CV_VERSION_ID  CVAUX_STR(CV_MAJOR_VERSION) CVAUX_STR(CV_MINOR_VERSION) CVAUX_STR(CV_SUBMINOR_VERSION)
#ifdef _DEBUG
#define cvLIB(name) "opencv_" name CV_VERSION_ID "d"
#else
#define cvLIB(name) "opencv_" name CV_VERSION_ID
#endif

#pragma comment( lib, cvLIB("calib3d") )
//#pragma comment( lib, cvLIB("contrib") )
#pragma comment( lib, cvLIB("core") )
#pragma comment( lib, cvLIB("features2d") )
#pragma comment( lib, cvLIB("flann") )
//#pragma comment( lib, cvLIB("gpu") )
#pragma comment( lib, cvLIB("highgui") )
#pragma comment( lib, cvLIB("imgproc") )
#pragma comment( lib, cvLIB("video") )
#pragma comment( lib, cvLIB("legacy") )
#pragma comment( lib, cvLIB("ml") )