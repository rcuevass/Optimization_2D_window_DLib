#pragma once
#include "Utilities.h"


class simple_optimizer
{
public:
	simple_optimizer(void);
	~simple_optimizer(void);


	//Point4D gradDescent(double xc, double yc, double r, double mu,double Lambda, int rPower);
	Point5D gradDescent(double xc, double yc, double rx, double ry, double alpha,
		double mu,double Lambda, int rPower);
	
	
	
	void SetImageData(cv::Mat mImage);


	cv::Mat m_mImageData;
};


