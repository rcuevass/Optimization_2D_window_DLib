// OpenCV includes
#include <opencv2/opencv.hpp>

#include "simple_optimizer.h"
#include "Utilities.h"
#include "functions.h"
#include <math.h>

simple_optimizer::simple_optimizer(void)
{
}


simple_optimizer::~simple_optimizer(void)
{
}

functions grad;

//Point4D simple_optimizer::gradDescent(double xc, double yc, double r, double mu,\
//	double Lambda, int rPower){

Point5D simple_optimizer::gradDescent(double xc, double yc, double rx, double ry,\
	double theta, double mu, double Lambda, int rPower){

	Point5D optPos;


	double xPos = xc;
	double yPos = yc;
	double rxPos = rx;
	double ryPos = ry;
	double thetaPos = theta;

	double muPos = mu;
	
	int maxCount = 0;
	double eps = 1e-3;
	//double alpha = 20.0;
	double alpha = 50.0;



	double epsAdaGrad = 1e-8;
	double SumGradX = 0.0;
	double SumGradY = 0.0;
	double SumGradRx = 0.0;
	double SumGradRy = 0.0;
	double SumGradTheta = 0.0;
	double SumGradMu = 0.0;


	double Gradthresh;




	Point5D gradVal;
	double functVal; 

	//gradVal = grad.ObjFunction4DGrad(xPos, yPos, rPos, muPos,Lambda,rPower);


	double incrAng = -30.0;

	do 
	{

		cv::Mat Temp = m_mImageData.clone();
		//cv::circle(Temp,cv::Point(xPos,yPos),rPos,1.0);

		double Pi = 3.141592654;

	

		//cout << " Test " << cos(2*Pi) << endl;
		//cin.get();

		// Reference for rotation:
		// http://docs.opencv.org/2.4/modules/imgproc/doc/geometric_transformations.html?highlight=warpaffine#getrotationmatrix2d
		//cv::Mat rotMat =  getRotationMatrix2D(cv::Point(xPos,yPos),45.0,1.0);
		

		//cv::circle(Temp,cv::Point(xPos,yPos),rPos,0.5,-1);
		
		//cv::rectangle(Temp,cv::Point(xPos-rPos/2,yPos - rPos/2),cv::Point(xPos+rPos/2,yPos + rPos/2),1.0);
		

		// Example of rotating a rectangle in OpenCV 
		// http://docs.opencv.org/2.4/modules/core/doc/basic_structures.html#RotatedRect
		

		incrAng -= 0.0; //10.0;
		//cv::RotatedRect rRect = cv::RotatedRect(cv::Point2f(xPos,yPos),cv::Size2f(rPos,rPos/2),incrAng);
		
		cv::RotatedRect rRect = cv::RotatedRect(cv::Point2f(xPos,yPos),
																						cv::Size2f(rxPos,ryPos),-thetaPos);
		cv::Point2f vertices[4];
		rRect.points(vertices);
		for (int i = 0; i < 4; i++){
			//line(Temp, vertices[i], vertices[(i+1)%4], cv::Scalar(0,255,0));
			line(Temp, vertices[i], vertices[(i+1)%4], 0.5);
		}
		//cv::Rect brect = rRect.boundingRect();
		//cv::rectangle(Temp, brect, cv::Scalar(255,0,0));
		//cv::circle(Temp,cv::Point(xPos,yPos),rPos,1.0);



		/*

		Vector3D endVectorX, endVectorY;
		Vector3D vVect;
		endVectorX[0] = xPos + rxPos - xPos;
		endVectorX[1] = yPos+0 - yPos;
		endVectorY[0] = xPos+0 - xPos;
		endVectorY[1] = yPos - ryPos - yPos;
		vVect[2] = 1.0;


		functions rotatePoint;
		endVectorX = rotatePoint.RotateVector(endVectorX,-thetaPos,vVect);
		endVectorY = rotatePoint.RotateVector(endVectorY,-thetaPos,vVect);


		cv::arrowedLine(Temp,cv::Point(xPos,yPos),cv::Point(endVectorX[0]+xPos,endVectorX[1]+yPos),
			cv::Scalar(0.5,0.5,0.5),2);

		cv::arrowedLine(Temp,cv::Point(xPos,yPos),cv::Point(endVectorY[0]+xPos,endVectorY[1]+yPos),
			cv::Scalar(0.5,0.5,0.5),2);

    */   

		cv::imshow("Partial_result" , Temp);
		cv::waitKey(10);



		Vector3D uVect;
		uVect[2] = 1.0; 
		//functVal = grad.ObjFunction4D(xPos, yPos, rPos, muPos,Lambda,rPower);
		functVal = grad.ObjFunction2DRot(xPos,yPos,rxPos,ryPos,muPos,Lambda,rPower,thetaPos,uVect);
			
			

		//gradVal = grad.ObjFunction4DGrad(xPos, yPos, rPos, muPos,Lambda,rPower);
		gradVal = grad.ObjFunction2DRotGrad(xPos,yPos,rxPos,ryPos,muPos,Lambda,rPower,thetaPos,uVect);



		SumGradX += pow(gradVal[0],2);
		SumGradY += pow(gradVal[1],2);
		SumGradRx += pow(gradVal[2],2);
		SumGradRy += pow(gradVal[3],2);
		SumGradTheta += pow(gradVal[4],2);

		//GradX += pow(gradVal[0],2); 


		xPos += alpha*gradVal[0]/sqrt(SumGradX + epsAdaGrad);
		yPos += alpha*gradVal[1]/sqrt(SumGradY + epsAdaGrad);
		//rPos += 1.0*alpha*gradVal[2]/sqrt(SumGradR + epsAdaGrad);
		
		// Fernando's idea 
		//rPos += 0.05*gradVal[2];
		
		
		//rxPos += alpha*gradVal[2]/sqrt(SumGradRx + epsAdaGrad);
		//ryPos += alpha*gradVal[3]/sqrt(SumGradRy + epsAdaGrad);

		//thetaPos += alpha*gradVal[4]; // /sqrt(SumGradTheta + epsAdaGrad);
		//thetaPos += 1.5*gradVal[4];
		//thetaPos = thetaPos - 360.0*floor(thetaPos/360.0);




		//rPos += alpha*gradVal[2]/sqrt(SumGradR + epsAdaGrad);
		
		if (rxPos < 1 )
		{
			rxPos = 1;
		}


		if (ryPos < 1 )
		{
			ryPos = 1;
		}
		//rPos += alpha*gradVal[2];
		//muPos *=  1.1; // alpha*gradVal[3];
		//muPos += alpha*gradVal[3];

		//muPos *= 1.01;
		maxCount += 1;


		Gradthresh = sqrt( pow(gradVal[0],2) + pow(gradVal[1],2) ) ;

		printf("%5i\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%5.5f\t%5.5f\t||%5.5f\t%5.5f\t%5.5f\n",
			maxCount,xPos, yPos,rxPos,ryPos,thetaPos,
			gradVal[0],gradVal[1],Gradthresh);

		





	} while (maxCount<50  && (Gradthresh > eps) );

	




	
	optPos[0] = xPos;
	optPos[1] = yPos;
	//optPos[2] = rPos;
	optPos[2] = rxPos;
	optPos[3] = ryPos;
	optPos[4] = muPos;

	return optPos;

}


void simple_optimizer::SetImageData(cv::Mat mImage)
{
	m_mImageData = mImage;

	grad.SetImage(mImage);
}