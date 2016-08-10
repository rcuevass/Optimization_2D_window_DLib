// The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt
/*

    This is an example illustrating the use the general purpose non-linear 
    optimization routines from the dlib C++ Library.

    The library provides implementations of the conjugate gradient,  BFGS,
    L-BFGS, and BOBYQA optimization algorithms.  These algorithms allow you to
    find the minimum of a function of many input variables.  This example walks
    though a few of the ways you might put these routines to use.


		This is a useful video to setup dlib in visual studio:

		https://www.youtube.com/watch?v=lhfIbaygA-s

		and a website related with the video:

		http://quantlabs.net/blog/2014/04/how-to-get-great-open-source-dlib-c-bayesian-network-library-working-on-visual-studio/#sthash.9OJNS8p6.dpuf

*/
// OpenCV includes
#include <opencv2/opencv.hpp>



//#include <dlib/optimization.h>
#include <iostream>


#include "functions.h"
#include "Utilities.h"
//#include <vector>

#include "simple_optimizer.h"
#include <iostream>
#include <math.h>



using namespace std;
//using namespace cv;
//using namespace std;
//using namespace dlib;




// ----------------------------------------------------------------------------------------

// In dlib, the general purpose solvers optimize functions that take a column
// vector as input and return a double.  So here we make a typedef for a
// variable length column vector of doubles.  This is the type we will use to
// represent the input to our objective functions which we will be minimizing.

//typedef matrix<double,0,1> column_vector;

typedef matrix<double,0,1> column_vector;


// Parameters of our sigmoid-like 
//double mu = 0.5;
// Class for functions
functions getFunct;

//double xCoord = double(35);
//double yCoord = double(35);
//double r = double(15);
//double mu = 0.5;
 
//int xMin = 30;
//int yMin = 30;
//int rMin = 5;
//double muMin = 0.1;


//int xMax = 70;
//int yMax = 70;
//int rMax = 20;
//double muMax = 1.0;


// ----------------------------------------------------------------------------------------

int main()
{
	


	//cout << "Let's test dlib idiot!";
	//cout << endl;
	//cin.get();

	//imshow("Culero", Test);
	//cv::waitKey();
	//exit(0);


	//for (int i=0;i<1000;i++)
	//	pTest[i+50000] = 1.0;


	//for (int j = 250000; j<750000;j++)
	//	for (int i = 250; i < 750; i++)
	//		pTest[i + j] = 1.0;




	//cv::circle(Test,cv::Point(500,500),100,0.5,-1);
	//cv::rectangle(Test,cv::Rect(200,200,100,100),0.5,1);
	//cv::imshow("puto",Test);
	//cv::waitKey();
	//exit(0);

	simple_optimizer myGrad;

	// Data parameters
	int xInt = 220;
	int yInt = 220;  
	int radiusIntX = 60;
	int radiusIntY = 60;
	double muInit = 0.8;
	double beta = 78.46; 


	// Synthetic image
	cv::Mat Test(1000,1000,CV_64F,0.0);
	double * pTest = &Test.at<double>(0,0);
	
	// RECTANGLE
	cv::rectangle(Test,cv::Point(xInt-radiusIntX/1.0,yInt - radiusIntY/1.0),
		cv::Point(xInt+radiusIntX/1.0,yInt + radiusIntY/1.0),1.0,-1);
	
	// Reference for rotations:
	// http://felix.abecassis.me/2011/10/opencv-rotation-deskewing/
	// http://docs.opencv.org/2.4/modules/imgproc/doc/geometric_transformations.html?highlight=warpaffine#warpaffine
	cv::Mat rot_mat =  getRotationMatrix2D(cv::Point(xInt,yInt),beta,1.0);
	//cv::Mat rotated; 
	
	// Apply rotation...
	cv::warpAffine(Test, Test, rot_mat, Test.size());
	//cv::imshow("Rotated?" , Test);
	//cv::waitKey(10);
	//cin.get();


	// CIRCLE
	//cv::circle(Test,cv::Point(xInt,yInt),radiusInt,0.5,-1);


	// ROTATED RECTANGLE

	//cv::RotatedRect rRect = cv::RotatedRect(cv::Point2f(xInt,yInt),
	//	cv::Size2f(radiusInt,2.0*radiusInt),25.0);
	//cv::Point2f vertices[4];
	//rRect.points(vertices);
	//for (int i = 0; i < 4; i++){
		//line(Temp, vertices[i], vertices[(i+1)%4], cv::Scalar(0,255,0));
	//	line(Test, vertices[i], vertices[(i+1)%4], 0.5);
	//}


  ////////////////////////////////////////////////////////
	// 
	//				PCA - Principal component analysis
	//
	////////////////////////////////////////////////////////


	cv::Mat TestCopy;
	Test(cv::Rect(100,100,200,200)).copyTo(TestCopy);
	//Test(cv::Rect(100,100,400,400)).copyTo(TestCopy);

	//cv::imshow("Sub-image" , Test);

	/*
	Vector3D endVectorX, endVectorY;
	Vector3D uVect;
	endVectorX[0] = xInt + radiusIntX - xInt;
	endVectorX[1] = yInt+0 - yInt;
	endVectorY[0] = xInt+0 - xInt;
	endVectorY[1] = yInt - radiusIntY - yInt;
	uVect[2] = 1.0;

	
	functions rotatePoint;
	endVectorX = rotatePoint.RotateVector(endVectorX,-beta,uVect);
	endVectorY = rotatePoint.RotateVector(endVectorY,-beta,uVect);

	*/

	//cv::arrowedLine(Test,cv::Point(xInt,yInt),cv::Point(endVectorX[0]+xInt,endVectorX[1]+yInt),
	//	cv::Scalar(0.5,0.5,0.5),2);

	//cv::arrowedLine(Test,cv::Point(xInt,yInt),cv::Point(endVectorY[0]+xInt,endVectorY[1]+yInt),
	//	cv::Scalar(0.5,0.5,0.5),2);
	



	cv::PCA PCAFeatures(TestCopy,cv::Mat(),CV_PCA_DATA_AS_ROW);
	Point2D m_pCenter = Point2D(PCAFeatures.mean.at<double>(0,0),
														  PCAFeatures.mean.at<double>(0,1));
															//PCAFeatures.mean.at<double>(0,2));
	// Reference: http://stackoverflow.com/questions/9360375/getting-error-for-ambiguous-symbol-and-need-help-to-remove-it
	// std indicates we are using OpenCV vector as opposed to Dlib's
	std::vector<Point2D> m_vpEigenVectors;
	std::vector<double> m_vdEigenValues;
	// if we wanted to use dlib's one:
	//dlib::vector<Point2D>  m_vpEigenvectors;




	//Store the center of the object
	cv::Point cntr = cv::Point((PCAFeatures.mean.at<double>(0, 0)),
		(PCAFeatures.mean.at<double>(0, 1)));

	cout << "Center = " << cntr << endl;


	// Store eigenvalues and eigenvectors:
  m_vpEigenVectors.resize(2);
	m_vdEigenValues.resize(2);




	for (int i = 0; i < 2 ; i++){
		m_vpEigenVectors[i] = Point2D(PCAFeatures.eigenvectors.at<double>(i,0),
																	PCAFeatures.eigenvectors.at<double>(i,1));
															 //,PCAFeatures.eigenvectors.at<double>(i,2));
		m_vdEigenValues[i] = PCAFeatures.eigenvalues.at<double>(0,1);

		cout << "Element " << i << endl;
		cout << "Eigenvalue " << m_vdEigenValues[i] << endl; 
		cout << "Eigenvector " << m_vpEigenVectors[i][0] << " " << 
			m_vpEigenVectors[i][2] << endl; // " " << m_vpEigenVectors[i][3] << endl;
	}



	// Draw the principal components
	//cv::circle(img, cntr, 3, cv::Scalar(255, 0, 255), 2);
	//cv::Point p1 = cntr + 0.02 * cv::Point(static_cast<int>(eigen_vecs[0].x * eigen_val[0]),
	//	static_cast<int>(eigen_vecs[0].y * eigen_val[0]));
	//cv::Point p2 = cntr - 0.02 * cv::Point(static_cast<int>(eigen_vecs[1].x * eigen_val[1]),
	//	static_cast<int>(eigen_vecs[1].y * eigen_val[1]));
	//drawAxis(img, cntr, p1, cv::Scalar(0, 255, 0), 1);
	//drawAxis(img, cntr, p2, cv::Scalar(255, 255, 0), 5);
	//double angle = atan2(eigen_vecs[0].y, eigen_vecs[0].x); // orientation in radians



//cin.get();






	// Set the synthetic image as the input for the gradient descent fn.
	myGrad.SetImageData(Test);


	// Initial guess
	double xCoord;
	double yCoord;
	double rx;
	double ry;
	double r;
	double theta;
	double mu;

	// for regularization...
	int rPower = 2;
	double lambdaVal = 0.0; //2.0;

	xCoord = 300;
	yCoord = 300;
	r = 60;
	rx = 60;
	ry = 60;
	theta = 0.0;


	//mu = 0.5;
	mu = 0.02; // 0.05;


	cout << " Initial coordinates ..." << endl;
	cout << xInt << " " << " " << yInt << " " <<
		radiusIntX << " " << radiusIntY << " " << muInit << endl;

	//myGrad.gradDescent(xCoord, yCoord, r, mu, lambdaVal,rPower);

	Point5D finalPositions;

	finalPositions = myGrad.gradDescent(xCoord,yCoord,rx,ry,theta,mu,lambdaVal,rPower);



	cv::Mat FinalIm = Test.clone();


	// Arrows on target image

	
	
	Vector3D tarVectorX, tarVectorY;
	Vector3D uVect;
	tarVectorX[0] = xInt + radiusIntX - xInt;
	tarVectorX[1] = yInt+0 - yInt;
	tarVectorY[0] = xInt+0 - xInt;
	tarVectorY[1] = yInt - radiusIntY - yInt;
	uVect[2] = 1.0;

	
	functions rotatePoint;
	tarVectorX = rotatePoint.RotateVector(tarVectorX,-beta,uVect);
	tarVectorY = rotatePoint.RotateVector(tarVectorY,-beta,uVect);

	
	cv::arrowedLine(Test,cv::Point(xInt,yInt),cv::Point(tarVectorX[0]+xInt,tarVectorX[1]+yInt),
		0.0,2);

	cv::arrowedLine(Test,cv::Point(xInt,yInt),cv::Point(tarVectorY[0]+xInt,tarVectorY[1]+yInt),
		0.0,2);
	


	// Arrows on learned image


	double xPos = finalPositions[0];
	double yPos = finalPositions[1];
	double rxPos = finalPositions[2];
	double ryPos = finalPositions[3];


	Vector3D endVectorX, endVectorY;
	Vector3D vVect;
	endVectorX[0] = xPos + rxPos - xPos;
	endVectorX[1] = yPos+0 - yPos;
	endVectorY[0] = xPos+0 - xPos;
	endVectorY[1] = yPos - ryPos - yPos;
	vVect[2] = 1.0;


	//functions rotatePoint;
	endVectorX = rotatePoint.RotateVector(endVectorX,-theta,vVect);
	endVectorY = rotatePoint.RotateVector(endVectorY,-theta,vVect);


	cv::arrowedLine(Test,cv::Point(xPos,yPos),cv::Point(endVectorX[0]+xPos,endVectorX[1]+yPos),
		cv::Scalar(0.5,0.5,0.5),2);

	cv::arrowedLine(Test,cv::Point(xPos,yPos),cv::Point(endVectorY[0]+xPos,endVectorY[1]+yPos),
		cv::Scalar(0.5,0.5,0.5),2);


	cv::RotatedRect rRect = cv::RotatedRect(cv::Point2f(xPos,yPos),
		cv::Size2f(rxPos,ryPos),-theta);
	cv::Point2f vertices[4];
	rRect.points(vertices);
	for (int i = 0; i < 4; i++){
		//line(Temp, vertices[i], vertices[(i+1)%4], cv::Scalar(0,255,0));
		line(Test, vertices[i], vertices[(i+1)%4], 0.5);
	}


	cv::imshow("Partial_rotation" , Test);
	cv::waitKey(10);
	double dotVal = endVectorX.dot(tarVectorX);
	double len1 = endVectorX.length();
	double len2 = tarVectorX.length();
	Vector3D crossVect = tarVectorX.cross(endVectorX);
	double lenCross  = crossVect.length();

	for (int i = 0 ; i <=2 ; i++)
	{
		crossVect[i] /= lenCross;
	}

	cout << "Cross and dot products " << endl;
	cout << dotVal << endl;
	
	cout << crossVect[0] << " " << crossVect[1] << " " << crossVect[2] << endl;
	double angle = acos(dotVal/(len1*len2))*180/3.141592654;
	
	
	cout << " Angle = " << angle;


	cin.get();

	// Draw final rectangle

	cout << "Drawing final image..." << endl;

	endVectorX = rotatePoint.RotateVector(endVectorX,-angle,crossVect);
	endVectorY = rotatePoint.RotateVector(endVectorY,-angle,crossVect);



	cv::arrowedLine(FinalIm,cv::Point(xPos,yPos),cv::Point(endVectorX[0]+xPos,endVectorX[1]+yPos),
		cv::Scalar(0.5,0.5,0.5),2);

	cv::arrowedLine(FinalIm,cv::Point(xPos,yPos),cv::Point(endVectorY[0]+xPos,endVectorY[1]+yPos),
		cv::Scalar(0.5,0.5,0.5),2);


	cv::RotatedRect rRectFin = cv::RotatedRect(cv::Point2f(xPos,yPos),
		cv::Size2f(rxPos,ryPos),-angle);
	cv::Point2f verticesFin[4];
	rRectFin.points(verticesFin);
	for (int i = 0; i < 4; i++){
		line(FinalIm, verticesFin[i], verticesFin[(i+1)%4], 0.5);
	}


	cv::imshow("Final_result" , FinalIm);
	cv::waitKey(10);



	cin.get();

	//exit(0);




	// Rectangle

	int uppXint = 100;
	int uppYint = 100;
	int lowXint = 200;
	int lowYint = 200;



	// For constrained search...


	int xMin = 60;
	int yMin = 60;
	int rMin = 40;
	double muMin = 0.004;


	int xMax = 1000;
	int yMax = 1000;
	int rMax = 500;
	double muMax = 2.0;

		


}

