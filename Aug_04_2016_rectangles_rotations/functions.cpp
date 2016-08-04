
// OpenCV includes
#include <opencv2/opencv.hpp>

#include "functions.h"
#include <math.h>
#include "Utilities.h"
//#include <vector>

#include <dlib/matrix.h>
using namespace dlib;






functions::functions(void)
{
}

functions::~functions(void)
{
}


int kmin = 0;
int kmax = 999;
double valImage = 255; //1e2;


// sFunction and derivatives...
double functions::sFunction(double x, double xc, double r, double mu){
	double valFunt = exp(-mu * ( x - xc + r) );      
	return 1e0/(1e0 + valFunt);
}

double functions::sFunctionDerXc(double x, double xc, double r, double mu){
	double valFunct = sFunction(x,xc,r,mu) - 1.0; 
	return valFunct*mu*sFunction(x,xc,r,mu);
}


double functions::sFunctionDerR(double x, double xc, double r, double mu){
	double valFunct = 1.0 - sFunction(x,xc,r,mu); 
	return valFunct*mu*sFunction(x,xc,r,mu);
}

double functions::sFunctionDerMu(double x, double xc, double r, double mu){
	double valFunct = ( x - (xc - r) );
	return valFunct*sFunction(x,xc,r,mu)*( 1.0 - sFunction(x,xc,r,mu) );
}


// Sigma (sigmoid-like) function and derivatives...
double functions::Sigma(double x, double xc, double r, double mu){
	return sFunction(x,xc,r,mu)*sFunction(xc,x,r,mu);
}

double functions::SigmaDerXc(double x, double xc, double r, double mu){
	double valFunct = sFunction(x,xc,r,mu) - sFunction(xc,x,r,mu);
	return valFunct*mu*Sigma(x,xc,r,mu);
}


double functions::SigmaDerR(double x, double xc, double r, double mu){
	double valFunct = 2.0 - sFunction(x,xc,r,mu) - sFunction(xc,x,r,mu);
	return valFunct*mu*Sigma(x,xc,r,mu);
}

double functions::SigmaDerMu(double x, double xc, double r, double mu){
	double valFunct1 = r*( 2.0 - sFunction(x,xc,r,mu) - sFunction(xc,x,r,mu) );
	double valFunct2 = (x - xc)*( sFunction(xc,x,r,mu) - sFunction(x,xc,r,mu) );
	return Sigma(x,xc,r,mu)*(valFunct1 + valFunct2);
}

double functions::SigmaDer2Xc(double x, double xc, double r, double mu){
	double valFunct1 = sFunction(x,xc,r,mu) - sFunction(xc,x,r,mu);
	valFunct1 = valFunct1*SigmaDerXc(x,xc,r,mu);
	
	double valFunct2 = sFunctionDerXc(x,xc,r,mu) + sFunctionDerXc(xc,x,r,mu);
	valFunct2 = valFunct2*Sigma(x,xc,r,mu);
	
	return mu*(valFunct1 + valFunct2);
}


double functions::SigmaDer2R(double x, double xc, double r, double mu){
	double valFunct1 = 2.0 - sFunction(x,xc,r,mu) - sFunction(xc,x,r,mu);
	valFunct1 = valFunct1*SigmaDerR(x,xc,r,mu);

	double valFunct2 = -sFunctionDerR(x,xc,r,mu) - sFunctionDerR(xc,x,r,mu);
	valFunct2 = valFunct2*Sigma(x,xc,r,mu);

	return mu*(valFunct1 + valFunct2);
}

double functions::SigmaDer2Mu(double x, double xc, double r, double mu){
	double valFunct1 = r*( 2.0 - sFunction(x,xc,r,mu) - sFunction(xc,x,r,mu) );
	double valFunct2 = (x - xc)*( sFunction(xc,x,r,mu) - sFunction(x,xc,r,mu) );
	double valFin1 = SigmaDerMu(x,xc,r,mu)*(valFunct1 + valFunct2);
	double valFin2 = ( ( (x - xc - r)*SigmaDerMu(xc,x,r,mu) ) - \
		( (x - xc + r)*SigmaDerMu(x,xc,r,mu) ) );
	valFin2 = valFin2*Sigma(x,xc,r,mu);

	return (valFin1 + valFin2); 
}


double functions::SigmaDer2XcR(double x, double xc, double r, double mu){
	double valFunct1 = Sigma(x,xc,r,mu)*
										 ( sFunctionDerR(x,xc,r,mu) - sFunction(xc,x,r,mu) );
	double valFunct2 = SigmaDerR(x,xc,r,mu)*
										 ( sFunction(x,xc,r,mu) - sFunction(xc,x,r,mu) );
	return mu*(valFunct1 + valFunct2);
}

double functions::SigmaDer2XcMu(double x, double xc, double r, double mu){
	double valFunct1 = ( Sigma(x,xc,r,mu) + (mu*SigmaDerMu(x,xc,r,mu)) )*
										 ( sFunction(x,xc,r,mu) - sFunction(xc,x,r,mu) );
	double valFunct2 = mu*Sigma(x,xc,r,mu)*
		( sFunctionDerMu(x,xc,r,mu) - sFunctionDerMu(xc,x,r,mu) );
	return (valFunct1 + valFunct2);
}

double functions::SigmaDer2RMu(double x, double xc, double r, double mu){
	double valFunct1 = ( Sigma(x,xc,r,mu) + (mu*SigmaDerMu(x,xc,r,mu)) )*
		( 2.0 - sFunction(x,xc,r,mu) - sFunction(xc,x,r,mu) );
	double valFunct2 = mu*Sigma(x,xc,r,mu)*
		( sFunctionDerMu(x,xc,r,mu) + sFunctionDerMu(xc,x,r,mu) );
	return (valFunct1 - valFunct2);
	
}




// 2D Objective function (xC and yC) and derivatives ...
double functions::ObjFunction2D(double xc, double yc, double r, double mu){

	double valFunct = 0.0;

	for (int i = kmin; i<= kmin; i++)
	{
		for (int j = kmin; j<= kmax; j++)
		{
			valFunct += valImage*Sigma(double (i),xc,r,mu)*Sigma(double(j),yc,r,mu);
		}
	}

	return valFunct;
}

Point2D functions::ObjFunction2DGrad(double xc, double yc, double r, double mu){

	//int kmin = 30;
	//int kmax = 60;
	Point2D valFunct;

	for (int i = kmin; i<= kmin; i++)
	{
		for (int j = kmin; j<= kmax; j++)
		{
			valFunct[0] += valImage*SigmaDerXc(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu);
			valFunct[1] += valImage*Sigma(double(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu);
		}
	}

	return valFunct;
}

Point3D functions::ObjFunction2DHessian(double xc, double yc, double r, double mu){

	//int kmin = 30;
	//int kmax = 60;
	Point3D valFunct;

	for (int i = kmin; i<= kmin; i++)
	{
		for (int j = kmin; j<= kmax; j++)
		{
			// Second partial derivative with respect to Xc
			valFunct[0] += valImage*SigmaDer2Xc(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu);

			// Second partial derivative with respect to Xc and Yc
			valFunct[1] += valImage*SigmaDerXc(double(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu);


			// Second partial derivative with respect to Yc
			valFunct[2] += valImage*Sigma(double(i),xc,r,mu)*SigmaDer2Xc(double(j),yc,r,mu);

		}
	}

	return valFunct;
}


// 3D Objective function (xC, yC and r) and derivatives

double functions::ObjFunction3D(double xc, double yc, double r, double mu){

	double valFunct = 0.0;

	for (int i = kmin; i<= kmin; i++)
	{
		for (int j = kmin; j<= kmax; j++)
		{
			valFunct += valImage*Sigma(double (i),xc,r,mu)*Sigma(double(j),yc,r,mu);
		}
	}

	return valFunct;
}



Point3D functions::ObjFunction3DGrad(double xc, double yc, double r, double mu){

	//int kmin = 30;
	//int kmax = 60;
	Point3D valFunct;

	for (int i = kmin; i<= kmin; i++)
	{
		for (int j = kmin; j<= kmax; j++)
		{
			valFunct[0] += valImage*SigmaDerXc(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu);
			valFunct[1] += valImage*Sigma(double(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu);
			valFunct[2] += valImage*( (SigmaDerR(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu)) + 
										 (Sigma(double(i),xc,r,mu)*SigmaDerR(double(j),yc,r,mu)) ); 

		// for j ends
		}

	// for i ends
	}

	return valFunct;
}


Point3D functions::ObjFunction3DHessianDiag(double xc, double yc, double r, double mu){
	
	Point3D HessVal;
	

	for (int i = kmin; i<= kmin; i++)
	{
		for (int j = kmin; j<= kmax; j++)
		{
			// Second partial derivative with respect to Xc: XX
			HessVal[0] += valImage*SigmaDer2Xc(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu);

			// Second partial derivative with respect to Yc: YY
			HessVal[1] += valImage*Sigma(double(i),xc,r,mu)*SigmaDer2Xc(double(j),yc,r,mu);

			// Second partial derivative with respect to R: RR
			HessVal[2] += valImage*
				( 
				(SigmaDer2R(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu)) +
				(Sigma(double(i),xc,r,mu)*SigmaDer2R(double(j),yc,r,mu)) +
				(2.0*SigmaDerR(double(i),xc,r,mu)*SigmaDerR(double(j),yc,r,mu))
				);
		
			
		// for j ends
		}
	// for i ends
	}

	return HessVal;
}



Point3D functions::ObjFunction3DHessianCross(double xc, double yc, double r, double mu){

	Point3D HessVal;


	for (int i = kmin; i<= kmin; i++)
	{
		for (int j = kmin; j<= kmax; j++)
		{

			// Second partial derivative with respect to Xc and Yc: XY
			HessVal[0] += valImage*SigmaDerXc(double(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu);


			// Second partial derivative with respect to Xc and R: XR
			HessVal[1] += valImage*
				( ( SigmaDer2XcR(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu) ) +
				( SigmaDerXc(int(i),xc,r,mu)*SigmaDerR(double(j),yc,r,mu) ) );

			// Second partial derivative with respect to Yc and R: YR
			HessVal[2] += valImage*
				( ( Sigma(double(i),xc,r,mu)*SigmaDer2Xc(double(j),yc,r,mu) ) +
				( SigmaDerR(int(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu) ) );


			// for j ends
		}
		// for i ends
	}

	return HessVal;
}




matrix<double,3,3> functions::ObjFunction3DHessian(double xc, double yc,
	double r, double mu){
	
	matrix<double,3,3> HessVal;

	for (int i = kmin; i<= kmin; i++)
	{
		for (int j = kmin; j<= kmax; j++)
		{
			// Second partial derivative with respect to Xc: XX
			HessVal(0,0) += valImage*SigmaDer2Xc(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu);

			// Second partial derivative with respect to Yc: YY
			HessVal(1,1) += valImage*Sigma(double(i),xc,r,mu)*SigmaDer2Xc(double(j),yc,r,mu);

			// Second partial derivative with respect to R: RR
			HessVal(2,2) += valImage*
				( 
				(SigmaDer2R(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu)) +
				(Sigma(double(i),xc,r,mu)*SigmaDer2R(double(j),yc,r,mu)) +
				(2.0*SigmaDerR(double(i),xc,r,mu)*SigmaDerR(double(j),yc,r,mu))
				);

			// Second partial derivative with respect to Xc and Yc: XY
			HessVal(0,1) += valImage*SigmaDerXc(double(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu);
			
			HessVal(1,0) = HessVal(0,1);

			// Second partial derivative with respect to Xc and R: XR
			HessVal(0,2) += valImage*
				( ( SigmaDer2XcR(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu) ) +
				( SigmaDerXc(int(i),xc,r,mu)*SigmaDerR(double(j),yc,r,mu) ) );

			HessVal(2,0) = HessVal(0,2);


			// Second partial derivative with respect to Yc and R: YR
			HessVal(1,2) += valImage*
				( ( Sigma(double(i),xc,r,mu)*SigmaDer2Xc(double(j),yc,r,mu) ) +
				( SigmaDerR(int(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu) ) );

			HessVal(2,1) = HessVal(1,2);
		}
	}

	return HessVal;
}


// 4D Objective function (xC, yC, r and mu) and derivatives

double functions::ObjFunction4D(double xc, double yc, double r, double mu,double Lambda, int rPow){

	double valFunct = 0.0;
	double mData = 0;
	double functVal = 0;

	for (int i = kmin; i< kmax; i++)
	{
		for (int j = kmin; j< kmax; j++)
		{
			double valAtPixel = m_mImage.at<double>(i,j);



			// Image TIMES window
			//valFunct += valAtPixel*Sigma(double (i),xc,r,mu)*Sigma(double(j),yc,r,mu);


			// Image MINUS window
			valFunct += valAtPixel - (Sigma(double (i),xc,r,mu)*Sigma(double(j),yc,r,mu) );
			
			// Log-loss; image MINUS window
			//functVal = (Sigma(double (i),xc,r,mu)*Sigma(double(j),yc,r,mu) );
			//valFunct +=  ( valAtPixel*log(functVal) ) + ( (1 - valAtPixel)*log(1 - functVal) );  
				
			//mData += 1;

		}
	}

	//valFunct -= Lambda*pow(r,rPow);

	//cout << " val function " << " " <<  valFunct << endl; 

	//valFunct *= (1.0/mData);
	return valFunct; 
}



Point4D functions::ObjFunction4DGrad(double xc, double yc, double r, double mu,double Lambda, int rPow){

	//int kmin = 30;
	//int kmax = 60;
	Point4D valFunct;

	double mData = 0;

	double gradFactor = 0.0;

	for (int i = kmin; i<= kmax; i++)
	{
		for (int j = kmin; j<= kmax; j++)
		{
			double ThisPixelValue = m_mImage.at<double>(i,j);


			// Image MINUS window
			gradFactor = 2.0*
				( ThisPixelValue - ( Sigma(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu) ) );

			// Log-loss image MINUS window:

			//gradFactor = mu*
			//	( ThisPixelValue - ( Sigma(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu) ) );



			// with respect to Xc:
			
			// Image TIMES window
			//valFunct[0] += ThisPixelValue*SigmaDerXc(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu);

			// Image MINUS window

			valFunct[0] += gradFactor*SigmaDerXc(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu);
			

			// Log-loss image MINUS window

			//valFunct[0] += gradFactor*(sFunction(double(i),xc,r,mu) -  sFunction(xc,double(i),r,mu))/
			//	(1.0 - Sigma(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu));


			
			// with respect to Yc:

			
			// Image TIMES window

			//valFunct[1] += ThisPixelValue*Sigma(double(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu);
			
			
			// Image MINUS window

			valFunct[1] +=  gradFactor*Sigma(double(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu);
			

			// Log-loss image MINUS window

			//valFunct[1] += gradFactor*(sFunction(double(j),yc,r,mu) -  sFunction(yc,double(j),r,mu))/
			//	(1.0 - Sigma(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu));

			
			// with respect to R


			// Image TIMES window
			//valFunct[2] += ThisPixelValue*( (SigmaDerR(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu)) + 
			//	(Sigma(double(i),xc,r,mu)*SigmaDerR(double(j),yc,r,mu)) ); 

			// Image MINUS window

			valFunct[2] += gradFactor*
				( (SigmaDerR(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu)) + 
				(Sigma(double(i),xc,r,mu)*SigmaDerR(double(j),yc,r,mu)) ); 



			// Log-loss image MINUS window

			//valFunct[2] += gradFactor*
			//							(4.0 - sFunction(double(i),xc,r,mu) -
			//										 sFunction(xc,double(i),r,mu) - 
			//										 sFunction(double(j),yc,r,mu) -
			//										 sFunction(yc,double(j),r,mu)	)/
			//											(1.0 - ( Sigma(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu) ) );

			// HACK
			//valFunct[2] = 0;


			// with respect to mu
			valFunct[3] += ThisPixelValue*( (SigmaDerMu(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu)) + 
																			(SigmaDerMu(double(j),yc,r,mu)*Sigma(double(i),xc,r,mu) )); 

			// HACK
			valFunct[3] = 0;



			//mData += 1;

			// for j ends
		}

		// for i ends
	}

	// adding term coming from regularization...
	valFunct[2] -= rPow*Lambda*pow(r , rPow - 1);

	// HACK
	//valFunct[2] = 0.0;

	//cout << " " << " Value function in gradient  " << sum_ << endl;


	//for (int k = 0; k<=3; k++)
	//{
	//	valFunct[k] *= (1.0/mData);
	//}

	return valFunct;
}


matrix<double,4,4> functions::ObjFunction4DHessian(double xc, double yc,
	double r, double mu, double Lambda, int rPow){

		matrix<double,4,4> HessVal;

		for (int i=0; i<=3 ; i++)
		{
			HessVal(i,i) = 0.0;
			for (int j=i+1; j<=3; j++)
			{
				HessVal(i,j) = HessVal(i,j) = 0.0;
			}
		}

		for (int i = kmin; i<= kmin; i++)
		{
			for (int j = kmin; j<= kmax; j++)
			{
				// Second partial derivative with respect to Xc: XX
				HessVal(0,0) += valImage*
					SigmaDer2Xc(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu);

				// Second partial derivative with respect to Yc: YY
				HessVal(1,1) += valImage*
					Sigma(double(i),xc,r,mu)*SigmaDer2Xc(double(j),yc,r,mu);

				// Second partial derivative with respect to R: RR
				HessVal(2,2) += valImage*
					( 
					(SigmaDer2R(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu)) +
					(Sigma(double(i),xc,r,mu)*SigmaDer2R(double(j),yc,r,mu)) +
					(2.0*SigmaDerR(double(i),xc,r,mu)*SigmaDerR(double(j),yc,r,mu))
					);


				// Second partial derivative with respect to Mu: MuMu
				HessVal(3,3) += valImage*
					(
					(SigmaDer2Mu(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu)) +
					(SigmaDerMu(double(i),xc,r,mu)*SigmaDerMu(double(j),yc,r,mu)) +
					(Sigma(double(i),xc,r,mu)*SigmaDer2Mu(double(j),yc,r,mu)) 
					);

				// Second partial derivative with respect to Xc and Yc: XY
				HessVal(0,1) += valImage*
					SigmaDerXc(double(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu);

				HessVal(1,0) = HessVal(0,1);

				// Second partial derivative with respect to Xc and R: XR
				HessVal(0,2) += valImage*
					( ( SigmaDer2XcR(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu) ) +
					( SigmaDerXc(int(i),xc,r,mu)*SigmaDerR(double(j),yc,r,mu) ) );

				HessVal(2,0) = HessVal(0,2);



				// Second partial derivative with respect to Xc and Mu: XMu
				HessVal(0,3) += valImage*
					(
					(SigmaDer2XcMu(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu)) +
					(SigmaDerXc(double(i),xc,r,mu)*SigmaDerMu(double(j),yc,r,mu))
					);

				HessVal(3,0) = HessVal(0,3);



				// Second partial derivative with respect to Yc and R: YR
				HessVal(1,2) += valImage*
					( ( Sigma(double(i),xc,r,mu)*SigmaDer2Xc(double(j),yc,r,mu) ) +
					( SigmaDerR(int(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu) ) );

				HessVal(2,1) = HessVal(1,2);


				// Second partial derivative with respect to Yc and Mu: YMu
				HessVal(1,3) += valImage*
					(
					(Sigma(double(i),xc,r,mu)*SigmaDer2XcMu(double(j),yc,r,mu)) +
					(SigmaDerMu(double(i),xc,r,mu)*SigmaDerXc(double(j),yc,r,mu))
					);

				HessVal(3,1) = HessVal(1,3);


				// Second partial derivative with respect to R and Mu: RMu
				HessVal(2,3) += valImage*
					(
					(SigmaDer2RMu(double(i),xc,r,mu)*Sigma(double(j),yc,r,mu)) +
					(SigmaDerR(double(i),xc,r,mu)*SigmaDerMu(double(j),yc,r,mu)) +
					(Sigma(double(i),xc,r,mu)*SigmaDer2RMu(double(j),yc,r,mu)) +
					(SigmaDerMu(double(i),xc,r,mu)*SigmaDerR(double(j),yc,r,mu))
					);

				HessVal(3,2) = HessVal(2,3);

			}
		}

		HessVal(2,2) -= rPow*(rPow - 1)*Lambda*pow(r,rPow-2);

		return HessVal;
}

void functions::SetImage(cv::Mat mImage)
{
	m_mImage = mImage;
}

