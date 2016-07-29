#pragma once
#include "Utilities.h"
#include <vector>
//#include <dlib/matrix/matrix.h>

#include <dlib/matrix.h>
using namespace dlib;

class functions
{
public:
	functions(void);
	~functions(void);

	// sFunction and derivatives
	double sFunction(double x, double xc, double r, double mu); // done
	double sFunctionDerXc(double x, double xc, double r, double mu); // done
	double sFunctionDerR(double x, double xc, double r, double mu); // done
	double sFunctionDerMu(double x, double xc, double r, double mu); // done


	// Sigma (sigmoid-like) function and derivatives:
	// First derivatives //
	double Sigma(double x, double xc, double r, double mu); // done
	double SigmaDerXc(double x, double xc, double r, double mu); // done
	double SigmaDerR(double x, double xc, double r, double mu); // done
	double SigmaDerMu(double x, double xc, double r, double mu); // done

	// Second derivatives //
	double SigmaDer2Xc(double x, double xc, double r, double mu); // done 
	double SigmaDer2R(double x, double xc, double r, double mu); // done
	double SigmaDer2XcR(double x, double xc, double r, double mu); // done
	double SigmaDer2Mu(double x, double xc, double r, double mu); // done
	double SigmaDer2XcMu(double x, double xc, double r, double mu); // done
	double SigmaDer2RMu(double x, double xc, double r, double mu); // done
	


	// 2D Objective function (xC and yC), gradient and hessian 
	double ObjFunction2D(double xc, double yc, double r, double mu); // done
	Point2D ObjFunction2DGrad(double xc, double yc, double r, double mu); // done
	Point3D ObjFunction2DHessian(double xc, double yc, double r, double mu); // done


	// 3D Objective function (xC, yC and r), gradient and hessian
	double ObjFunction3D(double xc, double yc, double r, double mu); // done
	Point3D ObjFunction3DGrad(double xc, double yc, double r, double mu); // done
	
	
	// The next two functions were a useful but expensive attempt...
	Point3D ObjFunction3DHessianDiag(double xc, double yc, double r, double mu); // done
	Point3D ObjFunction3DHessianCross(double xc, double yc, double r, double mu); // done

	// ... this one is the definite and more efficient way!
	matrix<double,3,3> ObjFunction3DHessian(double xc, double yc, double r, double mu); // done



	// 4D objective function (xC, yC, r, mu), gradient and hessian
	double ObjFunction4D(double xc, double yc, double r, double mu); // done
	Point4D ObjFunction4DGrad(double xc, double yc, double r, double mu); // done
	matrix<double,4,4> ObjFunction4DHessian(double xc, double yc, double r, double mu); 

};




