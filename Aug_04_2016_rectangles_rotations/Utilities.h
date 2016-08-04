
#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <set>
#include <iomanip>
#include <sstream>
//#include <opencv2/core/core.hpp>
using namespace std;


class Utilities
{
public:
	Utilities(void);
	~Utilities(void);
};


//////////////////////////////////
//			Point2D			Begin
class Point2D {
public:
	Point2D(); 
	Point2D(double x, double y);  
	Point2D(const Point2D& other); 
	double SqrDistanceTo(const Point2D& other);

	double Abs();

	//	Point2D operator +(const Point2D& p, const Point2D& q); 
	Point2D& operator =(const Point2D& other); 
	Point2D& operator +=(const Point2D& other);
	Point2D& operator /=(const double & other);




	void normalize();

	double& operator[](int i); 
	double operator[](int i) const; 

private:
	double m_data[2];
};
//			Point2D			End
//////////////////////////////////


//////////////////////////////////
//			Point3D			Begin
class Point3D {
public:
	Point3D(); 
	Point3D(double x, double y, double z);  
	Point3D(const Point3D& other); 
	Point3D(const Point2D& other);
	Point3D(const Point2D& other, double z);

	Point3D& operator =(const Point3D& other); 
	Point3D& operator +=(const Point3D& other); 

	double& operator[](int i); 
	double operator[](int i) const; 

	double length() const; 
	double dot(const Point3D& other) const;
	double SqrDistanceTo(const Point3D& other) const;
private:
	double m_data[3];
};



// Added by Rogelio Cuevas
//////////////////////////////////
//			Point4D			Begin
class Point4D {
public:
	Point4D(); 
	Point4D(double x, double y, double z, double w);  
	Point4D(const Point4D& other); 
	Point4D(const Point3D& other);
	Point4D(const Point3D& other, double w);

	Point4D& operator =(const Point4D& other); 
	Point4D& operator +=(const Point4D& other); 

	double& operator[](int i); 
	double operator[](int i) const; 

	double length() const; 
	double dot(const Point4D& other) const;
	double SqrDistanceTo(const Point4D& other) const;
private:
	double m_data[4];
};



//////////////////////////////////
//			Vector3D		Begin
class Vector3D {
public:
	Vector3D(); 
	Vector3D(double x, double y, double z); 
	Vector3D(const Vector3D& other); 
	Vector3D(const Point3D& other);

	Vector3D& operator =(const Vector3D& other); 
	Vector3D& operator +=(const Vector3D& other); 
	Vector3D& operator /=(const double other);
	Vector3D& operator *=(const double other);

	double& operator[](int i);  
	double operator[](int i) const;  

	double length() const; 
	double lengthSquare() const;
	double normalize();
	Vector3D normalized() const;
	double dot(const Vector3D& other) const; 
	Vector3D cross(const Vector3D& other) const; 
	Vector3D elementwise_product(const Vector3D& other) const; 
	Point3D AsPoint();
	void EpsilonCorrect(const Vector3D& v);		//!< if the vector is almost equal to the origin substitute it with v

private:
	double m_data[3];
};
//			Vector3D		End
//////////////////////////////////



