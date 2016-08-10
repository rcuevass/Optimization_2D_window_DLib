
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



// Added by RCS
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






// Added by RCS
//////////////////////////////////
//			Point5D			Begin
class Point5D {
public:
	Point5D(); 
	Point5D(double x, double y, double z, double w, double u);  
	Point5D(const Point5D& other); 
	Point5D(const Point4D& other);
	Point5D(const Point4D& other, double w);

	Point5D& operator =(const Point5D& other); 
	Point5D& operator +=(const Point5D& other); 

	double& operator[](int i); 
	double operator[](int i) const; 

	double length() const; 
	double dot(const Point5D& other) const;
	double SqrDistanceTo(const Point5D& other) const;
private:
	double m_data[5];
};







// Added by RCS
// QuaternionVecs!!!
template<class T = double>
class QuaternionVec
{
public:
	T w, x, y, z;

	// Numerical constructor
	QuaternionVec(const T &w, const T &x, const T &y, const T &z): w(w), x(x), y(y), z(z) {};
	QuaternionVec(const T &x, const T &y, const T &z): w(T()), x(x), y(y), z(z) {}; // For 3-rotations
	QuaternionVec(const T &r): w(r), x(T()), y(T()), z(T()) {};
	QuaternionVec(): w(T()), x(T()), y(T()), z(T()) {};

	// Copy constructor and assignment
	QuaternionVec(const QuaternionVec &q): w(q.w), x(q.x), y(q.y), z(q.z) {};
	QuaternionVec& operator=(const QuaternionVec &q) { w=q.w; x=q.x; y=q.y; z=q.z; return *this; }

	// Unary operators
	QuaternionVec operator-() const { return QuaternionVec(-w, -x, -y, -z); }
	QuaternionVec operator~() const { return QuaternionVec(w, -x, -y, -z); } // Conjugate
	





	// Norm-squared. SQRT would have to be made generic to be used here
	T normSquared() const { return w*w + x*x + y*y + z*z; }


	// In-place operators
	QuaternionVec& operator+=(const T &r) 
	{ w += r; return *this; }
	QuaternionVec& operator+=(const QuaternionVec &q) 
	{ w += q.w; x += q.x; y += q.y; z += q.z; return *this; }

	QuaternionVec& operator-=(const T &r) 
	{ w -= r; return *this; }
	QuaternionVec& operator-=(const QuaternionVec &q) 
	{ w -= q.w; x -= q.x; y -= q.y; z -= q.z; return *this; }

	QuaternionVec& operator*=(const T &r) 
	{ w *= r; x *= r; y *= r; z *= r; return *this; }
	QuaternionVec& operator*=(const QuaternionVec &q) 
	{ 
		T oldW(w), oldX(x), oldY(y), oldZ(z);
		w = oldW*q.w - oldX*q.x - oldY*q.y - oldZ*q.z;
		x = oldW*q.x + oldX*q.w + oldY*q.z - oldZ*q.y;
		y = oldW*q.y + oldY*q.w + oldZ*q.x - oldX*q.z;
		z = oldW*q.z + oldZ*q.w + oldX*q.y - oldY*q.x;
		return *this;
	}

	QuaternionVec& operator/=(const T &r) 
	{ w /= r; x /= r; y /= r; z /= r; return *this; }
	QuaternionVec& operator/=(const QuaternionVec &q) 
	{ 
		T oldW(w), oldX(x), oldY(y), oldZ(z), n(q.normSquared());
		w = (oldW*q.w + oldX*q.x + oldY*q.y + oldZ*q.z) / n;
		x = (oldX*q.w - oldW*q.x + oldY*q.z - oldZ*q.y) / n;
		y = (oldY*q.w - oldW*q.y + oldZ*q.x - oldX*q.z) / n;
		z = (oldZ*q.w - oldW*q.z + oldX*q.y - oldY*q.x) / n;
		return *this;
	}

	// Binary operators based on in-place operators
	QuaternionVec operator+(const T &r) const { return QuaternionVec(*this) += r; }
	QuaternionVec operator+(const QuaternionVec &q) const { return QuaternionVec(*this) += q; }
	QuaternionVec operator-(const T &r) const { return QuaternionVec(*this) -= r; }
	QuaternionVec operator-(const QuaternionVec &q) const { return QuaternionVec(*this) -= q; }
	QuaternionVec operator*(const T &r) const { return QuaternionVec(*this) *= r; }
	QuaternionVec operator*(const QuaternionVec &q) const { return QuaternionVec(*this) *= q; }
	QuaternionVec operator/(const T &r) const { return QuaternionVec(*this) /= r; }
	QuaternionVec operator/(const QuaternionVec &q) const { return QuaternionVec(*this) /= q; }

	// Comparison operators, as much as they make sense
	bool operator==(const QuaternionVec &q) const 
	{ return (w == q.w) && (x == q.x) && (y == q.y) && (z == q.z); }
	bool operator!=(const QuaternionVec &q) const { return !operator==(q); }

	// The operators above allow QuaternionVec op real. These allow real op QuaternionVec.
	// Uses the above where appropriate.
	template<class T> friend QuaternionVec<T> operator+(const T &r, const QuaternionVec<T> &q);
	template<class T> friend QuaternionVec<T> operator-(const T &r, const QuaternionVec<T> &q);
	template<class T> friend QuaternionVec<T> operator*(const T &r, const QuaternionVec<T> &q);
	template<class T> friend QuaternionVec<T> operator/(const T &r, const QuaternionVec<T> &q);

	// Allows cout << q 
	template<class T> friend ostream& operator<<(ostream &io, const QuaternionVec<T> &q);
};

// Friend functions need to be outside the actual class definition
template<class T>
QuaternionVec<T> operator+(const T &r, const QuaternionVec<T> &q) 
{ return q+r; }

template<class T>
QuaternionVec<T> operator-(const T &r, const QuaternionVec<T> &q)
{ return QuaternionVec<T>(r-q.w, q.x, q.y, q.z); }

template<class T>
QuaternionVec<T> operator*(const T &r, const QuaternionVec<T> &q) 
{ return q*r; }

template<class T>
QuaternionVec<T> operator/(const T &r, const QuaternionVec<T> &q)
{
	T n(q.normSquared());
	return QuaternionVec(r*q.w/n, -r*q.x/n, -r*q.y/n, -r*q.z/n);
}

template<class T>
ostream& operator<<(ostream &io, const QuaternionVec<T> &q)
{ 
	io << q.w;
	(q.x < T()) ? (io << " - " << (-q.x) << "i") : (io << " + " << q.x << "i");
	(q.y < T()) ? (io << " - " << (-q.y) << "j") : (io << " + " << q.y << "j");
	(q.z < T()) ? (io << " - " << (-q.z) << "k") : (io << " + " << q.z << "k");
	return io;

}






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

	double angle(const Vector3D& other) const; 
	Vector3D cross(const Vector3D& other) const; 
	Vector3D elementwise_product(const Vector3D& other) const; 
	Point3D AsPoint();
	void EpsilonCorrect(const Vector3D& v);		//!< if the vector is almost equal to the origin substitute it with v

private:
	double m_data[3];
};
//			Vector3D		End
//////////////////////////////////



