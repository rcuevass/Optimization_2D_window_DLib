#include "Utilities.h"
#include <math.h>
#include <iostream>
using namespace std;


Utilities::Utilities(void)
{
}


Utilities::~Utilities(void)
{
}




//////////////////////////////////
//			Point2D			Begin
Point2D::Point2D() {
	m_data[0] = 0.0;
	m_data[1] = 0.0;
}

Point2D::Point2D(double x, double y) { 
	m_data[0] = x;
	m_data[1] = y;
}

Point2D::Point2D(const Point2D& other) {
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
}


Point2D& Point2D::operator +=(const Point2D& other) {
	m_data[0] += other.m_data[0];
	m_data[1] += other.m_data[1];
	return *this;
}






Point2D& Point2D::operator /=(const double & other) {
	m_data[0] /= other;
	m_data[1] /= other;
	return *this;
}

Point2D& Point2D::operator =(const Point2D& other) {
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
	return *this;
}

double Point2D::SqrDistanceTo(const Point2D& other){
	return (m_data[0]-other.m_data[0])*(m_data[0]-other.m_data[0])+
		(m_data[1]-other.m_data[1])*(m_data[1]-other.m_data[1]);
}

// added by RCS


double Point2D::Abs(){
	return sqrt((m_data[0]*m_data[0]) + (m_data[1]*m_data[1]));
}


double& Point2D::operator[](int i) {
	return m_data[i];
}

double Point2D::operator[](int i) const {
	return m_data[i];
}

Point2D operator +(const Point2D & u, const Point2D & v)
{
	return Point2D (u[0]+v[0], u[1]+v[1]);
}

Point2D operator -(const Point2D & u, const Point2D & v)
{
	return Point2D (u[0]-v[0], u[1]-v[1]);
}

void Point2D::normalize() 
{
	double denom = sqrt(m_data[0]*m_data[0]+m_data[1]*m_data[1]);
	if (denom<1e-10) return;
	m_data[0]/=denom;
	m_data[1]/=denom;
}




//////////////////////////////////
//			Point3D			Begin

Point3D::Point3D() {
	m_data[0] = 0.0;
	m_data[1] = 0.0;
	m_data[2] = 0.0;
}

Point3D::Point3D(double x, double y, double z) { 
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}

Point3D::Point3D(const Point3D& other) {
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
	m_data[2] = other.m_data[2];
}

Point3D::Point3D(const Point2D& other) {
	m_data[0] = other[0];
	m_data[1] = other[1];
	m_data[2] = 1.0;
}

Point3D::Point3D(const Point2D& other, double z) {
	m_data[0] = other[0];
	m_data[1] = other[1];
	m_data[2] = z;
}

Point3D& Point3D::operator =(const Point3D& other) {
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
	m_data[2] = other.m_data[2];
	return *this;
}

Point3D& Point3D::operator +=(const Point3D& other) {
	m_data[0] += other.m_data[0];
	m_data[1] += other.m_data[1];
	m_data[2] += other.m_data[2];
	return *this;
}

//Point3D& Point3D::operator *=(const double other) {
//	m_data[0] *= other;
//	m_data[1] *= other;
//	m_data[2] *= other;
//	return *this;
//}

double& Point3D::operator[](int i) {
	return m_data[i];
}

double Point3D::operator[](int i) const {
	return m_data[i];
}

double Point3D::dot(const Point3D& other) const
{
	return m_data[0]*other.m_data[0] + 
		m_data[1]*other.m_data[1] + 
		m_data[2]*other.m_data[2];
}

double Point3D::SqrDistanceTo(const Point3D& other) const{
	return (m_data[0]-other.m_data[0])*(m_data[0]-other.m_data[0])+
		(m_data[1]-other.m_data[1])*(m_data[1]-other.m_data[1])+
		(m_data[2]-other.m_data[2])*(m_data[2]-other.m_data[2]);
}


double Point3D::length() const
{
	return sqrt(dot(*this));
}

//Point3D Point3D::operator +(const Point3D& u) const
//{
//	return Point3D(m_data[0]+v[0], m_data[1]+v[1], m_data[2]+v[2]);
//}




//////////////////////////////////
//			Point4D			Begin

Point4D::Point4D() {
	m_data[0] = 0.0;
	m_data[1] = 0.0;
	m_data[2] = 0.0;
	m_data[3] = 0.0;
}

Point4D::Point4D(double x, double y, double z, double w) { 
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
	m_data[3] = w;
}

Point4D::Point4D(const Point4D& other) {
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
	m_data[2] = other.m_data[2];
	m_data[3] = other.m_data[3];
}

Point4D::Point4D(const Point3D& other) {
	m_data[0] = other[0];
	m_data[1] = other[1];
	m_data[2] = other[2];
	m_data[3] = 1.0;
}

Point4D::Point4D(const Point3D& other, double w) {
	m_data[0] = other[0];
	m_data[1] = other[1];
	m_data[2] = other[2];
	m_data[3] = w;
}

Point4D& Point4D::operator =(const Point4D& other) {
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
	m_data[2] = other.m_data[2];
	m_data[3] = other.m_data[3];
	return *this;
}

Point4D& Point4D::operator +=(const Point4D& other) {
	m_data[0] += other.m_data[0];
	m_data[1] += other.m_data[1];
	m_data[2] += other.m_data[2];
	m_data[3] += other.m_data[3];
	return *this;
}

//Point3D& Point3D::operator *=(const double other) {
//	m_data[0] *= other;
//	m_data[1] *= other;
//	m_data[2] *= other;
//	return *this;
//}

double& Point4D::operator[](int i) {
	return m_data[i];
}

double Point4D::operator[](int i) const {
	return m_data[i];
}

double Point4D::dot(const Point4D& other) const
{
	return m_data[0]*other.m_data[0] + 
		m_data[1]*other.m_data[1] + 
		m_data[2]*other.m_data[2] +
		m_data[3]*other.m_data[3];
}

double Point4D::SqrDistanceTo(const Point4D& other) const{
	return (m_data[0]-other.m_data[0])*(m_data[0]-other.m_data[0])+
		(m_data[1]-other.m_data[1])*(m_data[1]-other.m_data[1])+
		(m_data[2]-other.m_data[2])*(m_data[2]-other.m_data[2])+
		(m_data[3]-other.m_data[3])*(m_data[3]-other.m_data[3]);
}


double Point4D::length() const
{
	return sqrt(dot(*this));
}



// Vector3D begin


Vector3D::Vector3D() {
	m_data[0] = 0.0;
	m_data[1] = 0.0;
	m_data[2] = 0.0;
}

Vector3D::Vector3D(double x, double y, double z) { 
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}


Vector3D::Vector3D(const Vector3D& other) {
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
	m_data[2] = other.m_data[2];
}


Vector3D::Vector3D(const Point3D& other) {
	m_data[0] = other[0];
	m_data[1] = other[1];
	m_data[2] = other[2];
}

Vector3D& Vector3D::operator =(const Vector3D& other) {
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
	m_data[2] = other.m_data[2];
	return *this;
}

Vector3D& Vector3D::operator +=(const Vector3D& other) {
	m_data[0] += other.m_data[0];
	m_data[1] += other.m_data[1];
	m_data[2] += other.m_data[2];
	return *this;
}

Vector3D& Vector3D::operator /=(const double other) {
	m_data[0] /= other;
	m_data[1] /= other;
	m_data[2] /= other;
	return *this;
}

Vector3D& Vector3D::operator *=(const double other) {
	m_data[0] *= other;
	m_data[1] *= other;
	m_data[2] *= other;
	return *this;
}

Point3D Vector3D::AsPoint(){
	return Point3D(m_data[0],m_data[1],m_data[2]);
}

double& Vector3D::operator[](int i) {
	return m_data[i];
}
double Vector3D::operator[](int i) const {
	return m_data[i];
}

double Vector3D::length() const
{
	return sqrt(dot(*this));
}

double Vector3D::lengthSquare() const
{
	return (dot(*this));
}


//double Vector3D::normalize() 
//{
//	double len=length();
//	(*this)=(1.0/len)*(*this);
//	return len;
//}


//Vector3D Vector3D::normalized() const
//{
//	double len=length();
//	return Vector3D((1.0/len)*(*this));
//}


double Vector3D::dot(const Vector3D& other) const
{
	return m_data[0]*other.m_data[0] + 
		m_data[1]*other.m_data[1] + 
		m_data[2]*other.m_data[2];
}

Vector3D Vector3D::elementwise_product(const Vector3D& other) const
{
	return Vector3D(
		m_data[0]*other[0],
		m_data[1]*other[1],
		m_data[2]*other[2]);
}


Vector3D Vector3D::cross(const Vector3D& other) const
{
	return Vector3D(
		m_data[1]*other[2] - m_data[2]*other[1],
		m_data[2]*other[0] - m_data[0]*other[2],
		m_data[0]*other[1] - m_data[1]*other[0]);
}





// Operators on points and vectors
Point2D operator *(double s, const Point2D& v); 
Point2D operator +(const Point2D& u, const Point2D& v); 
Point2D operator -(const Point2D& u, const Point2D& v); 
Vector3D operator *(double s, const Vector3D& v); 
Point3D operator *(double s, const Point3D& v); 
Vector3D operator +(const Vector3D& u, const Vector3D& v); 
Point3D operator +(const Point3D& u, const Point3D& v); 
Point3D operator -(const Point3D& u, const Point3D& v); 
Vector3D operator -(const Vector3D& u, const Vector3D& v); 
Vector3D operator -(const Vector3D& u); 
Point3D operator -(const Point3D& u, const Vector3D& v); 
Vector3D cross(const Vector3D& u, const Vector3D& v); 
std::ostream& operator <<(std::ostream& o, const Point3D& p); 
std::ostream& operator <<(std::ostream& o, const Vector3D& v); 
