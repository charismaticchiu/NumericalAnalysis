/*
Numerical Analysis HW10
VEC.h
u100060011      陳增鴻

*/
#ifndef VEC_H_
#define VEC_H_
#include "Complex.h"
class VEC{
private:
	int dim; // length of the vector
	Complex *val; // array of the value
public:
	VEC(int n);
	VEC(const VEC &v1);
	~VEC();
	int len();
	VEC operator-();
	VEC &operator=(const VEC v1);
	VEC &operator+=(const VEC v1);
	VEC &operator-=(const VEC v1);
	VEC &operator*=(double d1);
	VEC &operator*=(Complex c1);
	VEC &operator/=(double d1);
	VEC &operator/=(Complex c1);
	VEC operator+(const VEC v1);
	VEC operator-(const VEC v1);
	//double operator*(VEC v1); // inner product
	VEC operator*(double d1);
	VEC operator*(Complex d1);
	VEC operator/(double d1);
	VEC operator/(Complex d1);
	Complex &operator[](int n);
	friend VEC operator*(double a,const VEC v1);
	friend VEC operator*(Complex a,const VEC v1);
	friend VEC *newVEC(int n);
};
VEC operator *(double a, const VEC v1);
VEC operator *(Complex a, const VEC v1);
VEC *newVEC(int n);
#endif
