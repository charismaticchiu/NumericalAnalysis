/*
Numerical Analysis HW10
VEC.cpp
u100060011      陳增鴻

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "VEC.h"
#include "Complex.h"
VEC::VEC(int n){
	dim = n;
	val = (Complex *)calloc(n,sizeof(Complex));
}
VEC::VEC(const VEC &v1)
{
	dim = v1.dim;
	
	double x,y;
	val = (Complex *)calloc(dim,sizeof(Complex));
	for (int i = 0; i < dim ; i++){
		x = v1.val[i].r();
		y = v1.val[i].i();
		Complex r(x,y);
		val[i] = r;
	}
}

VEC::~VEC()
{
	free(val);
}
int VEC::len()
{
	return dim;
}
VEC VEC::operator -()
{
	double x,y;
	for (int i = 0 ; i<dim ; i++){
		x = -val[i].r();
		y = -val[i].i();
		Complex r(x,y);
		val[i] = r;
	}
	return *this;
}
VEC &VEC::operator=(const VEC v1){
	dim = v1.dim;
	for(int i = 0; i < dim;i++){
		val[i] = v1.val[i];
	}
	return *this;
}
VEC &VEC::operator+=(const VEC v1)
{
	double x,y;
	for (int i = 0; i <dim ; i++){
		x += v1.val[i].r();
		y += v1.val[i].i();
		Complex r(x,y);
		val[i] = r;
	}
	return *this;
}
VEC &VEC::operator-=(const VEC v1)
{
	double x,y;
	for (int i = 0; i < dim ; i++){
		x -= v1.val[i].r();
		y -= v1.val[i].i();
		Complex r(x,y);
		val[i] = r;
	}
	return *this;
}
VEC &VEC::operator*=(double d1){
	for (int i = 0 ; i < dim ; i++){
		val[i]*=d1;
	}
	return *this;
}
VEC &VEC::operator*=(Complex c1){
	for (int i = 0 ; i < dim ; i++){
		val[i]*=c1;
	}
	return *this;
}
VEC &VEC::operator/=(double d1){
	for (int i = 0 ; i < dim ; i++){
		val[i]/=d1;
	}
	return *this;
}
VEC &VEC::operator/=(Complex c1){
	for (int i = 0 ; i < dim ; i++){
		val[i]/=c1;
	}
	return *this;
}
VEC VEC::operator+(const VEC v1)
{
	VEC s(*this);
	for(int i = 0 ; i < dim ; i++) s.val[i]+=v1.val[i];
	return s;
}
VEC VEC::operator-(const VEC v1){
	VEC s(*this);
	for (int i = 0; i < dim ; i++) s.val[i]-=v1.val[i];
	return s;
}

VEC VEC::operator*(double d1){
	VEC s(*this);
	for(int i = 0 ; i < dim ; i++){
		s.val[i] *= d1;
	}
	return s;
}
VEC VEC::operator*(Complex c1){
	VEC s(*this);
	for(int i = 0 ; i < dim ; i++){
		s.val[i] *= c1;
	}
	return s;
}
VEC VEC::operator/(double d1){
	VEC s(*this);
	for (int i = 0; i < dim ; i++){
		s.val[i] /= d1;
	}
	return s;
}
VEC VEC::operator/(Complex c1){
	VEC s(*this);
	for (int i = 0; i < dim ; i++){
		s.val[i] /= c1;
	}
	return s;
}
Complex &VEC::operator[](int n){
	if (n<0) n = 0;
	else if (n>=dim) n = dim -1;
	return val[n];
}
VEC *newVEC(int n){
	VEC *vptr;
	vptr=(VEC *)malloc(sizeof(VEC));
	vptr->dim = n;
	vptr->val=(Complex *)calloc(n,sizeof(Complex));
	return vptr;
}
VEC operator*(double a ,const VEC v1){
	VEC s(v1);
	for(int i = 0 ; i < v1.dim ; i++ ) s[i]*=a;
	return s;
}
VEC operator*(Complex a ,const VEC v1){
	VEC s(v1);
	for(int i = 0 ; i < v1.dim ; i++ ) s[i]*=a;
	return s;
}


