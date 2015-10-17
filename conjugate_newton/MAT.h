/*
Numerical Analysis HW10
MAT.h
u100060011      陳增鴻

*/
#ifndef MAT_H
#define MAT_H
#include "VEC.h"
#include "Complex.h"
class MAT{
private:
	int n;
	VEC **va;
public:
	MAT(int dim);
	MAT(const MAT &m1);
	MAT(int dim,Complex *v);
	~MAT();
	int dim();
	MAT tpose();
	MAT &operator-();
	MAT &operator=(MAT m1);
	MAT &operator+=(MAT &m1);
	MAT &operator-=(MAT &m1);
	MAT &operator*=(double a);
	MAT &operator/=(double a);
	MAT operator+(MAT m1);
	MAT operator-(MAT m1);
	MAT operator*(MAT m1);
	VEC &operator[](int m);
	VEC operator*(VEC v1);
	MAT operator*(double a);
	MAT operator/(double a);
	friend MAT operator*(double a,MAT &m1);
	friend VEC operator*(VEC &v1,MAT &m1);
};
/*
MAT &luFact(MAT &m1);
VEC fwdSubs(MAT &m1, VEC b);	//for lu
VEC fwdSubs1(MAT &m1, VEC b);	//for cholesky
VEC bckSubs(MAT &m1, VEC b);
MAT operator*(double a,const MAT &m1);
VEC operator*(VEC &v1,MAT &m1);
MAT &cholesky(MAT &A);
VEC choSolve(MAT &L, VEC b);
int jacobi(MAT &A,VEC b,VEC &x,int maxIter,double tol);
int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol);
int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol);
int cg(MAT &A,VEC b,VEC &x,int maxIter,double tol);
double norm2(VEC x);
double power(MAT A,VEC q,int maxiter,double eps,double &v);
double invpower(MAT A,VEC q,int maxiter,double eps,double &v,int aw);
void QRDecomp(MAT &A,MAT &Q, MAT &R);
double ConvErr(MAT A);
int EVqr(MAT &A,double tol,int maxiter);
int EVqrShifted(MAT &A,double mu, double tol,int maxiter);
double Lagrange(double x,VEC &XDATA,VEC &YDATA);
void splineM(int N,VEC &X,VEC &Y,VEC &M);
double spline(double x,int N,VEC &X,VEC &Y,VEC &M);
*/
VEC polyRoots(Complex x,VEC &A,int maxiter, double eps);
#endif
