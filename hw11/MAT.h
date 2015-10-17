//hw11 Numerical Integrations
//100060007 Ming-Chang Chiu
//Date: 2015.5.15
#ifndef MAT_h
#define MAT_h
#include "VEC.h"
class MAT {
    private:
    
    int n; // define nxn matrix
    VEC **va;// array of n pointers to vectors

    public:

    MAT(int dim); // uninit constructor
    MAT(const MAT &m1); // copy constructor
    MAT(int dim,double *v); // init constructor
    ~MAT();    // destructor
    int dim();  // return dimension of the matrix
    MAT tpose(); // transpose
    MAT &operator-(); // unary operator, negative value
    MAT &operator=(MAT m1); // assignment
    MAT &operator+=(MAT &m1);// m += m1;
    MAT &operator-=(MAT &m1); // m -= m1;
    MAT &operator*=(double a); // m *= dbl;
    MAT &operator/=(double a); //m /= dbl
    MAT operator+(MAT m1); //m1+m2
    MAT operator-(MAT m1);//m1-m2
    MAT operator*(MAT m1);//m1*m2
    VEC & operator[](int m); // m'th row
    VEC operator*(VEC v1); //mxv1
    MAT operator*(double a);//m*dbl
    MAT operator/(double a); //m/dbl
    friend MAT operator*(double a, MAT &m1); // dbl x m
    friend VEC operator*(VEC &v1,MAT &m1);  // vT x m  ￼

    friend MAT &luFact(MAT &m1);
    friend VEC fwdSubs(MAT &m1,VEC b);
    friend VEC bckSubs(MAT &m1,VEC b);
    friend MAT &cholesky(MAT &A);     // Cholesky decomposition
    friend VEC choSolve(MAT &L,VEC b);// forward and backward substitutions
    friend int jacobi(MAT &A,VEC b,VEC &x,int maxIter,double tol);
    friend int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol);
    friend int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol);
    friend double l2norm(VEC x);
    friend double l1norm(VEC x);
    friend double linfnorm(VEC x);
    friend int cg(MAT &A,VEC b,VEC &x,int maxIter, double tol);
    friend MAT inverse(MAT &A);
    friend double invPowerShift(MAT &A, double w, VEC q,int maxIter);
    friend double powerMethod(MAT &A,VEC q, int maxIter);
    friend int EVqr(MAT &A,double tol,int maxiter);
    friend int EVqrShifted(MAT &A,double mu,double tol,int maxiter);
    
};
MAT operator*(double a, MAT &m1);      // dbl x m
VEC operator*(VEC &v1,MAT &m1);             // vT x m
VEC fwdSubs(MAT &m1,VEC b);
VEC bckSubs(MAT &m1,VEC b);
MAT &cholesky(MAT &A);     // Cholesky decomposition
VEC choSolve(MAT &L,VEC b);// forward and backward substitutions
int jacobi(MAT &A,VEC b,VEC &x,int maxIter,double tol);
int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol);
int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol); // symmetric gauss-seide
double l2norm(VEC x);
double l1norm(VEC x);
double linfnorm(VEC x);
int cg(MAT &A,VEC b,VEC &x,int maxIter, double tol);
MAT inverse(MAT &A);
double invPowerShift(MAT &A, double w, VEC q, int maxIter);
double powerMethod(MAT &A, VEC q, int maxIter);
void QRDecomp(MAT &A,MAT &Q, MAT &R);
int EVqr(MAT &A,double tol,int maxiter);
int EVqrShifted(MAT &A,double mu,double tol,int maxiter);
double QRerror(MAT &A);
double Lagrange(double x,VEC &XDATA,VEC &YDATA);
void splineM(int N,VEC &X,VEC &Y,VEC &M);
 // generate spline momentum M
double spline(double x,int N,VEC &X,VEC &Y,VEC &M); // spline interp at x
double integ(VEC &X,VEC &Y,int n); // composite nth order Newton-Cotes integral
#endif


