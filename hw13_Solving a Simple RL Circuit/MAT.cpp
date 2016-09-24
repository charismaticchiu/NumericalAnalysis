//hw11 Numerical Integrations
//100060007 Ming-Chang Chiu
//Date: 2015.5.15
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MAT.h"
#include <algorithm>
#include <cmath>
MAT::MAT(int dim)
{
    n=dim;
    va=(VEC **)malloc(n*sizeof(VEC*));
    for (int i=0;i<n;i++)
	va[i]=newVEC(n);
}
MAT::MAT(const MAT &m1)
{
    VEC **vsrc=m1.va; //to get around not indexing const MAT
    n=m1.n;
    va=(VEC **)malloc(n*sizeof(VEC*));
    for(int i=0;i<n;i++){
	va[i]=newVEC(n);
	(*va[i])=(*vsrc[i]);
    }
}
MAT::MAT(int dim,double *v)
{
    n=dim;
    va=(VEC **)malloc(n*sizeof(VEC*));
    for(int i=0;i<n;i++){
	va[i]=newVEC(n);
	for(int j=0;j<n;j++){
	    (*va[i])[j]= *(v++); //array indexing + VEC indexing
	}
    }
}
MAT::~MAT()
{
    for (int i=n-1;i>=0;i--) free(va[i]);
    free(va);
}
int MAT::dim()
{
    return n;
}
MAT MAT::tpose()
{
    MAT mnew(n);
    for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
	    mnew[i][j]=(*va[j])[i];
	}
    }
    return mnew;
}
MAT & MAT::operator-() {// Need no returning?
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
	    (*va[i])[j]=-(*va[i])[j];
}
MAT &MAT::operator=(MAT m1)
{
    for (int i=0; i<n; i++)
        (*va[i])=m1[i];
    return *this;
}
MAT &MAT::operator+=(MAT &m1)
{
    for (int i=0; i<n; i++)
        (*va[i])+=m1[i];
    return *this;
}
MAT &MAT::operator-=(MAT &m1)
{
    for (int i=0; i<n; i++)
        (*va[i])-=m1[i];
    return *this;
}
MAT &MAT::operator*=(double a)
{
    for (int i=0; i<n; i++)
        (*va[i])*=a;
    return *this;
}
MAT &MAT::operator/=(double a)
{
    for (int i=0; i<n; i++)
        (*va[i])/=a;
    return *this;
}
MAT MAT::operator+(MAT m1)
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])+m1[i];
    return s;
}
MAT MAT::operator-(MAT m1)
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])-m1[i];
    return s;
}
MAT MAT::operator*(MAT m1)
{
    MAT z(n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            z[i][j]=0;
	    for (int k=0; k<n; k++)
		z[i][j]+=((*va[i])[k]*m1[k][j]);
	}
    }
    return z;
}
VEC &MAT::operator[](int m)
{
    return *va[m];
}
VEC MAT::operator*(VEC v1)
{
    VEC s(n);
    for (int i=0; i<n; i++) {
        s[i]=(*va[i])*v1;
    }
    return s;
}
MAT MAT::operator*(double a)
{
    MAT z(n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            z[i][j]=(*va[i])[j] * a;//why use (*va[i])[j] rather than va[i][j];
        }
    }
    return z;    
}
MAT MAT::operator/(double a)
{
    MAT z(n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            z[i][j]=(*va[i])[j] / a;
        }
    }
    return z;
}
MAT operator*(double a, MAT &m1)
{
    MAT v2(m1.n);
    for (int i=0; i<m1.n; i++) {
        for (int j=0; j<m1.n; j++) {
            v2[i][j]= a * m1[i][j];
        }
    }
    return v2;
}
VEC operator*(VEC &v1,MAT &m1)
{
    VEC v2(m1.n);
    for (int i=0; i<m1.n; i++) {
        v2[i]=0;
        for (int j=0; j<m1.n; j++) {
            v2[i] += v1[j]*m1[j][i];
        }
    }
    return v2;
}
MAT &luFact(MAT &m1)
{
    int i,j,k;
    for(i=0;i<m1.dim();i++){
	//copy m1[i][j] to u[i][j] needs no action due to in-place LU
	for(j=i+1; j<m1.dim(); j++)// form l[j][i]
	    m1[j][i] /= m1[i][i];
	for(j=i+1; j<m1.dim(); j++){//update lower submatrix
	    for(k=i+1; k<m1.dim(); k++){
		m1[j][k] -= m1[j][i]*m1[i][k];
	    }
	}
    }
    return m1;
}
VEC fwdSubs(MAT &m1,VEC b)//fwdSubs must be modified from the one used by LU Decomposition
{
    VEC y = b;
    //y[0]=b[0]/m1[0][0];
    for(int i=0; i<m1.dim();i++){
      //y[i]/=m1[i][i];//Add this step for modification about Cholesky method
		for(int j=i+1; j<m1.dim(); j++)
			y[j] -= m1[j][i]*y[i];
    }
    return y;
}
VEC bckSubs(MAT &m1,VEC b)//b as y
{
    int i,j,k;
    VEC x=b;
    for(i=m1.dim()-1 ;i>=0 ;i--){
		x[i] /= m1[i][i];
		for(j=i-1; j>=0; j--)
			x[j] -= m1[j][i] * x[i];
    }
    return x;
}
MAT &cholesky(MAT &A)
{ 
    int i,j,k;
    for(i=0; i<A.dim(); i++){
		A[i][i] = sqrt(A[i][i]);
		for(j=i+1; j<A.dim(); j++){
			A[j][i] /= A[i][i];
		}
		for(j=i+1;j<A.dim();j++){
			for(k=i+1; k<=j;k++)
			A[j][k] -= A[j][i]*A[k][i];
		}
		for(k=i+1;k<A.dim();k++) A[i][k]=0; // since our Cholesky just gives the lower triangle and doesn't change the val of upper triangle, we have to clear that by hand   
    }
    return A;
}
VEC choSolve(MAT &L,VEC b)
{
    VEC x(b);

	//forward substitutions
	for (int i=0; i<x.len(); i++)
	{	
		x[i] /= L[i][i];
		for (int j=i+1; j<L.dim(); j++)
			x[j] -= L[j][i]*x[i];
		
	}
	
	//backward substitutions
	for (int i=(x.len()-1); i>=0; i--)
	{
		x[i] /= L[i][i];
		for (int j=i-1; j>=0; j--)
			x[j] -= L[i][j]*x[i];
	}	

	return x;
} 
double l2norm(VEC x){
    double sum = 0.0;
    for(int i = 0; i< x.len(); i++){
	sum += (x[i]*x[i]);
    }
    
    return sqrt(sum);
}
double l1norm(VEC x){
    double sum = 0.0;
    for (int i = 0; i < x.len(); i++){
	sum += std::abs(x[i]);
    }
    return sum;
}
double linfnorm(VEC x){
    double max=0.0;
    for (int i = 0 ;i< x.len();i++)
      max = std::max(std::abs(x[i]),max);
    return max;
    }
int jacobi(MAT &A,VEC b,VEC &x,int maxIter,double tol)
{
    int k = 0;
    double theta=0.0;
    MAT LU(A.dim()),a = A;
    VEC newx(A.dim()),ans(A.dim()),temp(A.dim());

    LU = luFact(a);       //LU in-place decomposition
    temp=fwdSubs(LU,b);   //forward substitution
    ans= bckSubs(LU,temp);//backward substitution
   
    while( k < maxIter)
    {      
		x = newx;
			for(int i=0; i < A.dim(); i++)
		{
			theta = 0.0;
				for(int j=0; j<A.dim();j++){
			if(i!=j){
				theta += (A[i][j]*x[j]);
			}
			}
			newx[i]= (b[i]-theta)/A[i][i];
			}
			k++;
		
		if (linfnorm(newx - x) < tol){
				break;
		}
    }    
    printf("Diffrence w.r.t. hw4: %E\n",linfnorm(newx-ans));
   
    return k;
}
int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol)
{
    int k = 0;
    double theta=0.0;
    VEC newx(A.dim()),ans(A.dim()),temp(A.dim());
    MAT LU(A.dim()),a = A;

    LU = luFact(a);       //LU in-place decomposition
    temp=fwdSubs(LU,b);   //forward substitution
    ans= bckSubs(LU,temp);//backward substitution
   
    while( k < maxIter)
    {
		x = newx;
		for(int i=0; i< A.dim(); i++)
		{
			theta = 0.0;
			for(int j=0; j < A.dim();j++){
				if(i!=j){
					theta += (A[i][j]*newx[j]);
				}
			}
			newx[i]= (b[i]-theta)/A[i][i];
		}
		k++;        
		if (linfnorm(newx - x) < tol)
			break;
    }
    printf("Diffrence w.r.t. hw4: %E\n",linfnorm(newx-ans));
    
    return k;
}
int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol)
{
    int k = 0;
    double theta=0.0;    
    MAT LU(A.dim()),a=A;
    VEC newx(A.dim()),ans(A.dim()),temp(A.dim());
    
    LU = luFact(a);       //LU in-place decomposition
    temp=fwdSubs(LU,b);   //forward substitution
    ans= bckSubs(LU,temp);//backward substitution
    
    while(k < maxIter)
    {
		x = newx;
		//forward
        for(int i=0; i< A.dim(); i++)
		{
			theta = 0.0;
			for(int j=0; j<A.dim();j++){
				if(i!=j){
					theta += (A[i][j]*newx[j]);
				}
			}
			newx[i]= (b[i]-theta)/A[i][i];
		}
		//backward
		for(int i=(A.dim()-1); i>= 0; i--)
		{	
			theta = 0.0;
			for(int j = (A.dim()-1); j >= 0; j--){
				if(i!=j){
					theta += (A[i][j]*newx[j]);
				}
			}
			newx[i]= (b[i]-theta)/A[i][i];
		}
        k++;        
        if (linfnorm(newx - x) < tol )
            break;
    }
    printf("Difference w.r.t. hw4: %E\n",linfnorm(newx-ans));
        
    return k;
}
int cg(MAT &A,VEC b,VEC &x,int maxIter, double tol)
{	
	//A_p: A * p
	VEC  A_p(A.dim());
    int node = std::sqrt(A.dim());
    //p: search direction; 
    VEC p(A.dim()), r(A.dim()), newr(A.dim()), newx(A.dim());//,ans(A.dim()),temp(A.dim());
    //MAT LU(A.dim()),a = A;
    //r2: rT * r
    double alpha, beta, r2, newr2, err;//,g;
    //g = (double)(node-1)/2000.0;
    int iter = 0;//
    /*
    LU = luFact(a);       //LU in-place decomposition
    temp=fwdSubs(LU,b);   //forward substitution
    ans= bckSubs(LU,temp);//backward substitution
    */

    
    //Initial condition for CGM
    p = r = b - (A * x);
    r2 = r * r;
   
    while(iter < maxIter)
    {
	A_p = A * p;
	alpha = r2 / (p * A_p);
	newx = x + (alpha * p);
	err = std::sqrt(r2/A.dim());
	if ( err < tol )
            break;

	newr = r - (alpha * A_p);
	newr2 = newr * newr;
	beta = newr2 / r2;
	p = newr + beta * p;


	//Re-initialization for next iteration
	x = newx;
	r = newr;
	r2 = newr2;
	//////////////////////////////////////
	iter++;		
    }
    //printf("cg:Vne: %lf; Vsw: %lf; Vse: %lf; R:%lf\n",newx[node-1],newx[(node-1)*node],newx[node*node-1], std::abs( 1/(g*(newx[0]-newx[1])+g*(newx[0]-newx[node] )) ) );
    //printf("LU:Vne: %lf; Vsw: %lf; Vse: %lf; R:%lf\n",ans[node-1],ans[(node-1)*node],ans[node*node-1], std::abs( 1/(g*(ans[0]-ans[1])+g*(ans[0]-ans[node] )) )  );
    //printf("Difference w.r.t. hw4: %E\n",linfnorm(ans - newx));
    return iter;
}
MAT inverse(MAT &A)
{
    VEC temp(A.dim());
    // I: identity matrix; inv: inverse of A; LU: LU in-place
    MAT LU(A.dim()), a = A, I(A.dim()), inv(A.dim());
    LU = luFact(a);       //LU in-place decomposition
    // Form identity matrix
    for (int i =0; i < A.dim(); i++){
	I[i][i] = 1;
    }
    // solve inverse of A. Answer is stored row by row
    for (int i = 0; i< A.dim(); i++){	
	temp=fwdSubs(LU,I[i]);   //forward substitution
	inv[i]= bckSubs(LU,temp);//backward substitution
    }
    // we need the result of inverse() column by column, so return transpose matrix of inv
    return inv.tpose();
}
double invPowerShift(MAT &A, double omega, VEC q, int maxIter){
    int k = 1;
    //I: identity matrix; LU: LU in-place; 
    MAT I(A.dim()), LU(A.dim()), shiftedA(A.dim());
    //q: 
    VEC z(A.dim()), temp(A.dim());// q(A.dim())
    
    double lambdaN = 100.0, newlambdaN=0.0;
    //Form identity matrix
    for (int i =0; i < A.dim(); i++){
	I[i][i] = 1;
    }
    shiftedA = A-(omega*I);
    LU = luFact(shiftedA);       //LU in-place decomposition
    
    while( (k < maxIter) && (std::abs(newlambdaN - lambdaN) >= 0.000000001) ){
		lambdaN = newlambdaN;	
		temp = fwdSubs(LU,q);   //forward substitution
		z = bckSubs(LU,temp);//backward substitution
		q = z/l2norm(z);
		newlambdaN = q * A * q;

    }
    printf("Iteration: %d\n",k);
    return newlambdaN;
}
double powerMethod(MAT &A, VEC q, int maxIter){
	VEC z(A.dim()), r(A.dim()), u(A.dim()), w(A.dim());
	MAT A_k(A.dim());
	double lambda1, tol= 1.0;
	//double new_lambda
	int iter = 0;
	
	while( (iter < maxIter) && (tol >= 0.000000001) ){
		z = A * q;
		iter++;
		q = z/l2norm(z);
		lambda1 = q * A * q;
		r = A * q - lambda1 * q;
		u = q * A;
		w = u/l2norm(u);
		tol = l2norm(r)/std::abs(w*q);	
	}
	return lambda1;
}
void QRDecomp(MAT &A,MAT &Q, MAT &R){
    int n=A.dim();
    MAT AT(A.tpose());
    MAT QT(Q.tpose());
    VEC S(n);                   //sum vector
    R[0][0]=sqrt(AT[0]*AT[0]);  //initial R
    QT[0]=AT[0]/R[0][0];                //initial Q

    for(int j=1;j<n;j++){
        for(int i=0;i<n;i++){
            S[i] = 0;
        }       //initialization of sum vector
        for(int i=0;i<=j;i++){
                R[i][j] = QT[i]*AT[j];
        }
        for(int i=0;i<=j-1;i++){
                S = S + R[i][j]*QT[i];          //do the summation
        }
        S = AT[j] - S;
        R[j][j]=sqrt(S*S);
        QT[j] = S/R[j][j];
    }
    Q=QT.tpose();

	
	
	/*
    for(int j=1;j<n;j++){
        
        for(int i=0;i<=j;i++){
                R[i][j] = QT[i]*AT[j];
        }
        for(int i=0;i<=j-1;i++){
                AT[j] = AT[j] - (R[i][j]*QT[i]);          //do the summation
        }
        
        R[j][j]=std::sqrt(AT[j]*AT[j]);
        QT[j] = AT[j]/R[j][j];
    }
    */

}	

int EVqr(MAT &A,double tol,int maxiter){
    int iter = 0 ;
    MAT T(A);
    MAT R(A.dim()), Q(A.dim());
    //since we are not able to access column of a matrix easily, so we access row of matrices as substitute.
    
    do{
	//QR decomposition
	QRDecomp(T,Q,R);	
	T = R * Q;
	
	iter++;
	
    }while(  (iter < maxiter) && (QRerror(T) > tol)  );
    
    
    printf("Largest 3 eigenvalues: ");
    for (int i = 0; i < 3; i++)
	printf(" %lf;",T[i][i]);
    printf("\nSmallest 3 eigenvalues: ");
    for (int i = A.dim()-1; i > A.dim()-4; i--)
	printf(" %lf;",T[i][i]);
    //printf("\niter %d\nerror: %E\n",iter,QRerror(A_k));
    return iter;

}
int EVqrShifted(MAT &A,double mu,double tol,int maxiter){
    int iter = 0 ;
    MAT T(A);
    MAT R(A.dim()), Q(A.dim()),muI(A.dim());
    
    //since we are not able to access column of a matrix easily, so we access row of matrices as substitute.
    for (int i =0 ; i < A.dim(); i++)
		muI[i][i] = mu;

    do{
	for(int i=0;i<A.dim();i++){           //subtract mu*I before QR factorization
		T[i][i] = T[i][i] - mu;
	}
	
	//QR decomposition
	QRDecomp(T,Q,R);	
	T = R * Q;
	//QR decomposition ends		
	for(int i=0;i<A.dim();i++){           //add mu*I back
		T[i][i] = T[i][i] + mu;
    }

	iter++;	
	
    }while(  (iter < maxiter) && (QRerror(T) > tol)  );
        
    printf("Largest 3 eigenvalues: ");
    for (int i = 0; i < 3; i++)
	printf(" %lf;",T[i][i]);
    printf("\nSmallest 3 eigenvalues: ");
    for (int i = A.dim()-1; i > A.dim()-4; i--)
	printf(" %lf;",T[i][i]);
    //printf("\niter %d\nerror: %E\n",iter,QRerror(A_k));
    return iter;

}
double QRerror(MAT &A){
    double max = std::abs(A[1][0]);
    for (int i = 2; i< A.dim(); i++)
	max = std::max(std::abs(A[i][i-1]),max);
    return max;
}
double Lagrange(double x,VEC &XDATA,VEC &YDATA){
    int n = XDATA.len();
    double NS[n];
    int i,j,k;
    for (i=0; i<n; i++) NS[i]=YDATA[i];
    for (k=1; k<n; k++) {
	for (j=0; j<n-k; j++) {
	    NS[j]=((x-XDATA[j])*NS[j+1]-(x-XDATA[k+j])*NS[j])
		/(XDATA[j+k]-XDATA[j]);
	}
    }
    return NS[0];
}
void splineM(int N,VEC &X,VEC &Y,VEC &M){// generate spline momentum M
    	//Zero boundary moment
    	MAT s(N+1),LU(N+1);
    	VEC d(N+1),temp(N+1);
	//li: lambda_i,hi1:h_(i+1)
	double ui,li,hi,hi1;
	for(int i = 1; i< N; i++){// i = 1~n-1
		hi = X[i] - X[i-1];
		hi1= X[i+1]-X[i];
		ui = hi/(hi+hi1);
		li = hi1/(hi+hi1);
		d[i] = (6/(hi+hi1))*( ((Y[i+1]-Y[i])/hi1) - ((Y[i]-Y[i-1])/hi));
		s[i][i-1] = ui;
		s[i][i] = 2.0;
		s[i][i+1]= li;
	}
	s[0][0] = s[N][N] = 2.0;
	LU = luFact(s);       //LU in-place decomposition
    	temp=fwdSubs(LU,d);   //forward substitution
    	M = bckSubs(LU,temp);//backward substitution
}

double spline(double x,int N,VEC &X,VEC &Y,VEC &M){ // spline interp at x
    int I=0;//Interval number
    double h,ans; //X_i - X_i_1
    for (int i = 1; i < N+1; i++){
		if ((X[i-1]<=x) && (x <= X[i])){
			I = i;
			break;
		}	    
    }
    h = X[I] - X[I-1];
    splineM(N, X, Y, M);
    ans =( M[I-1]*pow(X[I]-x,3)/(6*h) ) + (M[I]*pow(x-X[I-1],3)/(6*h))+ (((Y[I] - Y[I-1])/h) - (h*(M[I]-M[I-1])/6) ) * (x-X[I-1])+ Y[I-1] - (M[I-1]*pow(h,2)/6);
    return ans;
    
}
double integ(VEC &X,VEC &Y,int n){ // composite nth order Newton-Cotes integral
    //iter: m in the lecture notes
    //len()-1 is the interval numbers
    // m*n should be exactly the number of intervals
    int iter = (X.len()-1)/n;
    double space = X[X.len()] - X[0];
    double h = space/(X.len()-1);
    double sum = 0.0;
    //w can acutally be pre-calculated 
    VEC w(n+1);
    if(n == 1)
	w[0] = w[1] = 0.5;
    else if (n == 2)
	w[0] = w[2] =  1./3. , w[1] = 4./3.;
    else if (n == 3)
	w[0] = w[3] =  3./8. , w[1] = w[2] = 9./8.;
    else if (n == 4)
	w[0] = w[4] =  14./45. , w[1] = w[3] = 64./45. , w[2] = 24./45.;
    else if (n == 5)
	w[0] = w[5] =  95./288. , w[1] = w[4] = 375./288. , w[2] = w[3] = 250./288.;
    else if (n == 6)
	w[0] = w[6] =  41./140. , w[1] = w[5] = 216./140. , w[2] = w[4] = 27./140. , w[3] = 272./140.;
    
    for(int i = 0; i< iter; i++){
	for(int j = 0; j<= n; j++){
	    sum += (w[j] * Y[ n*i + j ]);
	    //printf("%d\n",n*i+j);
	}
    }
    return sum * h;
  
}  

double newton(double ini,double limit, double eps){
    int k = 0;
    double err = 1. + eps;
    double x = ini;
    double newx = 0.0;
    while(err > eps){
	newx = x - ( (x-pow(log(x),2)-0.9) / (1- 2* log(x)/x)) ;// modify this line as f(x) is chaging
	printf("newx: %lf\n",newx);
	if(newx > x + limit) newx = x + limit;
	else if(newx<x - limit)newx = x -limit;
	printf("nodified newx: %lf\n",newx);
	k++;
	x = newx;
	err = fabs((x-pow(log(x),2)-0.9));
	printf("err: %lf\n\n",err);
    }
    return x;
  
} 
MAT genJ1(VEC x){
    MAT m(x.len());
    double r ;//r is a dependent var of voltage, so it is specific to each resistor
    int N = 3;//3 nodes per side
    
    int i,j,k = N * N;
    double alpha = 0.1;
    for(int p = 0; p < N; p++){
	for(int q = 0; q < N; q++){
	    //If the node we are traversing is the upper left sub-network
	    if(p < N -1 && q < N-1){
		//indexing i,j
		i = N * p + q;// node we are traversing
		j = i + 1;// node to the right of i
		r = 1 + alpha * fabs(x[i] - x[j]);
		m[i][i] += 1/r;
		m[j][j] += 1/r;
		m[i][j] -= 1/r;
		m[j][i] -= 1/r;
		
		j = (p+1)*N + q;// node below i
		r = 1 + alpha * fabs(x[i] - x[j]);
		m[i][i] += 1/r;
		m[j][j] += 1/r;
		m[i][j] -= 1/r;
		m[j][i] -= 1/r;
		
	    }
	    //If the node we are traversing is the right side boundary node
	    //, excluding bottom one
	    else if(p < N-1 && q == N-1){
		i = N * p + q;
		j = (p+1)*N + q;// node below i
		r = 1 + alpha * fabs(x[i] - x[j]);
		m[i][i] += 1/r;
		m[j][j] += 1/r;
		m[i][j] -= 1/r;
		m[j][i] -= 1/r;
		
	    }
	    //If the node is at bottom boundary side,
	    //excluding lowest right node
	    if(p ==N - 1 && q < N - 1){
		i = N*p+q;
		j = i + 1;// node to the right of i
		r = 1 + alpha * fabs(x[i] - x[j]);
		m[i][i] += 1/r;
		m[j][j] += 1/r;
		m[i][j] -= 1/r;
		m[j][i] -= 1/r;		
	    }		
	}
    }
    //Make row of Jacobian of ground node and voltage node be correct
    for(int p = 0; p < x.len(); p++){
	m[(N-1)*N + (N-1)/2 ][p] = m[0][p] = 0.0; 
    }
    m[(N-1)*N + (N-1)/2][(N-1)*N + (N-1)/2] = m[0][0] = 1.;

    
    return m;
}

// index of f follow the order in hw12.pdf
VEC genF1(VEC x, double v){
    VEC f(x.len());//f: the system values
    VEC r(x.len()+12);//r the resistors' values
    r[9] = 1 + 0.1 * fabs(x[0] - x[1]);
    r[10] = 1 + 0.1 * fabs(x[0] - x[3]);
    r[11] = 1 + 0.1 * fabs(x[1] - x[2]);
    r[12] = 1 + 0.1 * fabs(x[1] - x[4]);
    r[13] = 1 + 0.1 * fabs(x[2] - x[5]);
    r[14] = 1 + 0.1 * fabs(x[3] - x[4]);
    r[15] = 1 + 0.1 * fabs(x[3] - x[6]);
    r[16] = 1 + 0.1 * fabs(x[4] - x[5]);
    r[17] = 1 + 0.1 * fabs(x[4] - x[7]);
    r[18] = 1 + 0.1 * fabs(x[5] - x[8]);
    r[19] = 1 + 0.1 * fabs(x[6] - x[7]);
    r[20] = 1 + 0.1 * fabs(x[7] - x[8]);
    //Kichkoff's current equations
    f[0] = x[0] - v;
    f[1] = (x[1] - x[0]) /r[9]  + (x[1] - x[4]) /r[12]+ (x[1] - x[2]) /r[11];
    f[2] = (x[2] - x[1]) /r[11] + (x[2] - x[5]) /r[13];
    f[3] = (x[3] - x[0]) /r[10] + (x[3] - x[4]) /r[14]+ (x[3] - x[6]) /r[15];
    f[4] = (x[4] - x[1]) /r[12] + (x[4] - x[3]) /r[14]+ (x[4] - x[5]) /r[16] + (x[4] - x[7]) /r[17];
    f[5] = (x[5] - x[2]) /r[13] + (x[5] - x[4]) /r[16]+ (x[5] - x[8]) /r[18];
    f[6] = (x[6] - x[3]) /r[15] + (x[6] - x[7]) /r[19];
    //    f[7] = (x[7] - x[6]) /r[19] + (x[7] - x[4]) /r[17]+ (x[7] - x[8]) /r[20];
    f[7] = 0.0;//ground node
    f[8] = (x[8] - x[7]) /r[20] + (x[8] - x[5]) /r[18];
    return f;
}


//Index of nodes follow the order described in hw12a.pdf
MAT genJ2(VEC x){
    MAT m(x.len());//create Jacobian matrix
    int N = 3;//3 nodes per side
    
    int i,j,k = N * N;//k: row number for each resistor
    double b= 1.;
    double r ;
    for(int p = 0; p < N; p++){
	for(int q = 0; q < N; q++){
	    //Enter this condition if traversing upper left nodes
	    if(p < N -1 && q < N-1){
		i = N * p + q;//indexing node i, the traversed node
		j = i + 1;    //indexing node j, the right hand node
		r = 1 + x[k]; //update resistance----r = r0 + kappa*T
	       	m[i][i] += 1/r;
		m[j][j] += 1/r;
		m[i][j] -= 1/r;
		m[j][i] -= 1/r;
		m[i][k] = (x[i] - x[j])/(r*r);
		m[k][i] = -2*b*(x[i] - x[j])/r;
		m[k][j] = 2*b*(x[i] - x[j])/r;
		m[k][k] = 1 + b*(x[i]-x[j])*(x[i] - x[j])/(r*r);
		k++;

		j = (p+1)*N + q;//indexing j, the node below i
		r = 1 + x[k];
		m[i][i] += 1/r;
		m[j][j] += 1/r;
		m[i][j] -= 1/r;
		m[j][i] -= 1/r;
		m[i][k] = (x[i] - x[j])/(r*r);
		m[k][i] = -2*b*(x[i] - x[j])/r;
		m[k][j] = 2*b*(x[i] - x[j])/r;
		m[k][k] = 1 + b*(x[i]-x[j])*(x[i] - x[j])/(r*r);
		k++;
	    }
	    //Enter if traversing right side boundary nodes
	    //, excluding the bottom one
	    else if(p < N-1 && q == N-1){
		i = N * p + q;
		j = (p+1)*N + q;//indexing j, the node below i
		r = 1 + x[k];
		m[i][i] += 1/r;
		m[j][j] += 1/r;
		m[i][j] -= 1/r;
		m[j][i] -= 1/r;
		m[i][k] = (x[i] - x[j])/(r*r);
		m[k][i] = -2*b*(x[i] - x[j])/r;
		m[k][j] = 2*b*(x[i] - x[j])/r;
		m[k][k] = 1 + b*(x[i]-x[j])*(x[i] - x[j])/(r*r);
		k++;
	    }
	    //Enter if traversing bottom boundary nodes
	    //, excluding the lowest right node.
	    if(p ==N - 1 && q < N - 1){
		i = N*p+q;//
		j = i + 1;//indexing j, the node to the right of i
		r = 1 + x[k];
		m[i][i] += 1/r;
		m[j][j] += 1/r;
		m[i][j] -= 1/r;
		m[j][i] -= 1/r;
		m[i][k] = (x[i] - x[j])/(r*r);
		m[k][i] = -2*b*(x[i] - x[j])/r;
		m[k][j] = 2*b*(x[i] - x[j])/r;
		m[k][k] = 1 + b*(x[i]-x[j])*(x[i] - x[j])/(r*r);
		k++;
	    }
		
	}
    }
    
    //Make rows of Jacobian of ground node and voltage node be correct
    for(int p = 0; p < x.len(); p++){
	m[(N-1)*N + (N-1)/2 ][p] = m[0][p] = 0.0; 
    }
    m[(N-1)*N + (N-1)/2][(N-1)*N + (N-1)/2] = m[0][0] = 1.;

    
    return m;
}

//Index of f follows the order described in hw12a.pdf
VEC genF2(VEC x, double v){
    VEC f(x.len());
    VEC r(x.len());
    int resisNum = 12;
    for(int i = 0; i < resisNum; i++){
	r[i+9] = 1 + x[i+9];
    }
    //Kichkoff's current equations
    f[0] = x[0] - v;
    f[1] = (x[1] - x[0]) /r[9] + (x[1] - x[4]) /r[12]+ (x[1] - x[2]) /r[11];
    f[2] = (x[2] - x[1]) /r[11] + (x[2] - x[5]) /r[13];
    f[3] = (x[3] - x[0]) /r[10] + (x[3] - x[4]) /r[14]+ (x[3] - x[6]) /r[15];
    f[4] = (x[4] - x[1]) /r[12] + (x[4] - x[3]) /r[14]+ (x[4] - x[5]) /r[16] + + (x[4] - x[7]) /r[17];
    f[5] = (x[5] - x[2]) /r[13] + (x[5] - x[4]) /r[16]+ (x[5] - x[8]) /r[18];
    f[6] = (x[6] - x[3]) /r[15] + (x[6] - x[7]) /r[19];
    //    f[7] = (x[7] - x[6]) /r[19] + (x[7] - x[4]) /r[17]+ (x[7] - x[8]) /r[20];
    //Node 7's equation is not needed for it is gnd node.
    f[7] = 0.0;
    f[8] = (x[8] - x[7]) /r[20] + (x[8] - x[5]) /r[18];
    //Temperature's equation
    f[9]  = x[9]  - (x[0] - x[1])*(x[0] - x[1])/r[9];
    f[10] = x[10] - (x[0] - x[3])*(x[0] - x[3])/r[10];
    f[11] = x[11] - (x[1] - x[2])*(x[1] - x[2])/r[11];
    f[12] = x[12] - (x[1] - x[4])*(x[1] - x[4])/r[12];
    f[13] = x[13] - (x[2] - x[5])*(x[2] - x[5])/r[13];
    f[14] = x[14] - (x[3] - x[4])*(x[3] - x[4])/r[14];
    f[15] = x[15] - (x[3] - x[6])*(x[3] - x[6])/r[15];
    f[16] = x[16] - (x[4] - x[5])*(x[4] - x[5])/r[16];
    f[17] = x[17] - (x[4] - x[7])*(x[4] - x[7])/r[17];
    f[18] = x[18] - (x[5] - x[8])*(x[5] - x[8])/r[18];
    f[19] = x[19] - (x[6] - x[7])*(x[6] - x[7])/r[19];
    f[20] = x[20] - (x[7] - x[8])*(x[7] - x[8])/r[20];
    return f;
}
