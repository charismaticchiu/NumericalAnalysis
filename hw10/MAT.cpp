//hw10 Spline Interpolations
//100060007 Ming-Chang Chiu
//Date: 2015.4.12
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
