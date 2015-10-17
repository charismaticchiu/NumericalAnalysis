/*
Numerical Analysis HW10
MAT.cpp
u100060011      陳增鴻

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MAT.h"
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )		//compare a and b, if a>b return a, else return b
#endif
MAT::MAT(int dim)
{
	n = dim;
	va=(VEC **)malloc(n*sizeof(VEC*));
	for(int i = 0 ; i < n ; i++){
		va[i]=newVEC(n);
	}
}
MAT::MAT(const MAT &m1)
{
	VEC **vsrc=m1.va;
	n=m1.n;
	va=(VEC **)malloc(n*sizeof(VEC*));
	for(int i = 0 ; i < n; i++){
		va[i] = newVEC(n);
		(*va[i])=(*vsrc[i]);
	}	
}
MAT::MAT(int dim,double *v)
{
	n=dim;
	va = (VEC **)malloc(n*sizeof(VEC*));
	for(int i = 0 ; i < n ; i++){
		va[i] = newVEC(n);
		for(int j = 0; j < n ; j++){
			(*va[i])[j]=*(v++);
		}
	}
}
MAT::~MAT()
{
	for(int i = n-1;i>=0;i--) free(va[i]);
	free(va);
}
int MAT::dim()
{
	return n;
}
MAT MAT::tpose()
{
	MAT mnew(n);
	for(int i = 0; i < n ; i++){
		for(int j = 0; j<n ; j++){
			mnew[i][j]=(*va[j])[i];
		}
	}
	return mnew;
}
MAT & MAT::operator-()
{
	for (int i=0; i < n ; i++)
		for(int j=0; j < n ; j++)
			(*va[i])[j]=-(*va[i])[j];
		return *this;
}
MAT & MAT::operator=(MAT m1)
{
	for (int i = 0; i<n;i++)
		(*va[i])=m1[i];
	return *this;
}
MAT & MAT::operator+=(MAT &m1)
{
	for(int i=0;i<n;i++){
		(*va[i])+=m1[i];
	}
	return *this;
}
MAT & MAT::operator-=(MAT &m1)
{
	for(int i=0;i<n;i++)
		(*va[i])-=m1[i];
	return *this;
}
MAT & MAT::operator*=(double a)
{
	for(int i = 0 ; i < n ; i++)
		(*va[i])*=a;
	return *this;
}
MAT & MAT:: operator/=(double a)
{
	for(int i = 0; i < n; i++)
		(*va[i])/=a;
	return *this;
}
MAT MAT::operator+(MAT m1)
{
	MAT s(n);
	for(int i = 0;i<n;i++)
		s[i]=(*va[i])+m1[i];
	return s;
}
MAT MAT::operator-(MAT m1)
{
	MAT s(n);
	for(int i = 0;i<n;i++)
		s[i]=(*va[i])-m1[i];
	return s;
}
MAT MAT::operator*(MAT m1)
{
	MAT z(n);
	for(int i = 0; i < n;i++){
		for(int j = 0; j < n;j++){
			z[i][j]=0;
			for(int k=0;k<n;k++)
				z[i][j]+=((*va[i])[k]*m1[k][j]);
		}
	}
	return z;
}

VEC &MAT::operator[](int m)
{
	if(m<0)m=0;
	else if (m>n)m=n-1;
	return *va[m];
}
VEC MAT::operator*(VEC v1)
{
	VEC s(n);
	for(int i = 0; i<n ; i++){
		s[i]=(*va[i])*v1;
	}
	return s;
}
MAT MAT::operator*(double a)
{
	MAT s(*this);
	for(int i = 0 ; i < s.n ; i++)
		for(int j = 0 ; j < s.n;j++)
			s[i][j]*=a;
	return s;
}
MAT MAT::operator/(double a)
{
	MAT s(*this);
	for(int i = 0; i < s.n ; i++)
		for(int j = 0; j < s.n;j++)
			s[i][j]/=a;
	return s;
}
VEC operator*(VEC &v1,MAT &m1)
{
	VEC v2(m1.n);
	for(int i = 0; i <m1.dim();i++){
		v2[i]=0;
		for(int j = 0; j<m1.dim();j++){
			v2[i]+=v1[j]*m1[j][i];
		}
	}
	return v2;
}
MAT operator*(double a,const MAT &m1)
{
	MAT m2(m1);
	for(int i = 0 ; i < m2.dim();i++)	
		for(int j = 0; j < m2.dim() ; j++)
			m2[i][j]*=a;
	return m2;
}
MAT & luFact(MAT &m1)
{
	for(int i = 0 ; i < m1.dim() ; i++){	//form lj column
		for(int j = i+1 ; j < m1.dim() ; j++){
			m1[j][i] /= m1[i][i];
		}
		for(int j = i + 1; j < m1.dim() ; j++){ //update lower submatrix
			for(int k = i+1; k<m1.dim(); k++){
				m1[j][k] -= m1[j][i] * m1[i][k];
			}
		}
	}
	return m1;
}
VEC fwdSubs(MAT &m1, VEC b) {
        VEC y(b);
        for(int i = 0 ; i < m1.dim() ; i++) y[i] = b[i];
        for(int i = 0 ; i < m1.dim() ; i++){
		
                for(int j = i+1 ; j<m1.dim() ; j++){
                        y[j] -= m1[j][i] * y[i];
                }
        }
        return y;
}
VEC fwdSubs1(MAT &m1, VEC b) {
        VEC y(b);
        for(int i = 0 ; i < m1.dim() ; i++) y[i] = b[i];
        for(int i = 0 ; i < m1.dim() ; i++){
		y[i] /= m1[i][i];		//because the diagonal elements are not 1's, modify by adding this line
                for(int j = i+1 ; j<m1.dim() ; j++){
                        y[j] -= m1[j][i] * y[i];
                }
        }
        return y;
}
VEC bckSubs(MAT &m1,VEC y){
        VEC x(y);
        for(int i = 0 ; i < m1.dim() ; i++) x[i] = y[i];
        for(int i = m1.dim()-1 ; i >= 0 ; i--){
                x[i] /= m1[i][i];
                for(int j = i-1 ; j >= 0 ; j--){
                        x[j] -= m1[j][i] *x[i];
                }
        }
        return x;
}
MAT & cholesky(MAT &A){					//Cholesky Decomposition
	for(int i=0; i<A.dim(); i++){
		A[i][i] = sqrt(A[i][i]);
		for(int j=i+1; j<A.dim(); j++){
			A[j][i] /= A[i][i];
		}
		for(int j=i+1; j<A.dim(); j++){
			for(int k=i+1; k<=j; k++){
				A[j][k] -= A[j][i]*A[k][i];
			}
		}
		for(int j=i+1; j<A.dim(); j++){		//make A a lower triangular matrix
			A[i][j] = 0;
		}
	}
	return A;
}
VEC choSolve(MAT &L, VEC b){
	VEC y(b);
	y = fwdSubs1(L,b);				//do the forward substitution to solve y
	L = L.tpose();					//compute the transpose of L
	VEC x(y);
	x = bckSubs(L,y);				//do the backward substitution to solve x
	return x;
}

int jacobi(MAT &A,VEC b,VEC &x, int maxIter, double tol){
	int n = A.dim();
	VEC S(n);			//summation vector
	VEC G(n);			//diagonal of A
	VEC tmp(n);			//save the x before iterative method
	MAT L(n);			
	VEC cpr(n);			//the result of hw4
	//iterative
	int k=0;
	double err1;
	double err2;
	double errinf;
	double diff1,diff2,diffinf;
	while(k<maxIter){
		for(int i=0;i<n;i++){
			S[i] = 0;		//initialize the sum vector to 0
		}
		for(int i=0;i<n;i++){
			G[i] = A[i][i];		
			//save the previous vector x to tmp(for comparison between iteration k and k+1)
			tmp[i] = x[i]; 
			for(int j=0;j<n;j++){
				if(j!=i){
					S[i] = S[i] + A[i][j]*x[j];	//sum the Aij*xj, where i is not equal to j
				}
			}
		}
		for(int i=0;i<n;i++){
			x[i] = (1/G[i])*(b[i]-S[i]);		//compute the new solution of x at iteration k
		}
		k++;						//after each iteration, k=k+1
		err1 = 0;	//1-norm
		err2 = 0;	//2-norm
		errinf = 0;	//infinity-norm
		
		
		for(int i=0;i<n;i++){
			err1 += fabs(x[i] - tmp[i]);		//sum of absolute error
			err2 += pow(x[i] - tmp[i], 2);		//sum of square error
			errinf = max(errinf,fabs(x[i]-tmp[i]));	//take the maximum absolute difference as error
		}						//max(a,b) defined on the top
		err2 = pow(err2,0.5);				//find the square root of err2 and get the 2-norm
//you can change err1 to err2 or errinf to get different number of iterations k based on which p-norm you preferred
        		
		if(err1<tol){
			break;
		}
	}
	diff1 = 0;
        diff2 = 0;
        diffinf = 0;
	//the answer from hw4

        L = cholesky(A);                        //in-place Cholesky Decomposition
        cpr = choSolve(L,b);                      //solve the linear system by forward and backward substitution
        //use the same method to compute the error between hw4 and hw5, see if < 10^-7. if not, decrease the tol
        for(int i=0;i<n;i++){
        	diff1 += fabs(x[i] - cpr[i]);
                diff2 += pow(x[i] - cpr[i], 2);
                diffinf = max(diffinf,fabs(x[i]-cpr[i]));
        }
        diff2 = pow(diff2,0.5);
        printf("error between hw4 and hw5(1-norm 2-norm infinity-norm): %e %e %e\n",diff1,diff2,diffinf);
	return k;
}

int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol)
{
	int n = A.dim();
	VEC S(n);
	VEC G(n);
	VEC tmp(n);
	MAT L(n);
	VEC cpr(n);
	//iterative
	int k=0;
	double err1;
	double err2;
	double errinf;
	double diff1,diff2,diffinf;
	while(k<maxIter){
		for(int i=0;i<n;i++){
			S[i] = 0;
		}
		for(int i=0;i<n;i++){
			G[i] = A[i][i];
			tmp[i] = x[i];
			for(int j=0;j<n;j++){
				if(j!=i){
					S[i] = S[i] + A[i][j]*x[j];
				}
			}
			x[i] = (1/G[i])*(b[i]-S[i]);	
		//both x[i] and S[i] are computed in the same for loop, so some updated values are used in the current iteration
		}
		k++;
		err1 = 0;
		err2 = 0;
		errinf = 0;
		for(int i=0;i<n;i++){
			err1 += fabs(x[i] - tmp[i]);
			err2 += pow(x[i] - tmp[i], 2);
			errinf = max(errinf,fabs(x[i]-tmp[i]));
		}
		err2 = pow(err2,0.5);
	
		if(err1<tol){
			break;
		}
	}
	diff1 = 0;
        diff2 = 0;
        diffinf = 0;
	//the answer from hw4

        L = cholesky(A);                        //in-place Cholesky Decomposition
        cpr = choSolve(L,b);                      //solve the linear system by forward and backward substitution
        //use the same method to compute the error between hw4 and hw5, see if < 10^-7. if not, decrease the tol
        for(int i=0;i<n;i++){
        	diff1 += fabs(x[i] - cpr[i]);
                diff2 += pow(x[i] - cpr[i], 2);
                diffinf = max(diffinf,fabs(x[i]-cpr[i]));
        }
        diff2 = pow(diff2,0.5);
        printf("error between hw4 and hw5(1-norm 2-norm infinity-norm): %e %e %e\n",diff1,diff2,diffinf);

	return k;	
}

int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol)
{
        int n = A.dim();
        VEC S(n);
        VEC G(n);
        VEC tmp(n);
	MAT L(n);
	VEC cpr(n);
        //iterative
        int k=0;
        double err1;
        double err2;
        double errinf;
	double diff1,diff2,diffinf;
        while(k<maxIter){
                for(int i=0;i<n;i++){
                        S[i] = 0;
                }
		//forward gauss-seidel 
		//compute vector x from x[0] to x[n-1]
                for(int i=0;i<n;i++){
                        G[i] = A[i][i];
                        tmp[i] = x[i];
                        for(int j=0;j<n;j++){
                                if(j!=i){
                                        S[i] = S[i] + A[i][j]*x[j];
                                }
                        }
                        x[i] = (1/G[i])*(b[i]-S[i]);
              	}
		//backward gauss-seidel
		//use the reslults of forward gauss-seidel to compute vector x
		//direction : from x[n-1] to x[0]
		for(int i=0;i<n;i++){
                        S[i] = 0;	//we also need to initialize the sum vector.
                }
		for(int i=n-1;i>=0;i--){
                        for(int j=n-1;j>=0;j--){
                                if(j!=i){
                                        S[i] = S[i] + A[i][j]*x[j];
                                }
                        }
                        x[i] = (1/G[i])*(b[i]-S[i]);
                }
                k++;
                err1 = 0;
                err2 = 0;
                errinf = 0;
                for(int i=0;i<n;i++){
                        err1 += fabs(x[i] - tmp[i]);
                        err2 += pow(x[i] - tmp[i], 2);
                        errinf = max(errinf,fabs(x[i]-tmp[i]));
                }
                err2 = pow(err2,0.5);
               
		if(err1<tol){
			break;
		}

	}
	diff1 = 0;
        diff2 = 0;
        diffinf = 0;
	//the answer from hw4

        L = cholesky(A);                        //in-place Cholesky Decomposition
        cpr = choSolve(L,b);                      //solve the linear system by forward and backward substitution
        //use the same method to compute the error between hw4 and hw5, see if < 10^-7. if not, decrease the tol
        for(int i=0;i<n;i++){
        	diff1 += fabs(x[i] - cpr[i]);
                diff2 += pow(x[i] - cpr[i], 2);
                diffinf = max(diffinf,fabs(x[i]-cpr[i]));
        }
        diff2 = pow(diff2,0.5);
        printf("error between hw4 and hw5(1-norm 2-norm infinity-norm): %e %e %e\n",diff1,diff2,diffinf);
        return k;
}
int cg(MAT &A,VEC b,VEC &x,int maxIter,double tol){
	int n = A.dim();
	int k = 0;
	VEC x_next(n);
	VEC r_next(n);
	VEC p_next(n);
	VEC r(n);
	VEC p(n);
	MAT L(n);
    VEC cpr(n);
	double alpha,beta;
	double err;
	double diff;
	r = b - A*x;					//initial condition
	p = r;
	while(k<maxIter){				//conjugate gradient decent
		alpha = (r*r)/(p*(A*p));
		x_next = x + alpha*p;
		r_next = r - alpha*(A*p);
		beta = (r_next*r_next)/(r*r);
		p_next = r_next + beta*p;
		k++;
		x = x_next;					//assign to the x,r,p of the next iteration 
		r = r_next;
		p = p_next;
		err = pow((r*r)/n,0.5);
		if(err<tol)					//see if the error is smaller than the defined tolerance. if so, break out of the loop 
			break;
	}
        diff = 0;
        //the answer from hw4

        L = cholesky(A);                        //in-place Cholesky Decomposition
        cpr = choSolve(L,b);                      //solve the linear system by forward and backward substitution
        //use the same method to compute the error between hw4 and hw5, see if < 10^-7. if not, decrease the tol
        for(int i=0;i<n;i++){
                diff = max(diff,fabs(x[i]-cpr[i]));		//use infinity norm to compute the error
        }
        printf("error between hw4 and hw6(infinity-norm): %e\n",diff);

	return k;	
}


double norm2(VEC x){  //compute 2-norm
    int n=x.leng();
    double norm=0;
    for(int i=0;i<n;i++){
        norm += pow(x[i], 2);
    }
    norm = pow(norm, 0.5);
    return norm;
}

double power(MAT A,VEC q,int maxiter,double eps,double &v){
        double tol = 1 + eps;
        int k = 0;
        int n = A.dim();
        VEC z(n);
        VEC r(n);
        VEC Aq(n);
        VEC u(n);
        VEC w(n);
        while((tol>=eps) && (k<=maxiter)){

                Aq = A*q;
                k = k+1;
                q = Aq/norm2(Aq);
                Aq = A*q;
                v = q*Aq;
                r = Aq-v*q;
                u = Aq;
                w = u/norm2(u);
                tol = norm2(r)/fabs(q*w);

        }
        printf("accuracy of lambda max = %e\n",norm2(r));		//set 2-norm error of residue as accuracy term
        return k;
}
double invpower(MAT A,VEC q,int maxiter,double eps,double &v,int aw){
        double tol = 1+eps;
        int k=0;
        int n=A.dim();
        VEC y(n);
        VEC z(n);
        VEC r(n);
        VEC Aq(n);
        VEC u(n);
        VEC w(n);
        MAT A1(n);
        MAT S(n);
        for(int i=0;i<n;i++){
        	S[i][i] = aw;			//define the shifting matrix
        }
        A1 = A-S;
        luFact(A1);					//use LU to solve the equation
        while((tol>=eps) && (k<=maxiter)){
            k = k+1;
            y = fwdSubs(A1,q);		//forward subsitution
            z = bckSubs(A1,y);		//backward substitution
            q = z/norm2(z);
            Aq = A*q;
            v = q*Aq;
            r = Aq - v*q;
            u = Aq;
            w = u/norm2(u);
            tol = norm2(r)/fabs(q*w);
        }
        printf("accuracy of lambda min = %e\n",norm2(r));
        return k;
}

double ConvErr(MAT A){		//find the maximum among the elements just below the diagonal 
	double error = 0;
	for(int i=1;i<A.dim();i++){
		error = max(error,fabs(A[i][i-1]));
	}
	return error;
}

// A = Q * R  QR Decomposition
void QRDecomp(MAT &A,MAT &Q, MAT &R){
    int n=A.dim();
    MAT AT(A.tpose());
    MAT QT(Q.tpose());
    VEC S(n);			//sum vector
    R[0][0]=sqrt(AT[0]*AT[0]);	//initial R
    QT[0]=AT[0]/R[0][0];		//initial Q

    for(int j=1;j<n;j++){
    	for(int i=0;i<n;i++){
            S[i] = 0;
        }       //initialization of sum vector
        for(int i=0;i<=j;i++){
          	R[i][j] = QT[i]*AT[j];
        } 
        for(int i=0;i<=j-1;i++){ 
          	S = S + R[i][j]*QT[i];		//do the summation
        }
        S = AT[j] - S;  
        R[j][j]=sqrt(S*S);
        QT[j] = S/R[j][j];
    }
    Q=QT.tpose(); 
        
}
int EVqr(MAT &A,double tol,int maxiter){
    int n=A.dim();
    int k=0;
    double err=10;
    MAT Q(n),R(n),T(A); 
    while(err>tol && k<maxiter){
        QRDecomp(T,Q,R);  			//QR factorization
        T=R*Q;						//matrix multiplication
        err = ConvErr(T);			//get the error term in each iteration
        k++;
    }
        
    printf("the three largest EV:\n");
	printf("%lf %lf %lf\n",T[0][0],T[1][1],T[2][2]);
	printf("the three smallest EV:\n");
	printf("%lf %lf %lf\n",T[n-1][n-1],T[n-2][n-2],T[n-3][n-3]);
	printf("iter = %d\n",k);
	printf("error = %e\n",err);
    A=T;
    return k;
}
int EVqrShifted(MAT &A,double mu, double tol,int maxiter){
	int n=A.dim();
    int k=0;
    double err=10;
    MAT Q(n),R(n),T(A); 
    while(err>tol && k<maxiter){
    	for(int i=0;i<n;i++){		//subtract mu*I before QR factorization
    			T[i][i] = T[i][i] - mu;
    	}
        QRDecomp(T,Q,R);  			//QR factorization
        T=R*Q;						//matrix multiplication
        for(int i=0;i<n;i++){		//add mu*I back
    			T[i][i] = T[i][i] + mu;
    	}
        err = ConvErr(T);
        k++;
    }
        
    printf("the three largest EV:\n");
	printf("%lf %lf %lf\n",T[0][0],T[1][1],T[2][2]);
	printf("the three smallest EV:\n");
	printf("%lf %lf %lf\n",T[n-1][n-1],T[n-2][n-2],T[n-3][n-3]);
	printf("iter = %d\n",k);
	printf("error = %e\n",err);
    A=T;
    return k;
}

double Lagrange(double x,VEC &XDATA,VEC &YDATA){	//use non-recursive Neville’s algorithm to calculate Lagrange Interpolation
	int n=XDATA.leng();			//get the length of the vector XDATA
	VEC NS(n);
	for(int i=0;i<n;i++){
		NS[i] = YDATA[i];
	}
	for(int k=1;k<n;k++){
		for(int j=0;j<n-k;j++){
			NS[j] = ((x-XDATA[j])*NS[j+1]-(x-XDATA[k+j])*NS[j])/(XDATA[j+k]-XDATA[j]);
		}
	}
	return NS[0];
}
void splineM(int N,VEC &X,VEC &Y,VEC &M){
	M[0]=0;
	M[N]=0;
	VEC h(N);
	h[0] = 0;
	for(int i=1;i<N;i++){			// h is the length of the subinterval 
		h[i] = X[i] - X[i-1];
	}
	VEC MU(N+1);
	VEC LAMBDA(N);
	VEC d(N+1);
	MAT W(N+1);
	MAT L(N+1);
	LAMBDA[0] = 0;
	d[0] = 0;
	MU[N] = 0;
	d[N] = 0;
	for(int i=1;i<N;i++){			//prepare the variable for the tridiagonal matrix
		MU[i] = h[i]/(h[i]+h[i+1]);
		LAMBDA[i] = h[i+1]/(h[i]+h[i+1]);
		d[i] = (6/(h[i]+h[i+1]))*(((Y[i+1]-Y[i])/h[i+1])-((Y[i]-Y[i-1])/h[i]));
	}
	for(int i=0;i<N+1;i++){		//initialization
		for(int j=0;j<N+1;j++){
			W[i][j] = 0;
		}
	}
	for(int i=0;i<N+1;i++){		
		W[i][i] = 2;			//diagonal
		
	}
	for(int i=0;i<N;i++){
		W[i][i+1] = LAMBDA[i];	//above diagonal
	}
	for(int i=1;i<N+1;i++){
		W[i][i-1] = MU[i];	//below diagonal
		
	}
	
	L = cholesky(W);                        //in-place Cholesky Decomposition
    M = choSolve(L,d);                      //solve the linear system by forward and backward substitution
}
double spline(double x,int N,VEC &X,VEC &Y,VEC &M){
	double S;
	VEC h(N);
	h[0] = 0;
	
	for(int i=1;i<N;i++){				// h is the length of the subinterval 
		h[i] = X[i] - X[i-1];
	}
	for(int i=1;i<N;i++){
		if(x<=X[i]&&x>=X[i-1])				
			S = Y[i-1]+((Y[i]-Y[i-1])/h[i]-(h[i]/6)*(M[i]+2*M[i-1]))*(x-X[i-1])+(M[i-1]/2)*pow(x-X[i-1],2)+(M[i]-M[i-1])/(6*h[i])*pow(x-X[i-1],3);
	}
	return S;
}

VEC polyRoots(double x,VEC &A,int maxiter, double eps){
	int n = A.leng()-1;
	double err,f,df,B_i,C_i;
	VEC Z(n);
	int k;
	VEC B(n+1);
	VEC C(n+1);
	while(n>=1){
		err = 1+eps;
		k = 0;
		while((err>=eps)&&(k<maxiter)){
			B[n-1] = A[n];
			C[n-1] = B[n-1];
			for(int j=n-2;j>=0;j--)
				B[j] = A[j+1]+x*B[j+1];
			for(int j=n-3;j>=0;j--)
				C[j] = B[j+1]+x*C[j+1];
			B_i = A[0]+x*B[0];
			C_i = B[0]+x*C[0];
			f = B_i;
			df = C_i;
			x = x - f/df;
			err = fabs(f);
			k++;
		}
		Z[n-1] = x;
		for(int j=0;j<n;j++)
			A[j] = B[j];
		x = Z[n-1];
		n--;
	}
	return Z;
}
