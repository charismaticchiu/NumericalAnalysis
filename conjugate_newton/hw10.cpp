#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "VEC.h"
#include "POLY.h"
#include "Complex.h"
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )		//compare a and b, if a>b return a, else return b
#endif
int main()
{
	int n = 3;			//order+1
	double h = 0;
	int maxiter = 100000;
	double eps = pow(10,-10);
	VEC A(n);
	VEC R1(n-1);
/*VEC R2(n-1);
	VEC R3(n-1);
	VEC R4(n-1);*/
	/*
	Complex a0(4,0);
	Complex a1(1,0);
	Complex a2(5,0);
	Complex a3(7,0);
	Complex a4(3,0);
	Complex a5(6,0);
	Complex a6(0,0);
	Complex a7(1,0);
	A[0] = a0;		//x0
	A[1] = a1;		//x1
	A[2] = a2;		//x2
	A[3] = a3;
	A[4] = a4;
	A[5] = a5;
	A[6] = a6;
	A[7] = a7;
	*/
	Complex a0(-1,0);
	Complex a1(0,0);
        Complex a2(1,0);
//Complex a3(1,0);
	A[0] = a0;		//x0
	A[1] = a1;		//x1
	A[2] = a2;		//x2
// A[3] = a3;
	for(int i=0;i<n;i++){
		h = max(h,fabs(A[i]/A[n-1]));
	}
	Complex x1(1+h,0);
	
	R1 = polyRoots(x1,A,maxiter,eps);
	
	for(int i=0;i<n-1;i++){
		printf("%lf+j%lf\n",R1[i].r(),R1[i].i());
		
		
	}
	return 0;
}
