#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MAT.h"
#include "Complex.h"
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )		//compare a and b, if a>b return a, else return b
#endif

VEC polyRoots(Complex x,VEC &A,int maxiter, double eps){
	int n = A.len()-1;
	double err;
	Complex f(x);
	Complex df(x);
	Complex B_i(x);
	Complex C_i(x);
	VEC Z(n);
	int k;
	VEC B(n+1);
	VEC C(n+1);
	while(n>=1){
		err = 1+eps;
		k = 0;
		while((err>=eps)&&(k<maxiter)){
			B[n-1] = A[n];
			C[n-2] = B[n-1];
			for(int j=n-2;j>=0;j--)		//n-2
				B[j] = A[j+1]+x*B[j+1];
			for(int j=n-3;j>=0;j--)		//n-3
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
		for(int j=0;j<n+1;j++)
			A[j] = B[j];
		x = Z[n-1];
		n--;
	}
	return Z;
}
