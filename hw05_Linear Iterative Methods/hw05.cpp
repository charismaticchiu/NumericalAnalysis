//hw05 Linear Iterative Methods
//100060007 Ming-Chang Chiu
//Date: 2015.4.5
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MAT.h"
#include <cmath>
#include <ctime>
#include <climits>
using namespace std;
int main()
{
    //node: node per side; resisPerSide: resistance per side; iter: number of iteration
    int resisPerSide, node, iter;
    resisPerSide = 20;
    node = resisPerSide+1;

    //r:resistance; g:conductance; V:voltage of vlt source
    double r, g, V=1.0;

    //a: stamping matrix; LU: matrix containing L and U after matrix "a" perform LU decomposition
    MAT a(node*node), LU(node*node);

    //b:RHS; ans: answer to the voltage vector; temp: temporary answer after fwdSubs
    VEC b(node*node),ans(node*node),temp(node*node), x1(node*node),x2(node*node),x3(node*node);
    //1-norm, 2-norm, or inf-norm
    char norm = 'i';

    switch (resisPerSide){
	case  2:
	    r= 1000.0;
            break;
	case  4:
	    r = 500.0;
            break;
	case 10:
	    r = 200.0;
            break;
	case 20:
	    r = 100.0;
            break;
        default:
	    printf("Wrong input! Please execute the program again!\n");
	    exit(0);
    }
    g = 1.0/r;

    /**Ordering is like
       0 1 2 3 4  5
       6 7 8 9 10 11
       ...
     **/
     
    //Construct matrix using Kichhoff Method

    //indexing: row nultiply "node per side" plus column
    for (int i = 0; i<node; i++){
	for(int j = 0; j<node; j++){
	    if( (i==(node-1)) && (j==(resisPerSide/2))){//ground node
		a[i*node+j][i*node+j] = 1;
	    }
	    else if(j==0){// first column
		if(i==0){//first row
		    /*
		    a[i*node+j][i*node+j] += 2*g;
		    a[i*node+j][i*node+j+1] -= g;//right
		    a[i*node+j][(i+1)*node+j] -= g;//down
		    */
		    a[i*node+j][i*node+j] = 1;
		    //printf("%lf\n",a[i][i]);
		}
		else if(i==node-1){//last row
		    a[i*node+j][i*node+j] += 2*g;
		    a[i*node+j][(i-1)*node+j] -= g;//up
		    a[i*node+j][i*node+j+1] -= g;//right
		    
		}
		else{
		    a[i*node+j][i*node+j] += 3*g;
		    a[i*node+j][(i-1)*node+j] -= g;//up
		    a[i*node+j][i*node+j+1] -= g;//right
		    a[i*node+j][(i+1)*node+j] -= g;//down
		    //printf("a(%d,%d)= %lf",i*node+j,i*node+j, a[i*node+j][i*node+j]);
		}
	    }
	    else if(j==node-1){//last column
		if(i==0){
		    a[i*node+j][i*node+j] += 2*g;
		    a[i*node+j][(i+1)*node+j] -= g;//down
		    a[i*node+j][i*node+j-1] -= g;//left
		}
		else if(i==node-1){
		    a[i*node+j][i*node+j] += 2*g;
		    a[i*node+j][(i-1)*node+j] -= g;//up
		    a[i*node+j][(i)*node+j-1] -= g;//left
		}
		else{
		    a[i*node+j][i*node+j] += 3*g;
		    a[i*node+j][(i-1)*node+j] -= g;//up
		    a[i*node+j][(i+1)*node+j] -= g;//down
		    a[i*node+j][(i)*node+j-1] -= g;//left
		}
	    }
	    else{
		if(i==0){
		    a[i*node+j][i*node+j] += 3*g;
		    a[i*node+j][i*node+j+1] -= g;//right
		    a[i*node+j][(i+1)*node+j] -= g;//down
		    a[i*node+j][(i)*node+j-1] -= g;//left

		}
		else if(i==node-1){
		    a[i*node+j][i*node+j] += 3*g;
		    a[i*node+j][(i-1)*node+j] -= g;//up
		    a[i*node+j][i*node+j+1] -= g;//right
		    a[i*node+j][(i)*node+j-1] -= g;//left
		    //printf("a(%d,%d)= %lf",i*node+j,i*node+j, a[i*node+j][i*node+j]);
		}
		else{
		    a[i*node+j][i*node+j] += 4*g;
		    a[i*node+j][(i-1)*node+j] -= g;//up
		    a[i*node+j][i*node+j+1] -= g;//right
		    a[i*node+j][(i+1)*node+j] -= g;//down
		    a[i*node+j][(i)*node+j-1] -= g;//left
		    //printf("a(%d,%d)= %lf",i*node+j,i*node+j, a[i*node+j][i*node+j]);
		}
	    }
	}
    }   
     
    //Matrix formation END


    //set the Voltage of vlt source
    b[0]=V;
   
    clock_t t1, t2, t3, t4, t5, t6;
    
    t1 = clock();
    iter = jacobi(a,b,x1,50000,0.0000001);
    t2 = clock();
    printf("Time for jacobi: %f\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
    printf("Number of iter for jacobi: %d\n",iter);
    
    
    t3 = clock();
    iter = gaussSeidel(a,b,x2,50000,0.0000001);
    t4 = clock();
    printf("Time for Gauss-Seidel: %f\n", (double)(t4 - t3) / CLOCKS_PER_SEC);
    printf("Number of iter for Gauss-Seidel: %d\n",iter);
    
    
    t5 = clock();
    iter = sgs(a,b,x3,50000,0.0000001);
    t6 = clock();
    printf("Time for sgs: %f\n", (double)(t6 - t5) / CLOCKS_PER_SEC);
    printf("Number of iter for symmetric Gauss-Seidel: %d\n",iter);
    

    

}
