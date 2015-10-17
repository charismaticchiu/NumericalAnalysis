//hw06 Conjugate Gradient Methods
//100060007 Ming-Chang Chiu
//Date: 2015.4.12
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MAT.h"
#include <cmath>
#include <ctime>
#include <climits>
using namespace std;
int main(int argc, char * argv[])
{
    //node: node per side; resisPerSide: resistance per side; iter: number of iteration
    int resisPerSide, node, iter;
    resisPerSide = atoi(argv[1]);
    node = resisPerSide+1;

    //r:resistance; g:conductance; V:voltage of vlt source
    double r, g, V=1.0;

    //a: stamping matrix; LU: matrix containing L and U after matrix "a" perform LU decomposition
    MAT a(node*node), LU(node*node);

    //b:RHS; ans: answer to the voltage vector; temp: temporary answer after fwdSubs
    VEC b(node*node),ans(node*node),temp(node*node), x1(node*node),x2(node*node),x3(node*node);

    r = 2000.0 / ((double) resisPerSide);
    //printf("r:%g",r);
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
    b[node] = b[1] = g * V;
    a[1][0] = a[node][0] = 0;
    //In this case the Vlt Source is connected to node 1 so we start from i=1; otherwise
    //we shall modify the i accordingly
    for(int i = 1; (i<node*node);i++)
    {
	if ( i != node*(node-1)+(resisPerSide/2))
	    a[i][node*(node-1)+(resisPerSide/2)]=0;
	//printf("%d; ",i);
    }
	//Matrix formation END
    /*
    for(int i = 0 ; i< node*node; i++){
      for (int j = 0; j <node*node; j++)
	printf("%lf ",a[i][j]);
      printf("\n");
    }
    */
    iter = cg(a,b,x1,50000,0.0000001);
    //printf("Time for cg: %f\n", (double)(t6 - t5) / CLOCKS_PER_SEC);
    printf("Number of iter for Congugate Gradient: %d\n",iter);
    

}
