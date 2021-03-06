//hw04 Resistor Networks
//100060007 Ming-Chang Chiu
//Date: 2015.3.25
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MAT.h"
#include <cmath>
using namespace std;
int main(int argc,char * argv[])
{
    //node: node per side; resisPerSide: resistance per side
    int resisPerSide, node;
    resisPerSide=atoi(argv[1]);
    node = resisPerSide+1;

    //r:resistance; g:conductance; V:voltage of vlt source
    double r,g,V=1.0;

    //a: stamping matrix; LU: matrix containing L and U after matrix "a" perform LU decomposition
    MAT a(node*node), LU(node*node);

    //b:RHS; ans: answer to the voltage vector; temp: temporary answer after fwdSubs
    VEC b(node*node),ans(node*node),temp(node*node);

    r = 2000.0 / ((double) resisPerSide);
    g = 1.0/r;
    
    /**Ordering is like
       0 1 2 3 4  5
       6 7 8 9 10 11
       ...
     **/
    
    //Construct matrix using Stamping Method

    //indexing: row nultiply "node per side" plus column
    for (int i = 0; i<node; i++){
	for(int j = 0; j<node; j++){
	    if( (i==(node-1)) && (j==(resisPerSide/2))){//ground node
		a[i*node+j][0] = a[0][i*node+j] = 1;
	    }
	    else if(j==0){
		if(i==0){
		    a[i*node+j][i*node+j] += 2*g;
		    a[i*node+j][i*node+j+1] -= g;//right
		    a[i*node+j][(i+1)*node+j] -= g;//down
		    //printf("%lf\n",a[i][i]);
		}
		else if(i==node-1){
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
	    else if(j==node-1){
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
    //Make sure the column of Iv is correct
    //In this case the Vlt Source is connected to node 1 so we start from i=1; otherwise
    //we shall modify the i accordingly
    for(int i = 1;i<node*node;i++)
	a[i][node*(node-1)+(resisPerSide/2)]=0;
    //Matrix formation END


    //set the Voltage of vlt source
    b[node*(node-1)+(resisPerSide/2)]=V;
    /*
    for(int i = 0;i<node*node;i++){
		for(int j=0;j<node*node;j++)
			printf(" %4lf",a[i][j]);
		printf("\n");
    }
    */
    //Execute LU Dcomposition
    LU = luFact(a);       //LU in-place decomposition
    temp=fwdSubs(LU,b);   //forward substitution
    ans= bckSubs(LU,temp);//backward substitution
    //for(int i=0;i<node*node;i++) printf(" %lf\n",ans[i]);
    printf("Vne=%lf; Vsw=%lf; Vse=%lf; Req=%lf\n", ans[node-1], ans[node*(node-1)], ans[node*node-1], abs(V/ans[node*(node-1)+(resisPerSide/2)]));
}
