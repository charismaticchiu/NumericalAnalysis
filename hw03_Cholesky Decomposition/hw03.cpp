//hw03 Cholesky Decomposition
//100060007 Ming-Chang Chiu
//Date: 2015.3.18
#include <stdio.h>
#include <stdlib.h>
#include "MAT.h"
#include <math.h>
using namespace std;
int main()
{
    int dim;
    scanf("%d",&dim);
    MAT a(dim),L(dim),U(dim);
    VEC b(dim),ans(dim),temp(dim);
    for(int i=0;i<dim;i++)
	for(int j=0;j<dim;j++)
	    scanf("%lf",&a[i][j]);
    for(int i=0;i<dim;i++) scanf("%lf",&b[i]);
    a = cholesky(a);       //Cholesky in-place decomposition
    ans = choSolve(a,b);   //forward and backward substitution
    /*for(int i=0;i<dim;i++){
	for(int j=0;j<dim;j++){
	    printf(" %lf",a[i][j]);
	}
	printf("\n");
	}*/
    for(int i=0;i<dim;i++) printf(" %lf\n",ans[i]);
    //printf("%d",dim^2);
}
