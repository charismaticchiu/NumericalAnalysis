//hw02 LU Decomposition
//100060007 Ming-Chang Chiu
//Date: 2015.3.17
#include <stdio.h>
#include <stdlib.h>
#include "MAT.h"

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
    a = luFact(a);       //LU in-place decomposition
    temp=fwdSubs(a,b);   //forward substitution
    ans= bckSubs(a,temp);//backward substitution
    for(int i=0;i<dim;i++) printf(" %lf\n",ans[i]);

}
