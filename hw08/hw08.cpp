//hw08 Matrix Eigenvalues
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
    int dim, iter1, iter2;
    scanf("%d",&dim);
    MAT a(dim);
    for(int i=0;i<dim;i++)
	for(int j=0;j<dim;j++)
	    scanf("%lf",&a[i][j]);

    iter1 = EVqr(a,0.000000001, 50000);    
    //Uncomment when need to use Shifted form
    //iter1 = EVqrShifted(a,0.5,0.000000001, 4000);
    
    printf("\nNumber of Iteration is: %d\n",iter1);
    
    return 0;
}
