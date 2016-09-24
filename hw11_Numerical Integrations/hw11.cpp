//hw11 Numerical Integrations
//100060007 Ming-Chang Chiu
//Date: 2015.5.15
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MAT.h"
#include <cmath>
#include <ctime>
#include <climits>
#include <iostream>
using namespace std;
int main(int argc, char * argv[])
{
    int i = 0, k = 0, n = 0;
    i = atoi(argv[1]);
    n = atoi(argv[2]);
    VEC XDATA(i),YDATA(i);// i == n+1 in the spline()    
    double ans;
    char a,b;
    
    scanf("%c %c",&a,&b);    
	//XDATA is the xvalues; YDATA is the y values
    for(int j = 0;j < i;j++) scanf("%lf %lf",&XDATA[j],&YDATA[j]);
        
    // Newton-Cotes intergration    
    ans = integ(XDATA,YDATA,n);
    printf("Ans: %lf\n", ans);  
    return 0;
}
