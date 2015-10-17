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
    double limit = 0.5;
	double eps = pow(10,-8);
    double ans;
    
    ans = newton(2.0,limit,eps);
    printf("Ans: %lf\n", ans);  
    return 0;
}
