//hw10 Spline Interpolations
//100060007 Ming-Chang Chiu
//Date: 2015.5.8
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
    int i = 0, k = 0;
    //xdata, yadata will eat f*.dat, which is the given points
    //xinput,yinput will eat f301.dat, which is the original waveform
    //diff is the error for each X coordinate
    //youtput will store the interpolated value for each X
    i = atoi(argv[1]);
    VEC XDATA(i),YDATA(i),M(i);// i == n+1 in the spline()
    double xdata[21]={0.0},xinput[301]={0.0},diff[301] = {0.0}, maxerr = 0.0;
    double ydata[21]={0.0},yinput[301]={0.0},youtput[301] = {0.0},youtput2[301] = {0.0},ans;
    char a,b;
    double x = 475.0, y ;
    FILE *fp, *fp2;
    fp = fopen("f301.dat","r");
    fp2= fopen("data.txt","w");
    
    scanf("%c %c",&a,&b);  
    
    //feed f*.dat to arrays
    for(int j = 0;j < i;j++) scanf("%lf %lf",&XDATA[j],&YDATA[j]);
    //printf("Dimension = %d\n",i);
    

    //read the first line of f301
    fscanf(fp,"%c %c",&a,&b);
    //read the rest of f301
    while(fscanf(fp,"%lf %lf",&xinput[k],&yinput[k])!=EOF){	
	k++;
    }

    //convert eaten given points array to VEC type
    /* VEC XDATA(i),YDATA(i),M(i);// i == n+1 in the spline()
    for (int j = 0 ; j < i; j++){ 
	XDATA[j] = xdata[j]; 
	YDATA[j] = ydata[j];
	}*/

    // Do Spline interpolation    
    for (int j = 0; j < 301; j++){
	//N = i-1 for we have n+1 support points(i.e. i support points)
	youtput[j] = spline(xinput[j],i-1 ,XDATA, YDATA,M);
	diff[j] = fabs(yinput[j]-youtput[j] );
	fprintf(fp2,"%lf\t%lf\t%lf\t%lf\n",xinput[j],youtput[j],diff[j],yinput[j]);
	if (maxerr < diff[j])
	    maxerr = diff[j];
    }
    printf("maximum absolute error: %lf\n", maxerr);  
    return 0;
}
