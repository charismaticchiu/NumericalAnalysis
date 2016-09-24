//hw09 Polynomial Interpolations
//100060007 Ming-Chang Chiu
//Date: 2015.4.12
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
    double xdata[21]={0.0},xinput[301]={0.0},diff[301] = {0.0};
    double ydata[21]={0.0},yinput[301]={0.0},youtput[301] = {0.0};
    char a,b;
    double x = 475.0, y ;
    FILE *fp, *fp2;
    fp = fopen("f301.dat","r");
    fp2= fopen("data.txt","w");
    
    scanf("%c %c",&a,&b);  
    
    //feed f*.dat to arrays
    while(scanf("%lf %lf",&xdata[i],&ydata[i])!=EOF){
	i++;
    }// i will than be the dimension of VECTOR
    //printf("Dimension = %d\n",i);
    

    //read the first line
    fscanf(fp,"%c %c",&a,&b);
    
    while(fscanf(fp,"%lf %lf",&xinput[k],&yinput[k])!=EOF){	
	k++;
    }

    //convert eaten given points array to VEC type
    VEC XDATA(i),YDATA(i);
    for (int j = 0 ; j < i; j++){ 
	XDATA[j] = xdata[j]; 
	YDATA[j] = ydata[j];
    }

    // Do Lagrange interpolation
    for (int j = 0; j < 301; j++){
	youtput[j] = Lagrange(xinput[j], XDATA, YDATA);
	diff[j] = fabs(yinput[j]-youtput[j] );
	fprintf(fp2,"%lf\t%lf\t%lf\t%lf\n",xinput[j],youtput[j],diff[j],yinput[j]);
    }
    /*
    printf("Difference:\n");
    for (int j = 0; j < 301; j++){
	
	fprintf("%lf\n",);
    }
    */
    
    return 0;
}
