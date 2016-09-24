//hw12 Nonlinear Resistor Networks
//100060007 Ming-Chang Chiu
//Date: 2015.5.29
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MAT.h"
#include <cmath>
using namespace std;
int main(int argc,char * argv[])
{
    int ch = atoi(argv[1]);
    //h: step size
    //R: resistance
    //L: inductance
    //t: time, start from 0.0
    //i: current, initial value = 1.0
    double h,R = 1, L = pow(10,-9),t = 0.;
    double i = 1., V,analytical,error;
    FILE *fp;
    FILE *fp2;
    FILE *fp3;
           
    if (ch == 1){//forward
	fp  = fopen("data1.csv","w");
	h = 2 * pow(10, -11);//step size
	
	while(t <= 5.01 * pow(10, -9)){
	   
	    if(t == 0){
		V = 1;
		i = i + (h/L) * (V - R * i);		
		analytical = exp(-1*R*t/L);
		error = fabs(analytical - i);
		fprintf(fp,"%g,%lf,%lf,%lf\n",t,i,analytical, error);		
	    }
	    else{
		V = 0;
		i = i + (h/L) * (V - R * i);		
		analytical = exp(-1*R*t/L);
		error = fabs(analytical - i);
		fprintf(fp,"%g,%lf,%lf,%lf\n",t,i,analytical, error);		
	    }	    
	    t += h;
	}
	printf("t = %g\n",t);
	fclose(fp);
    }
    else if (ch == 2){//backward
	fp2 = fopen("data2.csv","w");
	h = 2* pow(10, -11);//step size
	while(t <= 5.01 * pow(10, -9)){
	    double y = (1 + (h*R/L));
	    if(t == 0){
		V = 1;
		i = (i + (h * V / L) ) / y;		
		analytical = exp(-1*R*t/L);
		error = fabs(analytical - i);
		fprintf(fp2,"%g,%lf,%lf,%lf\n",t,i,analytical, error);	
	    }
	    else{
		V = 0;
		i = (i + (h * V / L) ) / y;		
		analytical = exp(-1*R*t/L);
		error = fabs(analytical - i);
		fprintf(fp2,"%g,%lf,%lf,%lf\n",t,i,analytical, error);
	    }
	    
	    t += h;
	}
	fclose(fp2);
    }
    else if (ch == 3){//trapezoidal
	fp3 = fopen("data3.csv","w");
	h = 2 * pow(10, -10);//step size
	while(t <= 5.01 * pow(10, -9)){
	    double y = (1 + (h * R /(L* 2)));
	    if(t == 0){
		V = 1;
		i = (i + (h * V/L) - (h * R * i / (2.*L)) )/y;		
		analytical = exp(-1*R*t/L);
		error = fabs(analytical - i);
		fprintf(fp3,"%g,%lf,%lf,%lf\n",t,i,analytical, error);		
	    }
	    else{
		V = 0;
		i = (i + (h * V/L) - (h * R * i / (2.*L)) )/y;		
		analytical = exp(-1*R*t/L);
		error = fabs(analytical - i);
		fprintf(fp3,"%g,%lf,%lf,%lf\n",t,i,analytical, error);		
	    }
	    
	    t += h;
	}
	fclose(fp3);
    }
    else{
	printf("Wrong input!!! Please start over!!");
    }
    
    
    
    return 0;
}
