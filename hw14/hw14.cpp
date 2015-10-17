//hw14 Nonlinear Resistor Networks
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
    double h,R1 = 1, R2 = 1, L2 = pow(10,-9), L1 = pow(10,-12),t = 0.;
    double i1 = 0., i2 = 0., V = 1., V0 = 0., analytical,error;
    FILE *fp1;
    FILE *fp2;
    FILE *fp3;
    FILE *fp4;
    FILE *fp5;
    FILE *fp6;
    FILE *fp7;
    VEC temp(2);
    VEC x(2);
    
    if (ch == 1){//trapezoidal, h = 10^-13
	MAT A(2);
	h = pow(10, -13);//step size
	A[0][0] = 2*L1+h*R1;
	A[1][1] = h*R2+2*L2+h*R1;
	A[0][1] = -1*h*R1;
	A[1][0] = -1*h*R1;
	MAT LU = luFact(A);
	VEC b(2);
	fp1 = fopen("data1.csv","w");
	
	while(t <= 5.01 * pow(10, -8)){
	    
	    if(t == 0){		
		b[0] = 2*L1*i1+h*(V+V0-R1*(i1-i2));
		b[1] = 2*L2*i2+h*(R1*(i1-i2)-i2*R2);
		temp = fwdSubs(LU, b);
		x = bckSubs(LU,temp);		
		fprintf(fp1,"%g,%lf,%lf\n",t,i1,i2);		
	    }
	    else{
		
		b[0] = 2*L1*i1+h*(V+V-R1*(i1-i2));
		b[1] = 2*L2*i2+h*(R1*(i1-i2)-i2*R2); 
		temp = fwdSubs(LU, b);
		x = bckSubs(LU,temp);		
		fprintf(fp1,"%g,%lf,%lf\n",t,i1,i2);		
	    }
	    i1 = x[0];
	    i2 = x[1];
	    t += h;
	}
	fclose(fp1);
    }
    else if (ch == 2){//trapezoidal, h = 10^-11
	MAT A(2);
	h = pow(10, -11);//step size
	A[0][0] = 2*L1+h*R1;
	A[1][1] = h*R2+2*L2+h*R1;
	A[0][1] = -1*h*R1;
	A[1][0] = -1*h*R1;
	MAT LU = luFact(A);
	VEC b(2);
	fp2 = fopen("data2.csv","w");
	
	while(t <= 5.01 * pow(10, -8)){
	    
	    if(t == 0){		
		b[0] = 2*L1*i1+h*(V+V0-R1*(i1-i2));
		b[1] = 2*L2*i2+h*(R1*(i1-i2)-i2*R2);
		temp = fwdSubs(LU, b);
		x = bckSubs(LU,temp);		
		fprintf(fp2,"%g,%lf,%lf\n",t,i1,i2);		
	    }
	    else{
		
		b[0] = 2*L1*i1+h*(V+V-R1*(i1-i2));
		b[1] = 2*L2*i2+h*(R1*(i1-i2)-i2*R2);
		temp = fwdSubs(LU, b);
		x = bckSubs(LU,temp);		
		fprintf(fp2,"%g,%lf,%lf\n",t,i1,i2);		
	    }
	    i1 = x[0];
	    i2 = x[1];
	    t += h;
	}
	fclose(fp2);
    }
    else if (ch == 3){//trapezoidal, h = 10^-9
	MAT A(2);
	h = pow(10, -9);//step size
	A[0][0] = 2*L1+h*R1;
	A[1][1] = h*R2+2*L2+h*R1;
	A[0][1] = -1*h*R1;
	A[1][0] = -1*h*R1;
	MAT LU = luFact(A);
	VEC b(2);
	fp3 = fopen("data3.csv","w");
	
	while(t <= 5.01 * pow(10, -8)){
	    
	    if(t == 0){		
		b[0] = 2*L1*i1+h*(V+V0-R1*(i1-i2));
		b[1] = 2*L2*i2+h*(R1*(i1-i2)-i2*R2);
		temp = fwdSubs(LU, b);
		x = bckSubs(LU,temp);		
		fprintf(fp3,"%g,%lf,%lf\n",t,i1,i2);		
	    }
	    else{
		
		b[0] = 2*L1*i1+h*(V+V-R1*(i1-i2));
		b[1] = 2*L2*i2+h*(R1*(i1-i2)-i2*R2);
		temp = fwdSubs(LU, b);
		x = bckSubs(LU,temp);		
		fprintf(fp3,"%g,%lf,%lf\n",t,i1,i2);	
	    }
	    i1 = x[0];
	    i2 = x[1];
	    t += h;
	}
	fclose(fp3);
    }
    else if(ch == 4){//forward,h = 10^-11
	h = pow(10, -11);//step size
	fp4 = fopen("data4.csv","w");
	int count = 0;
	//i1n_1: i(n-1)
	//newi1: i(n+1)
	double i1n_1, i2n_1, newi1,newi2;
	while(t <= 5.01 * pow(10, -8)){
	    if(t == 0){//use first order fwd to generate x(h)
		
		newi1 = i1 + h * (V0 - (i1-i2) * R1)/L1;
		newi2 = i2 + h * ((i1-i2)*R1 - i2*R2)/L2;
		fprintf(fp4,"%g,%lf,%lf\n",t,i1,i2);
		count += 1;
	    }
	    else if(count == 1){
		newi1 = i1 + (h/(2*L1))*(3*(V - (i1-i2)*R1) - (V0 - (i1n_1-i2n_1)*R1));
		newi2 = i2 + (h/(2*L2))*((3*((i1-i2)*R1-i2*R2)) - ( (i1n_1-i2n_1)*R1-i2n_1*R2));
		fprintf(fp4,"%g,%lf,%lf\n",t,i1,i2);
	    }
	    else {
		newi1 = i1 + (h/(2*L1))*(3*(V - (i1-i2)*R1) - (V  - (i1n_1-i2n_1)*R1));
		newi2 = i2 + (h/(2*L2))*((3*((i1-i2)*R1-i2*R2)) - ( (i1n_1-i2n_1)*R1-i2n_1*R2));
		fprintf(fp4,"%g,%lf,%lf\n",t,i1,i2);
	    }
	    i1n_1 = i1;
	    i2n_1 = i2;
	    i1 = newi1;
	    i2 = newi2;
	    t+=h;
	}        
	fclose(fp4);
    }
    else if(ch == 5){//forward,h = 10^-9
	h = pow(10, -9);//step size
	fp5 = fopen("data5.csv","w");
	int count = 0;
	//i1n_1: i(n-1)
	//newi1: i(n+1)
	double i1n_1, i2n_1, newi1,newi2;
	i1n_1= i2n_1= newi1=newi2=0.;
	while(t <= 5.01 * pow(10, -8)){
	    
	    if(t == 0){//use first order fwd to generate x(h)
		
		newi1 = i1 + h * (V0 - (i1-i2) * R1);
		newi2 = i2 + h * ((i1-i2)*R1 - i2*R2);
		fprintf(fp5,"%g,%lf,%lf\n",t,i1,i2);
		count += 1;
	    }
	    else if(count == 1){//the second step, x(2h)
		
		newi1 = i1+ h/(2*L1)*(   3.*(V-(i1-i2)*R1)  -  (V0-(i1n_1-i2n_1)*R1)  );	
		newi2 = i2 + h/(2*L2)*( (3.*( (i1-i2)*R1-i2*R2) ) - ( (i1n_1-i2n_1)*R1-i2n_1*R2) );
		
		fprintf(fp5,"%g,%lf,%lf\n",t,i1,i2);
		count += 1;
	    }
	    else {
		newi1 =i1+ h/(2*L1)*(3.*(V-(i1-i2)*R1)-(V-(i1n_1-i2n_1)*R1));
		newi2 = i2 + h/(2*L2)*((3.*((i1-i2)*R1-i2*R2)) - ( (i1n_1-i2n_1)*R1-i2n_1*R2));
		printf("%g,%g,%g\n",t,i1*L1,i2*L2);
		fprintf(fp5,"%g,%lf,%lf\n",t,i1,i2);
	    }
	    i1n_1 = i1;
	    i2n_1 = i2;
	    i1 = newi1;
	    i2 = newi2;
	    t += h;
	}        
	fclose(fp5);
    }
    else if(ch == 6){//Gear 2, h = 10^-11
	h = pow(10, -11);//step size
	fp6 = fopen("data6.csv","w");
	MAT A(2);
	A[0][0] = (2*h*R1)/(3*L1) + 1;
	A[1][1] = 2*h*(R2+R1)/(3*L2) + 1;
	A[0][1] = (-2 * h*R1)/(3*L1);
	A[1][0] = (-2 * h*R1)/(3*L2);
	MAT LU = luFact(A);
	VEC b(2);
        double i1n_1,i2n_1;
	while(t <= 5.01 * pow(10, -8)){
	    
	    if(t == 0){//use first order fwd to generate x(h)
		
		x[0] = i1 + h * (V0 - (i1-i2) * R1)/L1;
		x[1] = i2 + h * ((i1-i2)*R1 - i2*R2)/L2;
		fprintf(fp6,"%g,%lf,%lf\n",t,i1,i2);
		
	    }
	    else{
		
		b[0] = 4.*i1/3. - 1.*i1n_1/3. + 2*h*V/(3*L1);
		b[1] = 4.*i2/3. - 1.*i2n_1/3.;
		temp = fwdSubs(LU, b);
		x = bckSubs(LU,temp);		
		fprintf(fp6,"%g,%lf,%lf\n",t,i1,i2);	
	    }
	    i1n_1 = i1;
	    i2n_1 = i2;
	    i1 = x[0];
	    i2 = x[1];
	    t += h;
	}
	fclose(fp6);
    }
    else if(ch == 7){//Gear 2, h = 10^-11
	h = pow(10, -9);//step size
	fp7 = fopen("data7.csv","w");
	MAT A(2);
	A[0][0] = (2*h*R1)/(3*L1) + 1;
	A[1][1] = 2*h*(R2+R1)/(3*L2) + 1;
	A[0][1] = (-2 * h*R1)/(3*L1);
	A[1][0] = (-2 * h*R1)/(3*L2);
	MAT LU = luFact(A);
	VEC b(2);
        double i1n_1,i2n_1;
	while(t <= 5.01 * pow(10, -8)){
	    
	    if(t == 0){//use first order fwd to generate x(h)
		
		x[0] = i1 + h * (V0 - (i1-i2) * R1)/L1;
		x[1] = i2 + h * ((i1-i2)*R1 - i2*R2)/L2;
		fprintf(fp7,"%g,%lf,%lf\n",t,i1,i2);
		
	    }
	    else{
		
		b[0] = 4.*i1/3. - 1.*i1n_1/3. + 2*h*V/(3*L1);
		b[1] = 4.*i2/3. - 1.*i2n_1/3.;
		temp = fwdSubs(LU, b);
		x = bckSubs(LU,temp);		
		fprintf(fp7,"%g,%lf,%lf\n",t,i1,i2);	
	    }
	    i1n_1 = i1;
	    i2n_1 = i2;
	    i1 = x[0];
	    i2 = x[1];
	    t += h;
	    //printf("%lf\n",t);
	}
	fclose(fp7);
    }
    else{
	printf("Wrong input!!! Please start over!!");
    }
    
    
    
    return 0;
}
