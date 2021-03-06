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
    int Q = atoi(argv[1]);//Q specifies the Question number
    double v ;    
    double precision = 0.000000001;
    double err1 = 1 + precision,err2 = 1 + precision;
    int k1 = 0, k2 = 0;

    
    if(Q == 1){
	FILE *fp;
	fp = fopen("data1.txt","w");
	v = 0.0;
	for(int l = 0; l < 51; l++){
	    VEC x1(9);
	    x1[0] = v;
	    VEC temp1(x1.len());
	    VEC deltaX1(x1.len());
	    printf("%lf\n", v);
	    err1 = 1 + precision;
	    VEC f = genF1(x1,v);
	    while(err1 > precision){
		
		MAT j = genJ1(x1);
		MAT LU = luFact(j);
		temp1 = fwdSubs(LU, -1. * f);
		deltaX1 = bckSubs(LU, temp1);
		x1 = x1 + deltaX1;
		k1++;
		f = genF1(x1,v);
		err1 = linfnorm(f);
	    }
	    // calculate the resistances after Newton's iteration
	    // let alone the ordering
	    VEC r(x1.len()+12);
	    r[9] = 1 + 0.1 * fabs(x1[0] - x1[1]);
	    r[10] = 1 + 0.1 * fabs(x1[0] - x1[3]);
	    r[11] = 1 + 0.1 * fabs(x1[1] - x1[2]);
	    r[12] = 1 + 0.1 * fabs(x1[1] - x1[4]);
	    r[13] = 1 + 0.1 * fabs(x1[2] - x1[5]);
	    r[14] = 1 + 0.1 * fabs(x1[3] - x1[4]);
	    r[15] = 1 + 0.1 * fabs(x1[3] - x1[6]);
	    r[16] = 1 + 0.1 * fabs(x1[4] - x1[5]);
	    r[17] = 1 + 0.1 * fabs(x1[4] - x1[7]);
	    r[18] = 1 + 0.1 * fabs(x1[5] - x1[8]);
	    r[19] = 1 + 0.1 * fabs(x1[6] - x1[7]);
	    r[20] = 1 + 0.1 * fabs(x1[7] - x1[8]);

	    //Print out informations
	    //printf("Iter:%d\n",k1);
	    printf("Ans for Q1:\n");
	    double curr= (x1[0]-x1[1])/r[9] + (x1[0]-x1[3])/r[10];
	    double i2= (x1[3]-x1[6])/r[15];
	    double i7= (x1[4]-x1[7])/r[17],i12 = (x1[5]-x1[8])/r[18];
	    printf("Tot current: %lf\n",curr);
	    printf("R2  current: %lf\n",i2);
	    printf("R7  current: %lf\n",i7);
	    printf("R12 current: %lf\n",i12);
	    fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",v,curr,i2,i7,i12);
	    v += 0.1;
	 }
	 fclose(fp);
    }
    else if (Q == 2){
	FILE *fp;
	fp = fopen("data2.txt","w");
	v = 0;
	for(int l = 0; l < 51; l++){
	    VEC x2(21);
	    x2[0] = v;
	    VEC temp2(x2.len());
	    VEC deltaX2(x2.len());
	    err2 = 1 + precision;
	    printf("%lf\n", v);
	    VEC f = genF2(x2, v);
	    while(err2 > precision){
		
		MAT j = genJ2(x2);
		MAT LU = luFact(j);
		temp2 = fwdSubs(LU, -1. * f);
		deltaX2 = bckSubs(LU, temp2);
		x2 = x2 + deltaX2;
		k2++;
		f = genF2(x2,v);
		err2 = linfnorm(f);
	    }

	    //Calculate the resistances
	    //let alone the ordering
	    VEC r(x2.len());
	    int resisNum = 12;
	    for(int i = 0; i < resisNum; i++){
		r[i+9] = 1 + x2[i+9];
	    }

	    //Print out info
	    //printf("Iter:%d\n",k2);
	    printf("Ans for Q2:\n");	
	    double curr = (x2[0]-x2[1])/r[9] + (x2[0]-x2[3])/r[10];
	    double t2 = (x2[3]-x2[6])*(x2[3]-x2[6])/r[15];
	    double t7 = (x2[4]-x2[7])*(x2[4]-x2[7])/r[17];
	    double t12 = (x2[5]-x2[8])*(x2[5]-x2[8])/r[18];
	    printf("Tot current: %lf\n",curr);
	    printf("T2  : %lf\n",t2);
	    printf("T7  : %lf\n",t7);
	    printf("T12 : %lf\n",t12);
	    fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",v,curr,t2,t7,t12);
	    v += 0.1;
	   
	}
	fclose(fp);
    }
    return 0;
}
