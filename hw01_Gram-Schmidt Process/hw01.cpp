#include <stdio.h>
#include <stdlib.h>
#include "MAT.h"
#include <math.h>
using namespace std;

void gramSchmidtI(MAT &gT, MAT &aT,int dim);
void gramSchmidtII(MAT &gT, MAT &aT,int dim);
void gramSchmidtIII(MAT &gT, MAT &aT,int dim);

int main(int argc, char *argv[]){
  int dim = 0;
  double sigma = 0.0;
  scanf("%d",&dim);
  MAT a(dim);
  MAT aT(dim);
  MAT g(dim);
  MAT gT(dim);
  MAT d(dim);
  for (int i=0;i<dim;i++)
    for (int j=0;j<dim;j++)
      scanf("%lf",&a[i][j]);
  aT=a.tpose();
  gT[0] = aT[0];
  //*******Gram-Schmidt
  gramSchmidtI(gT, aT, dim);
  //*******END
  g = gT.tpose();
  d = gT*g;
  /*
  printf("D=\n");
  for (int i=0;i<dim;i++){
    for (int j=0;j<dim;j++)
      printf(" %lf",d[i][j]);
    printf("\n");
  }
  */
  for (int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      if(i==j)
  	j++;
      if (j==dim) break;
      sigma += (d[i][j]*d[i][j]);
    }
  }
  sigma=sqrt(sigma);
  printf("sigma = %lf\n",sigma);
  return 0;
}

void gramSchmidtI(MAT &gT, MAT &aT,int dim){
  VEC zero(dim);
  for (int k=1;k<dim;k++){
    gT[k]=zero;
    for (int i=0;i<k;i++){
      gT[k] = gT[k] + (((aT[k]*gT[i])*gT[i])/(gT[i]*gT[i]));
    }
    gT[k] = aT[k] - gT[k];
  }
}
void gramSchmidtII(MAT &gT, MAT &aT,int dim){

  for(int k=1;k<dim;k++){
    gT[k]=aT[k];
    for(int i=0;i<k;i++){//Why it is i<k rather than k-1?
      gT[k] = gT[k]- (((gT[k] * gT[i]) * gT[i]) / (gT[i] * gT[i]));
     }
  }
}
void gramSchmidtIII(MAT &gT, MAT &aT,int dim){
  
  for(int k=1;k<dim;k++){
    gT[k]=aT[k];
    for(int i=0;i<k;i++){//Why it is i<k rather than k-1?
      gT[k] = gT[k]- (((gT[k] * gT[i])  / (gT[i] * gT[i]))*gT[i]);
    }
  }
}
