#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MAT.h"
MAT::MAT(int dim)
{
  n=dim;
  va=(VEC **)malloc(n*sizeof(VEC*));
  for (int i=0;i<n;i++)
    va[i]=newVEC(n);
}
MAT::MAT(const MAT &m1)
{
  VEC **vsrc=m1.va; //to get around not indexing const MAT
  n=m1.n;
  va=(VEC **)malloc(n*sizeof(VEC*));
  for(int i=0;i<n;i++){
    va[i]=newVEC(n);
    (*va[i])=(*vsrc[i]);
  }
}
MAT::MAT(int dim,double *v)
{
  n=dim;
  va=(VEC **)malloc(n*sizeof(VEC*));
  for(int i=0;i<n;i++){
    va[i]=newVEC(n);
    for(int j=0;j<n;j++){
      (*va[i])[j]= *(v++); //array indexing + VEC indexing
    }
  }
}
MAT::~MAT()
{
  for (int i=n-1;i>=0;i--) free(va[i]);
  free(va);
}
int MAT::dim()
{
  return n;
}
MAT MAT::tpose()
{
    MAT mnew(n);
    for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
	    mnew[i][j]=(*va[j])[i];
	}
    }
    return mnew;
}
MAT & MAT::operator-() {
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
	    (*va[i])[j]=-(*va[i])[j];
}
MAT &MAT::operator=(MAT m1)
{
    for (int i=0; i<n; i++)
        (*va[i])=m1[i];
    return *this;
}
MAT &MAT::operator+=(MAT &m1)
{
    for (int i=0; i<n; i++)
        (*va[i])+=m1[i];
    return *this;
}
MAT &MAT::operator-=(MAT &m1)
{
    for (int i=0; i<n; i++)
        (*va[i])-=m1[i];
    return *this;
}
MAT &MAT::operator*=(double a)
{
    for (int i=0; i<n; i++)
        (*va[i])*=a;
    return *this;
}
MAT &MAT::operator/=(double a)
{
    for (int i=0; i<n; i++)
        (*va[i])/=a;
    return *this;
}
MAT MAT::operator+(MAT m1)
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])+m1[i];
    return s;
}
MAT MAT::operator-(MAT m1)
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])-m1[i];
    return s;
}
MAT MAT::operator*(MAT m1)
{
    MAT z(n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            z[i][j]=0;
	    for (int k=0; k<n; k++)
		z[i][j]+=((*va[i])[k]*m1[k][j]);
	}
    }
    return z;
}
VEC &MAT::operator[](int m)
{
    return *va[m];
}
VEC MAT::operator*(VEC v1)
{
    VEC s(n);
    for (int i=0; i<n; i++) {
        s[i]=(*va[i])*v1;
    }
    return s;
}
MAT MAT::operator*(double a)
{
    MAT z(n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            z[i][j]=(*va[i])[j] * a;//why use (*va[i])[j] rather than va[i][j];
        }
    }
    return z;    
}
MAT MAT::operator/(double a)
{
    MAT z(n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            z[i][j]=(*va[i])[j] / a;
        }
    }
    return z;
}
MAT operator*(double a, MAT &m1)
{
    MAT v2(m1.n);
    for (int i=0; i<m1.n; i++) {
        for (int j=0; j<m1.n; j++) {
            v2[i][j]= a * m1[i][j];
        }
    }
    return v2;
}
VEC operator*(VEC &v1,MAT &m1)
{
    VEC v2(m1.n);
    for (int i=0; i<m1.n; i++) {
        v2[i]=0;
        for (int j=0; j<m1.n; j++) {
            v2[i] += v1[j]*m1[j][i];
        }
    }
    return v2;
}
