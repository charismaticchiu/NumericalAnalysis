//hw05 Linear Iterative Methods
//100060007 Ming-Chang Chiu
//Date: 2015.4.5
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "VEC.h"
//using namespace std;
VEC::VEC(int n)
{
  dim=n;
  val=(double *)calloc(n, sizeof(double));
}
VEC::VEC(const VEC &v1)
{
  dim=v1.dim;
  val=(double *)calloc(dim, sizeof(double));
  for (int i=0; i<dim; i++){
    val[i]=v1.val[i];
  }
}
VEC::VEC(int n, double *v)
{
  dim=n;
  val=(double *)calloc(n,sizeof(double));
  for (int i=0; i<n; i++) val[i]=v[i];
}
VEC::~VEC()
{
  free(val);
}
int VEC::len()
{
  return dim;
}
VEC VEC::operator-()
{
  for (int i=0; i<dim; i++) val[i]=-val[i];
  return *this;
}
VEC &VEC::operator=(const VEC v1)
{
  dim=v1.dim;
  for(int i=0; i<dim; i++){
    val[i]=v1.val[i];
  }
  return *this;
}
VEC &VEC::operator+=(const VEC v1)
{
  for (int i=0;i<dim;i++)
    val[i]+=v1.val[i];
  return *this;
}
VEC &VEC::operator-=(const VEC v1)
{
  for (int i=0;i<dim;i++)
    val[i]-=v1.val[i];
  return *this;
}
VEC &VEC::operator*=(double a)
{
  for (int i=0;i<dim;i++)
    val[i]*=a;
  return *this;
}
VEC &VEC::operator/=(double a)
{
  for (int i=0;i<dim;i++)
    val[i]/=a;
  return *this;
}
VEC VEC::operator+(const VEC v1)
{
  VEC s(*this);
  for (int i=0;i<dim;i++) s.val[i]+=v1.val[i];
  return s;
}
VEC VEC::operator-(const VEC v1)
{
  VEC s(*this);
  for (int i=0;i<dim;i++) s.val[i]-=v1.val[i];
  return s;
}
double VEC::operator*(VEC v1)
{
  double temp=0.0;
  for(int i=0;i<dim;i++)
    temp += (val[i]*v1.val[i]);
  return temp;
}
VEC VEC::operator*(double a)
{
  VEC s(*this);
  for(int i=0;i<dim;i++) s.val[i] *= a;
  return s;
}

VEC VEC::operator/(double a)
{
  VEC s(*this);
  for(int i=0;i<dim;i++) s.val[i] /= a;
  return s;
}
double &VEC::operator[](int n)
{
  if (n<0) n=0;
  else if (n>=dim) n = dim - 1;
  return val[n];
}
VEC operator*(double a, const VEC v1)
{
  VEC s(v1.dim);
  for (int i =0; i < v1.dim; i++) s.val[i]=a * v1.val[i];
  return s;
}
VEC *newVEC(int n)
{
  VEC *vptr;
  vptr=(VEC *)malloc(sizeof(VEC));
  vptr -> dim=n;
  vptr -> val=(double *)calloc(n,sizeof(double));
  return vptr;
}
