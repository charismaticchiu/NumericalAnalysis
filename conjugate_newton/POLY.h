#ifndef POLY_H
#define POLY_H
#include "VEC.h"
#include "Complex.h"

VEC polyRoots(Complex x,VEC &A,int maxiter, double eps);
#endif