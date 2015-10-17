// complex number class
#ifndef COMPLEX_H
#define COMPLEX_H
#include <stdio.h>
#include <stdlib.h>
class Complex {
  public:
    Complex(double,double);     // constructor
    Complex(const Complex&);          // copy constructor
    double r() const;                   // get real part
    double i() const;                   // get imaginary part
    Complex& operator+= (Complex&);     // C1 += C2;
    Complex& operator+= (double);       // C1 += dbl;
    Complex& operator-= (Complex&);     // C1 -= C2;
    Complex& operator-= (double);       // C1 -= dbl;
    Complex& operator*= (Complex&);     // C1 *= C2;
    Complex& operator*= (double);       // C1 *= dbl;
    Complex& operator/= (Complex&);     // C1 *= C2;
    Complex& operator/= (double);       // C1 *= dbl;
  private:
    double x,y;
};

Complex operator+(Complex);             // unary plus 
Complex operator+(Complex,Complex);     // C1 + C2
Complex operator+(Complex,double);      //C1+dbl 
Complex operator+(double,Complex);      //dbl+C1 
Complex operator-(Complex);             // unary minus 
Complex operator-(Complex,Complex);     //C1- C2 
Complex operator-(Complex,double);      //C1- dbl
Complex operator-(double,Complex);      //dbl-C1
Complex operator*(Complex,Complex);     //C1* C2 
Complex operator*(Complex,double);      //C1* dbl
Complex operator*(double,Complex);      //dbl*C1
Complex operator/(Complex,Complex);     //C1/ C2
Complex operator/(Complex,double);      //C1/ dbl
Complex operator/(double,Complex);      //dbl/C1 
double fabs(Complex z);                 // |C1|
Complex& operator++(Complex&);             // prefix increment
Complex operator++(Complex&,int);       // postfix increment
bool operator==(Complex,Complex);        // testing for equality
int operator!=(Complex,Complex); 
#endif
