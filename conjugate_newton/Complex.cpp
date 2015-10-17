// Function definitions
#include <math.h>
#include "Complex.h"

Complex::Complex(double r,double i)
{
    x=r; y=i; 
}
Complex::Complex(const Complex &z)
{
    x=z.x; y=z.y;
}
double Complex::r() const
{
return x; 
}
double Complex::i() const
{
return y;
}

Complex & Complex::operator+=(Complex &z)
{
    x+=z.x;
    y+=z.y;
    return *this;
}
Complex & Complex::operator+=(double d1)
{
	x+=d1;
    return *this;
}
Complex & Complex::operator-=(Complex &z)
{
    x-=z.x;
    y-=z.y;
    return *this;
}
Complex & Complex::operator-=(double d1)
{
    x-=d1;
    return *this;
}
Complex & Complex::operator*=(Complex &z)
{
    x=x*z.x-y*z.y;
    y=x*z.y+y*z.x;
    return *this;
}
Complex & Complex::operator*=(double d1)
{
    x=x*d1;
    y=y*d1;
    return *this;
}
Complex & Complex::operator/=(Complex &z)
{
    double p = z.x*z.x+z.y*z.y;
    x = (x*z.x-y*z.y)/p;
    y = (x*z.y+y*z.x)/p;
    return *this;
}

Complex & Complex::operator/=(double d1)
{
    x = x/d1;
    y = y/d1;
    return *this;
}
Complex operator+(Complex z)
{
    Complex z1(z);
	return z1; 
}
Complex operator+(Complex z1,Complex z2)
{
    Complex z(z1);
    return z+=z2;
}
Complex operator+(Complex z1,double d)
{
    Complex z(z1);
    return z+=d;
}
Complex operator+(double d,Complex z1)
{
    Complex z(z1);
    return z+=d;
}
Complex operator-(Complex z1,Complex z2)
{
    Complex z(z1);
    return z-=z2;
}
Complex operator-(double db1,Complex z2)
{
    Complex z(db1-z2.r(),z2.i());
    return z;
}
Complex operator-(Complex z1,double db2)
{
    Complex z(z1.r()-db2, z1.i());
    return z;
}
Complex operator*(Complex z1,Complex z2)
{
    double x,y;
    x=z1.r()*z2.r()-z1.i()*z2.i();
    y=z1.r()*z2.i()+z1.i()*z2.r();
    Complex z(x,y);
    return z;
}
Complex operator*(Complex z1,double d1)
{
    
    double x,y;
    x=z1.r()*d1;
    y=z1.i()*d1;
    Complex z(x,y);
    return z;
}
Complex operator*(double d1,Complex z1)
{
    
    double x,y;
    x=z1.r()*d1;
    y=z1.i()*d1;
    Complex z(x,y);
    return z;
}
Complex operator/(Complex z1,Complex z2)
{
    
    double x,y;
    double p = z2.r()*z2.r()+z2.i()*z2.i();
    x = (z1.r()*z2.r()+z1.i()*z2.i())/p;
    y = (-z1.r()*z2.i()+z1.i()*z2.r())/p;
    Complex z(x,y);
    return z;
}
Complex operator/(Complex z1,double db2)
{
    
    double x,y;
    //double p = z2.r()*z2.r()+z2.i()*z2.i();
    x = z1.r()/db2;
    y = z1.i()/db2;
    Complex z(x,y);
    return z;
}
Complex operator/(double db1,Complex z2)
{
    
    double x,y;
    //double p = z2.r()*z2.r()+z2.i()*z2.i();
    //x = z1.r()/db2;
    //y = z1.i()/db2;
    Complex temp(db1,0);
    
    return temp/z2;
}
double fabs(Complex z)
{
    return sqrt(z.r()*z.r()+z.i()*z.i());
}
Complex &operator++(Complex &z1)
{
    return z1+=1.0;
}
Complex operator++(Complex &z1,int)
{
    Complex z(z1);
    ++z1;
    return z;
}
bool operator==(Complex z1,Complex z2)
{
	if (z1.r()==z2.r() && z1.i()==z2.i()) return true;
	return false;
}
