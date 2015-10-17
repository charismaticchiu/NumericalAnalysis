#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "VEC.h"
//#include "MAT.h"
#include "Complex.h"
int main()
{
	Complex a(1,2);
	Complex b(3,4);
	Complex result(a);
	//result= 2./a   ;
	printf("%lf+j%lf\n",result.r(),result.i());
	//a=3.0;
	printf("%lf+j%lf\n",(1.+a).r(),(1.+a).i());
	return 0;
}
