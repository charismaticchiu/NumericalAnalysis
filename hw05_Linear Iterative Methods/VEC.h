//hw05 Linear Iterative Methods
//100060007 Ming-Chang Chiu
//Date: 2015.4.5
//vector class
#ifndef VEC_H
#define VEC_H
class VEC {
private:
	int dim;				   	//vector length
	double *val;			      	//array to store vector
public:
	VEC(int n);		      	//uninit constructor, val set to 0
	VEC(const VEC &v1);			//copy constructor
	VEC(int n,double *v);		      //init constructor
	~VEC();			      			//destructor
	int len();	      			//dimension of the vector
	VEC operator-();	      	//unary operator, negative value
	VEC &operator=(const VEC v1);	//assigment
	VEC &operator+=(const VEC v1);	//V += v1
	VEC &operator-=(const VEC v1);	//V -= v1
	VEC &operator*=(double a);		//V *= a
	VEC &operator/=(double a);		//V /= a
	VEC operator+(const VEC v1);	//V + v1
	VEC operator-(const VEC v1);	//V - v1
	double operator*(VEC v1);		//innter product
	VEC operator*(double a);		//V * a
	VEC operator/(double a);		//V / a
	double &operator[](int n);		//indexing
	friend VEC operator*(double a,const VEC v1);	//a * V
	friend VEC *newVEC(int n);		//alloc memory for VEC
};
VEC operator*(double a, const VEC v1);
VEC *newVEC(int n);
#endif
