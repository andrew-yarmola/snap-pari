#ifndef _complex_
#define _complex_

#ifdef __DECCXX
#define exception math_exeption
#include <math.h>
#undef exception
#include <stdcomp>
#else
#include <math.h>
#endif

#include <stdio.h>
#include <iostream>

#ifndef THINK_C
typedef unsigned char   Boolean;
#endif

typedef int FuncResult;
enum
{
    func_OK = 0,
    func_cancelled,
    func_failed,
    func_bad_input
};

class trace_init {
 public:
  trace_init(const char* name) { printf("global init in: %s\n", name); }
};

/* #define all our real constants here too. */

#define HALF              0.5
#define QUARTER           0.25
#define THREE_OVER_TWO    1.5  
#define ZERO              0.0 
#define ONE               1.0 
#define TWO               2.0
#define FOUR              4.0

#define FUDGE             0.6 

#ifndef  PI 
#define PI               3.14159265358979323846
#endif
#define TWO_PI           6.28318530717958647693
#define FOUR_PI         12.56637061435917295385
#define PI_OVER_2        1.57079632679489661923
#define PI_OVER_3        1.04719755119659774615
#define THREE_PI_OVER_2  4.71238898038468985769
#define ROOT_3_OVER_2    0.86602540378443864676
#define ROOT_3           1.73205080756887729352

#define TRUE            1
#define FALSE           0

class Complex
{
 public:
  double real,imag;

  Complex() {}
  Complex(double r, double i) : real(r), imag(i) {}

  double  operator [] (int i) const { return i?imag:real; }
  double& operator [] (int i) { return i?imag:real; }

  double arg() const { return atan2(imag,real); }
};

extern Complex complex_minus(Complex z0, Complex z1);
extern Complex complex_plus(Complex z0, Complex z1);
extern Complex complex_mult(Complex z0, Complex z1);
extern Complex complex_div(Complex z0, Complex z1);
extern Complex complex_sqrt(Complex z);
extern Complex complex_conjugate(Complex z);
extern Complex complex_negate(Complex z);
extern Complex complex_real_mult(double r, Complex z);
extern Complex complex_exp(Complex z);
extern Complex complex_log(Complex z, double approx_arg);

extern Complex complex_cosh(const Complex& z);
extern Complex complex_sinh(const Complex& z);
extern Complex complex_acosh(Complex z);

/* The following version of cross ratio gives cr(0, infinity; a, 1) = a. */
extern Complex cross_ratio(Complex p, Complex q, Complex r, Complex s);

extern double  complex_modulus(Complex z);
extern double  complex_modulus_squared(Complex z);

extern Boolean complex_nonzero(Complex z);
extern Boolean complex_infinite(Complex z);
extern Boolean complex_big(const Complex& z);
extern Boolean complex_small(const Complex& z, double eps = 1e-10);
extern Boolean complex_close(const Complex& a, const Complex& b, double eps = 1e-10);
extern Boolean same_point(const Complex& a, const Complex& b, double eps = 1e-10);

extern int uhs_compare(const Complex& a, const Complex& b, double eps = 1e-10);
extern Boolean uhs_less(const Complex& a, const Complex& b, double eps = 1e-10);

extern void print(const Complex& z, std::ostream& out = std::cout);
extern Boolean read(Complex& z, FILE *in = stdin); 
extern Boolean complex_from_string(const char* s, Complex& z);

extern const Complex Zero; 
extern const Complex One;
extern const Complex Two;
extern const Complex Four;
extern const Complex MinusOne;
extern const Complex I;
extern const Complex MinusI;
extern const Complex TwoPiI;
extern const Complex Infinity;

inline Complex operator + (const Complex& a, const Complex& b)
{
  return complex_plus(a,b); 
}
inline Complex operator - (const Complex& a, const Complex& b)
{
  return complex_minus(a,b); 
}
inline Complex operator * (const Complex& a, const Complex& b)
{
  return complex_mult(a,b); 
}
inline Complex operator * (double r, const Complex& a)
{
  return complex_real_mult(r,a); 
}
inline Complex operator / (const Complex& a, const Complex& b)
{
  return complex_div(a,b); 
}
inline Complex operator / (const Complex& a, double r)
{
  return Complex(a.real/r, a.imag/r); 
}
inline Complex operator - (const Complex& a)
{
  return complex_minus(Zero,a); 
}
inline Complex operator += (Complex& a, const Complex& b)
{
  return a = complex_plus(a,b); 
}
inline Complex operator -= (Complex& a, const Complex& b)
{
  return a = complex_minus(a,b); 
}
inline Complex operator *= (Complex& a, const Complex& b)
{
  return a = complex_mult(a,b); 
}
inline Complex operator /= (Complex& a, const Complex& b)
{
  return a = complex_div(a,b); 
}
inline Boolean operator == (const Complex& a, const Complex& b)
{
  return same_point(a,b); 
}
inline Boolean operator != (const Complex& a, const Complex& b)
{
  return !same_point(a,b); 
}
inline Complex cross_ratio(const Complex a[])
{ 
  return cross_ratio(a[0],a[1],a[2],a[3]); 
}

// Fixed width printing of doubles and Complexes. 

void fwprint(FILE* fp, double x, int w, char pad=' ');
void fwprint(FILE* fp, const Complex& z, int w);
void fwprint(std::ostream& out, double x, int w, char pad=' ');
void fwprint(std::ostream& out, const Complex& z, int w);

std::ostream& operator << (std::ostream& out, const Complex& z);

inline Boolean operator < (const Complex& a, const Complex& b)
{
  return (a.real < b.real - 1e-10) || 
    (a.real <= b.real + 1e-10 && a.imag < b.imag - 1e-10);
}

#endif

/* $Id: complex.h,v 1.1.1.1 2006/06/24 09:08:58 oag Exp $ */ 
