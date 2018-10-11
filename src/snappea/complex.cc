/*
 *  complex.c
 *
 *  This file provides the standard complex arithmetic functions:
 *
 *  Complex complex_minus           (Complex z0, Complex z1),
 *          complex_plus            (Complex z0, Complex z1),
 *          complex_mult            (Complex z0, Complex z1),
 *          complex_div             (Complex z0, Complex z1),
 *          complex_sqrt            (Complex z),
 *          complex_conjugate       (Complex z),
 *          complex_negate          (Complex z),
 *          complex_real_mult       (double r, Complex z),
 *          complex_exp             (Complex z),
 *          complex_log             (Complex z, double approx_arg);
 *  double  complex_modulus         (Complex z);
 *  double  complex_modulus_squared (Complex z);
 *  Boolean complex_nonzero         (Complex z);
 *  Boolean complex_infinite        (Complex z);
 */

#include "complex.h"
#include <float.h>
#include <stdio.h>
#include <string.h>

// const trace_init ti("complex.c"); 
using namespace std;

const Complex Zero      = Complex( 0.0, 0.0);
const Complex One       = Complex( 1.0, 0.0);
const Complex Two       = Complex( 2.0, 0.0);
const Complex Four      = Complex( 4.0, 0.0);
const Complex MinusOne  = Complex(-1.0, 0.0);
const Complex I         = Complex( 0.0, 1.0);
const Complex MinusI    = Complex( 0.0,-1.0); 
const Complex TwoPiI    = Complex( 0.0, TWO_PI);
const Complex Infinity  = Complex(1e34, 0.0);

#define BIG_MODULUS 1e10
#define SMALL_MODULUS 1e-10

Complex complex_plus(
    Complex z0,
    Complex z1)
{
    Complex sum;

    sum.real = z0.real + z1.real;
    sum.imag = z0.imag + z1.imag;

    return sum;
}


Complex complex_minus(
    Complex z0,
    Complex z1)
{
    Complex diff;

    diff.real = z0.real - z1.real;
    diff.imag = z0.imag - z1.imag;

    return diff;
}


Complex complex_div(
    Complex z0,
    Complex z1)
{
    double  mod_sq;
    Complex quotient;

    if (complex_nonzero(z1)) {
      mod_sq            =  z1.real * z1.real  +  z1.imag * z1.imag;
      quotient.real = (z0.real * z1.real  +  z0.imag * z1.imag)/mod_sq;
      quotient.imag = (z0.imag * z1.real  -  z0.real * z1.imag)/mod_sq;
    }
    else
      quotient = Infinity;

    return quotient;
}


Complex complex_mult(
    Complex z0,
    Complex z1)
{
    Complex product;

    product.real    = z0.real * z1.real  -  z0.imag * z1.imag;
    product.imag    = z0.real * z1.imag  +  z0.imag * z1.real;

    return product;
}


Complex complex_sqrt(
    Complex z)
{
    double  mod,
	    arg;
    Complex root;

    mod = sqrt(complex_modulus(z)); /* no need for safe_sqrt() */
    if (mod == 0.0)
	return Zero;
    arg = 0.5 * atan2(z.imag, z.real);
    root.real = mod * cos(arg);
    root.imag = mod * sin(arg);

    return root;
}


Complex complex_conjugate(
    Complex z)
{
    z.imag = - z.imag;

    return z;
}


Complex complex_negate(
    Complex z)
{
    z.real = - z.real;
    z.imag = - z.imag;

    return z;
}


Complex complex_real_mult(
    double  r,
    Complex z)
{
    Complex multiple;

    multiple.real   = r * z.real;
    multiple.imag   = r * z.imag;

    return multiple;
}

Complex complex_cosh(const Complex& z)
{
  double mod1,mod2;

  mod1 = exp(z.real); 
  mod2 = 1.0/mod1;

  return Complex(0.5 * (mod1 + mod2) * cos(z.imag), 0.5 * (mod1 - mod2) * sin(z.imag));
}

Complex complex_sinh(const Complex& z)
{
  double mod1,mod2;

  mod1 = exp(z.real); 
  mod2 = 1.0/mod1;

  return Complex(0.5 * (mod1 - mod2) * cos(z.imag), 0.5 * (mod1 + mod2) * sin(z.imag));
}
		     
Complex complex_exp(
    Complex z)
{
    double  modulus;
    Complex result;

    modulus = exp(z.real);
    result.real = modulus * cos(z.imag);
    result.imag = modulus * sin(z.imag);

    return result;
}


Complex complex_log(
    Complex z,
    double  approx_arg)
{
    Complex result;

    if (z.real == 0.0  &&  z.imag == 0.0)
    {
	printf("log(0 + 0i) encountered\n");
	result.real = - DBL_MAX;
	result.imag = approx_arg;
	return result;
    }

    result.real = 0.5 * log(z.real * z.real + z.imag * z.imag);

    result.imag = atan2(z.imag, z.real);
    while (result.imag - approx_arg > PI)
	result.imag -= TWO_PI;
    while (approx_arg - result.imag > PI)
	result.imag += TWO_PI;

    return result;
}

Complex cross_ratio(Complex p, Complex q, Complex r, Complex s)
{
  if (complex_small(complex_minus(p,r)) ||
      complex_small(complex_minus(q,s)))
    return Infinity; 

  if (complex_big(p)) 
    return complex_div(complex_minus(q,s),complex_minus(q,r));

  if (complex_big(q)) 
    return complex_div(complex_minus(p,r),complex_minus(p,s));

  if (complex_big(r)) 
    return complex_div(complex_minus(q,s),complex_minus(p,s));
  
  if (complex_big(s)) 
    return complex_div(complex_minus(p,r),complex_minus(q,r));

  return complex_div(
	      complex_mult(complex_minus(p,r),complex_minus(q,s)),
	      complex_mult(complex_minus(q,r),complex_minus(p,s))
	      );
}

double complex_modulus(
    Complex z)
{
    return sqrt(z.real * z.real + z.imag * z.imag); /* no need for safe_sqrt() */
}


double complex_modulus_squared(
    Complex z)
{
    return (z.real * z.real + z.imag * z.imag);
}


Boolean complex_nonzero(
    Complex z)
{
    return (z.real || z.imag);
}


Boolean complex_infinite(
    Complex z)
{
    return (z.real == Infinity.real && z.imag == Infinity.imag);
}

Boolean complex_big(const Complex& z)
{
  return (fabs(z.real) > BIG_MODULUS || fabs(z.imag) > BIG_MODULUS); 
}

Boolean complex_small(const Complex& z, double eps)
{
  return (fabs(z.real) < eps && fabs(z.imag) < eps); 
}

Boolean complex_close(const Complex& a, const Complex& b, double eps)
{
  return fabs(a.real-b.real) < eps && fabs(a.imag-b.imag) < eps;
}

Boolean same_point(const Complex& a, const Complex& b, double eps)
{
  if (complex_modulus(a) < 0.5) {
    if (complex_modulus(b) > 2.0) return FALSE; 
    return complex_small(b - a, eps); 
  }
  return complex_small(One/a - One/b, eps);
}

// Equivalent to conjugating a and b as necessary to map them 
// into the UHP then comparing lexicographically on components. 
// -1 if a' < b', 0 if a'==b', 1 if a' > b'. 

int uhs_compare(const Complex& a, const Complex& b, double eps)
{
  double rdiff = a.real - b.real; 
  if (rdiff > eps) return 1; 
  if (rdiff < -eps) return -1; 
  double idiff = fabs(a.imag) - fabs(b.imag); 
  if (idiff > eps) return 1; 
  if (idiff < -eps) return -1; 
  return 0; 
}

extern Boolean uhs_less(const Complex& a, const Complex& b, double eps)
{
  double rdiff = a.real - b.real; 
  if (rdiff < -eps) return TRUE; 
  if (rdiff > eps) return FALSE; 
  return (fabs(a.imag) - fabs(b.imag) < -eps); 
}

void print(Complex const& z, ostream& out)
{
  out << z; 
}

static int ldig0(double x, int w=0)
{
  double cut=10.;
  if (w>2) cut = 10. - 5.*pow(10.,2-w);
  int n = 1;
  if (x < 0.) {
    x = -x; 
  }
  while (x >= cut) { n++; x /= 10.; }
  return n; 
}

// print a double in exactly w characters

void fwprint(FILE* fp, double x, int w, char pad)
{
  int ld = ldig0(x,w); 
  int prec = w - ld - 2; 

  if (x == 0.) {
    fprintf(fp, "%c%.*f", pad, prec, 0.); 
  } else {
    if (x > 0.) fprintf(fp, "%c", pad); 
    fprintf(fp, "%.*f", prec, x); 
  }
}

void fwprint(ostream& out, double x, int w, char pad)
{
  ios::fmtflags old = out.setf(ios::fixed); 
  if (!w) w = out.precision()+2; 
  int ld = ldig0(x,w); 
  int newprec = w - ld - 2; 
  if (newprec < 1) newprec = 1; 
  if (!pad && x >= 0.) newprec++; 
  int oldprec = out.precision(newprec); 
  if (!pad) {
    if (x==0.) {
      out << 0.;
    } else {
      out << x; 
    }
  } else if (x == 0.) {
    out << pad << 0.; 
  } else {
    if (x > 0.) out << pad; 
    out << x; 
  }
  out.precision(oldprec); 
  out.flags(old); 
}

// Print a Complex in exactly 2*w + 2 characters. 

void fwprint(ostream& out, const Complex& z, int w)
{
  fwprint(out, z.real, w); 
  fwprint(out, z.imag, w, '+'); 
  out << "*i"; 
}

void fwprint(FILE* fp, const Complex& z, int w)
{
  fwprint(fp, z.real, w); 
  fwprint(fp, z.imag, w, '+'); 
  fprintf(fp, "*i"); 
}

// Print "infinity" with padw dashes before & after. 

static void print_infinity(ostream& out, int padw)
{
  int i; 
  for (i=0; i<padw; i++) out << '-';
  out << "infinity";
  for (i=0; i<padw; i++) out << '-';
}  

// Print a Complex in exactly 2*out.precision() + 6 chars. 

ostream& operator << (ostream& out, const Complex& z)
{
  if (complex_big(z)) print_infinity(out, out.precision()-1); 
  else fwprint(out, z, out.precision()+2); 
  return out; 
}

Boolean read(Complex& z, FILE *in)
{
  return fscanf(in, " %lf + %lf *i ", &z.real, &z.imag) == 2;
}

Complex complex_acosh(Complex z)
{
  return complex_log(z + complex_sqrt(z * z - One), 0.0);
}

Boolean complex_from_string(const char* s, Complex& z)
{
  char rest[128];
  char *buf = rest; 
  if (strncmp("i",s,1)==0) {
    z = I; 
    return TRUE;
  }
  if (strncmp("-i",s,2)==0) {
    z = MinusI;
    return TRUE;
  }
  z.imag = 0.;
  int narg = sscanf(s, "%lf%s", &z.real, rest);  
  if (narg==0) {
    return FALSE;
  }
  if (narg==1 || strlen(rest)==0) {
    return TRUE;     // 1, -1, etc.
  }
  if (strncmp("i",rest,1)==0 || strncmp("*i",rest,2)==0) { // .5*i, .5i
    z.imag = z.real; z.real = 0.0; 
    return TRUE;
  }
  if (strncmp("+",buf,1)==0) ++buf;             // 1+-i, 1+-2i, etc
  if (strncmp("i",buf,1)==0) {                  // 1+i
    z.imag = 1.0; 
    return TRUE; 
  }
  if (strncmp("-i",buf,2)==0) {                 // 1-i
    z.imag =-1.0;
    return TRUE; 
  }
  char b2[4];
  narg = sscanf(buf, "%lf%3s", &z.imag, b2);     // 1+2.5i, 1+2.5*i
  if (narg == 1 || strlen(b2)==0 || 
      (strncmp("i",b2,1)!=0 && strncmp("*i",b2,2)!=0)) {
    return FALSE; 
  }
  return TRUE; 
}
