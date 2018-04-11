/*************************** wnchyppr.cpp **********************************
* Author:        Agner Fog
* Date created:  2002-10-20
* Last modified: 2013-11-06
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* Calculation of univariate and multivariate Wallenius noncentral 
* hypergeometric probability distribution.
*
* This file contains source code for the class CWalleniusNCHypergeometric 
* and CMultiWalleniusNCHypergeometricMoments defined in stocc.h.
*
* Documentation:
* ==============
* The file stocc.h contains class definitions.
* The file nchyp.pdf, available from www.agner.org/random/theory 
* describes the theory of the calculation methods.
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2002-2013 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <string.h>                    // memcpy function
#include "stocc.h"                     // class definition
#include "erfres.h"                    // table of error function residues (Don't precompile this header)

/***********************************************************************
constants
***********************************************************************/
static const double LN2 = 0.693147180559945309417; // log(2)


/***********************************************************************
Log and Exp functions with special care for small x
***********************************************************************/
// These are functions that involve expressions of the types log(1+x)
// and exp(x)-1. These functions need special care when x is small to
// avoid loss of precision. There are three versions of these functions:
// (1) Assembly version in library randomaXX.lib
// (2) Use library functions log1p and expm1 if available
// (3) Use Taylor expansion if none of the above are available

#ifdef RANDOMA_H  
// (1)
// Assembly library randomaXX.lib is used.
// Nothing to include here.

#elif defined(__GNUC__) || defined (__INTEL_COMPILER) || defined(HAVE_EXPM1)
// (2) 
// Functions log1p(x) = log(1+x) and expm1(x) = exp(x)-1 are available
// in the math libraries of Gnu and Intel compilers 
// and in R.DLL (www.r-project.org).

double pow2_1(double q, double * y0 = 0) {
   // calculate 2^q and (1-2^q) without loss of precision.
   // return value is (1-2^q). 2^q is returned in *y0
   double y, y1;
   q *= LN2;
   if (fabs(q) > 0.1) {
      y = exp(q);                      // 2^q
      y1 = 1. - y;                     // 1-2^q
   }
   else { // Use expm1
      y1 = expm1(q);                   // 2^q-1
      y = y1 + 1;                      // 2^q
      y1 = -y1;                        // 1-2^q
   }
   if (y0) *y0 = y;                    // Return y if not void pointer
   return y1;                          // Return y1
}

double log1mx(double x, double x1) {
   // Calculate log(1-x) without loss of precision when x is small.
   // Parameter x1 must be = 1-x.
   if (fabs(x) > 0.03) {
      return log(x1);
   }
   else { // use log1p(x) = log(1+x)
      return log1p(-x);
   }
}

double log1pow(double q, double x) {
   // calculate log((1-e^q)^x) without loss of precision.
   // Combines the methods of the above two functions.
   double y, y1;

   if (fabs(q) > 0.1) {
      y = exp(q);                      // e^q
      y1 = 1. - y;                     // 1-e^q
   }
   else { // Use expm1
      y1 = expm1(q);                   // e^q-1
      y = y1 + 1;                      // e^q
      y1 = -y1;                        // 1-e^q
   }

   if (y > 0.1) { // (1-y)^x calculated without problem
      return x * log(y1);
   }
   else { // Use log1p
      return x * log1p(-y);
   }
}

#else
// (3)
// Functions log1p and expm1 are not available in MS and Borland compiler
// libraries. Use explicit Taylor expansion when needed.

double pow2_1(double q, double * y0 = 0) {
   // calculate 2^q and (1-2^q) without loss of precision.
   // return value is (1-2^q). 2^q is returned in *y0
   double y, y1, y2, qn, i, ifac;
   q *= LN2;
   if (fabs(q) > 0.1) {
      y = exp(q);
      y1 = 1. - y;
   }
   else { // expand 1-e^q = -summa(q^n/n!) to avoid loss of precision
      y1 = 0;  qn = i = ifac = 1;
      do {
         y2 = y1;
         qn *= q;  ifac *= i++;
         y1 -= qn / ifac;
      }
      while (y1 != y2);
      y = 1.-y1;
   }
   if (y0) *y0 = y;
   return y1;
}

double log1mx(double x, double x1) {
   // Calculate log(1-x) without loss of precision when x is small.
   // Parameter x1 must be = 1-x.
   if (fabs(x) > 0.03) {
      return log(x1);
   }
   else { // expand ln(1-x) = -summa(x^n/n)
      double y, z1, z2, i;
      y = i = 1.;  z1 = 0;
      do {
         z2 = z1;
         y *= x;
         z1 -= y / i++;
      }
      while (z1 != z2);
      return z1;
   }
}

double log1pow(double q, double x) {
   // calculate log((1-e^q)^x) without loss of precision
   // Uses various Taylor expansions to avoid loss of precision
   double y, y1, y2, z1, z2, qn, i, ifac;

   if (fabs(q) > 0.1) {
      y = exp(q);  y1 = 1. - y;
   }
   else { // expand 1-e^q = -summa(q^n/n!) to avoid loss of precision
      y1 = 0;  qn = i = ifac = 1;
      do {
         y2 = y1;
         qn *= q;  ifac *= i++;
         y1 -= qn / ifac;
      }
      while (y1 != y2);
      y = 1. - y1;
   }
   if (y > 0.1) { // (1-y)^x calculated without problem
      return x * log(y1);
   }
   else { // expand ln(1-y) = -summa(y^n/n)
      y1 = i = 1.;  z1 = 0;
      do {
         z2 = z1;
         y1 *= y;
         z1 -= y1 / i++;
      }
      while (z1 != z2);
      return x * z1;
   }
}

#endif

/***********************************************************************
Other shared functions
***********************************************************************/

double LnFacr(double x) {
   // log factorial of non-integer x
   int32 ix = (int32)(x);
   if (x == ix) return LnFac(ix);      // x is integer
   double r, r2, D = 1., f;
   static const double             
      C0 =  0.918938533204672722,      // ln(sqrt(2*pi))
      C1 =  1./12.,
      C3 = -1./360.,
      C5 =  1./1260.,
      C7 = -1./1680.;
   if (x < 6.) {
      if (x == 0 || x == 1) return 0;
      while (x < 6) D *= ++x;
   }
   r  = 1. / x;  r2 = r*r;
   f = (x + 0.5)*log(x) - x + C0 + r*(C1 + r2*(C3 + r2*(C5 + r2*C7)));
   if (D != 1.) f -= log(D);
   return f;
}


double FallingFactorial(double a, double b) {
   // calculates ln(a*(a-1)*(a-2)* ... * (a-b+1))

   if (b < 30 && int(b) == b && a < 1E10) {
      // direct calculation
      double f = 1.;
      for (int i = 0; i < b; i++) f *= a--;
      return log(f);
   }

   if (a > 100.*b && b > 1.) {
      // combine Stirling formulas for a and (a-b) to avoid loss of precision
      double ar = 1./a;
      double cr = 1./(a-b);
      // calculate -log(1-b/a) by Taylor expansion
      double s = 0., lasts, n = 1., ba = b*ar, f = ba;
      do {
         lasts = s;
         s += f/n;
         f *= ba;
         n++;
      }
      while (s != lasts);
      return (a+0.5)*s + b*log(a-b) - b + (1./12.)*(ar-cr)    
         //- (1./360.)*(ar*ar*ar-cr*cr*cr)
         ;
   }
   // use LnFacr function
   return LnFacr(a)-LnFacr(a-b);
}

double Erf (double x) {
   // Calculates the error function erf(x) as a series expansion or
   // continued fraction expansion.
   // This function may be available in math libraries as erf(x)
   static const double rsqrtpi  = 0.564189583547756286948; // 1/sqrt(pi)
   static const double rsqrtpi2 = 1.12837916709551257390;  // 2/sqrt(pi)
   if (x < 0.) return -Erf(-x);
   if (x > 6.) return 1.;
   if (x < 2.4) {
      // use series expansion
      double term;                     // term in summation
      double j21;                      // 2j+1
      double sum = 0;                  // summation
      double xx2 = x*x*2.;             // 2x^2
      int j;  
      term = x;  j21 = 1.;
      for (j=0; j<=50; j++) {          // summation loop
         sum += term;
         if (term <= 1.E-13) break;
         j21 += 2.;
         term *= xx2 / j21;
      }
      return exp(-x*x) * sum * rsqrtpi2;
   }
   else {
      // use continued fraction expansion
      double a, f;
      int n = int(2.25f*x*x - 23.4f*x + 60.84f); // predict expansion degree
      if (n < 1) n = 1;
      a = 0.5 * n;  f = x;
      for (; n > 0; n--) {             // continued fraction loop
         f = x + a / f;
         a -= 0.5;
      }
      return 1. - exp(-x*x) * rsqrtpi / f;
   }
}


int32 FloorLog2(float x) {
   // This function calculates floor(log2(x)) for positive x.
   // The return value is <= -127 for x <= 0.

   union UfloatInt {  // Union for extracting bits from a float
      float f;
      int32 i;
      UfloatInt(float ff) {f = ff;}  // constructor
   };

#if defined(_M_IX86) || defined(__INTEL__) || defined(_M_X64) || defined(__IA64__) || defined(__POWERPC__)
   // Running on a platform known to use IEEE-754 floating point format
   //int32 n = *(int32*)&x;
   int32 n = UfloatInt(x).i;
   return (n >> 23) - 0x7F;
#else
   // Check if floating point format is IEEE-754
   static const UfloatInt check(1.0f);
   if (check.i == 0x3F800000) {
      // We have the standard IEEE floating point format
      int32 n = UfloatInt(x).i;
      return (n >> 23) - 0x7F;
   }
   else {
      // Unknown floating point format
      if (x <= 0.f) return -127;
      return (int32)floor(log(x)*(1./LN2));
   }
#endif
}


int NumSD (double accuracy) {
   // Gives the length of the integration interval necessary to achieve
   // the desired accuracy when integrating/summating a probability 
   // function, relative to the standard deviation
   // Returns an integer approximation to 2*NormalDistrFractile(accuracy/2)
   static const double fract[] = {
      2.699796e-03, 4.652582e-04, 6.334248e-05, 6.795346e-06, 5.733031e-07,
      3.797912e-08, 1.973175e-09, 8.032001e-11, 2.559625e-12, 6.381783e-14};
   int i;
   for (i = 0; i < (int)(sizeof(fract)/sizeof(*fract)); i++) {
      if (accuracy >= fract[i]) break;
   }
   return i + 6;
}

