/*************************** stoc3.cpp **********************************
* Author:        Agner Fog
* Date created:  2002-10-02
* Last modified: 2008-11-21
* Project:       stocc.zip
* Source URL:    www.agner.org/random
*
* Description:
* Non-uniform random number generator functions.
*
* This file contains source code for the class StochasticLib3 derived
* from StochasticLib1 or StochasticLib2, defined in stocc.h.
*
* This class implements methods for sampling from the noncentral and extended
* hypergeometric distributions, as well as the multivariate versions of these.
*
* Documentation:
* ==============
* The file stocc.h contains class definitions.
* Further documentation at www.agner.org/random
*
* Copyright 2002-2008 by Agner Fog.
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*****************************************************************************/

#include <string.h>                    // memcpy function
#include "stocc.h"                     // class definitions
//#include "wnchyppr.cpp"              // calculate Wallenius noncentral hypergeometric probability
//#include "fnchyppr.cpp"              // calculate Fisher's noncentral hypergeometric probability


/******************************************************************************
         Methods for class StochasticLib3
******************************************************************************/


/***********************************************************************
             Constructor
***********************************************************************/
StochasticLib3::StochasticLib3(int seed) : StochasticLib1(seed) {
  SetAccuracy(1.E-8);                  // set default accuracy
}


/***********************************************************************
             SetAccuracy
***********************************************************************/
void StochasticLib3::SetAccuracy(double accur) {
  // define accuracy of calculations for
  // WalleniusNCHyp and MultiWalleniusNCHyp
  if (accur < 0.) accur = 0.;
  if (accur > 0.01) accur = 0.01;
  accuracy = accur;
}



/******************************************************************************
            Fisher's noncentral hypergeometric distribution
******************************************************************************/
int32 StochasticLib3::FishersNCHyp (int32 n, int32 m, int32 N, double odds) {
/*
   This function generates a random variate with Fisher's noncentral
   hypergeometric distribution.

   This distribution resembles Wallenius noncentral hypergeometric distribution
   and the two distributions are sometimes confused. A more detailed
   explanation of this distribution is given below under the multivariate
   Fisher's noncentral hypergeometric distribution (MultiFishersNCHyp).
   For further documentation see nchyp.pdf, awailable from www.agner.org/random

   This function uses inversion by chop-down search from zero when parameters
   are small, and the ratio-of-uniforms rejection method when the former
   method would be too slow or would give overflow.
*/
  int32 fak, addd;     // used for undoing transformations
  int32 x;             // result

  // check if parameters are valid
  if (n > N || m > N || n < 0 || m < 0 || odds <= 0.) {
    if (odds == 0.) {
      if (n > N-m) FatalError("Not enough items with nonzero weight in function FishersNCHyp");
      return 0;}
    FatalError("Parameter out of range in function FishersNCHyp");}

  if (odds == 1.) {
    // use hypergeometric function if odds == 1
    return Hypergeometric(n, m, N);}

  // symmetry transformations
  fak = 1;  addd = 0;
  if (m > N/2) {
    // invert m
    m = N - m;
    fak = -1;  addd = n;}

  if (n > N/2) {
    // invert n
    n = N - n;
    addd += fak * m;  fak = - fak;}

  if (n > m) {
    // swap n and m
    x = n;  n = m;  m = x;}

  // cases with only one possible result end here
  if (n == 0 || odds == 0.) return addd;

  if (fak == -1) {
    // reciprocal odds if inverting
    odds = 1. / odds;}

  // choose method
  if (n < 30 && N < 1024 && odds > 1.E-5 && odds < 1.E5) {
    // use inversion by chop down method
    x = FishersNCHypInversion (n, m, N, odds);}

  else {
    // use ratio-of-uniforms method
    x = FishersNCHypRatioOfUnifoms (n, m, N, odds);}

  // undo symmetry transformations
  return x * fak + addd;}


/***********************************************************************
             Subfunctions used by FishersNCHyp
***********************************************************************/

int32 StochasticLib3::FishersNCHypInversion
(int32 n, int32 m, int32 N, double odds) {
/*
  Subfunction for FishersNCHyp distribution.
  Implements Fisher's noncentral hypergeometric distribution by inversion
  method, using chop-down search starting at zero.

  Valid only for 0 <= n <= m <= N/2.
  Without overflow check the parameters must be limited to n < 30, N < 1024,
  and 1.E-5 < odds < 1.E5. This limitation is acceptable because this method
  is slow for higher n.

  The execution time of this function grows with n.

  See the file nchyp.pdf for theoretical explanation.
*/
  static int32 fnc_n_last = -1, fnc_m_last = -1, fnc_N_last = -1;
  static double   fnc_o_last = -1, fnc_f0, fnc_scale;

   int32 x;                       // x value
   int32 L;                       // derived parameter
   double f;                      // scaled function value
   double sum;                    // scaled sum of function values
   double a1, a2, b1, b2, f1, f2; // factors in recursive calculation
   double u;                      // uniform random variate

  L = N-m-n;

  if (n != fnc_n_last || m != fnc_m_last || N != fnc_N_last || odds != fnc_o_last) {
    // parameters have changed. set-up
    fnc_n_last = n; fnc_m_last = m; fnc_N_last = N; fnc_o_last = odds;

    // f(0) is set to an arbitrary value because it cancels out.
    // A low value is chosen to avoid overflow.
    fnc_f0 = 1.E-100;

    // calculate summation of e(x), using the formula:
    // f(x) = f(x-1) * (m-x+1)*(n-x+1)*odds / (x*(L+x))
    // All divisions are avoided by scaling the parameters
    sum = f = fnc_f0;  fnc_scale = 1.;
    a1 = m;  a2 = n;  b1 = 1;  b2 = L + 1;
    for (x = 1; x <= n; x++) {
      f1 = a1 * a2 * odds;
      f2 = b1 * b2;
      a1--;  a2--;  b1++;  b2++;
      f *= f1;
      sum *= f2;
      fnc_scale *= f2;
      sum += f;
      // overflow check. not needed if parameters are limited:
      // if (sum > 1E100) {sum *= 1E-100; f *= 1E-100; fnc_scale *= 1E-100;}
      }
    fnc_f0 *= fnc_scale;
    fnc_scale = sum;
    // now f(0) = fnc_f0 / fnc_scale.
    // We are still avoiding all divisions by saving the scale factor
    }

  // uniform random
  u = Random() * fnc_scale;

  // recursive calculation:
  // f(x) = f(x-1) * (m-x+1)*(n-x+1)*odds / (x*(L+x))
  f = fnc_f0;  x = 0;  a1 = m;  a2 = n;  b1 = 0;  b2 = L;
  do {
    u -= f;
    if (u <= 0) break;
    x++;  b1++;  b2++;
    f *= a1 * a2 * odds;
    u *= b1 * b2;
    // overflow check. not needed if parameters are limited:
    // if (u > 1.E100) {u *= 1E-100;  f *= 1E-100;}
    a1--;  a2--;}
  while (x < n);
  return x;}


int32 StochasticLib3::FishersNCHypRatioOfUnifoms
(int32 n, int32 m, int32 N, double odds) {
/*
  Subfunction for FishersNCHyp distribution.
  Valid for 0 <= n <= m <= N/2, odds != 1

  Fisher's noncentral hypergeometric distribution by ratio-of-uniforms
  rejection method.

  The execution time of this function is almost independent of the parameters.
*/
  static int32 fnc_n_last = -1, fnc_m_last = -1, fnc_N_last = -1; // previous parameters
  static double fnc_o_last = -1;
  static int32 fnc_bound;              // upper bound
  static double fnc_a;                 // hat center
  static double fnc_h;                 // hat width
  static double fnc_lfm;               // ln(f(mode))
  static double fnc_logb;              // ln(odds)
  int32 L;                             // N-m-n
  int32 mode;                          // mode
  double mean;                         // mean
  double variance;                     // variance
  double x;                            // real sample
  int32 k;                             // integer sample
  double u;                            // uniform random
  double lf;                           // ln(f(x))
  double AA, BB, g1, g2;               // temporary

  L = N-m-n;

  if (n != fnc_n_last || m != fnc_m_last || N != fnc_N_last || odds != fnc_o_last) {
    // parameters have changed. set-up
    fnc_n_last = n;  fnc_m_last = m;  fnc_N_last = N;  fnc_o_last = odds;

    // find approximate mean
    AA = (m+n)*odds+L; BB = sqrt(AA*AA - 4*odds*(odds-1)*m*n);
    mean = (AA-BB)/(2*(odds-1));

    // find approximate variance
    AA = mean * (m-mean); BB = (n-mean)*(mean+L);
    variance = N*AA*BB/((N-1)*(m*BB+(n+L)*AA));

    // compute log(odds)
    fnc_logb = log(odds);

    // find center and width of hat function
    fnc_a = mean + 0.5;
    fnc_h = 1.028 + 1.717*sqrt(variance+0.5) + 0.032*fabs(fnc_logb);

    // find safety bound
    fnc_bound = (int32)(mean + 4.0 * fnc_h);
    if (fnc_bound > n) fnc_bound = n;

    // find mode
    mode = (int32)(mean);
    g1 =(double)(m-mode)*(n-mode)*odds;
    g2 =(double)(mode+1)*(L+mode+1);
    if (g1 > g2 && mode < n) mode++;

    // value at mode to scale with:
    fnc_lfm = mode * fnc_logb - fc_lnpk(mode, L, m, n);}

  while(1) {
    u = Random();
    if (u == 0) continue;                       // avoid divide by 0
    x = fnc_a + fnc_h * (Random()-0.5)/u;
    if (x < 0. || x > 2E9) continue;            // reject, avoid overflow
    k = (int32)(x);                             // truncate
    if (k > fnc_bound) continue;                // reject if outside safety bound
    lf = k*fnc_logb - fc_lnpk(k,L,m,n) - fnc_lfm; // compute function value
    if (u * (4.0 - u) - 3.0 <= lf) break;       // lower squeeze accept
    if (u * (u-lf) > 1.0) continue;             // upper squeeze reject
    if (2.0 * log(u) <= lf) break;}             // final acceptance

  return k;}


/***********************************************************************
     Multivariate Fisher's noncentral hypergeometric distribution
***********************************************************************/
void StochasticLib3::MultiFishersNCHyp (int32 * destination,
int32 * source, double * weights, int32 n, int colors) {
/*
   This function generates a vector of random variates with the
   multivariate Fisher's noncentral hypergeometric distribution.

   This distribution is defined as the conditional distribution of 'colors'
   independent binomial variates
      x[i] = binomial(source[i], p[i])
   on the condition that the sum of all x[i] is n.
   p[i] = r * weights[i] / (1 + r * weights[i]),
   r is an arbitrary scale factor.

   Parameters:
   destination:    An output array to receive the number of balls of each
                   color. Must have space for at least 'colors' elements.
   source:         An input array containing the number of balls of each
                   color in the urn. Must have 'colors' elements.
                   All elements must be non-negative.
   weights:        The odds of each color. Must have 'colors' elements.
                   All elements must be non-negative.
   n:              The number of balls drawn from the urn.
                   Can't exceed the total number of balls with nonzero weight
                   in the urn.
   colors:         The number of possible colors.

   Method: The conditional method is used for generating a sample with the
   approximate distribution. This sample is used as a starting point for
   a Gibbs sampler. The accuracy depends on the number of scans with the
   Gibbs sampler.

   The function will reduce the number of colors, if possible, by eliminating
   colors with zero weight or zero number and pooling together colors with the
   same weight. A symmetry transformation is used if more than half the balls
   are taken. The problem thus reduced is handled in the arrays osource,
   oweights and osample of dimension colors2.
*/
  int order1[MAXCOLORS];       // sort order, index into source and destination
  int order2[MAXCOLORS];       // corresponding index into osource when equal weights pooled together
  int order3[MAXCOLORS];       // secondary index for sorting by variance
  int32 osource[MAXCOLORS];    // contents of source, sorted by weight with equal weights pooled together
  int32 osample[MAXCOLORS];    // balls sampled, sorted by weight
  double oweights[MAXCOLORS];  // sorted list of weights
  double var[MAXCOLORS];       // sorted list of variance
  int32 x = 0;                 // univariate sample
  int32 m;                     // number of items of one color
  int32 m1, m2;                // number of items in each weight group
  int32 msum;                  // total number of items of several or all colors
  int32 n0;                    // remaining balls to sample
  int32 n1, n2;                // sample size for each weight group
  double w = 0.;               // weight or variance of items of one color
  double w1, w2;               // mean weight of each weight group
  double wsum;                 // total weight of all items of several or all colors
  double odds;                 // weight ratio
  int i, j, k;                 // loop counters
  int a, b;                    // limits for weight group
  int c, c1, c2;               // color index
  int colors2;                 // reduced number of colors, number of entries in osource
  int ngibbs;                  // number of scans in Gibbs sampler
  int invert = 0;              // 1 if symmetry transformation used

  // check validity of parameters
  if (n < 0 || colors < 0 || colors > MAXCOLORS) FatalError("Parameter out of range in function MultiFishersNCHyp");
  if (colors == 0) return;
  if (n == 0) {for (i=0; i<colors; i++) destination[i] = 0; return;}

  // check validity of array parameters
  for (i=0, msum=0; i < colors; i++) {
    m = source[i];  w = weights[i];
    if (m < 0 || w < 0) FatalError("Parameter negative in function MultiFishersNCHyp");
    if (w) msum += m;}

  // sort by weight, heaviest first
  for (i=0; i < colors; i++) order1[i] = order3[i] = i;
  for (i=0; i < colors-1; i++) {
    c = order1[i];  k = i;
    w = weights[c];  if (source[c]==0) w = 0;
    for (j=i+1; j < colors; j++) {
      c2 = order1[j];
      if (weights[c2] > w && source[c2]) {
        w = weights[c2];  k = j;}}
      order1[i] = order1[k];  order1[k] = c;}

  // Skip any items with zero weight
  // this solves all problems with zero weights
  while (colors && (weights[c=order1[colors-1]]==0 || source[c]==0)) {
    colors--;  destination[c] = 0;}

  // check if we are taking all, or too many, balls
  if (n >= msum) {
    if (n > msum) FatalError("Taking more items than there are in function MultiFishersNCHyp");
    for (i = 0; i < colors; i++) {c = order1[i];  destination[c] = source[c];}
    return;}

  if (n > msum / 2) {
    // improve accuracy by symmetry transformation
    for (i=0, j=colors-1; i < j; i++, j--) { // reverse order list
      c = order1[i];  order1[i] = order1[j];  order1[j] = c;}
    n = msum - n;  invert = 1;}

  // copy source and weights into ordered lists and pool together colors with same weight
  for (i=0, c2=-1; i < colors; i++) {
    c = order1[i];
    if (i==0 || weights[c] != w) {
      c2++;
      x = source[c];
      oweights[c2] = w = invert ? 1./weights[c] : weights[c];}
    else {
      x += source[c];}
    osource[c2] = x;
    order2[i] = c2;
    osample[c2] = 0;}
  colors2 = c2 + 1;

  // simple cases
  if (colors2 == 1) osample[0] = n;
  if (colors2 == 2) {
    x = FishersNCHyp(n, osource[0], msum, oweights[0]/oweights[1]);
    osample[0] = x;  osample[1] = n - x;}

  if (colors2 > 2) {

    // divide weights into two groups, heavy and light
    a = 0;  b = colors2-1;
    w = sqrt(oweights[0] * oweights[colors2-1]);
    do {
      c = (a + b) / 2;
      if (oweights[c] > w) a = c; else b = c;}
    while (b > a + 1);
    a = 0; // heavy group goes from a to b-1, light group goes from b to colors2-1

    // calculate mean weight for heavy group
    for (i=a, m1=0, wsum=0; i < b; i++) {
      m1 += osource[i];  wsum += oweights[i] * osource[i];}
    w1 = wsum / m1;

    // calculate mean weight for light group
    for (i=b, m2=0, wsum=0; i < colors2; i++) {
      m2 += osource[i];  wsum += oweights[i] * osource[i];}
    w2 = wsum / m2;

    // split sample n into heavy (n1) and light (n2) groups
    n1 = FishersNCHyp(n, m1, m1+m2, w1/w2);
    n2 = n - n1;
    n0 = n1;

    // loop twice, for the two groops
    for (k=0; k < 2; k++) {

      // split group into single colors by calling FishersNCHyp b-a-1 times
      for (i = a; i < b-1; i++) {
        m = osource[i];  w = oweights[i];

        // calculate mean weight of remaining colors
        for (j=i+1, msum=0, wsum=0; j < b; j++) {
          m1 = osource[j];  w1 = oweights[j];
          msum += m1;  wsum += m1 * w1;}

        // split out color i
        if (w == w1) {
          x = Hypergeometric(n0, m, msum + m);}
        else {
          if (wsum == 0) {
            x = n0;}
          else {
            odds = w * msum / wsum;
            x = FishersNCHyp(n0, m, msum + m, odds);}}
        osample[i] += x;
        n0 -= x;}

      // get the last color in the group
      osample[i] += n0;

      // set parameters for second group
      a = b;  b = colors2;  n0 = n2;}

    // calculate variance
    CMultiFishersNCHypergeometric(n, osource, oweights, colors2).variance(var);

    // sort again, this time by variance
    for (i=0; i < colors2-1; i++) {
      c = order3[i];  k = i;
      w = var[c];
      for (j=i+1; j < colors2; j++) {
        c2 = order3[j];
        if (var[c2] > w) {
          w = var[c2];  k = j;}}
      order3[i] = order3[k];  order3[k] = c;}

    // determine number of scans (not fine-tuned):
    ngibbs = 4;  if (accuracy < 1E-6) ngibbs = 6;  if (colors2 > 5) ngibbs++;

    // Gibbs sampler
    for (k = 0; k < ngibbs; k++) {
      for (i = 0; i < colors2; i++) {
        c1 = order3[i];
        j = i + 1;  if (j == colors2) j = 0;
        c2 = order3[j];
        n1 = osample[c1] + osample[c2];
        x = FishersNCHyp(n1, osource[c1], osource[c1]+osource[c2], oweights[c1]/oweights[c2]);
        osample[c1] = x;
        osample[c2] = n1 - x;}}}

  if (invert) {
    // reverse symmetry transformation on result
    for (i=0; i < colors2; i++) {
      osample[i] = osource[i] - osample[i];}}

  // un-sort sample into destination
  for (i=0; i < colors; i++) {
    c1 = order1[i];  c2 = order2[i];
    if (source[c1] == osource[c2]) {
      destination[c1] = osample[c2];}
    else {
      x = Hypergeometric(osample[c2], source[c1], osource[c2]);
      destination[c1] = x;
      osample[c2] -= x;
      osource[c2] -= source[c1];}}}
