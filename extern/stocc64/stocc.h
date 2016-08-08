/*****************************   stocc.h   **********************************
* Author:        Agner Fog
* Date created:  2004-01-08
* Last modified: 2013-09-20
* Project:       randomc.h
* Source URL:    www.agner.org/random
*
* Description:
* This file contains function prototypes and class declarations for the C++ 
* library of non-uniform random number generators. Most functions are fast and 
* accurate, even for extreme values of the parameters.
*
*
* functions without classes:
* ==========================
*
* void EndOfProgram(void);
* System-specific exit code. You may modify this to make it fit your
* user interface.
*
* void FatalError(const char * ErrorText);
* Used for outputting error messages from the other functions and classes.
* You may have to modify this function to make it fit your user interface.
*
* double Erf (double x);
* Calculates the error function, which is the integral of the normal distribution.
*
* double LnFac(int64_t n);
* Calculates the natural logarithm of the factorial of n.
*
*
* class StochasticLib1:
* ====================
* This class can be derived from any of the uniform random number generators
* defined in randomc.h. StochasticLib1 provides the following non-uniform random 
* variate generators:
*
* int64_t Bernoulli(double p);
* Bernoulli distribution. Gives 0 or 1 with probability 1-p and p.
*
* double Normal(double m, double s);
* Normal distribution with mean m and standard deviation s.
*
* double NormalTrunc(double m, double s, double limit);
* Truncated normal distribution with tails cut off at m +/- limit
*
* int64_t Poisson (double L);
* Poisson distribution with mean L.
*
* int64_t Binomial (int64_t n, double p);
* Binomial distribution. n trials with probability p.
*
* int64_t Hypergeometric (int64_t n, int64_t m, int64_t N);
* Hypergeometric distribution. Taking n items out N, m of which are colored.
*
* void Multinomial (int64_t * destination, double * source, int64_t n, int64_t colors);
* void Multinomial (int64_t * destination, int64_t * source, int64_t n, int64_t colors);
* Multivariate binomial distribution.
*
* void MultiHypergeometric (int64_t * destination, int64_t * source, int64_t n, int64_t colors);
* Multivariate hypergeometric distribution.
*
* void Shuffle(int64_t * list, int64_t min, int64_t n);
* Shuffle a list of integers.
*
*
* class StochasticLib2:
* =====================
* This class is derived from class StochasticLib1. It redefines the functions
* Poisson, Binomial and HyperGeometric.
* In StochasticLib1, these functions are optimized for being called with 
* parameters that vary. In StochasticLib2, the same functions are optimized
* for being called repeatedly with the same parameters. If your parameters
* seldom vary, then StochasticLib2 is faster. The two classes use different
* calculation methods, both of which are accurate.
*
*
* class StochasticLib3:
* =====================
* This class can be derived from either StochasticLib1 or StochasticLib2, 
* whichever is preferred. It contains functions for generating variates with
* the univariate and multivariate Wallenius' and Fisher's noncentral
* hypergeometric distributions.
*
* int64_t WalleniusNCHyp (int64_t n, int64_t m, int64_t N, double odds);
* Sampling from Wallenius' noncentral hypergeometric distribution, which is 
* what you get when taking n items out N, m of which are colored, without 
* replacement, with bias.
*
* int64_t FishersNCHyp (int64_t n, int64_t m, int64_t N, double odds);
* Sampling from Fisher's noncentral hypergeometric distribution which is the
* conditional distribution of independent binomial variates given their sum n.
*
* void MultiWalleniusNCHyp (int64_t * destination, int64_t * source, double * weights, int64_t n, int64_t colors);
* Sampling from multivariate Wallenius' noncentral hypergeometric distribution.
*
* void MultiFishersNCHyp (int64_t * destination, int64_t * source, double * weights, int64_t n, int64_t colors);
* Sampling from multivariate Fisher's noncentral hypergeometric distribution.
*
*
* Uniform random number generators (integer and float) are also available, as
* these are inherited from the random number generator class that is the base
* class of StochasticLib1.
*
*
* class CWalleniusNCHypergeometric
* ================================
* This class implements various methods for calculating the probability 
* function and the mean and variance of the univariate Wallenius' noncentral 
* hypergeometric distribution. It is used by StochasticLib3 and can also be 
* used independently.
*
*
* class CMultiWalleniusNCHypergeometric
* =====================================
* This class implements various methods for calculating the probability func-
* tion and the mean of the multivariate Wallenius' noncentral hypergeometric
* distribution. It is used by StochasticLib3 and can also be used independently.
*
*
* class CMultiWalleniusNCHypergeometricMoments
* ============================================
* This class calculates the exact mean and variance of the multivariate
* Wallenius' noncentral hypergeometric probability distribution.
*
*
* class CFishersNCHypergeometric
* ==============================
* This class calculates the probability function and the mean and variance 
* of Fisher's noncentral hypergeometric distribution.
*
*
* class CMultiFishersNCHypergeometric
* ===================================
* This class calculates the probability function and the mean and variance 
* of the multivariate Fisher's noncentral hypergeometric distribution.
*
*
* source code:
* ============
* The code for EndOfProgram and FatalError is found in the file userintf.cpp.
* The code for the functions in StochasticLib1 is found in the file stoc1.cpp.
* The code for the functions in StochasticLib2 is found in the file stoc2.cpp.
* The code for the functions in StochasticLib3 is found in the file stoc3.cpp.
* The code for the functions in CWalleniusNCHypergeometric, 
* CMultiWalleniusNCHypergeometric and CMultiWalleniusNCHypergeometricMoments
* is found in the file wnchyppr.cpp.
* The code for the functions in CFishersNCHypergeometric and 
* CMultiFishersNCHypergeometric is found in the file fnchyppr.cpp
* LnFac is found in stoc1.cpp.
* Erf is found in wnchyppr.cpp.
*
*
* Examples:
* =========
* The file ex-stoc.cpp contains an example of how to use this class library.
*
* The file ex-cards.cpp contains an example of how to shuffle a list of items.
*
* The file ex-lotto.cpp contains an example of how to generate a sequence of
* random integers where no number can occur more than once.
*
* The file testbino.cpp contains an example of sampling from the binomial distribution.
*
* The file testhype.cpp contains an example of sampling from the hypergeometric distribution.
*
* The file testpois.cpp contains an example of sampling from the poisson distribution.
*
* The file testwnch.cpp contains an example of sampling from Wallenius noncentral hypergeometric distribution.
*
* The file testfnch.cpp contains an example of sampling from Fisher's noncentral hypergeometric distribution.
*
* The file testmwnc.cpp contains an example of sampling from the multivariate Wallenius noncentral hypergeometric distribution.
*
* The file testmfnc.cpp contains an example of sampling from the multivariate Fisher's noncentral hypergeometric distribution.
*
* The file evolc.zip contains examples of how to simulate biological evolution using this class library.
*
*
* Documentation:
* ==============
* The file ran-instructions.pdf contains further documentation and 
* instructions for these random number generators.
*
* The file distrib.pdf contains definitions of the standard statistic distributions:
* Bernoulli, Normal, Poisson, Binomial, Hypergeometric, Multinomial, MultiHypergeometric.
*
* The file sampmet.pdf contains theoretical descriptions of the methods used
* for sampling from these distributions.
*
* The file nchyp.pdf, available from www.agner.org/random/, contains
* definitions of the univariate and multivariate Wallenius and Fisher's 
* noncentral hypergeometric distributions and theoretical explanations of 
* the methods for calculating and sampling from these.
*
* Copyright 2004-2013 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*******************************************************************************/

#ifndef STOCC_H
#define STOCC_H

#include <math.h>
#include "randomc.h"

#ifdef R_BUILD
   #include "stocR.h"           // Include this when building R-language interface
#endif


/***********************************************************************
 Choose which uniform random number generator to base these classes on
***********************************************************************/

// STOC_BASE defines which base class to use for the non-uniform
// random number generator classes StochasticLib1, 2, and 3.
#ifndef STOC_BASE
   #ifdef R_BUILD
      // Inherit from StocRBase when building for R-language interface
      #define STOC_BASE StocRBase
   #else
      // #define STOC_BASE CRandomMersenne     // C++ Mersenne Twister
      #include "mt_wrapper.h"
      #define STOC_BASE MTWrapper     // C++ Mersenne Twister
      // Or choose any other random number generator base class, for example:
      //#include "randoma.h"
      //#define STOC_BASE CRandomSFMTA      // Binary library SFMT generator
   #endif
#endif

/***********************************************************************
         Other simple functions
***********************************************************************/

double LnFac(int64_t n);               // log factorial (stoc1.cpp)
double LnFacr(double x);               // log factorial of non-integer (wnchyppr.cpp)
double FallingFactorial(double a, double b); // Falling factorial (wnchyppr.cpp)
double Erf (double x);                 // error function (wnchyppr.cpp)
int64_t FloorLog2(float x);            // floor(log2(x)) for x > 0 (wnchyppr.cpp)
int64_t NumSD (double accuracy);           // used internally for determining summation interval


/***********************************************************************
         Constants and tables
***********************************************************************/

// Maximum number of colors in the multivariate distributions
#ifndef MAXCOLORS
   #define MAXCOLORS 32                // You may change this value
#endif

// constant for LnFac function:
static const int64_t FAK_LEN = 1024;       // length of factorial table

// The following tables are tables of residues of a certain expansion
// of the error function. These tables are used in the Laplace method
// for calculating Wallenius' noncentral hypergeometric distribution.
// There are ERFRES_N tables covering desired precisions from
// 2^(-ERFRES_B) to 2^(-ERFRES_E). Only the table that matches the
// desired precision is used. The tables are defined in erfres.h which
// is included in wnchyppr.cpp.

// constants for ErfRes tables:
static const int64_t ERFRES_B = 16;        // begin: -log2 of lowest precision
static const int64_t ERFRES_E = 40;        // end:   -log2 of highest precision
static const int64_t ERFRES_S =  2;        // step size from begin to end
static const int64_t ERFRES_N = (ERFRES_E-ERFRES_B)/ERFRES_S+1; // number of tables
static const int64_t ERFRES_L = 48;        // length of each table

// tables of error function residues:
extern "C" double ErfRes [ERFRES_N][ERFRES_L];

// number of std. deviations to include in integral to obtain desired precision:
extern "C" double NumSDev[ERFRES_N];


/***********************************************************************
         Class StochasticLib1
***********************************************************************/

class StochasticLib1 : public STOC_BASE {
   // This class encapsulates the random variate generating functions.
   // May be derived from any of the random number generators.
public:
   StochasticLib1 (int64_t seed);          // Constructor
   int64_t Bernoulli(double p);            // Bernoulli distribution
   double Normal(double m, double s);  // Normal distribution
   double NormalTrunc(double m, double s, double limit); // Truncated normal distribution
   int64_t Poisson (double L);         // Poisson distribution
   int64_t Binomial (int64_t n, double p); // Binomial distribution
   int64_t Hypergeometric (int64_t n, int64_t m, int64_t N); // Hypergeometric distribution
   void Multinomial (int64_t * destination, double * source, int64_t n, int64_t colors); // Multinomial distribution
   void Multinomial (int64_t * destination, int64_t * source, int64_t n, int64_t colors);// Multinomial distribution
   void MultiHypergeometric (int64_t * destination, int64_t * source, int64_t n, int64_t colors); // Multivariate hypergeometric distribution
   void Shuffle(int64_t * list, int64_t min, int64_t n); // Shuffle integers

   // functions used internally
protected:
   static double fc_lnpk(int64_t k, int64_t N_Mn, int64_t M, int64_t n); // used by Hypergeometric

   // subfunctions for each approximation method
   int64_t PoissonInver(double L);                         // poisson by inversion
   int64_t PoissonRatioUniforms(double L);                 // poisson by ratio of uniforms
   int64_t PoissonLow(double L);                           // poisson for extremely low L
   int64_t BinomialInver (int64_t n, double p);            // binomial by inversion
   int64_t BinomialRatioOfUniforms (int64_t n, double p);  // binomial by ratio of uniforms
   int64_t HypInversionMod (int64_t n, int64_t M, int64_t N);  // hypergeometric by inversion searching from mode
   int64_t HypRatioOfUnifoms (int64_t n, int64_t M, int64_t N);// hypergeometric by ratio of uniforms method

   // Variables specific to each distribution:
   // Variables used by Normal distribution
   double normal_x2;  int64_t normal_x2_valid;

   // Variables used by Hypergeometric distribution
   int64_t  hyp_n_last, hyp_m_last, hyp_N_last;            // Last values of parameters
   int64_t  hyp_mode, hyp_mp;                              // Mode, mode+1
   int64_t  hyp_bound;                                     // Safety upper bound
   double hyp_a;                                           // hat center
   double hyp_h;                                           // hat width
   double hyp_fm;                                          // Value at mode

   // Variables used by Poisson distribution
   double pois_L_last;                                     // previous value of L
   double pois_f0;                                         // value at x=0 or at mode
   double pois_a;                                          // hat center
   double pois_h;                                          // hat width
   double pois_g;                                          // ln(L)
   int64_t  pois_bound;                                    // upper bound

   // Variables used by Binomial distribution
   int64_t bino_n_last;                                    // last n
   double bino_p_last;                                     // last p
   int64_t bino_mode;                                      // mode
   int64_t bino_bound;                                     // upper bound
   double bino_a;                                          // hat center
   double bino_h;                                          // hat width
   double bino_g;                                          // value at mode
   double bino_r1;                                         // p/(1-p) or ln(p/(1-p))
};


/***********************************************************************
Class StochasticLib2
***********************************************************************/

class StochasticLib2 : public StochasticLib1 {
   // derived class, redefining some functions
public:
   int64_t Poisson (double L);                             // Poisson distribution
   int64_t Binomial (int64_t n, double p);                 // Binomial distribution
   int64_t Hypergeometric(int64_t n, int64_t M, int64_t N);// Hypergeometric distribution
   StochasticLib2(int64_t seed):StochasticLib1(seed){};        // Constructor  

   // subfunctions for each approximation method:
protected:
   int64_t PoissonModeSearch(double L);                    // poisson by search from mode
   int64_t PoissonPatchwork(double L);                     // poisson by patchwork rejection
   static double PoissonF(int64_t k, double l_nu, double c_pm); // used by PoissonPatchwork
   int64_t BinomialModeSearch(int64_t n, double p);        // binomial by search from mode
   int64_t BinomialPatchwork(int64_t n, double p);         // binomial by patchwork rejection
   double BinomialF(int64_t k, int64_t n, double l_pq, double c_pm); // used by BinomialPatchwork
   int64_t HypPatchwork (int64_t n, int64_t M, int64_t N); // hypergeometric by patchwork rejection

   // Variables used by Binomial distribution
   int64_t  Bino_k1, Bino_k2, Bino_k4, Bino_k5;
   double Bino_dl, Bino_dr, Bino_r1, Bino_r2, Bino_r4, Bino_r5, 
      Bino_ll, Bino_lr, Bino_l_pq, Bino_c_pm,
      Bino_f1, Bino_f2, Bino_f4, Bino_f5, 
      Bino_p1, Bino_p2, Bino_p3, Bino_p4, Bino_p5, Bino_p6;

   // Variables used by Poisson distribution
   int64_t  Pois_k1, Pois_k2, Pois_k4, Pois_k5;
   double Pois_dl, Pois_dr, Pois_r1, Pois_r2, Pois_r4, Pois_r5, 
      Pois_ll, Pois_lr, Pois_l_my, Pois_c_pm,
      Pois_f1, Pois_f2, Pois_f4, Pois_f5, 
      Pois_p1, Pois_p2, Pois_p3, Pois_p4, Pois_p5, Pois_p6;

   // Variables used by Hypergeometric distribution
   int64_t  Hyp_L, Hyp_k1, Hyp_k2, Hyp_k4, Hyp_k5;
   double Hyp_dl, Hyp_dr, 
      Hyp_r1, Hyp_r2, Hyp_r4, Hyp_r5, 
      Hyp_ll, Hyp_lr, Hyp_c_pm, 
      Hyp_f1, Hyp_f2, Hyp_f4, Hyp_f5, 
      Hyp_p1, Hyp_p2, Hyp_p3, Hyp_p4, Hyp_p5, Hyp_p6;
};


/***********************************************************************
Class StochasticLib3
***********************************************************************/

class StochasticLib3 : public StochasticLib1 {
   // This class can be derived from either StochasticLib1 or StochasticLib2.
   // Adds more probability distributions
public:
   StochasticLib3(int64_t seed);           // Constructor
   void SetAccuracy(double accur);     // Define accuracy of calculations
   int64_t WalleniusNCHyp (int64_t n, int64_t m, int64_t N, double odds); // Wallenius noncentral hypergeometric distribution
   int64_t FishersNCHyp (int64_t n, int64_t m, int64_t N, double odds); // Fisher's noncentral hypergeometric distribution
   void MultiWalleniusNCHyp (int64_t * destination, int64_t * source, double * weights, int64_t n, int64_t colors); // Multivariate Wallenius noncentral hypergeometric distribution
   void MultiComplWalleniusNCHyp (int64_t * destination, int64_t * source, double * weights, int64_t n, int64_t colors); // Multivariate complementary Wallenius noncentral hypergeometric distribution
   void MultiFishersNCHyp (int64_t * destination, int64_t * source, double * weights, int64_t n, int64_t colors); // Multivariate Fisher's noncentral hypergeometric distribution
   // subfunctions for each approximation method
protected:
   int64_t WalleniusNCHypUrn (int64_t n, int64_t m, int64_t N, double odds); // WalleniusNCHyp by urn model
   int64_t WalleniusNCHypInversion (int64_t n, int64_t m, int64_t N, double odds); // WalleniusNCHyp by inversion method
   int64_t WalleniusNCHypTable (int64_t n, int64_t m, int64_t N, double odds); // WalleniusNCHyp by table method
   int64_t WalleniusNCHypRatioOfUnifoms (int64_t n, int64_t m, int64_t N, double odds); // WalleniusNCHyp by ratio-of-uniforms
   int64_t FishersNCHypInversion (int64_t n, int64_t m, int64_t N, double odds); // FishersNCHyp by inversion
   int64_t FishersNCHypRatioOfUnifoms (int64_t n, int64_t m, int64_t N, double odds); // FishersNCHyp by ratio-of-uniforms

   // variables
   double accuracy;                                        // desired accuracy of calculations

   // Variables for Fisher
   int64_t fnc_n_last, fnc_m_last, fnc_N_last;             // last values of parameters
   int64_t fnc_bound;                                      // upper bound
   double fnc_o_last;
   double fnc_f0, fnc_scale;
   double fnc_a;                                           // hat center
   double fnc_h;                                           // hat width
   double fnc_lfm;                                         // ln(f(mode))
   double fnc_logb;                                        // ln(odds)

   // variables for Wallenius
   int64_t wnc_n_last, wnc_m_last, wnc_N_last;             // previous parameters
   double wnc_o_last;
   int64_t wnc_bound1, wnc_bound2;                         // lower and upper bound
   int64_t wnc_mode;                                       // mode
   double wnc_a;                                           // hat center
   double wnc_h;                                           // hat width
   double wnc_k;                                           // probability value at mode
   int64_t UseChopDown;                                        // use chop down inversion instead
   #define WALL_TABLELENGTH  512                           // max length of table
   double wall_ytable[WALL_TABLELENGTH];                   // table of probability values
   int64_t wall_tablen;                                    // length of table
   int64_t wall_x1;                                        // lower x limit for table
};


/***********************************************************************
Class CWalleniusNCHypergeometric
***********************************************************************/

class CWalleniusNCHypergeometric {
   // This class contains methods for calculating the univariate
   // Wallenius' noncentral hypergeometric probability function
public:
   CWalleniusNCHypergeometric(int64_t n, int64_t m, int64_t N, double odds, double accuracy=1.E-8); // constructor
   void SetParameters(int64_t n, int64_t m, int64_t N, double odds); // change parameters
   double probability(int64_t x);                          // calculate probability function
   int64_t MakeTable(double * table, int64_t MaxLength, int64_t * xfirst, int64_t * xlast, double cutoff = 0.); // make table of probabilities
   double mean(void);                                      // approximate mean
   double variance(void);                                  // approximate variance (poor approximation)
   int64_t mode(void);                                     // calculate mode
   double moments(double * mean, double * var);            // calculate exact mean and variance
   int64_t BernouilliH(int64_t x, double h, double rh, StochasticLib1 *sto); // used by rejection method

   // implementations of different calculation methods
protected:
   double recursive(void);                                 // recursive calculation
   double binoexpand(void);                                // binomial expansion of integrand
   double laplace(void);                                   // Laplace's method with narrow integration interval
   double integrate(void);                                 // numerical integration

   // other subfunctions
   double lnbico(void);                                    // natural log of binomial coefficients
   void findpars(void);                                    // calculate r, w, E
   double integrate_step(double a, double b);              // used by integrate()
   double search_inflect(double t_from, double t_to);      // used by integrate()

   // parameters
   double omega;                                           // Odds
   int64_t n, m, N, x;                                     // Parameters
   int64_t xmin, xmax;                                     // Minimum and maximum x
   double accuracy;                                        // Desired precision
   // parameters used by lnbico
   int64_t xLastBico;
   double bico, mFac, xFac;
   // parameters generated by findpars and used by probability, laplace, integrate:
   double r, rd, w, wr, E, phi2d;
   int64_t xLastFindpars;
};


/***********************************************************************
Class CMultiWalleniusNCHypergeometric
***********************************************************************/

class CMultiWalleniusNCHypergeometric {
   // This class encapsulates the different methods for calculating the
   // multivariate Wallenius noncentral hypergeometric probability function
public:
   CMultiWalleniusNCHypergeometric(int64_t n, int64_t * m, double * odds, int64_t colors, double accuracy=1.E-8); // constructor
   void SetParameters(int64_t n, int64_t * m, double * odds, int64_t colors); // change parameters
   double probability(int64_t * x);                        // calculate probability function
   void mean(double * mu);                                 // calculate approximate mean

      // implementations of different calculation methods
protected:
   double binoexpand(void);                                // binomial expansion of integrand
   double laplace(void);                                   // Laplace's method with narrow integration interval
   double integrate(void);                                 // numerical integration

   // other subfunctions
   double lnbico(void);                                    // natural log of binomial coefficients
   void findpars(void);                                    // calculate r, w, E
   double integrate_step(double a, double b);              // used by integrate()
   double search_inflect(double t_from, double t_to);      // used by integrate()

   // parameters
   double * omega;                                         // odds
   double accuracy;                                        // desired accuracy
   int64_t n;                                              // sample size
   int64_t N;                                              // total items in urn
   int64_t * m;                                            // items of each color in urn
   int64_t * x;                                            // items of each color sampled
   int64_t colors;                                             // number of different colors
   int64_t Dummy_align;                                        // filler
   // parameters generated by findpars and used by probability, laplace, integrate:
   double r, rd, w, wr, E, phi2d;
   // generated by lnbico
   double bico;
};


/***********************************************************************
Class CMultiWalleniusNCHypergeometricMoments
***********************************************************************/

class CMultiWalleniusNCHypergeometricMoments: public CMultiWalleniusNCHypergeometric {
   // This class calculates the exact mean and variance of the multivariate
   // Wallenius noncentral hypergeometric distribution by calculating all the 
   // possible x-combinations with probability < accuracy
public:
   CMultiWalleniusNCHypergeometricMoments(int64_t n, int64_t * m, double * odds, int64_t colors, double accuracy=1.E-8) 
      : CMultiWalleniusNCHypergeometric(n, m, odds, colors, accuracy) {};
   double moments(double * mean, double * stddev, int64_t * combinations = 0);

protected:
   // functions used internally
   double loop(int64_t n, int64_t c);                          // recursive loops
   // data
   int64_t xi[MAXCOLORS];                                  // x vector to calculate probability of
   int64_t xm[MAXCOLORS];                                  // rounded approximate mean of x[i]
   int64_t remaining[MAXCOLORS];                           // number of balls of color > c in urn
   double sx[MAXCOLORS];                                   // sum of x*f(x)
   double sxx[MAXCOLORS];                                  // sum of x^2*f(x)
   int64_t sn;                                             // number of combinations
};


/***********************************************************************
Class CFishersNCHypergeometric
***********************************************************************/

class CFishersNCHypergeometric {
   // This class contains methods for calculating the univariate Fisher's
   // noncentral hypergeometric probability function
public:
   CFishersNCHypergeometric(int64_t n, int64_t m, int64_t N, double odds, double accuracy = 1E-8); // constructor
   double probability(int64_t x);                          // calculate probability function
   double probabilityRatio(int64_t x, int64_t x0);         // calculate probability f(x)/f(x0)
   double MakeTable(double * table, int64_t MaxLength, int64_t * xfirst, int64_t * xlast, double cutoff = 0.); // make table of probabilities
   double mean(void);                                      // calculate approximate mean
   double variance(void);                                  // approximate variance
   int64_t mode(void);                                     // calculate mode (exact)
   double moments(double * mean, double * var);            // calculate exact mean and variance

protected:
   double lng(int64_t x);                                  // natural log of proportional function

   // parameters
   double odds;                                            // odds ratio
   double logodds;                                         // ln odds ratio
   double accuracy;                                        // accuracy
   int64_t n, m, N;                                        // Parameters
   int64_t xmin, xmax;                                     // minimum and maximum of x

   // parameters used by subfunctions
   int64_t xLast;
   double mFac, xFac;                                      // log factorials
   double scale;                                           // scale to apply to lng function
   double rsum;                                            // reciprocal sum of proportional function
   int64_t ParametersChanged;
};


/***********************************************************************
Class CMultiFishersNCHypergeometric
***********************************************************************/

class CMultiFishersNCHypergeometric {
   // This class contains functions for calculating the multivariate
   // Fisher's noncentral hypergeometric probability function and its mean and 
   // variance. Warning: the time consumption for first call to 
   // probability or moments is proportional to the total number of
   // possible x combinations, which may be extreme!
public:
   CMultiFishersNCHypergeometric(int64_t n, int64_t * m, double * odds, int64_t colors, double accuracy = 1E-9); // constructor
   double probability(int64_t * x);                        // calculate probability function
   void mean(double * mu);                                 // calculate approximate mean
   void variance(double * var);                            // calculate approximate variance
   double moments(double * mean, double * stddev, int64_t * combinations = 0); // calculate exact mean and variance

protected:
   double lng(int64_t * x);                                // natural log of proportional function
   void SumOfAll(void);                                    // calculates sum of proportional function for all x combinations
   double loop(int64_t n, int64_t c);                          // recursive loops used by SumOfAll
   int64_t n, N;                                           // copy of parameters
   int64_t * m;
   double * odds;
   int64_t colors;
   double logodds[MAXCOLORS];                              // log odds
   double mFac;                                            // sum of log m[i]!
   double scale;                                           // scale to apply to lng function
   double rsum;                                            // reciprocal sum of proportional function
   double accuracy;                                        // accuracy of calculation

   // data used by used by SumOfAll
   int64_t xi[MAXCOLORS];                                  // x vector to calculate probability of
   int64_t xm[MAXCOLORS];                                  // rounded approximate mean of x[i]
   int64_t remaining[MAXCOLORS];                           // number of balls of color > c in urn
   double sx[MAXCOLORS];                                   // sum of x*f(x) or mean
   double sxx[MAXCOLORS];                                  // sum of x^2*f(x) or variance
   int64_t sn;                                             // number of possible combinations of x
};

#endif
