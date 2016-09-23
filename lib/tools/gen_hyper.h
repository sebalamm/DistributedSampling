#include <mpfr.h>
#include <gmp.h>

#include "definitions.h"

#ifndef _GEN_HYPER_H_
#define _GEN_HYPER_H_

class GenHyper {
    public:
        GenHyper(ULONG seed) {
            mpfr_init2(hyp_n_last, 200);
            mpfr_set_d(hyp_n_last, -1.0, MPFR_RNDN);
            mpfr_init2(hyp_m_last, 200);
            mpfr_set_d(hyp_m_last, -1.0, MPFR_RNDN);
            mpfr_init2(hyp_N_last, 200);
            mpfr_set_d(hyp_N_last, -1.0, MPFR_RNDN);

            mpfr_init2(hyp_h, 200);
            mpfr_set_d(hyp_h, -1.0, MPFR_RNDN);
            mpfr_init2(hyp_a, 200);
            mpfr_set_d(hyp_a, -1.0, MPFR_RNDN);
            mpfr_init2(hyp_fm, 200);
            mpfr_set_d(hyp_fm, -1.0, MPFR_RNDN);
            mpfr_init2(hyp_bound, 200);
            mpfr_set_d(hyp_bound, -1.0, MPFR_RNDN);

            gmp_randinit_mt(rng);
            gmp_randseed_ui(rng, seed);
        }

        void RandomInit(ULONG seed) {
            gmp_randseed_ui(rng, seed);
        }

        ULONG Hypergeometric(ULONG n, ULONG m, ULONG N) {
            ULONG fak, addd;
            ULONG x;

            fak = 1;  addd = 0;
            if (m > N/2) {
               // invert m
               m = N - m;
               fak = -1;  addd = n;
            }    
            if (n > N/2) {
               // invert n
               n = N - n;
               addd += fak * m;  fak = - fak;
            }    
            if (n > m) {
               // swap n and m
               x = n;  n = m;  m = x;
            }    

            mpfr_t nf;
            mpfr_init2(nf, 200);
            mpfr_set_ui(nf, n, MPFR_RNDN);

            mpfr_t mf;
            mpfr_init2(mf, 200);
            mpfr_set_ui(mf, m, MPFR_RNDN);

            mpfr_t Nf;
            mpfr_init2(Nf, 200);
            mpfr_set_ui(Nf, N, MPFR_RNDN);

            mpfr_t hyp;
            mpfr_init2(hyp, 200);
            HypRatioOfUnifoms(nf, mf, Nf, hyp);
            x = mpfr_get_ui(hyp, MPFR_RNDN);

            return x * fak + addd;
        }

    private:
        // ULONG orig_hyp_n_last, orig_hyp_m_last, orig_hyp_N_last;
        // float128 orig_hyp_h, orig_hyp_a, orig_hyp_fm, orig_hyp_bound;

        mpfr_t hyp_n_last, hyp_m_last, hyp_N_last;
        mpfr_t hyp_h, hyp_a, hyp_fm, hyp_bound;
        gmp_randstate_t rng;

        const double SHAT1 = 2.943035529371538573;    // 8/e
        const double SHAT2 = 0.8989161620588987408;   // 3-sqrt(12/e)

        // float128 LnFac(int64_t n) {
        //    // log factorial function. gives natural logarithm of n!

        //    // define constants
        //    static const float128                 // coefficients in Stirling approximation     
        //       C0 =  0.918938533204672722,      // ln(sqrt(2*pi))
        //       C1 =  1./12., 
        //       C3 = -1./360.;
        //    // C5 =  1./1260.,                  // use r^5 term if FAK_LEN < 50
        //    // C7 = -1./1680.;                  // use r^7 term if FAK_LEN < 20
        //    // static variables
        //    // static float128 fac_table[FAK_LEN];   // table of ln(n!):
        //    // static int64_t initialized = 0;         // remember if fac_table has been initialized

        //    // if (n < FAK_LEN) {
        //    //    if (n <= 1) {
        //    //       if (n < 0) FatalError("Parameter negative in LnFac function");  
        //    //       return 0;
        //    //    }
        //    //    if (!initialized) {              // first time. Must initialize table
        //    //       // make table of ln(n!)
        //    //       float128 sum = fac_table[0] = 0.;
        //    //       for (int64_t i=1; i<FAK_LEN; i++) {
        //    //          sum += log(float128(i));
        //    //          fac_table[i] = sum;
        //    //       }
        //    //       initialized = 1;
        //    //    }
        //    //    return fac_table[n];
        //    // }
        //    // not found in table. use Stirling approximation
        //    float128  n1, r;
        //    n1 = n;  r  = 1. / n1;
        //    return (n1 + 0.5)*log(n1) - n1 + C0 + r*(C1 + r*r*C3);
        // }

        void LnFac(mpfr_t n, mpfr_t result) {
           // log factorial function. gives natural logarithm of n!

           // define constants
           static const double                 // coefficients in Stirling approximation     
              C0 =  0.918938533204672722,      // ln(sqrt(2*pi))
              C1 =  1./12., 
              C3 = -1./360.;
           // C5 =  1./1260.,                  // use r^5 term if FAK_LEN < 50
           // C7 = -1./1680.;                  // use r^7 term if FAK_LEN < 20
           // static variables
           // static float128 fac_table[FAK_LEN];   // table of ln(n!):
           // static int64_t initialized = 0;         // remember if fac_table has been initialized

           // if (n < FAK_LEN) {
           //    if (n <= 1) {
           //       if (n < 0) FatalError("Parameter negative in LnFac function");  
           //       return 0;
           //    }
           //    if (!initialized) {              // first time. Must initialize table
           //       // make table of ln(n!)
           //       float128 sum = fac_table[0] = 0.;
           //       for (int64_t i=1; i<FAK_LEN; i++) {
           //          sum += log(float128(i));
           //          fac_table[i] = sum;
           //       }
           //       initialized = 1;
           //    }
           //    return fac_table[n];
           // }
           // not found in table. use Stirling approximation
           // float128  n1, r;
           // n1 = n;  r  = 1. / n1;
           // return (n1 + 0.5)*log(n1) - n1 + C0 + r*(C1 + r*r*C3);

            mpfr_t r;
            mpfr_init2(r, 200);
            mpfr_d_div(r, 1.0, n, MPFR_RNDN);

            mpfr_t logn;
            mpfr_init2(logn, 200);
            mpfr_log(logn, n, MPFR_RNDN);

            mpfr_t sub;
            mpfr_init2(sub, 200);
            mpfr_mul(sub, r, r, MPFR_RNDN);
            mpfr_mul_d(sub, sub, C3, MPFR_RNDN);
            mpfr_add_d(sub, sub, C1, MPFR_RNDN);
            mpfr_mul(sub, sub, r, MPFR_RNDN);
            mpfr_add_d(sub, sub, C0, MPFR_RNDN);
            mpfr_add(sub, sub, n, MPFR_RNDN);

            mpfr_add_d(result, n, 0.5, MPFR_RNDN);
            mpfr_mul(result, result, logn, MPFR_RNDN);
            mpfr_sub(result, result, sub, MPFR_RNDN);
        }

        // float128 fc_lnpk(ULONG k, ULONG L, ULONG m, ULONG n) {
        //     return(LnFac(k) + LnFac(m - k) + LnFac(n - k) + LnFac(L + k));
        // }

        void fc_lnpk(mpfr_t k, mpfr_t L, mpfr_t m, mpfr_t n, mpfr_t result) {
            // return(LnFac(k) + LnFac(m - k) + LnFac(n - k) + LnFac(L + k));
            
            mpfr_t tmp;
            mpfr_init2(tmp, 200);

            mpfr_t lnfac; 
            mpfr_init2(lnfac, 200);

            LnFac(k, lnfac);
            mpfr_set(result, lnfac, MPFR_RNDN);

            mpfr_sub(tmp, m, k, MPFR_RNDN);
            LnFac(tmp, lnfac);
            mpfr_add(result, result, lnfac, MPFR_RNDN);

            mpfr_sub(tmp, n, k, MPFR_RNDN);
            LnFac(tmp, lnfac);
            mpfr_add(result, result, lnfac, MPFR_RNDN);

            mpfr_add(tmp, L, k, MPFR_RNDN);
            LnFac(tmp, lnfac);
            mpfr_add(result, result, lnfac, MPFR_RNDN);
        }

        void HypRatioOfUnifoms (mpfr_t n, mpfr_t m, mpfr_t N, mpfr_t result) {
            /*
            Subfunction for Hypergeometric distribution using the ratio-of-uniforms
            rejection method.

            This code is valid for 0 < n <= m <= N/2.

            The computation time hardly depends on the parameters, except that it matters
            a lot whether parameters are within the range where the LnFac function is
            tabulated.

            Reference: E. Stadlober: "The ratio of uniforms approach for generating
            discrete random variates". Journal of Computational and Applied Mathematics,
            vol. 31, no. 1, 1990, pp. 181-189.
            */

            // Sanity
            // ULONG orig_N = mpfr_get_ui(N, MPFR_RNDN);
            // ULONG orig_m = mpfr_get_ui(m, MPFR_RNDN);
            // ULONG orig_n = mpfr_get_ui(n, MPFR_RNDN);
            // printf("N=%llu ", orig_N);
            // mpfr_out_str (stdout, 10, 0, N, MPFR_RNDD);
            // putchar ('\n');

            // ULONG orig_L;                          // N-m-n
            // ULONG orig_mode;                       // mode
            // ULONG orig_k;                          // integer sample
            // double orig_x;                           // real sample
            // double orig_rNN;                         // 1/(N*(N+2))
            // double orig_my;                          // mean
            // double orig_var;                         // variance
            // double orig_u;                           // uniform random
            // double orig_lf;                          // ln(f(x))
             
            // convert inputs 
            mpfr_t L;                          // N-m-n
            mpfr_t mode;                       // mode
            mpfr_t k;                          // integer sample
            mpfr_t x;                           // real sample
            mpfr_t rNN;                         // 1/(N*(N+2))
            mpfr_t my;                          // mean
            mpfr_t var;                         // variance
            mpfr_t u;                           // uniform random
            mpfr_t lf;                          // ln(f(x))

            mpfr_init2(L, 200);
            mpfr_init2(mode, 200);
            mpfr_init2(k, 200);
            mpfr_init2(x, 200);
            mpfr_init2(rNN, 200);
            mpfr_init2(my, 200);
            mpfr_init2(var, 200);
            mpfr_init2(u, 200);
            mpfr_init2(lf, 200);

            mpfr_t lnfac; 
            mpfr_init2(lnfac, 200);


            // L = N - m - n;
            mpfr_sub(L, N, m, MPFR_RNDN);
            mpfr_sub(L, L, n, MPFR_RNDN);

            // orig_L = orig_N - orig_m - orig_n;
            // printf("L=%llu ", orig_L);
            // mpfr_out_str (stdout, 10, 0, L, MPFR_RNDD);
            // putchar ('\n');

            // if (hyp_N_last != N || hyp_m_last != m || hyp_n_last != n) {
            if (mpfr_cmp(hyp_N_last, N) != 0
                    || mpfr_cmp(hyp_m_last, m) != 0
                    || mpfr_cmp(hyp_n_last, n) != 0) {

                // hyp_N_last = N;  hyp_m_last = m;  hyp_n_last = n;         // Set-up
                mpfr_set(hyp_N_last, N, MPFR_RNDN);
                mpfr_set(hyp_m_last, m, MPFR_RNDN);
                mpfr_set(hyp_n_last, n, MPFR_RNDN);

                // orig_hyp_N_last = orig_N; orig_hyp_m_last = orig_m; orig_hyp_n_last = orig_n;
                // printf("hypN=%llu ", orig_hyp_N_last);
                // mpfr_out_str (stdout, 10, 0, hyp_N_last, MPFR_RNDD);
                // putchar ('\n');

                // rNN = 1. / ((float128)N*(N+2));                             // make two divisions in one
                // N*(N+2)
                mpfr_t denom;
                mpfr_init2(denom, 200);
                mpfr_add_d(denom, N, 2.0, MPFR_RNDN);
                mpfr_mul(denom, denom, N, MPFR_RNDN);
                // 1.0 / denom
                mpfr_d_div(rNN, 1.0, denom, MPFR_RNDN);

                // orig_rNN = 1. / ((float128)orig_N*(orig_N+2));                             // make two divisions in one
                // printf("rNN=%.20f ", orig_rNN);
                // mpfr_out_str (stdout, 10, 0, rNN, MPFR_RNDD);
                // putchar ('\n');

                // my = (float128)n * m * rNN * (N+2);                         // mean = n*m/N
                // n * m * rNN * (N+2)
                mpfr_add_d(my, N, 2.0, MPFR_RNDN);
                mpfr_mul(my, my, rNN, MPFR_RNDN);
                mpfr_mul(my, my, m, MPFR_RNDN);
                mpfr_mul(my, my, n, MPFR_RNDN);

                // orig_my = (float128)orig_n * orig_m * orig_rNN * (orig_N+2);
                // printf("my=%.20f ", orig_my);
                // mpfr_out_str (stdout, 10, 0, my, MPFR_RNDD);
                // putchar ('\n');

                // mode = (int64_t)(float128(n+1) * float128(m+1) * rNN * N);    // mode = floor((n+1)*(m+1)/(N+2))
                // n+1
                mpfr_t nplusone;
                mpfr_init2(nplusone, 200);
                mpfr_set_d(nplusone, 1.0, MPFR_RNDN);
                mpfr_add(nplusone, nplusone, n, MPFR_RNDN);

                // m+1
                mpfr_t mplusone;
                mpfr_init2(mplusone, 200);
                mpfr_set_d(mplusone, 1.0, MPFR_RNDN);
                mpfr_add(mplusone, mplusone, m, MPFR_RNDN);

                // rNN * N * (m+1) * (n+1)
                mpfr_mul(mode, N, rNN, MPFR_RNDN);
                mpfr_mul(mode, mode, mplusone, MPFR_RNDN);
                mpfr_mul(mode, mode, nplusone, MPFR_RNDN);
                mpfr_floor(mode, mode);

                // orig_mode = (int64_t)(float128(orig_n+1) * float128(orig_m+1) * orig_rNN * orig_N);    // mode = floor((n+1)*(m+1)/(N+2))
                // printf("mode=%llu ", orig_mode);
                // mpfr_out_str (stdout, 10, 0, mode, MPFR_RNDD);
                // putchar ('\n');

                // var = (float128)n * m * (N-m) * (N-n) / ((float128)N*N*(N-1));// variance
                // N*N*(N-1)
                mpfr_sub_d(denom, N, 1.0, MPFR_RNDN);
                mpfr_mul(denom, denom, N, MPFR_RNDN);
                mpfr_mul(denom, denom, N, MPFR_RNDN);

                // N-m
                mpfr_t Nminusm;
                mpfr_init2(Nminusm, 200);
                mpfr_sub(Nminusm, N, m, MPFR_RNDN);

                // N-n
                mpfr_t Nminusn;
                mpfr_init2(Nminusn, 200);
                mpfr_sub(Nminusn, N, n, MPFR_RNDN);

                // n * m * (N-m) * (N-n) / denom
                mpfr_t var;
                mpfr_init2(var, 200);
                mpfr_mul(var, n, m, MPFR_RNDN);
                mpfr_mul(var, var, Nminusm, MPFR_RNDN);
                mpfr_mul(var, var, Nminusn, MPFR_RNDN);
                mpfr_div(var, var, denom, MPFR_RNDN);

                // orig_var = (float128)orig_n * orig_m * (orig_N-orig_m) * (orig_N-orig_n) / ((float128)orig_N*orig_N*(orig_N-1));// variance
                // printf("var=%.20f ", orig_var);
                // mpfr_out_str (stdout, 10, 0, var, MPFR_RNDD);
                // putchar ('\n');

                // hyp_h = sqrt(SHAT1 * (var+0.5)) + SHAT2;                  // hat width
                mpfr_t hat;
                mpfr_init2(hat, 200);
                mpfr_add_d(hat, var, 0.5, MPFR_RNDN);
                mpfr_mul_d(hat, hat, SHAT1, MPFR_RNDN);
                mpfr_sqrt(hyp_h, hat, MPFR_RNDN);
                mpfr_add_d(hyp_h, hyp_h, SHAT2, MPFR_RNDN);

                // orig_hyp_h = sqrt(SHAT1 * (orig_var+0.5)) + SHAT2;                  // hat width
                // printf("hyph=%.20f ", orig_hyp_h);
                // mpfr_out_str (stdout, 10, 0, hyp_h, MPFR_RNDD);
                // putchar ('\n');

                // hyp_a = my + 0.5;                                         // hat center
                mpfr_add_d(hyp_a, my, 0.5, MPFR_RNDN);

                // orig_hyp_a = orig_my + 0.5;                                         // hat center
                // printf("hypa=%.20f ", orig_hyp_a);
                // mpfr_out_str (stdout, 10, 0, hyp_a, MPFR_RNDD);
                // putchar ('\n');

                // hyp_fm = fc_lnpk(mode, L, m, n);                          // maximum
                fc_lnpk(mode, L, m, n, lnfac);
                mpfr_set(hyp_fm, lnfac, MPFR_RNDN);

                // orig_hyp_fm = fc_lnpk(orig_mode, orig_L, orig_m, orig_n);                          // maximum
                // printf("hypfm=%.20f ", orig_hyp_fm);
                // mpfr_out_str (stdout, 10, 0, hyp_fm, MPFR_RNDD);
                // putchar ('\n');

                // hyp_bound = (int64_t)(hyp_a + 4.0 * hyp_h);               // safety-bound
                mpfr_mul_d(hyp_bound, hyp_h, 4.0, MPFR_RNDN);
                mpfr_add(hyp_bound, hyp_bound, hyp_a, MPFR_RNDN);
                mpfr_floor(hyp_bound, hyp_bound);

                // orig_hyp_bound = (int64_t)(orig_hyp_a + 4.0 * orig_hyp_h);               // safety-bound
                // printf("hypbound=%.20f ", orig_hyp_bound);
                // mpfr_out_str (stdout, 10, 0, hyp_bound, MPFR_RNDD);
                // putchar ('\n');

                // if (hyp_bound > n) hyp_bound = n;
                if (mpfr_cmp(hyp_bound, n) > 0) mpfr_set(hyp_bound, n, MPFR_RNDN);

                // if (orig_hyp_bound > orig_n) orig_hyp_bound = orig_n;
                // printf("hypbound=%.20f ", orig_hyp_bound);
                // mpfr_out_str (stdout, 10, 0, hyp_bound, MPFR_RNDD);
                // putchar ('\n');
            }

            while(1) {
                // std::cout << "mpfr" << std::endl;
                // u = Random();                              // uniform random number
                mpfr_urandomb(u, rng);

                // double orig_u = mpfr_get_d(u, MPFR_RNDN);
                // printf("u=%.20f ", orig_u);
                // mpfr_out_str (stdout, 10, 0, u, MPFR_RNDD);
                // putchar ('\n');

                // if (u == 0) continue;                      // avoid division by 0
                if (mpfr_zero_p(u) != 0) continue;

                // x = hyp_a + hyp_h * (Random()-0.5) / u;    // generate hat distribution
                mpfr_urandomb(x, rng);
                // double rand = mpfr_get_d(x, MPFR_RNDN);
                mpfr_sub_d(x, x, 0.5, MPFR_RNDN);
                mpfr_div(x, x, u, MPFR_RNDN);
                mpfr_mul(x, x, hyp_h, MPFR_RNDN);
                mpfr_add(x, x, hyp_a, MPFR_RNDN);

                // if (x < 0. || x > 2E9) continue;           // reject, avoid overflow
                if (mpfr_sgn(x) < 0) continue;

                // k = (int64_t)x;
                mpfr_floor(k, x);

                // orig_x = orig_hyp_a + orig_hyp_h * (rand-0.5) / orig_u;    // generate hat distribution
                // orig_k = (int64_t)orig_x;
                // printf("k=%llu ", orig_k);
                // mpfr_out_str (stdout, 10, 0, k, MPFR_RNDD);
                // putchar ('\n');

                // if (k > hyp_bound) continue;               // reject if outside range
                // if (orig_k > orig_hyp_bound) printf("cont "); 
                if (mpfr_cmp(k, hyp_bound) > 0) {
                    // printf("cont out\n");
                    continue;
                }

                // lf = hyp_fm - fc_lnpk(k,L,m,n);            // ln(f(k))
                fc_lnpk(k, L, m, n, lnfac);
                mpfr_sub(lf, hyp_fm, lnfac, MPFR_RNDN);

                // orig_lf = orig_hyp_fm - fc_lnpk(orig_k,orig_L,orig_m,orig_n);            // ln(f(k))
                // printf("lf=%.20f ", orig_lf);
                // mpfr_out_str (stdout, 10, 0, lf, MPFR_RNDD);
                // putchar ('\n');

                // if (u * (4.0 - u) - 3.0 <= lf) break;      // lower squeeze accept
                mpfr_t lower;
                mpfr_init2(lower, 200);
                mpfr_d_sub(lower, 4.0, u, MPFR_RNDN);
                mpfr_mul(lower, lower, u, MPFR_RNDN);
                mpfr_sub_d(lower, lower, 3.0, MPFR_RNDN);

                // if (orig_u * (4.0 - orig_u) - 3.0 <= orig_lf) printf("break ");
                if (mpfr_cmp(lower, lf) < 0 || mpfr_cmp(lower, lf) == 0) {
                    // printf("break\n");
                    break;
                }

                // if (u * (u-lf) > 1.0) continue;            // upper squeeze reject
                mpfr_t upper;
                mpfr_init2(upper, 200);
                mpfr_sub(upper, u, lf, MPFR_RNDN);
                mpfr_mul(upper, upper, u, MPFR_RNDN);

                // if (orig_u * (orig_u-orig_lf) > 1.0) printf("cont ");
                if (mpfr_cmp_d(upper, 1.0) > 0) {
                    // printf("cont up\n");
                    continue;
                }

                // if (2.0 * log(u) <= lf) break;             // final acceptance
                mpfr_t accept;
                mpfr_init2(accept, 200);
                mpfr_log(accept, u, MPFR_RNDN);
                mpfr_mul_d(accept, accept, 2.0, MPFR_RNDN);
                if (mpfr_cmp(accept, lf) < 0 || mpfr_cmp(accept, lf) == 0) break;
            }

            // return k;
            mpfr_set(result, k, MPFR_RNDN);
        }
};

#endif _GEN_HYPER_H_
