/******************************************************************************
 * rhyper.h
 *
 * Source of the sampling routine
 ******************************************************************************
 * Copyright (C) 2016 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/


#ifndef _RHYPER_H_
#define _RHYPER_H_

#include <algorithm>
#include <iostream>
#include <random>

#include "definitions.h"
#include "mt_wrapper.h"

#ifdef HAVE_NEARYINT
#define R_forceint(x)   nearbyint()
#else
#define R_forceint(x)   round(x)
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi)) */
#endif

// Re-implementation of R "rhyper"
class HyperGen {
    public:
        HyperGen(ULONG seed) : mt(seed) { }
        virtual ~HyperGen() { }

        void RandomInit(ULONG seed) {
            mt.RandomInit(seed);
        }

        __float128 unif_rand() {
            return mt.Random();
        }

        __float128 afc(__int128 i) {
            __float128 di = i;
            return (di + 0.5) * log(di) - di + M_LN_SQRT_2PI + (1./12. - 1./360. / (di * di)) / di;
        }

        __float128 generate(__float128 nn1in, __float128 nn2in, __float128 kkin) {
             /* extern __float128 afc(__int128); */

            __int128 nn1, nn2, kk;
            __int128 ix; // return value (coerced to __float128 at the very end)
            bool setup1, setup2;

            /* These should become 'thread_local globals' : */
            static __int128 ks = -1, n1s = -1, n2s = -1;
            static __int128 m, minjx, maxjx;
            static __int128 k, n1, n2; // <- not allowing larger integer par
            static __float128 tn;

            // II :
            static __float128 w;
            // III:
            static __float128 a, d, s, xl, xr, kl, kr, lamdl, lamdr, p1, p2, p3;

            nn1in = R_forceint(nn1in);
            nn2in = R_forceint(nn2in);
            kkin  = R_forceint(kkin);

            nn1 = (__int128)nn1in;
            nn2 = (__int128)nn2in;
            kk  = (__int128)kkin;

            /* if new parameter values, initialize */
            if (nn1 != n1s || nn2 != n2s) {
                setup1 = true;	setup2 = true;
            } else if (kk != ks) {
                setup1 = false;	setup2 = true;
            } else {
                setup1 = false;	setup2 = false;
            }
            if (setup1) {
                n1s = nn1;
                n2s = nn2;
                tn = nn1 + nn2;
                if (nn1 <= nn2) {
                    n1 = nn1;
                    n2 = nn2;
                } else {
                    n1 = nn2;
                    n2 = nn1;
                }
            }
            if (setup2) {
                ks = kk;
                if (kk + kk >= tn) {
                    k = (__int128)(tn - kk);
                } else {
                    k = kk;
                }
            }
            if (setup1 || setup2) {
                m = (__int128) ((k + 1.) * (n1 + 1.) / (tn + 2.));
                minjx = std::max((__int128)0, k - n2);
                maxjx = std::min(n1, k);
#ifdef DEBUG
                printf("rhyper(nn1=%lld, nn2=%lld, kk=%lld), setup: floor(mean)= m=%lld, jx in (%lld..%lld)\n",
                     nn1, nn2, kk, m, minjx, maxjx);
#endif
            }
            /* generate random variate --- Three basic cases */

            if (minjx == maxjx) { /* I: degenerate distribution ---------------- */
#ifdef DEBUG
                printf("rhyper(), branch I (degenerate)\n");
#endif
                ix = maxjx;
                goto L_finis; // return appropriate variate

            } else if (m - minjx < 10) { // II: (Scaled) algorithm HIN (inverse transformation) ----
                const static __float128 scale = 1e25; // scaling factor against (early) underflow
                const static __float128 con = 57.5646273248511421;
                                  // 25*log(10) = log(scale) { <==> exp(con) == scale }
                if (setup1 || setup2) {
                    __float128 lw; // log(w);  w = exp(lw) * scale = exp(lw + log(scale)) = exp(lw + con)
                    if (k < n2) {
                        lw = afc(n2) + afc(n1 + n2 - k) - afc(n2 - k) - afc(n1 + n2);
                    } else {
                        lw = afc(n1) + afc(     k     ) - afc(k - n2) - afc(n1 + n2);
                    }
                    w = exp(lw + con);
                }
                __float128 p, u;
#ifdef DEBUG
                printf("rhyper(), branch II; w = %g > 0\n", w);
#endif
                  L10:
                p = w;
                ix = minjx;
                u = unif_rand() * scale;
#ifdef DEBUG
                printf("  _new_ u = %g\n", u);
#endif
                while (u > p) {
                    u -= p;
                    p *= ((__float128) n1 - ix) * (k - ix);
                    ix++;
                    p = p / ix / (n2 - k + ix);
#ifdef DEBUG
                    printf("       ix=%3d, u=%11g, p=%20.14g (u-p=%g)\n", ix, u, p, u-p);
#endif
                    if (ix > maxjx)
                        goto L10;
                    // FIXME  if(p == 0.)  we also "have lost"  => goto L10
                }
            } else { /* III : H2PE Algorithm --------------------------------------- */

            __float128 u,v;

            if (setup1 || setup2) {
                s = sqrt((tn - k) * k * n1 * n2 / (tn - 1) / tn / tn);
                /* remark: d is defined in reference without __int128. */
                /* the truncation centers the cell boundaries at 0.5 */

                d = (__int128) (1.5 * s) + .5;
                xl = m - d + .5;
                xr = m + d + .5;
                a = afc(m) + afc(n1 - m) + afc(k - m) + afc(n2 - k + m);
                kl = exp(a - afc((__int128) (xl)) - afc((__int128) (n1 - xl))
                     - afc((__int128) (k - xl))
                     - afc((__int128) (n2 - k + xl)));
                kr = exp(a - afc((__int128) (xr - 1))
                     - afc((__int128) (n1 - xr + 1))
                     - afc((__int128) (k - xr + 1))
                     - afc((__int128) (n2 - k + xr - 1)));
                lamdl = -log(xl * (n2 - k + xl) / (n1 - xl + 1) / (k - xl + 1));
                lamdr = -log((n1 - xr + 1) * (k - xr + 1) / xr / (n2 - k + xr));
                p1 = d + d;
                p2 = p1 + kl / lamdl;
                p3 = p2 + kr / lamdr;
            }
#ifdef DEBUG
            printf("rhyper(), branch III {accept/reject}: (xl,xr)= (%g,%g); (lamdl,lamdr)= (%g,%g)\n",
                 xl, xr, lamdl,lamdr);
            printf("-------- p123= c(%g,%g,%g)\n", p1,p2, p3);
#endif
            __int128 n_uv = 0;
              L30:
            u = unif_rand() * p3;
            v = unif_rand();
            n_uv++;
            if(n_uv >= 10000) {
                printf("rhyper() branch III: giving up after %lld rejections\n", n_uv);
                exit(1);
            }
#ifdef DEBUG
            printf(" ... L30: new (u=%g, v ~ U[0,1])[%lld]\n", u, n_uv);
            // std::cin.get();
#endif

            if (u < p1) {		/* rectangular region */
                ix = (__int128) (xl + u);
            } else if (u <= p2) {	/* left tail */
                ix = (__int128) (xl + log(v) / lamdl);
                if (ix < minjx)
                    goto L30;
                v = v * (u - p1) * lamdl;
            } else {		/* right tail */
                ix = (__int128) (xr - log(v) / lamdr);
                if (ix > maxjx)
                    goto L30;
                v = v * (u - p2) * lamdr;
            }

            /* acceptance/rejection test */
            bool reject = true;

            if (m < 100 || ix <= 50) {
                /* explicit evaluation */
                /* The original algorithm (and TOMS 668) have
                   f = f * i * (n2 - k + i) / (n1 - i) / (k - i);
                   in the (m > ix) case, but the definition of the
                   recurrence relation on p134 shows that the +1 is
                   needed. */
                __int128 i;
                __float128 f = 1.0;
                if (m < ix) {
                    for (i = m + 1; i <= ix; i++)
                        f = f * (n1 - i + 1) * (k - i + 1) / (n2 - k + i) / i;
                } else if (m > ix) {
                    for (i = ix + 1; i <= m; i++)
                        f = f * i * (n2 - k + i) / (n1 - i + 1) / (k - i + 1);
                }
                if (v <= f) {
                    reject = false;
                }
            } else {

                const static __float128 deltal = 0.0078;
                const static __float128 deltau = 0.0034;

                __float128 e, g, r, t, y;
                __float128 de, dg, dr, ds, dt, gl, gu, nk, nm, ub;
                __float128 xk, xm, xn, y1, ym, yn, yk, alv;

#ifdef DEBUG
                printf(" ... accept/reject 'large' case v=%g\n", v);
#endif
                /* squeeze using upper and lower bounds */
                y = ix;
                y1 = y + 1.0;
                ym = y - m;
                yn = n1 - y + 1.0;
                yk = k - y + 1.0;
                nk = n2 - k + y1;
                r = -ym / y1;
                s = ym / yn;
                t = ym / yk;
                e = -ym / nk;
                g = yn * yk / (y1 * nk) - 1.0;
                dg = 1.0;
                if (g < 0.0)
                dg = 1.0 + g;
                gu = g * (1.0 + g * (-0.5 + g / 3.0));
                gl = gu - .25 * (g * g * g * g) / dg;
                xm = m + 0.5;
                xn = n1 - m + 0.5;
                xk = k - m + 0.5;
                nm = n2 - k + xm;
                ub = y * gu - m * gl + deltau
                + xm * r * (1. + r * (-0.5 + r / 3.0))
                + xn * s * (1. + s * (-0.5 + s / 3.0))
                + xk * t * (1. + t * (-0.5 + t / 3.0))
                + nm * e * (1. + e * (-0.5 + e / 3.0));
                /* test against upper bound */
                alv = log(v);
                if (alv > ub) {
                    reject = true;
                } else {
                /* test against lower bound */
                dr = xm * (r * r * r * r);
                if (r < 0.0)
                    dr /= (1.0 + r);
                ds = xn * (s * s * s * s);
                if (s < 0.0)
                    ds /= (1.0 + s);
                dt = xk * (t * t * t * t);
                if (t < 0.0)
                    dt /= (1.0 + t);
                de = nm * (e * e * e * e);
                if (e < 0.0)
                    de /= (1.0 + e);
                if (alv < ub - 0.25 * (dr + ds + dt + de)
                    + (y + m) * (gl - gu) - deltal) {
                    reject = false;
                }
                else {
                    /* * Stirling's formula to machine accuracy
                     */
                    if (alv <= (a - afc(ix) - afc(n1 - ix)
                        - afc(k - ix) - afc(n2 - k + ix))) {
                    reject = false;
                    } else {
                    reject = true;
                    }
                }
                }
            } // else
            if (reject)
                goto L30;
            }


        L_finis:
            /* return appropriate variate */

            if (kk + kk >= tn) {
                if (nn1 > nn2) {
                    ix = kk - nn2 + ix;
                } else {
                    ix = nn1 - ix;
                }
            } else {
                if (nn1 > nn2)
                    ix = kk - ix;
            }
            return ix;
        }

    private:
        MTWrapper mt;
};

#endif 

