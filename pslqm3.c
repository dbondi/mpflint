/*
 * pslqm3.c — Production PSLQ integer relation detection engine
 *
 * Three-level multipair PSLQ with optional intermediate precision (IP) layer.
 * C implementation using FLINT/arb for multi-precision arithmetic.
 * Based on Bailey's MPFUN2020 Fortran implementation.
 *
 * Precision levels:
 *   DP   (53 bits)     — fast inner loop
 *   IP   (IP_BITS)     — optional intermediate layer (absorbs DP flushes cheaply)
 *   MPM  (ndpm digits) — medium precision updates
 *   Full (digits)      — high-precision sync
 *
 * Strategies:
 *   0 = "standard"       — no nudge (Bailey's original pair selection)
 *   1 = "predicted_swap" — Givens-aware bidirectional insertion sort
 *
 * Environment variables:
 *   STRATEGY     — "standard" or "predicted_swap" (default: predicted_swap)
 *   THREADS      — FLINT thread count (default: 1)
 *   IP_BITS      — intermediate precision in bits (0=disabled, default: 0)
 *   IPM_OVERRIDE — DP iterations per MPM batch (default: 10)
 *
 * Build:
 *   cc -O2 -march=native -shared -fPIC -o pslqm3.dylib pslqm3.c \
 *      -I/opt/homebrew/include/flint -I/opt/homebrew/include \
 *      -L/opt/homebrew/lib -lflint -lgmp -lmpfr -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <pthread.h>
#include "flint/flint.h"
#include "flint/arb.h"
#include "flint/arb_mat.h"
#include "flint/arf.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

#include "pslq_arb.h"

/* ================================================================
 * Global state
 * ================================================================ */

static int g_strategy = 1;        /* 0=standard, 1=predicted_swap */
static int g_ipm_override = 0;
static long g_predicted_swap_cnt = 0;
static long g_pairsel_swap_cnt = 0;

int g_colpar_threads = 1;

#include "pslq_sort.c"

/* ================================================================
 * Profiling counters
 * ================================================================ */

/* Profiling counters */
double g_mxmdm_sec = 0.0;
int    g_mxmdm_calls = 0;
double g_mxm_sec = 0.0;
int    g_mxm_calls = 0;
static double g_iterdp_sec = 0.0;
static int    g_iterdp_calls = 0;
static double g_itermpm_sec = 0.0;
static int    g_itermpm_calls = 0;
static double g_lqmpm_sec = 0.0;
static int    g_lqmpm_calls = 0;
static double g_updtmpm_sec = 0.0;
static int    g_updtmpm_calls = 0;
static double g_updtmp_sec = 0.0;
static int    g_updtmp_calls = 0;
static double g_updtmp_wby_sec = 0.0;
static double g_updtmp_wbb_sec = 0.0;
static double g_updtmp_wah_sec = 0.0;
static int    g_cnt_izd1_deps = 0, g_cnt_izd1_tmx1 = 0;
static int    g_last_fmp_it = 0;
static int    g_izmm_wy_small = 0;
static int    g_izmm_wy_ratio = 0;
static int    g_izmm_wab_large = 0;

static void mxmdm_fmpz(int n1, int n2, slong mpm_prec,
                        const fmpz_mat_t fA, arb_ptr b);
static void arb_to_fmpz_vec(fmpz_mat_t fB, slong *b_min_exp_out,
                             arb_ptr data, int rows, int cols);
static void fmpz_result_to_arb(arb_ptr b, const fmpz_mat_t fC,
                                slong b_min_exp, int n1, int n2, slong prec);

/* ================================================================
 * initmp — initialize H, B, y from input x
 * ================================================================ */

static void initmp(int idb, int n, slong full_prec,
                   arb_ptr b, arb_ptr h, arb_ptr x, arb_ptr y)
{
    int n1 = n - 1;
    arb_ptr s = _arb_vec_init(n);
    arb_t t1;
    arb_init(t1);

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++)
            arb_zero(AM(b, i, j, n));
        arb_one(AM(b, j, j, n));
    }

    arb_zero(t1);
    for (int i = n - 1; i >= 0; i--) {
        arb_t sq;
        arb_init(sq);
        mp_mul(sq, x+(i), x+(i), full_prec);
        mp_add(t1, t1, sq, full_prec);
        mp_sqrt(s+(i), t1, full_prec);
        arb_clear(sq);
    }

    mp_inv(t1, s+(0), full_prec);

    for (int i = 0; i < n; i++) {
        mp_mul(y+(i), t1, x+(i), full_prec);
        mp_mul(s+(i), t1, s+(i), full_prec);
    }

    for (int j = 0; j < n1; j++) {
        for (int i = 0; i < j; i++)
            arb_zero(AM(h, i, j, n1));

        mp_div(AM(h, j, j, n1), s+(j+1), s+(j), full_prec);

        arb_t denom;
        arb_init(denom);
        mp_mul(denom, s+(j), s+(j+1), full_prec);
        mp_div(t1, y+(j), denom, full_prec);
        arb_clear(denom);

        for (int i = j + 1; i < n; i++) {
            mp_mul(AM(h, i, j, n1), y+(i), t1, full_prec);
            mp_neg(AM(h, i, j, n1), AM(h, i, j, n1));
        }
    }

    arb_clear(t1);
    _arb_vec_clear(s, n);
}

/* ================================================================
 * initmpm — initialize wa, wb, wh, wy from H, y
 * ================================================================ */

static void initmpm(int idb, int n, int nsq, slong mpm_prec,
                    arb_ptr wa, arb_ptr wb, arb_ptr wh, arb_ptr wy,
                    arb_ptr wsyq, arb_srcptr h, arb_srcptr y)
{
    int n1 = n - 1;

    arf_t af1, afmax, afone, afscale;
    arf_init(af1); arf_init(afmax); arf_init(afone); arf_init(afscale);
    arf_one(afone);
    arf_zero(afmax);
    for (int i = 0; i < n; i++) {
        arf_abs(af1, arb_midref(y+(i)));
        if (arf_cmp(af1, afmax) > 0)
            arf_set(afmax, af1);
    }

    arf_div(afscale, afone, afmax, mpm_prec, ARF_RND_NEAR);

    arb_t t1;
    arb_init(t1);
    arb_zero(t1);
    arf_set(arb_midref(t1), afscale);
    for (int i = 0; i < n; i++) {
        mp_set_round(wy+(i), y+(i), mpm_prec);
        mp_mul(wy+(i), t1, wy+(i), mpm_prec);
    }
    arf_clear(af1); arf_clear(afmax); arf_clear(afone); arf_clear(afscale);

    for (int j = 0; j < n1; j++)
        for (int i = 0; i < n; i++)
            mp_set_round(AM(wh, i, j, n1), h+(i * n1 + j), mpm_prec);

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            arb_zero(AM(wa, i, j, n));
            arb_zero(AM(wb, i, j, n));
        }
        arb_one(AM(wa, j, j, n));
        arb_one(AM(wb, j, j, n));
    }

    for (int j = 0; j < nsq; j++)
        for (int i = 0; i < n; i++)
            arb_zero(AM(wsyq, i, j, nsq));

    arb_clear(t1);
}

/* ================================================================
 * initdp — initialize DP arrays from MPM arrays
 * ================================================================ */

static void initdp(int idb, int n, int nsq, slong mpm_prec,
                   double *da, double *db, double *dh, double *dy,
                   double *dsyq, arb_srcptr wh, arb_srcptr wy)
{
    int n1 = n - 1;
    arf_t af1, afmax, afone, afscale;
    arf_init(af1); arf_init(afmax); arf_init(afone); arf_init(afscale);
    arf_one(afone);

    arf_zero(afmax);
    for (int i = 0; i < n; i++) {
        arf_abs(af1, arb_midref(wy+(i)));
        if (arf_cmp(af1, afmax) > 0)
            arf_set(afmax, af1);
    }

    arf_div(afscale, afone, afmax, mpm_prec, ARF_RND_NEAR);

    for (int i = 0; i < n; i++) {
        arf_mul(af1, afscale, arb_midref(wy+(i)), mpm_prec, ARF_RND_NEAR);
        dy[i] = arf_get_d(af1, ARF_RND_NEAR);
    }

    arf_zero(afmax);
    for (int j = 0; j < n1; j++) {
        arf_abs(af1, arb_midref(AM(wh, j, j, n1)));
        if (arf_cmp(af1, afmax) > 0)
            arf_set(afmax, af1);
    }

    arf_div(afscale, afone, afmax, mpm_prec, ARF_RND_NEAR);

    for (int j = 0; j < n1; j++)
        for (int i = 0; i < n; i++) {
            arf_mul(af1, afscale, arb_midref(AM(wh, i, j, n1)), mpm_prec, ARF_RND_NEAR);
            DM(dh, i, j, n) = arf_get_d(af1, ARF_RND_NEAR);
        }

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            DM(da, i, j, n) = 0.0;
            DM(db, i, j, n) = 0.0;
        }
        DM(da, j, j, n) = 1.0;
        DM(db, j, j, n) = 1.0;
    }

    for (int j = 0; j < nsq; j++)
        for (int i = 0; i < n; i++)
            DM(dsyq, i, j, nsq) = 0.0;

    arf_clear(af1); arf_clear(afmax); arf_clear(afone); arf_clear(afscale);
}

/* ================================================================
 * savedp — save/restore DP arrays
 * ================================================================ */

static void savedp(int n, const double *da, const double *db,
                   const double *dh, const double *dy,
                   double *dsa, double *dsb, double *dsh, double *dsy)
{
    for (int i = 0; i < n; i++)
        dsy[i] = dy[i];

    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            DM(dsa, i, j, n) = DM(da, i, j, n);
            DM(dsb, i, j, n) = DM(db, i, j, n);
            DM(dsh, i, j, n) = DM(dh, i, j, n);
        }
}

/* ================================================================
 * lqdp — LQ decomposition (DP)
 * ================================================================ */

static void lqdp(int n, int m, int stride, double *dh)
{
    int lup = (m < n) ? m : n;

    for (int l = 0; l < lup; l++) {
        if (l == m - 1) continue;

        int ml = m - 1 - l;
        double t = 0.0;

        for (int i = 0; i <= ml; i++)
            t += DM(dh, l, l + i, stride) * DM(dh, l, l + i, stride);

        double nrmxl = sqrt(t);
        if (nrmxl == 0.0) continue;

        if (DM(dh, l, l, stride) != 0.0)
            nrmxl = copysign(nrmxl, DM(dh, l, l, stride));
        t = 1.0 / nrmxl;

        for (int i = 0; i <= ml; i++)
            DM(dh, l, l + i, stride) = t * DM(dh, l, l + i, stride);

        DM(dh, l, l, stride) = 1.0 + DM(dh, l, l, stride);

        for (int j = l + 1; j < n; j++) {
            t = 0.0;
            for (int i = 0; i <= ml; i++)
                t += DM(dh, l, l + i, stride) * DM(dh, j, l + i, stride);
            t = -t / DM(dh, l, l, stride);
            for (int i = 0; i <= ml; i++)
                DM(dh, j, l + i, stride) = DM(dh, j, l + i, stride) + t * DM(dh, l, l + i, stride);
        }

        DM(dh, l, l, stride) = -nrmxl;
    }

    for (int j = 0; j < m; j++)
        for (int i = 0; i < j; i++)
            DM(dh, i, j, stride) = 0.0;
}

/* ================================================================
 * lqmpm — LQ decomposition (MPM)
 * ================================================================ */

static void lqmpm(int n, int m, slong mpm_prec, arb_ptr h)
{
    struct timespec _t0, _t1;
    clock_gettime(CLOCK_MONOTONIC, &_t0);
    int lup = (m < n) ? m : n;

    arb_t t, nrmxl;
    arb_init(t); arb_init(nrmxl);

    for (int l = 0; l < lup; l++) {
        if (l == m - 1) continue;
        int ml = m - 1 - l;

        arb_zero(t);
        for (int i = 0; i <= ml; i++)
            mp_addmul(t, AM(h, l, l+i, m), AM(h, l, l+i, m), mpm_prec);

        mp_sqrt(nrmxl, t, mpm_prec);
        if (arf_is_zero(arb_midref(nrmxl))) continue;

        if (!arf_is_zero(arb_midref(AM(h, l, l, m)))) {
            if (arf_sgn(arb_midref(AM(h, l, l, m))) < 0)
                mp_neg(nrmxl, nrmxl);
        }

        mp_inv(t, nrmxl, mpm_prec);
        for (int i = 0; i <= ml; i++)
            mp_mul(AM(h, l, l+i, m), t, AM(h, l, l+i, m), mpm_prec);

        arb_t one; arb_init(one); arb_one(one);
        mp_add(AM(h, l, l, m), AM(h, l, l, m), one, mpm_prec);
        arb_clear(one);

        for (int j = l + 1; j < n; j++) {
            arb_zero(t);
            for (int i = 0; i <= ml; i++)
                mp_addmul(t, AM(h, l, l+i, m), AM(h, j, l+i, m), mpm_prec);
            mp_div(t, t, AM(h, l, l, m), mpm_prec);
            mp_neg(t, t);
            for (int i = 0; i <= ml; i++)
                mp_addmul(AM(h, j, l+i, m), t, AM(h, l, l+i, m), mpm_prec);
        }

        mp_neg(AM(h, l, l, m), nrmxl);
    }

    for (int j = 0; j < m; j++)
        for (int i = 0; i < j; i++)
            arb_zero(AM(h, i, j, m));

    arb_clear(t); arb_clear(nrmxl);
    clock_gettime(CLOCK_MONOTONIC, &_t1);
    g_lqmpm_sec += (double)(_t1.tv_sec - _t0.tv_sec) + (double)(_t1.tv_nsec - _t0.tv_nsec) * 1e-9;
    g_lqmpm_calls++;
}

/* ================================================================
 * iterdp — one DP iteration
 * ================================================================ */

static void iterdp(int idb, int it, int n, int nsq,
                   double *da, double *db, double *dh, double *dsyq,
                   double *dy, int *imq, int *izd)
{
    struct timespec _tdp0, _tdp1;
    clock_gettime(CLOCK_MONOTONIC, &_tdp0);
    int n1 = n - 1;
    double deps = 1.0e-14;
    double tmx1 = 1.0e13;
    static double tmx2 = 0.0;
    static double gam = 0.0;
    static double *gam_pow = NULL;
    static int gam_pow_n = 0;
    if (gam == 0.0) {
        tmx2 = pow(2.0, 52.0);
        gam = sqrt(4.0 / 3.0);
    }
    if (gam_pow_n != n1) {
        free(gam_pow);
        gam_pow = (double *)malloc(n1 * sizeof(double));
        for (int i = 0; i < n1; i++)
            gam_pow[i] = pow(gam, (double)(i + 1));
        gam_pow_n = n1;
    }

    static int *ip = NULL, *ir = NULL, *is_buf = NULL;
    static double *dq = NULL, *dt = NULL;
    static int alloc_n = 0;
    if (alloc_n != n) {
        free(ip); free(ir); free(is_buf); free(dq); free(dt);
        ip = (int *)malloc(n * sizeof(int));
        ir = (int *)malloc(n * sizeof(int));
        is_buf = (int *)malloc(n * sizeof(int));
        dq = (double *)malloc(n * sizeof(double));
        dt = (double *)malloc((size_t)n * n * sizeof(double));
        alloc_n = n;
    }
    memset(is_buf, 0, n * sizeof(int));
    memset(dt, 0, (size_t)n * n * sizeof(double));

    *izd = 0;

    int mpr = (int)(0.4 * n + 0.5);

    for (int i = 0; i < n1; i++)
        dq[i] = gam_pow[i] * fabs(DM(dh, i, i, n));

    qsortdp(n1, dq, ip);

    int mq;
    if (*imq == 0) {
        mq = mpr;
    } else {
        mq = 1;
        *imq = 0;
    }

    int actual_mq = 0;
    {
        int ii = n1;
        for (int i = 0; i < mq; i++) {
        full_reinit_dp:
            ii = ii - 1;
            if (ii < 0) {
                actual_mq = i;
                goto check_wy_dp;
            }
            int j1 = ip[ii];
            int j2 = j1 + 1;
            if (is_buf[j1] != 0 || is_buf[j2] != 0) goto full_reinit_dp;
            ir[i] = j1;
            is_buf[j1] = 1;
            is_buf[j2] = 1;
            actual_mq = i + 1;
        }
    }
check_wy_dp:
    mq = actual_mq;

    for (int j = 0; j < mq; j++) {
        int im = ir[j];
        int im1 = im + 1;
        double t1;

        g_pairsel_swap_cnt++;

        t1 = dy[im]; dy[im] = dy[im1]; dy[im1] = t1;

        for (int i = 0; i < n; i++) {
            t1 = DM(da, im, i, n); DM(da, im, i, n) = DM(da, im1, i, n); DM(da, im1, i, n) = t1;
            t1 = DM(db, im, i, n); DM(db, im, i, n) = DM(db, im1, i, n); DM(db, im1, i, n) = t1;
        }

        for (int i = 0; i < n1; i++) {
            t1 = DM(dh, im, i, n); DM(dh, im, i, n) = DM(dh, im1, i, n); DM(dh, im1, i, n) = t1;
        }
    }

    for (int j = 0; j < mq; j++) {
        int im = ir[j];
        int im1 = im + 1;
        if (im <= n - 3) {
            double t1 = DM(dh, im, im, n);
            double t2 = DM(dh, im, im1, n);
            double t3 = sqrt(t1 * t1 + t2 * t2);
            if (t3 == 0.0) continue;
            t1 = t1 / t3;
            t2 = t2 / t3;

            for (int i = im; i < n; i++) {
                double t3v = DM(dh, i, im, n);
                double t4 = DM(dh, i, im1, n);
                DM(dh, i, im, n) = t1 * t3v + t2 * t4;
                DM(dh, i, im1, n) = -t2 * t3v + t1 * t4;
            }
        }
    }

    /* Full Hermite reduction */
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n - i; j++) {
            int ij = i + j;
            for (int k = j + 1; k < ij; k++)
                DM(dh, ij, j, n) -= DM(dt, ij, k, n) * DM(dh, k, j, n);
            DM(dt, ij, j, n) = round(DM(dh, ij, j, n) / DM(dh, j, j, n));
            DM(dh, ij, j, n) -= DM(dt, ij, j, n) * DM(dh, j, j, n);
        }
    }

    double t1 = fabs(dy[n - 1]);
    for (int j = 0; j < n1; j++) {
        for (int i = j + 1; i < n; i++)
            dy[j] += DM(dt, i, j, n) * dy[i];
        t1 = fmin(t1, fabs(dy[j]));
    }

    for (int j = 0; j < n1; j++) {
        for (int i = j + 1; i < n; i++) {
            double dtij = DM(dt, i, j, n);
            if (dtij == 0.0) continue;
            double *da_i = da + (size_t)i * n;
            const double *da_j = da + (size_t)j * n;
            for (int k = 0; k < n; k++)
                da_i[k] -= dtij * da_j[k];
        }
    }
    for (int j = 0; j < n1; j++) {
        for (int i = j + 1; i < n; i++) {
            double dtij = DM(dt, i, j, n);
            if (dtij == 0.0) continue;
            double *db_j = db + (size_t)j * n;
            const double *db_i = db + (size_t)i * n;
            for (int k = 0; k < n; k++)
                db_j[k] += dtij * db_i[k];
        }
    }

    double t2 = 0.0;
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++) {
            double v = fabs(DM(da, i, k, n));
            if (v > t2) t2 = v;
            v = fabs(DM(db, i, k, n));
            if (v > t2) t2 = v;
        }

    if (t1 <= deps) {
        *izd = 1;
        g_cnt_izd1_deps++;
    }

    if (t2 > tmx1 && t2 <= tmx2) {
        *izd = 1;
        g_cnt_izd1_tmx1++;
    } else if (t2 > tmx2) {
        *izd = 2;
        goto cleanup_iterdp;
    }

    if (isnan(t1) || isnan(t2)) {
        *izd = 2;
        goto cleanup_iterdp;
    }

    for (int j = 0; j < nsq; j++) {
        double t1v = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = fabs(dy[i] - DM(dsyq, i, j, nsq));
            if (diff > t1v) t1v = diff;
        }
        if (t1v <= deps) {
            *imq = 1;
            break;
        }
    }

    {
        int k = it % nsq;
        for (int i = 0; i < n; i++)
            DM(dsyq, i, k, nsq) = dy[i];
    }

cleanup_iterdp:
    clock_gettime(CLOCK_MONOTONIC, &_tdp1);
    g_iterdp_sec += (double)(_tdp1.tv_sec - _tdp0.tv_sec) + (double)(_tdp1.tv_nsec - _tdp0.tv_nsec) * 1e-9;
    g_iterdp_calls++;
}

/* ================================================================
 * itermpm — one MPM iteration
 * ================================================================ */

static void itermpm(int idb, int it, int n, int nsq, slong mpm_prec,
                    const arb_t epsm, arb_ptr wa, arb_ptr wb, arb_ptr wh,
                    arb_ptr wsyq, arb_ptr wy, int *imq, int *izmm)
{
    struct timespec _tmpm0, _tmpm1;
    clock_gettime(CLOCK_MONOTONIC, &_tmpm0);
    int n1 = n - 1;

    arb_t tmx1, tmx2, gam_arb, t1, t2, t3, t4;
    arb_init(tmx1); arb_init(tmx2); arb_init(gam_arb);
    arb_init(t1); arb_init(t2); arb_init(t3); arb_init(t4);

    mp_inv(tmx1, epsm, mpm_prec);
    arb_t tmp_1e20; arb_init(tmp_1e20);
    mp_set_d(tmp_1e20, 1e20);
    mp_div(tmx1, tmx1, tmp_1e20, mpm_prec);
    arb_clear(tmp_1e20);

    arb_one(tmx2);
    mp_mul_2exp_si(tmx2, tmx2, mpm_prec);

    *izmm = 0;
    int mpr = (int)(0.4 * n + 0.5);

    mp_set_d(t1, 4.0);
    mp_set_d(t2, 3.0);
    mp_div(gam_arb, t1, t2, mpm_prec);
    mp_sqrt(gam_arb, gam_arb, mpm_prec);

    arb_ptr wq = _arb_vec_init(n);
    int *ip = (int *)malloc(n * sizeof(int));
    int *ir = (int *)malloc(n * sizeof(int));
    int *is_arr = (int *)calloc(n, sizeof(int));

    for (int i = 0; i < n1; i++) {
        arb_pow_ui(t1, gam_arb, (ulong)(i + 1), mpm_prec);
        mp_abs(t2, AM(wh, i, i, n1));
        mp_set_round(t2, t2, mpm_prec);
        mp_mul(wq+(i), t1, t2, mpm_prec);
    }

    qsortmpm(n1, wq, ip);

    int mq;
    if (*imq == 0) {
        mq = mpr;
    } else {
        mq = 1;
        *imq = 0;
    }

    int ii = n1;
    int actual_mq = 0;
    for (int i = 0; i < mq; i++) {
    full_reinit_mpm:
        ii = ii - 1;
        if (ii < 0) {
            actual_mq = i;
            goto check_wy_mpm;
        }
        int j1 = ip[ii];
        int j2 = j1 + 1;
        if (is_arr[j1] != 0 || is_arr[j2] != 0) goto full_reinit_mpm;
        ir[i] = j1;
        is_arr[j1] = 1;
        is_arr[j2] = 1;
        actual_mq = i + 1;
    }
check_wy_mpm:
    mq = actual_mq;

    for (int j = 0; j < mq; j++) {
        int im = ir[j];
        int im1 = im + 1;

        arb_swap(wy+(im), wy+(im1));

        for (int i = 0; i < n; i++) {
            arb_swap(AM(wa, im, i, n), AM(wa, im1, i, n));
            arb_swap(AM(wb, im, i, n), AM(wb, im1, i, n));
        }

        for (int i = 0; i < n1; i++)
            arb_swap(AM(wh, im, i, n1), AM(wh, im1, i, n1));
    }

    for (int j = 0; j < mq; j++) {
        int im = ir[j];
        int im1 = im + 1;
        if (im <= n - 3) {
            mp_set(t1, AM(wh, im, im, n1));
            mp_set(t2, AM(wh, im, im1, n1));
            mp_mul(t3, t1, t1, mpm_prec);
            mp_addmul(t3, t2, t2, mpm_prec);
            mp_sqrt(t3, t3, mpm_prec);
            if (arf_is_zero(arb_midref(t3))) continue;
            mp_div(t1, t1, t3, mpm_prec);
            mp_div(t2, t2, t3, mpm_prec);

            for (int i = im; i < n; i++) {
                mp_set(t3, AM(wh, i, im, n1));
                mp_set(t4, AM(wh, i, im1, n1));
                mp_mul(AM(wh, i, im, n1), t1, t3, mpm_prec);
                mp_addmul(AM(wh, i, im, n1), t2, t4, mpm_prec);
                mp_mul(AM(wh, i, im1, n1), t1, t4, mpm_prec);
                mp_submul(AM(wh, i, im1, n1), t2, t3, mpm_prec);
            }
        }
    }

    arb_ptr wt = _arb_vec_init((size_t)n * n);
    arb_t _dot_tmp;
    arb_init(_dot_tmp);

    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n - i; j++) {
            int ij = i + j;
            int klen = ij - j - 1;

            if (klen > 0) {
                arb_approx_dot(_dot_tmp, AM(wh, ij, j, n1), 1,
                               AM(wt, ij, j + 1, n), 1,
                               AM(wh, j + 1, j, n1), n1,
                               klen, mpm_prec);
                arb_set(AM(wh, ij, j, n1), _dot_tmp);
            }

            mp_div(AM(wt, ij, j, n), AM(wh, ij, j, n1), AM(wh, j, j, n1), mpm_prec);
            arb_anint(AM(wt, ij, j, n), AM(wt, ij, j, n));
            mp_submul(AM(wh, ij, j, n1), AM(wt, ij, j, n), AM(wh, j, j, n1), mpm_prec);
        }
    }
    arb_clear(_dot_tmp);

    mp_abs(t1, wy+(n - 1));
    for (int j = 0; j < n1; j++) {
        for (int i = j + 1; i < n; i++)
            mp_addmul(wy+(j), AM(wt, i, j, n), wy+(i), mpm_prec);
        mp_abs(t3, wy+(j));
        arb_min_mid(t1, t1, t3);
    }

    for (int j = 0; j < n1; j++) {
        for (int i = j + 1; i < n; i++) {
            arb_srcptr wtij = AM(wt, i, j, n);
            if (arf_is_zero(arb_midref(wtij))) continue;
            for (int k = 0; k < n; k++) {
                mp_submul(AM(wa, i, k, n), wtij, AM(wa, j, k, n), mpm_prec);
                mp_addmul(AM(wb, j, k, n), wtij, AM(wb, i, k, n), mpm_prec);
            }
        }
    }

    arb_zero(t2);
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++) {
            mp_abs(t3, AM(wa, i, k, n));
            arb_max_mid(t2, t2, t3);
            mp_abs(t3, AM(wb, i, k, n));
            arb_max_mid(t2, t2, t3);
        }

    if (arb_cmp_mid(t1, epsm) <= 0) {
        *izmm = 1;
        g_izmm_wy_small++;
    }

    if (arb_cmp_mid(t2, tmx1) > 0 && arb_cmp_mid(t2, tmx2) <= 0) {
        *izmm = 1;
        g_izmm_wab_large++;
    } else if (arb_cmp_mid(t2, tmx2) > 0) {
        *izmm = 2;
        goto cleanup_itermpm;
    }

    for (int j = 0; j < nsq; j++) {
        arb_zero(t1);
        for (int i = 0; i < n; i++) {
            mp_sub(t3, wy+(i), AM(wsyq, i, j, nsq), mpm_prec);
            mp_abs(t3, t3);
            arb_max_mid(t1, t1, t3);
        }
        if (arb_cmp_mid(t1, epsm) <= 0) {
            *imq = 1;
            break;
        }
    }

    {
        int k = it % nsq;
        for (int i = 0; i < n; i++)
            mp_set(AM(wsyq, i, k, nsq), wy+(i));
    }

cleanup_itermpm:
    _arb_vec_clear(wq, n);
    _arb_vec_clear(wt, (size_t)n * n);
    free(ip); free(ir); free(is_arr);
    arb_clear(tmx1); arb_clear(tmx2); arb_clear(gam_arb);
    arb_clear(t1); arb_clear(t2); arb_clear(t3); arb_clear(t4);
    clock_gettime(CLOCK_MONOTONIC, &_tmpm1);
    g_itermpm_sec += (double)(_tmpm1.tv_sec - _tmpm0.tv_sec) + (double)(_tmpm1.tv_nsec - _tmpm0.tv_nsec) * 1e-9;
    g_itermpm_calls++;
}

#include "pslq_matmul.c"
#include "pslq_4level.c"

/* ================================================================
 * setprec — set working precision of arrays
 * ================================================================ */

static void setprec(slong mpm_prec, int n, int nsq,
                    arb_ptr wa, arb_ptr wb, arb_ptr wh, arb_ptr wsyq, arb_ptr wy)
{
    int n1 = n - 1;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            mp_set_round(AM(wa, i, j, n), AM(wa, i, j, n), mpm_prec);
            mp_set_round(AM(wb, i, j, n), AM(wb, i, j, n), mpm_prec);
            if (j < n1)
                mp_set_round(AM(wh, i, j, n1), AM(wh, i, j, n1), mpm_prec);
        }

    for (int j = 0; j < nsq; j++)
        for (int i = 0; i < n; i++)
            mp_set_round(AM(wsyq, i, j, nsq), AM(wsyq, i, j, nsq), mpm_prec);

    for (int i = 0; i < n; i++)
        mp_set_round(wy+(i), wy+(i), mpm_prec);
}

/* ================================================================
 * dynrange / dynrangem — min/max ratio
 * ================================================================ */

static void dynrange(int n, slong mpm_prec, arb_srcptr y, arb_t result)
{
    arb_t t1, t2, t3;
    arb_init(t1); arb_init(t2); arb_init(t3);
    mp_set_d(t1, 1e300);
    arb_zero(t2);
    for (int i = 0; i < n; i++) {
        mp_set_round(t3, y+(i), mpm_prec);
        mp_abs(t3, t3);
        arb_min_mid(t1, t1, t3);
        arb_max_mid(t2, t2, t3);
    }
    mp_div(result, t1, t2, mpm_prec);
    arb_clear(t1); arb_clear(t2); arb_clear(t3);
}

static void dynrangem(int n, slong mpm_prec, arb_srcptr wy, arb_t result)
{
    arb_t t1, t2, t3;
    arb_init(t1); arb_init(t2); arb_init(t3);
    mp_set_d(t1, 1e300);
    arb_zero(t2);
    for (int i = 0; i < n; i++) {
        mp_set_round(t3, wy+(i), mpm_prec);
        mp_abs(t3, t3);
        arb_min_mid(t1, t1, t3);
        arb_max_mid(t2, t2, t3);
    }
    mp_div(result, t1, t2, mpm_prec);
    arb_clear(t1); arb_clear(t2); arb_clear(t3);
}

/* ================================================================
 * bound_fn — compute norm bound from wh
 * ================================================================ */

static void bound_fn(int n, slong mpm_prec, arb_ptr wh, arb_t result)
{
    int n1 = n - 1;
    lqmpm(n, n1, mpm_prec, wh);

    arb_t t1;
    arb_init(t1);
    arb_zero(t1);
    for (int i = 0; i < n1; i++) {
        arb_t absv;
        arb_init(absv);
        mp_abs(absv, AM(wh, i, i, n1));
        arb_max_mid(t1, t1, absv);
        arb_clear(absv);
    }
    mp_inv(result, t1, mpm_prec);
    arb_clear(t1);
}

/* ================================================================
 * updtmpm — update MPM from DP
 * ================================================================ */

static void updtmpm(int idb, int it, int n, slong mpm_prec,
                    const arb_t epsm, const double *da, const double *db,
                    double dreps, arb_ptr wa, arb_ptr wb, arb_ptr wh,
                    arb_ptr wy, int *izmm)
{
    struct timespec _tup0, _tup1;
    clock_gettime(CLOCK_MONOTONIC, &_tup0);
    int n1 = n - 1;
    arb_t tmx1, tmx2, t1, t2, w1;
    arb_init(tmx1); arb_init(tmx2); arb_init(t1); arb_init(t2); arb_init(w1);

    mp_inv(tmx1, epsm, mpm_prec);

    arb_one(tmx2);
    mp_mul_2exp_si(tmx2, tmx2, mpm_prec);

    *izmm = 0;

    fmpz_mat_t fDA, fDB;
    fmpz_mat_init(fDA, n, n);
    fmpz_mat_init(fDB, n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            fmpz_set_d(fmpz_mat_entry(fDA, i, j), DM(da, i, j, n));
            fmpz_set_d(fmpz_mat_entry(fDB, i, j), DM(db, i, j, n));
        }

    mxmdm_vec(n, mpm_prec, db, wy);

    mp_set_d(t1, 1e300);
    arb_zero(t2);
    for (int i = 0; i < n; i++) {
        arb_t absv;
        arb_init(absv);
        mp_abs(absv, wy+(i));
        arb_min_mid(t1, t1, absv);
        arb_max_mid(t2, t2, absv);
        arb_clear(absv);
    }
    mp_div(w1, t1, t2, mpm_prec);

    if (arb_cmp_mid(t1, epsm) <= 0) {
        *izmm = 1;
        g_izmm_wy_small++;
    }

    mxmdm_fmpz(n, n, mpm_prec, fDA, wa);
    mxmdm_fmpz(n, n, mpm_prec, fDB, wb);

    arb_zero(t2);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            arb_t absv; arb_init(absv);
            mp_abs(absv, AM(wa, i, j, n));
            arb_max_mid(t2, t2, absv);
            mp_abs(absv, AM(wb, i, j, n));
            arb_max_mid(t2, t2, absv);
            arb_clear(absv);
        }
    if (arb_cmp_mid(t2, tmx1) > 0 && arb_cmp_mid(t2, tmx2) <= 0) {
        *izmm = 1;
    } else if (arb_cmp_mid(t2, tmx2) > 0) {
        *izmm = 2;
    }

    {
        arb_t dreps_arb;
        arb_init(dreps_arb);
        mp_set_d(dreps_arb, dreps);
        if (arb_cmp_mid(w1, dreps_arb) < 0) {
            *izmm = 1;
            g_izmm_wy_ratio++;
        }
        arb_clear(dreps_arb);
    }

    mxmdm_fmpz(n, n1, mpm_prec, fDA, wh);

    fmpz_mat_clear(fDA);
    fmpz_mat_clear(fDB);
    arb_clear(tmx1); arb_clear(tmx2); arb_clear(t1); arb_clear(t2); arb_clear(w1);
    clock_gettime(CLOCK_MONOTONIC, &_tup1);
    g_updtmpm_sec += (double)(_tup1.tv_sec - _tup0.tv_sec) + (double)(_tup1.tv_nsec - _tup0.tv_nsec) * 1e-9;
    g_updtmpm_calls++;
}

/* ================================================================
 * updtmp — update MP from MPM
 * ================================================================ */

static void updtmp(int idb, int it, int n, slong full_prec, slong mpm_prec,
                   arb_srcptr wa, arb_srcptr wb, const arb_t eps,
                   arb_ptr b, arb_ptr h, arb_ptr y, int *izm)
{
    struct timespec _tum0, _tum1;
    clock_gettime(CLOCK_MONOTONIC, &_tum0);
    int n1 = n - 1;
    arb_t t1, t2;
    arb_init(t1); arb_init(t2);

    *izm = 0;

    { struct timespec _ts0, _ts1; clock_gettime(CLOCK_MONOTONIC, &_ts0);
    mxm(n, 1, full_prec, wb, y);
    clock_gettime(CLOCK_MONOTONIC, &_ts1);
    g_updtmp_wby_sec += (double)(_ts1.tv_sec - _ts0.tv_sec) + (double)(_ts1.tv_nsec - _ts0.tv_nsec) * 1e-9; }

    int i1 = 0;
    arb_set_d(t1, 1e300);
    arb_zero(t2);
    for (int i = 0; i < n; i++) {
        arb_t absv;
        arb_init(absv);
        arb_abs(absv, y+(i));
        if (arb_cmp_mid(absv, t1) < 0) {
            i1 = i;
            arb_set(t1, absv);
        }
        arb_max_mid(t2, t2, absv);
        arb_clear(absv);
    }

    if (idb >= 2) {
        double d1, d2; int n1_v, n2_v;
        mpdecmd(t1, &d1, &n1_v);
        mpdecmd(t2, &d2, &n2_v);
        fprintf(stderr, "Iteration %8d  updtmp: Min, max of y = %11.6f e %6d %11.6f e %6d\n",
                it, d1, n1_v, d2, n2_v);
    }

    { struct timespec _ts0, _ts1; clock_gettime(CLOCK_MONOTONIC, &_ts0);
    mxm(n, n, full_prec, wb, b);
    clock_gettime(CLOCK_MONOTONIC, &_ts1);
    g_updtmp_wbb_sec += (double)(_ts1.tv_sec - _ts0.tv_sec) + (double)(_ts1.tv_nsec - _ts0.tv_nsec) * 1e-9; }

    arb_zero(t2);
    for (int i = 0; i < n; i++) {
        arb_t absv;
        arb_init(absv);
        arb_abs(absv, AM(b, i1, i, n));
        arb_max_mid(t2, t2, absv);
        arb_clear(absv);
    }

    arb_t threshold;
    arb_init(threshold);
    arb_mul(threshold, t2, eps, full_prec);

    if (arb_cmp_mid(t1, threshold) <= 0) {
        if (idb >= 2) {
            double d1; int n1_v;
            mpdecmd(t1, &d1, &n1_v);
            fprintf(stderr, "Iteration %8d  updtmp: Small value in y = %11.6f e %6d\n", it, d1, n1_v);
        }
        *izm = 1;
    }

    { struct timespec _ts0, _ts1; clock_gettime(CLOCK_MONOTONIC, &_ts0);
    mxm_colexp(n, n1, full_prec, wa, h);
    clock_gettime(CLOCK_MONOTONIC, &_ts1);
    g_updtmp_wah_sec += (double)(_ts1.tv_sec - _ts0.tv_sec) + (double)(_ts1.tv_nsec - _ts0.tv_nsec) * 1e-9; }

    arb_clear(threshold);
    arb_clear(t1); arb_clear(t2);
    clock_gettime(CLOCK_MONOTONIC, &_tum1);
    g_updtmp_sec += (double)(_tum1.tv_sec - _tum0.tv_sec) + (double)(_tum1.tv_nsec - _tum0.tv_nsec) * 1e-9;
    g_updtmp_calls++;
}

/* ================================================================
 * dp_swap_and_rotate — swap adjacent rows im,im+1 and apply Givens
 * ================================================================ */

static void dp_swap_and_rotate(int n, int im, double *da, double *db,
                               double *dh, double *dy)
{
    int im1 = im + 1;
    double tmp;

    tmp = dy[im]; dy[im] = dy[im1]; dy[im1] = tmp;

    for (int k = 0; k < n; k++) {
        tmp = DM(da, im, k, n); DM(da, im, k, n) = DM(da, im1, k, n); DM(da, im1, k, n) = tmp;
        tmp = DM(db, im, k, n); DM(db, im, k, n) = DM(db, im1, k, n); DM(db, im1, k, n) = tmp;
    }

    for (int k = 0; k < n - 1; k++) {
        tmp = DM(dh, im, k, n); DM(dh, im, k, n) = DM(dh, im1, k, n); DM(dh, im1, k, n) = tmp;
    }

    if (im <= n - 3) {
        double t1 = DM(dh, im, im, n), t2 = DM(dh, im, im1, n);
        double t3 = sqrt(t1 * t1 + t2 * t2);
        if (t3 > 0) {
            t1 /= t3; t2 /= t3;
            for (int ii = im; ii < n; ii++) {
                double a = DM(dh, ii, im, n), b = DM(dh, ii, im1, n);
                DM(dh, ii, im, n) = t1 * a + t2 * b;
                DM(dh, ii, im1, n) = -t2 * a + t1 * b;
            }
        }
    }

    g_predicted_swap_cnt++;
}

/* ================================================================
 * predicted_swap — Givens-aware bidirectional insertion sort
 * ================================================================ */

static void predicted_swap(int n, double *da, double *db, double *dh, double *dy)
{
    /* Forward pass: sort left to right */
    for (int i = 1; i < n - 1; i++) {
        int j = i;
        while (j > 0) {
            int r = j - 1;
            double os = fabs(DM(dh, r, r, n));
            double nr, nr1, dv = 0;
            if (r <= n - 3) {
                os += fabs(DM(dh, r+1, r+1, n));
                double a = DM(dh, r+1, r, n), b = DM(dh, r+1, r+1, n);
                dv = sqrt(a*a + b*b);
                nr = dv;
                nr1 = (dv > 1e-300) ? (-DM(dh, r, r, n)*b + DM(dh, r, r+1, n)*a)/dv : DM(dh, r, r+1, n);
            } else {
                nr = DM(dh, r+1, r, n);
                nr1 = 0.0;
            }
            if (fabs(nr) + fabs(nr1) >= os) break;

            int ok = 1;
            if (r > 0 && r <= n - 3 && dv > 1e-300) {
                double x = DM(dh, r+1, r-1, n);
                double d2 = sqrt(x*x + dv*dv);
                double A = fabs(DM(dh, r-1, r-1, n));
                if (d2 > 1e-300) {
                    if (d2 + dv * A / d2 - A - dv >= 0.0) ok = 0;
                }
            }
            if (ok) { dp_swap_and_rotate(n, r, da, db, dh, dy); j--; }
            else break;
        }
    }

    /* Backward pass: sort right to left */
    for (int i = n - 3; i >= 0; i--) {
        int j = i;
        while (j < n - 2) {
            int r = j;
            double os = fabs(DM(dh, r, r, n));
            double nr, nr1, dv = 0;
            double av = 0, bv = 0;
            if (r <= n - 3) {
                os += fabs(DM(dh, r+1, r+1, n));
                av = DM(dh, r+1, r, n); bv = DM(dh, r+1, r+1, n);
                dv = sqrt(av*av + bv*bv);
                nr = dv;
                nr1 = (dv > 1e-300) ? (-DM(dh, r, r, n)*bv + DM(dh, r, r+1, n)*av)/dv : DM(dh, r, r+1, n);
            } else {
                nr = DM(dh, r+1, r, n);
                nr1 = 0.0;
            }
            if (fabs(nr) + fabs(nr1) >= os) break;

            int ok = 1;
            if (r <= n - 4 && dv > 1e-300) {
                double y = (-bv * DM(dh, r+2, r, n) + av * DM(dh, r+2, r+1, n)) / dv;
                double B = fabs(DM(dh, r+2, r+2, n));
                double C = fabs(nr1);
                double d3 = sqrt(y*y + B*B);
                if (d3 > 1e-300) {
                    if (d3 + B * C / d3 - C - B >= 0.0) ok = 0;
                }
            }
            if (ok) { dp_swap_and_rotate(n, r, da, db, dh, dy); j++; }
            else break;
        }
    }
}

/* ================================================================
 * strip_radius — zero error balls on wa/wb (Fortran has no balls)
 * ================================================================ */

static void strip_radius(int n, arb_ptr wa, arb_ptr wb)
{
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            mag_zero(arb_radref(AM(wa, i, j, n)));
            mag_zero(arb_radref(AM(wb, i, j, n)));
        }
}

/* ================================================================
 * pslqm3 — main loop
 * ================================================================ */

static int pslqm3(int idb, int n, slong full_prec, slong mpm_prec,
                  int ndr, int nrb, int nep, int nepm, int itm_param,
                  arb_ptr x, arb_ptr r)
{
    int n1 = n - 1;
    int iq = 0, it = 0, its = 0, imq = 0, izm = 0, izmm = 0, izd = 0;
    int nsq = 8;
    int ipm = (g_ipm_override > 0) ? g_ipm_override : 10;
    int ipi = 500;
    double dreps = 1e-10;

    /* 4-level IP arrays */
    arb_ptr ih = NULL, iy =NULL;
    fmpz_mat_t ia_accum, ib_accum, ia_tmp;
    int ip_active = (g_ip_prec > 0);
    if (ip_active) {
        ih = _arb_vec_init((size_t)n * n1);
        iy =_arb_vec_init(n);
        fmpz_mat_init(ia_accum, n, n);
        fmpz_mat_init(ib_accum, n, n);
        fmpz_mat_init(ia_tmp, n, n);
        fprintf(stderr, "4-level PSLQ enabled: ip_prec=%ld bits (%ld digits)\n",
                g_ip_prec, (slong)(g_ip_prec * 0.30103));
    }

    int cnt_izd0 = 0, cnt_izd1 = 0, cnt_izd2 = 0, cnt_fullmp = 0, cnt_savedp = 0;
    int cnt_fullmp_izmm1 = 0, cnt_fullmp_izd2 = 0, cnt_fullmp_mpm = 0;
    g_predicted_swap_cnt = 0;
    g_pairsel_swap_cnt = 0;
    g_izmm_wy_small = 0;
    g_izmm_wy_ratio = 0;
    g_izmm_wab_large = 0;

    /* Allocate arrays */
    double *da = (double *)calloc((size_t)n * n, sizeof(double));
    double *db = (double *)calloc((size_t)n * n, sizeof(double));
    double *dh = (double *)calloc((size_t)n * n, sizeof(double));
    double *dsa = (double *)calloc((size_t)n * n, sizeof(double));
    double *dsb = (double *)calloc((size_t)n * n, sizeof(double));
    double *dsh = (double *)calloc((size_t)n * n, sizeof(double));
    double *dsyq = (double *)calloc((size_t)n * nsq, sizeof(double));
    double *dy = (double *)calloc(n, sizeof(double));
    double *dsy = (double *)calloc(n, sizeof(double));

    arb_ptr b = _arb_vec_init((size_t)n * n);
    arb_ptr h = _arb_vec_init((size_t)n * n1);
    arb_ptr y = _arb_vec_init(n);

    arb_ptr wa = _arb_vec_init((size_t)n * n);
    arb_ptr wb = _arb_vec_init((size_t)n * n);
    arb_ptr wh = _arb_vec_init((size_t)n * n1);
    arb_ptr wy = _arb_vec_init(n);
    arb_ptr wsyq = _arb_vec_init((size_t)n * nsq);

    arb_t eps, epsm, wn, w1, w2;
    arb_init(eps); arb_init(epsm); arb_init(wn); arb_init(w1); arb_init(w2);
    arb_zero(wn);

    {
        arb_t ten;
        arb_init(ten);
        arb_set_d(ten, 10.0);
        if (nep >= 0) {
            arb_pow_ui(eps, ten, (ulong)nep, full_prec);
        } else {
            arb_pow_ui(eps, ten, (ulong)(-nep), full_prec);
            arb_inv(eps, eps, full_prec);
        }
        arb_clear(ten);
    }

    {
        arb_t ten;
        arb_init(ten);
        mp_set_d(ten, 10.0);
        if (nepm >= 0) {
            arb_pow_ui(epsm, ten, (ulong)nepm, mpm_prec);
        } else {
            arb_pow_ui(epsm, ten, (ulong)(-nepm), mpm_prec);
            mp_inv(epsm, epsm, mpm_prec);
        }
        arb_clear(ten);
    }

    if (idb >= 2) fprintf(stderr, "PSLQM3 integer relation detection: n = %d\n", n);

    if (idb >= 2) fprintf(stderr, "Iteration %8d  MP initialization\n", it);
    initmp(idb, n, full_prec, b, h, x, y);

    /* State machine (Fortran tag origins):
     *   full_reinit=100  check_wy=110  dp_iter=120
     *   mpm_section=130  mpm_iter=140  detect=150  cleanup=160 */
full_reinit: /* tag100 */
    dynrange(n, mpm_prec, y, w1);
    if (idb >= 2) {
        double d1; int n1_v;
        mpdecmd(w1, &d1, &n1_v);
        fprintf(stderr, "Iteration %8d  Min/max ratio of y = %11.6f e %6d\n", it, d1, n1_v);
    }

    initmpm(idb, n, nsq, mpm_prec, wa, wb, wh, wy, wsyq, h, y);

    if (ip_active)
        initip(n, g_ip_prec, ih, iy, wh, wy, ia_accum, ib_accum);

    /* ======== tag 110 ======== */
check_wy:
    dynrangem(n, mpm_prec, wy, w1);

    {
        arb_t dreps_arb;
        arb_init(dreps_arb);
        mp_set_d(dreps_arb, dreps);
        if (arb_cmp_mid(w1, dreps_arb) < 0 || izd == 2) {
            arb_clear(dreps_arb);
            goto mpm_section;
        }
        arb_clear(dreps_arb);
    }

    /* Start DP iterations */
    if (ip_active)
        initdp(idb, n, nsq, g_ip_prec, da, db, dh, dy, dsyq, ih, iy);
    else
        initdp(idb, n, nsq, mpm_prec, da, db, dh, dy, dsyq, wh, wy);

    savedp(n, da, db, dh, dy, dsa, dsb, dsh, dsy);
    its = it;

    lqdp(n, n1, n, dh);

    /* ======== DP iteration loop: tag 120 ======== */
dp_iter:
    it = it + 1;
    if (idb >= 4 || (idb >= 2 && it % ipi == 0))
        fprintf(stderr, "Iteration %8d\n", it);

    iterdp(idb, it, n, nsq, da, db, dh, dsyq, dy, &imq, &izd);

    if (izd == 0) {
        cnt_izd0++;
        if ((it - its) % ipm == 0) {
            savedp(n, da, db, dh, dy, dsa, dsb, dsh, dsy);
            cnt_savedp++;
            its = it;
        }

        if (g_strategy == 1)
            predicted_swap(n, da, db, dh, dy);

        goto dp_iter;
    } else {
        /* izd == 1 or 2 */
        if (izd == 1) cnt_izd1++;
        if (izd == 2) {
            cnt_izd2++;
            it = its;
            if (!ip_active)
                savedp(n, dsa, dsb, dsh, dsy, da, db, dh, dy);
        }

        if (ip_active) {
            /* 4-level: DP -> IP -> MPM -> Full */
            int izip = 0;
            if (izd == 2) {
                izip = 1;
                updtip(it, n, g_ip_prec, da, db, dreps,
                       ia_accum, ib_accum, ia_tmp, ih, iy,
                       mpm_prec, wy, &izip);
                izip = 1;
            } else {
                updtip(it, n, g_ip_prec, da, db, dreps,
                       ia_accum, ib_accum, ia_tmp, ih, iy,
                       mpm_prec, wy, &izip);
                if (izip == 0)
                    goto check_wy;
            }
            flush_ip_to_mpm(idb, it, n, mpm_prec, epsm, dreps,
                            ia_accum, ib_accum, wa, wb, wh, wy, &izmm);
            initip(n, g_ip_prec, ih, iy, wh, wy, ia_accum, ib_accum);
            if (izd == 2 && izmm == 0) izd = 1;
        } else {
            updtmpm(idb, it, n, mpm_prec, epsm, da, db, dreps, wa, wb, wh, wy, &izmm);
        }

        if (izmm == 0 && izd != 2) {
            goto check_wy;
        } else if (izmm == 1 || izd == 2) {
            cnt_fullmp++;
            if (izmm == 1) cnt_fullmp_izmm1++;
            if (izd == 2) cnt_fullmp_izd2++;
            fprintf(stderr, "FMP_TRIGGER: it=%d fmp=%d izmm=%d izd=%d since_last=%d\n",
                it, cnt_fullmp, izmm, izd, it - g_last_fmp_it);
            g_last_fmp_it = it;

            strip_radius(n, wa, wb);

            if (idb >= 2) fprintf(stderr, "Iteration %8d  MP update\n", it);
            updtmp(idb, it, n, full_prec, mpm_prec, wa, wb, eps, b, h, y, &izm);

            bound_fn(n, mpm_prec, wh, w1);
            double d3; int n3;
            mpdecmd(w1, &d3, &n3);
            arb_max_mid(wn, wn, w1);
            double d4; int n4;
            mpdecmd(wn, &d4, &n4);
            if (idb >= 2)
                fprintf(stderr, "Iteration %8d  Norm bound = %11.6f e %5d   Max. bound = %11.6f e %5d\n",
                        it, d3, n3, d4, n4);

            if (dplog10_arb(wn) > nrb) {
                if (idb >= 1) fprintf(stderr, "Norm bound limit exceeded. %d\n", nrb);
                goto cleanup;
            }
            if (it > itm_param) {
                if (idb >= 1) fprintf(stderr, "Iteration limit exceeded %d\n", itm_param);
                goto cleanup;
            }

            if (izm == 0) goto full_reinit;
            else if (izm == 1) goto detect;
            else goto cleanup;
        } else if (izmm == 2) {
            goto cleanup;
        }
    }

    /* ======== MPM iteration section: tag 130 ======== */
mpm_section:
    izd = 0;
    its = it;

    dynrangem(n, mpm_prec, wy, w1);
    {
        double d1; int n1_v;
        mpdecmd(w1, &d1, &n1_v);
        int ndrm = 25;
        int ndpm2 = ndrm - n1_v;
        int nwdsm2 = (int)(ndpm2 / MPDPW) + 2;
        if (nwdsm2 < 4) nwdsm2 = 4;
        ndpm2 = (int)((nwdsm2 - 1) * MPDPW + 0.5);
        slong mpm_prec2 = (slong)nwdsm2 * MPNBT;

        if (idb >= 2)
            fprintf(stderr, "Iteration %8d  Start MPM iter: precision = %4d words, %5d digits\n",
                    it, nwdsm2, ndpm2);

        arb_t epsm2;
        arb_init(epsm2);
        arb_one(epsm2);
        mp_mul_2exp_si(epsm2, epsm2, 70 - (slong)nwdsm2 * MPNBT);

        setprec(mpm_prec2, n, nsq, wa, wb, wh, wsyq, wy);

        lqmpm(n, n1, mpm_prec2, wh);

    mpm_iter:
        it = it + 1;
        if (idb >= 2) fprintf(stderr, "Iteration %8d\n", it);

        itermpm(idb, it, n, nsq, mpm_prec2, epsm2, wa, wb, wh, wsyq, wy, &imq, &izmm);

        if (izmm == 0 && (it - its) % ipm == 0) {
            dynrangem(n, mpm_prec, wy, w1);
            izmm = 1;
        }

        if (izmm == 0) {
            goto mpm_iter;
        } else if (izmm == 1) {
            strip_radius(n, wa, wb);

            cnt_fullmp++;
            cnt_fullmp_mpm++;
            fprintf(stderr, "FMP_TRIGGER: it=%d fmp=%d izmm=%d izd=%d since_last=%d (mpm_path)\n",
                it, cnt_fullmp, izmm, izd, it - g_last_fmp_it);
            g_last_fmp_it = it;
            if (idb >= 2) fprintf(stderr, "Iteration %8d  MP update\n", it);
            updtmp(idb, it, n, full_prec, mpm_prec, wa, wb, eps, b, h, y, &izm);

            bound_fn(n, mpm_prec, wh, w1);
            {
                double d3_v; int n3_v;
                mpdecmd(w1, &d3_v, &n3_v);
                arb_max_mid(wn, wn, w1);
                double d4_v; int n4_v;
                mpdecmd(wn, &d4_v, &n4_v);
                if (idb >= 2)
                    fprintf(stderr, "Iteration %8d  Norm bound = %11.6f e %5d   Max. bound = %11.6f e %5d\n",
                            it, d3_v, n3_v, d4_v, n4_v);
            }
            if (dplog10_arb(wn) > nrb) {
                arb_clear(epsm2);
                goto cleanup;
            }
            if (it > itm_param) {
                arb_clear(epsm2);
                goto cleanup;
            }

            if (izm == 0) { arb_clear(epsm2); goto full_reinit; }
            else if (izm == 1) { arb_clear(epsm2); goto detect; }
            else { arb_clear(epsm2); goto cleanup; }
        } else if (izmm == 2) {
            arb_clear(epsm2);
            goto cleanup;
        }

        arb_clear(epsm2);
    }

    /* ======== Detection section: tag 150 ======== */
detect:
    {
        arb_t det_t1, det_t2;
        arb_init(det_t1); arb_init(det_t2);
        arb_set_d(det_t1, 1e300);
        arb_zero(det_t2);

        int j1 = 0;
        for (int j = 0; j < n; j++) {
            arb_t absv;
            arb_init(absv);
            arb_abs(absv, y+(j));
            if (arb_cmp_mid(absv, det_t1) < 0) {
                j1 = j;
                arb_set(det_t1, absv);
            }
            arb_max_mid(det_t2, det_t2, absv);
            arb_clear(absv);
        }

        for (int i = 0; i < n; i++)
            arb_set(r+(i), AM(b, j1, i, n));

        arb_zero(w1);
        for (int i = 0; i < n; i++) {
            arb_t ri_mpm;
            arb_init(ri_mpm);
            arb_set_round(ri_mpm, r+(i), mpm_prec);
            arb_addmul(w1, ri_mpm, ri_mpm, mpm_prec);
            arb_clear(ri_mpm);
        }
        arb_sqrt(w1, w1, mpm_prec);

        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
                if (j < n1)
                    arb_set_round(AM(wh, i, j, n1), h+(i * n1 + j), mpm_prec);
        bound_fn(n, mpm_prec, wh, w2);
        arb_max_mid(wn, wn, w2);

        if (idb >= 1) {
            double d1, d2, d3, d4;
            int n1_v, n2_v, n3_v, n4_v;
            mpdecmd(det_t1, &d1, &n1_v);
            mpdecmd(det_t2, &d2, &n2_v);
            mpdecmd(w1, &d3, &n3_v);
            mpdecmd(wn, &d4, &n4_v);
            fprintf(stderr, "Iteration %8d  Relation detected\n", it);
            fprintf(stderr, "Min, max of y = %11.6f e %6d %11.6f e %6d\n", d1, n1_v, d2, n2_v);
            fprintf(stderr, "Max. bound = %11.6f e %6d\n", d4, n4_v);
            fprintf(stderr, "Index of relation = %4d   Norm = %11.5f e %5d   Residual = %11.6f e %6d\n",
                    j1 + 1, d3, n3_v, d1, n1_v);
        }

        {
            double d1; int nn1;
            mpdecmd(det_t1, &d1, &nn1);
            if (d1 == 0.0) nn1 = nep;
            double d3; int nn3;
            mpdecmd(w1, &d3, &nn3);
            double d2; int nn2;
            mpdecmd(det_t2, &d2, &nn2);

            if (nn3 <= nrb && nn2 - nn1 >= ndr)
                iq = 1;
            else if (idb >= 2)
                fprintf(stderr, "Relation is too large or insufficient dynamic range in y at detection.\n");
        }

        if (iq == 1) {
            FILE *fp = fopen("/tmp/v29_rel.txt", "w");
            if (fp) {
                for (int i = 0; i < n; i++) {
                    slong nbits = arf_bits(arb_midref(r+(i)));
                    slong ndig = (nbits > 0) ? (slong)(nbits * 0.30103) + 10 : 50;
                    if (ndig < 50) ndig = 50;
                    char *s = arb_get_str(r+(i), ndig, 0);
                    fprintf(fp, "%s\n", s);
                    flint_free(s);
                }
                fclose(fp);
                fprintf(stderr, "Relation coefficients written to /tmp/v29_rel.txt\n");
            }
        }

        arb_clear(det_t1); arb_clear(det_t2);
    }

    /* ======== Final section: tag 160 ======== */
cleanup:
    free(da); free(db); free(dh); free(dsa); free(dsb); free(dsh);
    free(dsyq); free(dy); free(dsy);

    _arb_vec_clear(b, (size_t)n * n);
    _arb_vec_clear(h, (size_t)n * n1);
    _arb_vec_clear(y, n);
    _arb_vec_clear(wa, (size_t)n * n);
    _arb_vec_clear(wb, (size_t)n * n);
    _arb_vec_clear(wh, (size_t)n * n1);
    _arb_vec_clear(wy, n);
    if (ip_active) {
        _arb_vec_clear(ih, (size_t)n * n1);
        _arb_vec_clear(iy, n);
        fmpz_mat_clear(ia_accum);
        fmpz_mat_clear(ib_accum);
        fmpz_mat_clear(ia_tmp);
    }
    _arb_vec_clear(wsyq, (size_t)n * nsq);

    arb_clear(eps); arb_clear(epsm); arb_clear(wn); arb_clear(w1); arb_clear(w2);

    fprintf(stderr, "[PROFILE] mxm: %.1fs (%d calls)  mxmdm: %.1fs (%d calls)\n",
            g_mxm_sec, g_mxm_calls, g_mxmdm_sec, g_mxmdm_calls);
    fprintf(stderr, "[PROFILE] iterdp: %.1fs (%d calls)  itermpm: %.1fs (%d calls)  lqmpm: %.1fs (%d calls)\n",
            g_iterdp_sec, g_iterdp_calls, g_itermpm_sec, g_itermpm_calls, g_lqmpm_sec, g_lqmpm_calls);
    fprintf(stderr, "[PROFILE] updtmpm: %.1fs (%d calls)  updtmp: %.1fs (%d calls)\n",
            g_updtmpm_sec, g_updtmpm_calls, g_updtmp_sec, g_updtmp_calls);
    fprintf(stderr, "[PROFILE] updtmp sub: wby=%.1fs wbb=%.1fs wah=%.1fs\n",
            g_updtmp_wby_sec, g_updtmp_wbb_sec, g_updtmp_wah_sec);
    if (g_ip_prec > 0)
        fprintf(stderr, "[PROFILE] updtip: %.1fs (%d calls)  ip_flushes=%d\n",
                g_updtip_sec, g_updtip_calls, g_ip_flushes);
    fprintf(stderr, "[PROFILE] izd: 0=%d 1=%d(deps=%d,tmx1=%d) 2=%d  fullMP=%d  savedp=%d  ipm=%d\n",
            cnt_izd0, cnt_izd1, g_cnt_izd1_deps, g_cnt_izd1_tmx1, cnt_izd2, cnt_fullmp, cnt_savedp, ipm);

    fprintf(stderr, "PREDICTED_SWAPS: %ld\n", g_predicted_swap_cnt);

    {
        const char *_sf_path = getenv("V2_STATS_FILE");
        FILE *sf = fopen(_sf_path ? _sf_path : "/tmp/v2_stats.json", "w");
        if (sf) {
            fprintf(sf, "{\"iq\":%d,\"it\":%d,\"izd0\":%d,\"izd1\":%d,\"izd2\":%d,"
                        "\"fullmp\":%d,\"updtmpm\":%d,\"itermpm\":%d,\"mxmdm\":%d,"
                        "\"savedp\":%d,\"mxm_sec\":%.2f,\"mxmdm_sec\":%.2f,"
                        "\"updtmpm_sec\":%.2f,\"itermpm_sec\":%.2f,"
                        "\"iterdp_sec\":%.2f,\"iterdp_calls\":%d,"
                        "\"izd1_deps\":%d,\"izd1_tmx1\":%d,"
                        "\"ipm\":%d,\"strat\":%d,"
                        "\"predicted_swaps\":%ld,"
                        "\"fullmp_izmm1\":%d,\"fullmp_izd2\":%d,\"fullmp_mpm\":%d,"
                        "\"izmm_wy_small\":%d,\"izmm_wy_ratio\":%d,\"izmm_wab_large\":%d}\n",
                    iq, it, cnt_izd0, cnt_izd1, cnt_izd2,
                    cnt_fullmp, g_updtmpm_calls, g_itermpm_calls, g_mxmdm_calls,
                    cnt_savedp, g_mxm_sec, g_mxmdm_sec,
                    g_updtmpm_sec, g_itermpm_sec,
                    g_iterdp_sec, g_iterdp_calls,
                    g_cnt_izd1_deps, g_cnt_izd1_tmx1,
                    ipm, g_strategy,
                    g_predicted_swap_cnt,
                    cnt_fullmp_izmm1, cnt_fullmp_izd2, cnt_fullmp_mpm,
                    g_izmm_wy_small, g_izmm_wy_ratio, g_izmm_wab_large);
            fclose(sf);
        }
    }

    return iq;
}

/* ================================================================
 * C API entry point — ctypes compatible
 * ================================================================ */

int pslqm3_c(
    int n,
    const char **x_strs,
    int digits,
    int ndpm,
    int ndr,
    int nrb,
    int itm,
    double unused1,
    int64_t *result,
    int strat)
{
    (void)unused1;
    memset(result, 0, (size_t)n * sizeof(int64_t));

    /* Threading */
    {
        const char *thr = getenv("THREADS");
        if (!thr) thr = getenv("STRAT_THREADS");
        int nt = thr ? atoi(thr) : 1;
        if (nt > 1) { g_colpar_threads = nt; flint_set_num_threads(nt); }
    }

    /* Strategy: env var overrides function param */
    {
        const char *ns = getenv("STRATEGY");
        if (!ns) ns = getenv("NUDGE_STRATEGY");
        if (ns) {
            if (strcmp(ns, "standard") == 0 || strcmp(ns, "0") == 0)
                g_strategy = 0;
            else if (strcmp(ns, "predicted_swap") == 0 || strcmp(ns, "114") == 0 || strcmp(ns, "1") == 0)
                g_strategy = 1;
            else
                g_strategy = atoi(ns);
        } else {
            g_strategy = (strat == 0) ? 0 : 1;
        }
        fprintf(stderr, "STRATEGY=%s (%d)\n",
                g_strategy == 0 ? "standard" : "predicted_swap", g_strategy);
    }

    /* IPM override */
    {
        const char *ipm_env = getenv("IPM_OVERRIDE");
        g_ipm_override = ipm_env ? atoi(ipm_env) : 0;
        if (g_ipm_override > 0)
            fprintf(stderr, "IPM_OVERRIDE=%d\n", g_ipm_override);
    }

    /* IP precision */
    {
        const char *ip_env = getenv("IP_BITS");
        if (!ip_env) ip_env = getenv("QP_BITS");
        if (!ip_env) ip_env = getenv("QP_PREC");
        g_ip_prec = ip_env ? atol(ip_env) : 0;
        if (g_ip_prec > 0)
            fprintf(stderr, "IP_BITS=%ld\n", g_ip_prec);
    }

    /* Compute precision parameters */
    int nwds = (int)(digits / MPDPW + 2);
    int nwdsm = (int)(ndpm / MPDPW + 2);
    slong full_prec = (slong)nwds * MPNBT;
    slong mpm_prec = (slong)nwdsm * MPNBT;
    int nep = 30 - digits;
    int nepm = 20 - ndpm;

    fprintf(stderr, "pslqm3: n=%d digits=%d ndpm=%d full_prec=%ld mpm_prec=%ld\n",
            n, digits, ndpm, (long)full_prec, (long)mpm_prec);

    /* Parse input strings into arb_t */
    arb_ptr x = _arb_vec_init(n);
    for (int i = 0; i < n; i++)
        arb_set_str(x+(i), x_strs[i], full_prec);

    arb_ptr r = _arb_vec_init(n);

    int iq = pslqm3(1, n, full_prec, mpm_prec,
                    ndr, nrb, nep, nepm, itm, x, r);

    if (iq == 1) {
        for (int i = 0; i < n; i++) {
            const arf_struct *mid = arb_midref(r+(i));
            if (arf_is_zero(mid)) {
                result[i] = 0;
            } else {
                fmpz_t coeff;
                fmpz_init(coeff);
                arf_t rounded;
                arf_init(rounded);
                arf_nint(rounded, mid);
                arf_get_fmpz(coeff, rounded, ARF_RND_NEAR);

                if (fmpz_fits_si(coeff)) {
                    result[i] = fmpz_get_si(coeff);
                } else {
                    result[i] = (fmpz_sgn(coeff) > 0) ? INT64_MAX : INT64_MIN;
                }
                arf_clear(rounded);
                fmpz_clear(coeff);
            }
        }
    }

    _arb_vec_clear(x, n);
    _arb_vec_clear(r, n);

    return iq;
}
