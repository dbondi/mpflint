/*
 * pslq_arb.h — MPFUN-compatible arb wrappers and utility functions
 *
 * Provides round-to-nearest (ARF_RND_NEAR) wrappers matching MPFUN2020
 * semantics, plus Fortran-compatible helpers (mpdecmd, arb_anint, etc.).
 */

#ifndef PSLQ_ARB_H
#define PSLQ_ARB_H

#include <math.h>
#include "flint/arb.h"
#include "flint/arf.h"

/* ================================================================
 * Constants
 * ================================================================ */

static const double MPDPW = 18.061799739838871713;  /* log10(2^60) */
static const int    MPNBT = 60;                      /* bits per MPFUN word */

/* ================================================================
 * Flat 2D array accessors (row-major, 0-indexed)
 * ================================================================ */

#define DM(arr, i, j, nc)  ((arr)[(size_t)(i)*(nc) + (j)])
#define AM(arr, i, j, nc)  ((arr) + (size_t)(i)*(nc) + (j))

/* ================================================================
 * MPFUN-compatible rounding wrappers
 *
 * MPFUN2020 uses round-to-nearest in mproun(). FLINT's arb_* operations
 * use ARF_RND_DOWN (truncation toward zero) for the midpoint. Over 50K+
 * iterations, this systematic bias accumulates enough to cause different
 * integer rounding in the Hermite reduction, corrupting the B matrix.
 *
 * These wrappers operate directly on arb_t midpoints with ARF_RND_NEAR,
 * matching MPFUN2020 semantics. Radii are zeroed since we never use them.
 * ================================================================ */

#define MPFUN_RND ARF_RND_NEAR

static inline void mp_mul(arb_t c, const arb_t a, const arb_t b, slong prec) {
    arf_mul(arb_midref(c), arb_midref(a), arb_midref(b), prec, MPFUN_RND);
    mag_zero(arb_radref(c));
}

static inline void mp_add(arb_t c, const arb_t a, const arb_t b, slong prec) {
    arf_add(arb_midref(c), arb_midref(a), arb_midref(b), prec, MPFUN_RND);
    mag_zero(arb_radref(c));
}

static inline void mp_sub(arb_t c, const arb_t a, const arb_t b, slong prec) {
    arf_sub(arb_midref(c), arb_midref(a), arb_midref(b), prec, MPFUN_RND);
    mag_zero(arb_radref(c));
}

static inline void mp_div(arb_t c, const arb_t a, const arb_t b, slong prec) {
    arf_div(arb_midref(c), arb_midref(a), arb_midref(b), prec, MPFUN_RND);
    mag_zero(arb_radref(c));
}

static inline void mp_sqrt(arb_t c, const arb_t a, slong prec) {
    arf_sqrt(arb_midref(c), arb_midref(a), prec, MPFUN_RND);
    mag_zero(arb_radref(c));
}

static inline void mp_addmul(arb_t c, const arb_t a, const arb_t b, slong prec) {
    arf_addmul(arb_midref(c), arb_midref(a), arb_midref(b), prec, MPFUN_RND);
    mag_zero(arb_radref(c));
}

static inline void mp_submul(arb_t c, const arb_t a, const arb_t b, slong prec) {
    arf_submul(arb_midref(c), arb_midref(a), arb_midref(b), prec, MPFUN_RND);
    mag_zero(arb_radref(c));
}

static inline void mp_inv(arb_t c, const arb_t a, slong prec) {
    arf_t one;
    arf_init(one);
    arf_one(one);
    arf_div(arb_midref(c), one, arb_midref(a), prec, MPFUN_RND);
    mag_zero(arb_radref(c));
    arf_clear(one);
}

static inline void mp_neg(arb_t c, const arb_t a) {
    arf_neg(arb_midref(c), arb_midref(a));
    mag_zero(arb_radref(c));
}

static inline void mp_abs(arb_t c, const arb_t a) {
    arf_abs(arb_midref(c), arb_midref(a));
    mag_zero(arb_radref(c));
}

static inline void mp_set(arb_t c, const arb_t a) {
    arf_set(arb_midref(c), arb_midref(a));
    mag_zero(arb_radref(c));
}

static inline void mp_set_d(arb_t c, double d) {
    arf_set_d(arb_midref(c), d);
    mag_zero(arb_radref(c));
}

static inline void mp_mul_2exp_si(arb_t c, const arb_t a, slong e) {
    arf_mul_2exp_si(arb_midref(c), arb_midref(a), e);
    mag_zero(arb_radref(c));
}

static inline void mp_set_round(arb_t c, const arb_t a, slong prec) {
    arf_set_round(arb_midref(c), arb_midref(a), prec, MPFUN_RND);
    mag_zero(arb_radref(c));
}

/* ================================================================
 * Utility functions
 * ================================================================ */

static inline double gam_pow_fn(double base, int exp) {
    double result = 1.0;
    while (exp > 0) {
        if (exp & 1) result *= base;
        base *= base;
        exp >>= 1;
    }
    return result;
}

static void mpdecmd(const arb_t a, double *d, int *n)
{
    const arf_struct *mid = arb_midref(a);
    if (arf_is_zero(mid)) {
        *d = 0.0;
        *n = 0;
        return;
    }
    if (!arf_is_finite(mid)) {
        *d = 9.9;
        *n = 999999;
        return;
    }

    slong E = ARF_EXP(mid);

    arf_t man_arf;
    arf_init(man_arf);
    arf_mul_2exp_si(man_arf, mid, -E);
    double man_d = arf_get_d(man_arf, ARF_RND_NEAR);
    arf_clear(man_arf);

    double log10_v = (double)E * 0.30102999566398119 + log10(fabs(man_d));
    *n = (int)floor(log10_v);

    double a_d = arf_get_d(mid, ARF_RND_NEAR);
    if (fabs(a_d) > 1e-300 && fabs(a_d) < 1e300) {
        *d = a_d / pow(10.0, (double)(*n));
    } else {
        double frac = log10_v - (double)(*n);
        *d = pow(10.0, frac);
        if (arf_sgn(mid) < 0) *d = -(*d);
    }
}

static double dplog10_arb(const arb_t a)
{
    double da;
    int ia;
    mpdecmd(a, &da, &ia);
    if (da == 0.0) return -999999.0;
    return log10(fabs(da)) + (double)ia;
}

static void arb_anint(arb_t result, const arb_t x)
{
    arf_t half, tmp;
    arf_init(half); arf_init(tmp);
    arf_set_d(half, 0.5);
    if (arf_sgn(arb_midref(x)) >= 0) {
        arf_add(tmp, arb_midref(x), half, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_floor(arb_midref(result), tmp);
    } else {
        arf_sub(tmp, arb_midref(x), half, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_ceil(arb_midref(result), tmp);
    }
    mag_zero(arb_radref(result));
    arf_clear(half); arf_clear(tmp);
}

static int arb_cmp_mid(const arb_t a, const arb_t b)
{
    return arf_cmp(arb_midref(a), arb_midref(b));
}

static void arb_min_mid(arb_t result, const arb_t a, const arb_t b)
{
    if (arb_cmp_mid(a, b) <= 0)
        arb_set(result, a);
    else
        arb_set(result, b);
}

static void arb_max_mid(arb_t result, const arb_t a, const arb_t b)
{
    if (arb_cmp_mid(a, b) >= 0)
        arb_set(result, a);
    else
        arb_set(result, b);
}

#endif /* PSLQ_ARB_H */
