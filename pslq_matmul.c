/*
 * pslq_matmul.c — Matrix multiply infrastructure for PSLQ
 *
 * CRT-based and arb_mat_approx_mul paths for mxmdm (double×MPM)
 * and mxm (MPM×full-precision). Includes ColPar (column-parallel CRT)
 * with per-column exponents and pthread column splitting.
 */

#include <pthread.h>
#include "pslq_arb.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "flint/arb_mat.h"

extern double g_mxmdm_sec;
extern int    g_mxmdm_calls;
extern double g_mxm_sec;
extern int    g_mxm_calls;
extern int    g_colpar_threads;

static void arb_to_fmpz_vec(fmpz_mat_t fB, slong *b_min_exp_out,
                            arb_ptr b, int n1, int n2)
{
    slong b_min_exp = LONG_MAX;
    fmpz_t man, exp;
    fmpz_init(man); fmpz_init(exp);
    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++) {
            const arf_struct *m = arb_midref(b + (i * n2 + j));
            if (arf_is_zero(m) || !arf_is_finite(m)) continue;
            arf_get_fmpz_2exp(man, exp, m);
            slong e = fmpz_get_si(exp);
            if (e < b_min_exp) b_min_exp = e;
        }
    if (b_min_exp == LONG_MAX) b_min_exp = 0;
    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++) {
            const arf_struct *m = arb_midref(b + (i * n2 + j));
            if (arf_is_zero(m) || !arf_is_finite(m)) {
                fmpz_zero(fmpz_mat_entry(fB, i, j));
                continue;
            }
            arf_get_fmpz_2exp(man, exp, m);
            slong e = fmpz_get_si(exp);
            fmpz_set(fmpz_mat_entry(fB, i, j), man);
            slong up = e - b_min_exp;
            if (up > 0) fmpz_mul_2exp(fmpz_mat_entry(fB, i, j),
                                       fmpz_mat_entry(fB, i, j), up);
        }
    *b_min_exp_out = b_min_exp;
    fmpz_clear(man); fmpz_clear(exp);
}

static void fmpz_result_to_arb(arb_ptr b, const fmpz_mat_t fC,
                                slong b_min_exp, int n1, int n2, slong prec)
{
    fmpz_t exp_fmpz;
    fmpz_init(exp_fmpz);
    fmpz_set_si(exp_fmpz, b_min_exp);
    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++) {
            arf_set_fmpz_2exp(arb_midref(b + (i * n2 + j)),
                              fmpz_mat_entry(fC, i, j), exp_fmpz);
            arf_set_round(arb_midref(b + (i * n2 + j)),
                          arb_midref(b + (i * n2 + j)), prec, ARF_RND_NEAR);
            mag_zero(arb_radref(b + (i * n2 + j)));
        }
    fmpz_clear(exp_fmpz);
}

static void colpar_dispatch(int nt, const fmpz_mat_struct *fA, slong a_exp,
                            arb_ptr b, int n1, int n2, slong prec);

static void mxmdm_fmpz(int n1, int n2, slong mpm_prec,
                        const fmpz_mat_t fA, arb_ptr b)
{
    struct timespec _t0, _t1;
    clock_gettime(CLOCK_MONOTONIC, &_t0);

    int nt = g_colpar_threads;
    if (nt > 1 && n2 >= nt) {
        colpar_dispatch(nt, fA, 0, b, n1, n2, mpm_prec);
    } else {
        fmpz_mat_t fB, fC;
        slong b_min_exp;
        fmpz_mat_init(fB, n1, n2);
        arb_to_fmpz_vec(fB, &b_min_exp, b, n1, n2);
        fmpz_mat_init(fC, n1, n2);
        fmpz_mat_mul(fC, fA, fB);
        fmpz_result_to_arb(b, fC, b_min_exp, n1, n2, mpm_prec);
        fmpz_mat_clear(fB);
        fmpz_mat_clear(fC);
    }

    clock_gettime(CLOCK_MONOTONIC, &_t1);
    g_mxmdm_sec += (double)(_t1.tv_sec - _t0.tv_sec) + (double)(_t1.tv_nsec - _t0.tv_nsec) * 1e-9;
    g_mxmdm_calls++;
}

static void mxmdm_vec(int n, slong mpm_prec, const double *a, arb_ptr b)
{
    struct timespec _t0, _t1;
    clock_gettime(CLOCK_MONOTONIC, &_t0);

    arb_ptr c = _arb_vec_init(n);
    arb_ptr a_row = _arb_vec_init(n);
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++)
            arb_set_d(a_row + k, DM(a, i, k, n));
        arb_approx_dot(c + i, NULL, 0, a_row, 1, b, 1, n, mpm_prec);
    }
    for (int i = 0; i < n; i++)
        arb_swap(b + i, c + i);
    _arb_vec_clear(a_row, n);
    _arb_vec_clear(c, n);

    clock_gettime(CLOCK_MONOTONIC, &_t1);
    g_mxmdm_sec += (double)(_t1.tv_sec - _t0.tv_sec) + (double)(_t1.tv_nsec - _t0.tv_nsec) * 1e-9;
    g_mxmdm_calls++;
}

/* ================================================================
 * mxm — MPM x full-precision matrix multiply: b = a * b
 * ================================================================ */

static void arb_mat_wrap(arb_mat_t mat, arb_ptr data, int rows, int cols)
{
    arb_mat_init(mat, rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            arb_swap(arb_mat_entry(mat, i, j), data + (i * cols + j));
}

static void arb_mat_unwrap(arb_mat_t mat, arb_ptr data, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            arb_swap(arb_mat_entry(mat, i, j), data + (i * cols + j));
    arb_mat_clear(mat);
}

static void arb_to_fmpz_mat(fmpz_mat_t fm, slong *min_exp_out,
                            arb_srcptr data, int rows, int cols)
{
    slong min_exp = LONG_MAX;
    fmpz_t man, exp;
    fmpz_init(man); fmpz_init(exp);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) {
            const arf_struct *m = arb_midref(data + (i * cols + j));
            if (arf_is_zero(m) || !arf_is_finite(m)) continue;
            arf_get_fmpz_2exp(man, exp, m);
            slong e = fmpz_get_si(exp);
            if (e < min_exp) min_exp = e;
        }
    if (min_exp == LONG_MAX) min_exp = 0;

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) {
            const arf_struct *m = arb_midref(data + (i * cols + j));
            if (arf_is_zero(m) || !arf_is_finite(m)) {
                fmpz_zero(fmpz_mat_entry(fm, i, j));
                continue;
            }
            arf_get_fmpz_2exp(man, exp, m);
            slong e = fmpz_get_si(exp);
            fmpz_set(fmpz_mat_entry(fm, i, j), man);
            slong up = e - min_exp;
            if (up > 0) fmpz_mul_2exp(fmpz_mat_entry(fm, i, j),
                                       fmpz_mat_entry(fm, i, j), up);
        }

    *min_exp_out = min_exp;
    fmpz_clear(man); fmpz_clear(exp);
}

typedef struct {
    const fmpz_mat_struct *fA;
    slong a_exp;
    arb_ptr b;
    int n1, n2;
    int col_start, col_count, col_stride;
    slong full_prec;
    fmpz_mat_t fB_slice, fC_slice;
    slong *col_exp;
} colpar_arg_t;

static void *colpar_thread_fn(void *varg) {
    colpar_arg_t *a = (colpar_arg_t *)varg;
    int n1 = a->n1, n2 = a->n2;
    fmpz_t man, exp;
    fmpz_init(man); fmpz_init(exp);

    fmpz_mat_init(a->fB_slice, n1, a->col_count);
    a->col_exp = flint_malloc((size_t)a->col_count * sizeof(slong));

    for (int jj = 0; jj < a->col_count; jj++) {
        int j = a->col_start + jj * a->col_stride;
        slong min_e = LONG_MAX;
        for (int i = 0; i < n1; i++) {
            const arf_struct *m = arb_midref(a->b + (i * n2 + j));
            if (arf_is_zero(m) || !arf_is_finite(m)) continue;
            arf_get_fmpz_2exp(man, exp, m);
            slong e = fmpz_get_si(exp);
            if (e < min_e) min_e = e;
        }
        if (min_e == LONG_MAX) min_e = 0;
        a->col_exp[jj] = min_e;

        for (int i = 0; i < n1; i++) {
            const arf_struct *m = arb_midref(a->b + (i * n2 + j));
            if (arf_is_zero(m) || !arf_is_finite(m)) {
                fmpz_zero(fmpz_mat_entry(a->fB_slice, i, jj));
                continue;
            }
            arf_get_fmpz_2exp(man, exp, m);
            slong e = fmpz_get_si(exp);
            fmpz_set(fmpz_mat_entry(a->fB_slice, i, jj), man);
            slong up = e - min_e;
            if (up > 0) fmpz_mul_2exp(fmpz_mat_entry(a->fB_slice, i, jj),
                                       fmpz_mat_entry(a->fB_slice, i, jj), up);
        }
    }

    fmpz_mat_init(a->fC_slice, n1, a->col_count);
    fmpz_mat_mul(a->fC_slice, a->fA, a->fB_slice);

    fmpz_t exp_fmpz;
    fmpz_init(exp_fmpz);
    for (int jj = 0; jj < a->col_count; jj++) {
        int j = a->col_start + jj * a->col_stride;
        slong total_exp = a->a_exp + a->col_exp[jj];
        fmpz_set_si(exp_fmpz, total_exp);
        for (int i = 0; i < n1; i++) {
            arf_set_fmpz_2exp(arb_midref(a->b + (i * n2 + j)),
                              fmpz_mat_entry(a->fC_slice, i, jj), exp_fmpz);
            arf_set_round(arb_midref(a->b + (i * n2 + j)),
                          arb_midref(a->b + (i * n2 + j)), a->full_prec, ARF_RND_NEAR);
            mag_zero(arb_radref(a->b + (i * n2 + j)));
        }
    }
    fmpz_clear(exp_fmpz);
    fmpz_clear(man); fmpz_clear(exp);
    return NULL;
}

static void colpar_dispatch(int nt, const fmpz_mat_struct *fA, slong a_exp,
                            arb_ptr b, int n1, int n2, slong prec)
{
    flint_set_num_threads(1);
    pthread_t *threads = malloc((size_t)nt * sizeof(pthread_t));
    colpar_arg_t *args = malloc((size_t)nt * sizeof(colpar_arg_t));
    for (int t = 0; t < nt; t++) {
        args[t].fA = fA;
        args[t].a_exp = a_exp;
        args[t].b = b;
        args[t].n1 = n1;
        args[t].n2 = n2;
        args[t].col_start = t;
        args[t].col_stride = nt;
        args[t].col_count = (n2 - t + nt - 1) / nt;
        args[t].full_prec = prec;
        pthread_create(&threads[t], NULL, colpar_thread_fn, &args[t]);
    }
    for (int t = 0; t < nt; t++) {
        pthread_join(threads[t], NULL);
        fmpz_mat_clear(args[t].fB_slice);
        fmpz_mat_clear(args[t].fC_slice);
        flint_free(args[t].col_exp);
    }
    free(threads);
    free(args);
    flint_set_num_threads(g_colpar_threads);
}

static void mxm_colexp(int n1, int n2, slong full_prec, arb_srcptr a, arb_ptr b)
{
    fmpz_mat_t fA;
    slong a_exp;

    fmpz_mat_init(fA, n1, n1);
    arb_to_fmpz_mat(fA, &a_exp, a, n1, n1);

    int nt = g_colpar_threads;
    if (nt > 1 && n2 >= nt) {
        colpar_dispatch(nt, fA, a_exp, b, n1, n2, full_prec);
    } else {
        colpar_arg_t arg = { .fA = fA, .a_exp = a_exp, .b = b,
                             .n1 = n1, .n2 = n2, .col_start = 0,
                             .col_count = n2, .col_stride = 1, .full_prec = full_prec };
        colpar_thread_fn(&arg);
        fmpz_mat_clear(arg.fB_slice);
        fmpz_mat_clear(arg.fC_slice);
        flint_free(arg.col_exp);
    }

    fmpz_mat_clear(fA);
}

static void mxm(int n1, int n2, slong full_prec, arb_srcptr a, arb_ptr b)
{
    struct timespec _t0, _t1;
    clock_gettime(CLOCK_MONOTONIC, &_t0);

    slong max_bits = 0;
    for (int i = 0; i < n1 && max_bits <= 8000; i++)
        for (int j = 0; j < n2; j++) {
            slong bits = arf_bits(arb_midref(b + (i * n2 + j)));
            if (bits > max_bits) max_bits = bits;
        }

    if (max_bits > 8000 && n2 > 1) {
        fmpz_mat_t fA, fB, fC;
        slong a_exp, b_exp;

        fmpz_mat_init(fA, n1, n1);
        arb_to_fmpz_mat(fA, &a_exp, a, n1, n1);

        fmpz_mat_init(fB, n1, n2);
        arb_to_fmpz_mat(fB, &b_exp, b, n1, n2);

        fmpz_mat_init(fC, n1, n2);
        fmpz_mat_mul(fC, fA, fB);

        slong total_exp = a_exp + b_exp;
        fmpz_t exp_fmpz;
        fmpz_init(exp_fmpz);
        fmpz_set_si(exp_fmpz, total_exp);
        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2; j++) {
                arf_set_fmpz_2exp(arb_midref(b + (i * n2 + j)),
                                  fmpz_mat_entry(fC, i, j), exp_fmpz);
                arf_set_round(arb_midref(b + (i * n2 + j)),
                              arb_midref(b + (i * n2 + j)), full_prec, ARF_RND_NEAR);
                mag_zero(arb_radref(b + (i * n2 + j)));
            }
        fmpz_clear(exp_fmpz);
        fmpz_mat_clear(fA);
        fmpz_mat_clear(fB);
        fmpz_mat_clear(fC);
    } else {
        arb_mat_t A, B, C;
        arb_mat_wrap(A, (arb_ptr)a, n1, n1);
        arb_mat_wrap(B, (arb_ptr)b, n1, n2);
        arb_mat_init(C, n1, n2);
        arb_mat_approx_mul(C, A, B, full_prec);
        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2; j++)
                arb_swap(arb_mat_entry(C, i, j), b + (i * n2 + j));
        arb_mat_unwrap(A, (arb_ptr)a, n1, n1);
        arb_mat_clear(B);
        arb_mat_clear(C);
    }

    clock_gettime(CLOCK_MONOTONIC, &_t1);
    g_mxm_sec += (double)(_t1.tv_sec - _t0.tv_sec) + (double)(_t1.tv_nsec - _t0.tv_nsec) * 1e-9;
    g_mxm_calls++;
}
