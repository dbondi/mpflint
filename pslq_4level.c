/*
 * pslq_4level.c — Intermediate precision layer for 4-level PSLQ
 *
 * Adds an intermediate precision (IP) layer between DP (53 bits) and MPM (~10K bits):
 *   DP -> IP (~200-500 bits) -> MPM -> Full MP
 *
 * The IP layer absorbs DP flushes cheaply (small matrix multiply at ip_prec)
 * and only flushes to MPM when IP precision exhausts. This reduces the number
 * of expensive mpm_prec matrix multiplies.
 *
 * Controlled by IP_BITS env var (0 = disabled, default = 0).
 * Optimal value: H_col_spread / 6, typically 130-200 bits.
 * Requires: H_spread > 370 bits AND baseline izd=2 = 0.
 *
 * Works for Poisson psi2 problems (e.g. psi_s24) where H_spread is wide.
 * Does NOT work for Poisson phi2 problems (e.g. phi_s29) where izd=2 > 0.
 */

static slong g_ip_prec = 0;
static int g_ip_max_absorb = 0;
static double g_updtip_sec = 0.0;
static int g_updtip_calls = 0;
static int g_ip_flushes = 0;
static int g_ip_dp_absorbed = 0;

static void initip(int n, slong ip_prec, arb_ptr ih, arb_ptr iy,
                   arb_srcptr wh, arb_srcptr wy,
                   fmpz_mat_t ia_accum, fmpz_mat_t ib_accum)
{
    int n1 = n - 1;

    arf_t afmax, af1, afscale, afone;
    arf_init(afmax); arf_init(af1); arf_init(afscale); arf_init(afone);
    arf_one(afone); arf_zero(afmax);
    for (int i = 0; i < n; i++) {
        arf_abs(af1, arb_midref(wy + i));
        if (arf_cmp(af1, afmax) > 0) arf_set(afmax, af1);
    }
    arf_div(afscale, afone, afmax, ip_prec, ARF_RND_NEAR);

    arb_t scale;
    arb_init(scale);
    arf_set(arb_midref(scale), afscale);
    for (int i = 0; i < n; i++) {
        mp_set_round(iy + i, wy + i, ip_prec);
        mp_mul(iy + i, scale, iy + i, ip_prec);
    }
    arb_clear(scale);
    arf_clear(afmax); arf_clear(af1); arf_clear(afscale); arf_clear(afone);

    for (int j = 0; j < n1; j++)
        for (int i = 0; i < n; i++)
            mp_set_round(ih + (i * n1 + j), wh + (i * n1 + j), ip_prec);

    fmpz_mat_one(ia_accum);
    fmpz_mat_one(ib_accum);
}

static void updtip(int it, int n, slong ip_prec,
                   const double *da, const double *db,
                   double dreps,
                   fmpz_mat_t ia_accum, fmpz_mat_t ib_accum, fmpz_mat_t ia_tmp,
                   arb_ptr ih, arb_ptr iy,
                   slong mpm_prec, arb_ptr wy,
                   int *izip)
{
    struct timespec _t0, _t1;
    clock_gettime(CLOCK_MONOTONIC, &_t0);
    int n1 = n - 1;
    *izip = 0;

    fmpz_mat_t fDA, fDB;
    fmpz_mat_init(fDA, n, n);
    fmpz_mat_init(fDB, n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            fmpz_set_d(fmpz_mat_entry(fDA, i, j), DM(da, i, j, n));
            fmpz_set_d(fmpz_mat_entry(fDB, i, j), DM(db, i, j, n));
        }

    mxmdm_vec(n, ip_prec, db, iy);

    if (wy != NULL)
        mxmdm_vec(n, mpm_prec, db, wy);

    {
        arb_t t1, t2;
        arb_init(t1); arb_init(t2);
        mp_set_d(t1, 1e300);
        arb_zero(t2);
        for (int i = 0; i < n; i++) {
            arb_t absv; arb_init(absv);
            mp_abs(absv, iy + i);
            arb_min_mid(t1, t1, absv);
            arb_max_mid(t2, t2, absv);
            arb_clear(absv);
        }
        arb_t ratio; arb_init(ratio);
        mp_div(ratio, t1, t2, ip_prec);
        arb_t dreps_arb; arb_init(dreps_arb);
        mp_set_d(dreps_arb, dreps);
        if (arb_cmp_mid(ratio, dreps_arb) < 0)
            *izip = 1;

        arb_t epsi; arb_init(epsi);
        arb_one(epsi);
        mp_mul_2exp_si(epsi, epsi, 70 - ip_prec);
        if (arb_cmp_mid(t1, epsi) <= 0)
            *izip = 1;
        arb_clear(epsi);
        arb_clear(t1); arb_clear(t2); arb_clear(ratio); arb_clear(dreps_arb);
    }

    fmpz_mat_mul(ia_tmp, fDA, ia_accum);
    fmpz_mat_swap(ia_accum, ia_tmp);
    fmpz_mat_mul(ia_tmp, fDB, ib_accum);
    fmpz_mat_swap(ib_accum, ia_tmp);

    mxmdm_fmpz(n, n1, ip_prec, fDA, ih);

    if (*izip) {
        arb_t _imin, _imax;
        arb_init(_imin); arb_init(_imax);
        mp_set_d(_imin, 1e300); arb_zero(_imax);
        for (int i = 0; i < n; i++) {
            arb_t absv; arb_init(absv);
            mp_abs(absv, iy + i);
            arb_min_mid(_imin, _imin, absv);
            arb_max_mid(_imax, _imax, absv);
            arb_clear(absv);
        }
        slong ia_max_bits = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                slong b = fmpz_bits(fmpz_mat_entry(ia_accum, i, j));
                if (b > ia_max_bits) ia_max_bits = b;
            }
        double d1, d2; int n1v, n2v;
        mpdecmd(_imin, &d1, &n1v);
        mpdecmd(_imax, &d2, &n2v);
        fprintf(stderr, "IP_EXHAUST: it=%d iy_min=%.2fe%d iy_max=%.2fe%d ia_bits=%ld dp_absorbed=%d\n",
                it, d1, n1v, d2, n2v, ia_max_bits, g_ip_dp_absorbed);
        arb_clear(_imin); arb_clear(_imax);
    }

    {
        slong max_bits = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                slong b = fmpz_bits(fmpz_mat_entry(ia_accum, i, j));
                if (b > max_bits) max_bits = b;
                b = fmpz_bits(fmpz_mat_entry(ib_accum, i, j));
                if (b > max_bits) max_bits = b;
            }
        if (max_bits > ip_prec / 2)
            *izip = 1;
    }

    fmpz_mat_clear(fDA); fmpz_mat_clear(fDB);

    g_ip_dp_absorbed++;

    if (g_ip_max_absorb > 0 && g_ip_dp_absorbed >= g_ip_max_absorb)
        *izip = 1;

    g_updtip_calls++;
    clock_gettime(CLOCK_MONOTONIC, &_t1);
    g_updtip_sec += (double)(_t1.tv_sec - _t0.tv_sec) + (double)(_t1.tv_nsec - _t0.tv_nsec) * 1e-9;
}

static void flush_ip_to_mpm(int idb, int it, int n, slong mpm_prec,
                            const arb_t epsm, double dreps,
                            fmpz_mat_t ia_accum, fmpz_mat_t ib_accum,
                            arb_ptr wa, arb_ptr wb, arb_ptr wh, arb_ptr wy,
                            int *izmm)
{
    int n1 = n - 1;
    *izmm = 0;

    {
        arb_t t1, t2;
        arb_init(t1); arb_init(t2);
        mp_set_d(t1, 1e300); arb_zero(t2);
        for (int i = 0; i < n; i++) {
            arb_t absv; arb_init(absv);
            mp_abs(absv, wy + i);
            arb_min_mid(t1, t1, absv);
            arb_max_mid(t2, t2, absv);
            arb_clear(absv);
        }
        if (arb_cmp_mid(t1, epsm) <= 0)
            *izmm = 1;
        arb_t w1, dreps_arb;
        arb_init(w1); arb_init(dreps_arb);
        mp_div(w1, t1, t2, mpm_prec);
        mp_set_d(dreps_arb, dreps);
        if (arb_cmp_mid(w1, dreps_arb) < 0)
            *izmm = 1;
        arb_clear(t1); arb_clear(t2); arb_clear(w1); arb_clear(dreps_arb);
    }

    mxmdm_fmpz(n, n, mpm_prec, ia_accum, wa);
    mxmdm_fmpz(n, n, mpm_prec, ib_accum, wb);
    mxmdm_fmpz(n, n1, mpm_prec, ia_accum, wh);

    {
        arb_t _t2, _tmx1;
        arb_init(_t2); arb_init(_tmx1);
        arb_zero(_t2);
        mp_inv(_tmx1, epsm, mpm_prec);
        for (int _j = 0; _j < n; _j++)
            for (int _i = 0; _i < n; _i++) {
                arb_t absv; arb_init(absv);
                mp_abs(absv, AM(wa, _i, _j, n));
                arb_max_mid(_t2, _t2, absv);
                mp_abs(absv, AM(wb, _i, _j, n));
                arb_max_mid(_t2, _t2, absv);
                arb_clear(absv);
            }
        if (arb_cmp_mid(_t2, _tmx1) > 0)
            *izmm = 1;
        arb_clear(_t2); arb_clear(_tmx1);
    }

    fmpz_mat_one(ia_accum);
    fmpz_mat_one(ib_accum);

    g_ip_flushes++;

    {
        arb_t _wmin, _wmax;
        arb_init(_wmin); arb_init(_wmax);
        mp_set_d(_wmin, 1e300); arb_zero(_wmax);
        for (int _i = 0; _i < n; _i++) {
            arb_t absv; arb_init(absv);
            mp_abs(absv, wy + _i);
            arb_min_mid(_wmin, _wmin, absv);
            arb_max_mid(_wmax, _wmax, absv);
            arb_clear(absv);
        }
        double d1, d2; int n1v, n2v;
        mpdecmd(_wmin, &d1, &n1v);
        mpdecmd(_wmax, &d2, &n2v);
        fprintf(stderr, "IP_FLUSH: it=%d absorbed=%d wy_min=%.2fe%d wy_max=%.2fe%d ratio=e%d\n",
                it, g_ip_dp_absorbed, d1, n1v, d2, n2v, n1v - n2v);
        arb_clear(_wmin); arb_clear(_wmax);
    }
    g_ip_dp_absorbed = 0;
}
