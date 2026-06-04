/*
 * pslq_sort.c — Quicksort for doubles and arb_t (exact Fortran translation)
 *
 * Fortran arrays are 1-indexed; C arrays are 0-indexed.
 * ip[] returns 0-based indices.
 */

#include "pslq_arb.h"

static void qsortdp(int n, const double *a, int *ip)
{
    int i, iq, it, j, jq, jz, k, l;
    double s0;
    int ik[50], jk[50];

    for (i = 0; i < n; i++)
        ip[i] = i;

    if (n == 1) return;

    k = 1;
    ik[0] = 0;
    jk[0] = n - 1;

tag130:
    i = ik[k-1];
    j = jk[k-1];
    iq = i;
    jq = j;
    it = (i + j + 1) / 2;
    l = ip[j];
    ip[j] = ip[it];
    ip[it] = l;
    s0 = a[ip[j]];
    j = j - 1;

tag140:
    for (l = i; l <= j; l++) {
        if (s0 < a[ip[l]]) goto tag160;
    }
    i = j;
    goto tag190;

tag160:
    i = l;

    for (l = j; l >= i; l--) {
        if (s0 > a[ip[l]]) goto tag180;
    }
    j = i;
    goto tag190;

tag180:
    j = l;
    if (i >= j) goto tag190;
    l = ip[i];
    ip[i] = ip[j];
    ip[j] = l;
    goto tag140;

tag190:
    if (s0 >= a[ip[i]]) goto tag200;
    l = ip[jq];
    ip[jq] = ip[i];
    ip[i] = l;

tag200:
    k = k - 1;
    jz = 0;
    if (j == iq) goto tag210;
    k = k + 1;
    jk[k-1] = j;
    jz = 1;

tag210:
    i = i + 1;
    if (i == jq) goto tag220;
    k = k + 1;
    ik[k-1] = i;
    jk[k-1] = jq;
    if (jz == 0) goto tag220;
    if (j - iq >= jq - i) goto tag220;
    ik[k-2] = i;
    jk[k-2] = jq;
    ik[k-1] = iq;
    jk[k-1] = j;

tag220:
    if (k > 0) goto tag130;
}

static void qsortmpm(int n, arb_srcptr a, int *ip)
{
    int i, iq, it, j, jq, jz, k, l;
    int ik[50], jk[50];

    for (i = 0; i < n; i++)
        ip[i] = i;

    if (n == 1) return;

    k = 1;
    ik[0] = 0;
    jk[0] = n - 1;

tag130:
    i = ik[k-1];
    j = jk[k-1];
    iq = i;
    jq = j;
    it = (i + j + 1) / 2;
    l = ip[j];
    ip[j] = ip[it];
    ip[it] = l;
    int s0_idx = ip[j];
    j = j - 1;

tag140_m:
    for (l = i; l <= j; l++) {
        if (arb_cmp_mid(a + s0_idx, a + ip[l]) < 0) goto tag160_m;
    }
    i = j;
    goto tag190_m;

tag160_m:
    i = l;

    for (l = j; l >= i; l--) {
        if (arb_cmp_mid(a + s0_idx, a + ip[l]) > 0) goto tag180_m;
    }
    j = i;
    goto tag190_m;

tag180_m:
    j = l;
    if (i >= j) goto tag190_m;
    l = ip[i];
    ip[i] = ip[j];
    ip[j] = l;
    goto tag140_m;

tag190_m:
    if (arb_cmp_mid(a + s0_idx, a + ip[i]) >= 0) goto tag200_m;
    l = ip[jq];
    ip[jq] = ip[i];
    ip[i] = l;

tag200_m:
    k = k - 1;
    jz = 0;
    if (j == iq) goto tag210_m;
    k = k + 1;
    jk[k-1] = j;
    jz = 1;

tag210_m:
    i = i + 1;
    if (i == jq) goto tag220_m;
    k = k + 1;
    ik[k-1] = i;
    jk[k-1] = jq;
    if (jz == 0) goto tag220_m;
    if (j - iq >= jq - i) goto tag220_m;
    ik[k-2] = i;
    jk[k-2] = jq;
    ik[k-1] = iq;
    jk[k-1] = j;

tag220_m:
    if (k > 0) goto tag130;
}
