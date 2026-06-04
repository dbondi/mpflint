# PSLQ Integer Relation Detection — C/FLINT Implementation

A high-performance reimplementation of Bailey's multipair PSLQ algorithm in C using the [FLINT](https://flintlib.org/) arbitrary-precision library, achieving 4–14× speedups over the reference Fortran/MPFUN2020 implementation across Poisson summation problems.

## How PSLQ Works

Given a vector of real numbers **x** = (x₁, x₂, ..., xₙ), PSLQ finds integer coefficients **m** = (m₁, m₂, ..., mₙ) such that m₁x₁ + m₂x₂ + ... + mₙxₙ = 0, or proves that no such relation exists below a given norm bound. This is used to discover BBP-type formulas, minimal polynomials, and Poisson summation identities.

The multipair variant operates at three precision levels in a cycle:

1. **DP (double precision, 53 bits)** — Fast inner loop. Selects disjoint pairs of rows using γ-weighted diagonal magnitudes, swaps them, applies Givens rotations, and performs Hermite reduction. Runs up to `ipm` iterations per batch.

2. **MPM (medium precision, ndpm digits)** — Applies the accumulated DP row operations to the medium-precision H, y, wa, wb matrices. Detects when DP precision is exhausted (izd=1) or overflows (izd=2), triggering a flush upward.

3. **Full MP (full precision, digits)** — Periodic high-precision sync. Updates the full-precision H, B, y matrices and checks for relation detection. This is the most expensive step (n² × full_prec matrix multiply).

The algorithm converges when min|y| drops below ε relative to max|y| × max|B_row|, at which point the corresponding row of B gives the integer relation.

## What This Implementation Changes

### 1. FLINT's CRT Architecture — Better Baseline

The single largest source of speedup comes from switching the arithmetic library from MPFUN2020 to FLINT.

FLINT's `fmpz_mat_mul` uses **CRT (Chinese Remainder Theorem) multi-modular arithmetic**: it reduces a big-integer matrix multiply into `(ndpm_bits / 59)` independent word-level (64-bit) matrix multiplies modulo small primes, then reconstructs the result via CRT. Because these word-level multiplies are independent, FLINT can parallelize them across CPU cores — both `mxmdm` (MPM-level matrix multiplies) and `mxm` (fullMP-level matrix multiplies) scale with thread count. The ColPar path additionally splits fullMP multiplies across columns using pthreads with per-column exponent tracking for numerical stability.

MPFUN2020 uses direct big-integer multiplication with Karatsuba/FFT for element-wise operations within its matrix multiply loops. MPFUN is thread-safe but cannot parallelize a single PSLQ run — there are no OpenMP directives in the PSLQ driver code.

### 2. ndpm Tuning — CRT Makes Low Precision Cheap

The medium precision level (ndpm) controls the bit-width of the MPM matrices. Each `mxmdm` call does a CRT-based matrix multiply at ndpm precision — FLINT decomposes this into `ceil(ndpm_bits / 59)` independent word-level n×n multiplies. Lowering ndpm directly reduces the number of passes.

Bailey's Fortran code uses conservative ndpm defaults (e.g., ndpm=3000 for s=24, which is `ceil(3000 × 3.32 / 59) ≈ 169` CRT passes per mxmdm). At ndpm=1000, this drops to `ceil(1000 × 3.32 / 59) ≈ 57` passes — 3× fewer.

The tradeoff: lower ndpm exhausts MPM precision faster, triggering more fullMP updates. But at small n (e.g., n=65) each fullMP is cheap (~4s), so the mxmdm savings dominate.

**psi_s24 ndpm sweep (predicted_swap, 1 thread):**

| ndpm | Wall time | mxmdm time | fullMP count |
|------|-----------|-----------|--------------|
| 3000 (Bailey default) | 481s | 827s | 7 |
| 1500 | 222s | — | 13 |
| 1000 | 176s | 991s | 18 |
| 700 | 178s | 1042s | 26 |

The optimal ndpm balances mxmdm savings against fullMP cost. For n=65 problems, ndpm=700–1000 is optimal. For larger n (S29, n=197) where each fullMP costs ~34s, the optimal ndpm stays closer to Bailey's default.

This optimization is **free** — it requires no code changes, just setting a lower ndpm value.

### 3. Intermediate Precision Layer (IP/QP)

The standard 3-level algorithm flushes DP operations directly to MPM, which is expensive when ndpm is large. The IP layer adds a buffer between DP and MPM:

```
DP (53 bits) → IP (~200-500 bits) → MPM (~10K bits) → Full MP
```

When DP precision exhausts (izd=1), instead of flushing to MPM at full mpm_prec cost, the DP operations are applied to the IP-precision copies of H and y. The IP layer accumulates these updates as exact integer matrices (qa, qb) and only flushes to MPM when the IP precision itself exhausts.

Each IP-level matrix multiply is much cheaper than MPM-level (200 bits vs 10000 bits = 50× fewer CRT passes), so absorbing 5-10 DP flushes at IP cost before one MPM flush saves significant time.

**When IP works:** The IP layer requires two conditions:
- **H matrix column spread > 370 bits** — determined by |log₁₀(α)| × (n-1) × 3.32. When H entries span a wide range, IP has room to absorb DP flushes before exhausting.
- **Baseline izd=2 = 0** — izd=2 means DP fully overflows (da/db entries exceed 2^52). When this happens, the savedp restore mechanism kicks in, and IP can't help because the DP state needs full reconstruction.

**Optimal IP value:** H_col_spread / 6 (validated across 8 problems within ±10%). This could be computed automatically at runtime from the initial H matrix.

**When IP doesn't work:** Problems with narrow H_spread (e.g., Poisson φ₂ s=30, H_spread=241 bits) or nonzero izd=2 counts (e.g., Poisson φ₂ S29) cannot use IP. The Poisson ψ₂ family (s=24, s=22, s=17, etc.) generally has wide H_spread and works well with IP.

### 4. Predicted Swap Strategy

Bailey's original PSLQ selects swap pairs based on γ^i × |H[i,i]| ranking, with no additional sorting. The **predicted_swap** strategy adds a Givens-aware bidirectional insertion sort after each DP iteration:

For each position, before committing to a swap, it:
1. **Simulates the Givens rotation** that would follow the swap, predicting the new diagonal values
2. **Checks the diagsum criterion**: only swaps if |new_d[r]| + |new_d[r+1]| < |old_d[r]| + |old_d[r+1]|
3. **Checks neighbor impact**: verifies the swap won't worsen the adjacent position's diagonal, preventing oscillation

This produces more swaps per iteration (each guaranteed to improve diagsum) and better convergence:

| Problem | Standard iterations | Predicted_swap iterations | Ratio |
|---------|-------------------|-------------------------|-------|
| psi_s24 (n=65) | 94,272 | 27,594 | 3.4× |
| S29 (n=197) | ~166,000 | ~34,000 | 4.8× |

The iteration reduction doesn't always translate directly to wall-time improvement — the strategy adds O(n) work per DP iteration for the sorting passes. The real benefit comes when fewer iterations mean fewer of a specific expensive operation:

- **Fewer mxmdm calls**: Each MPM update batch costs O(n² × ndpm). Fewer DP iterations = fewer batches.
- **Fewer fullMP triggers**: Some fullMP events are triggered by wa/wb overflow (izmm=1). Fewer DP iterations means less wa/wb growth between fullMPs, sometimes eliminating overflow-triggered fullMPs entirely.

On psi_s24, predicted_swap gives 1.1× at ndpm=3000 (mxmdm dominates equally) but compounds with ndpm tuning for 3.0× total. On S29, the 4.8× iteration reduction translates to 1.9× wall time because the mxmdm savings from fewer iterations are significant at n=197.

### 5. Threading

The tradeoffs from lower ndpm (more fullMPs) and predicted_swap (more mxmdm calls per iteration) both become less costly when those operations are threaded — which FLINT supports and MPFUN does not.

**psi_s24 thread scaling (predicted_swap + ndpm=1000 + IP=300):**

| Threads | Time | Speedup |
|---------|------|---------|
| 1 | 118s | 1.0× |
| 4 (Apple M-series) | 38s | 3.1× |

Scaling is near-linear to 4 threads on Apple Silicon. Previous AWS benchmarks showed useful scaling to 8-16 threads on AMD EPYC (64 vCPU), with diminishing returns past 16 for n=65 matrices.

## Combined Results

All optimizations stack multiplicatively:

**Poisson ψ₂ s=24 (n=65) vs MPFUN2020 baseline (532s, ndpm=3000, 1 thread):**

| Optimization | Wall | Cumulative speedup |
|-------------|------|--------------------|
| MPFUN2020 ndpm=3000, 1 thread | 532s | baseline |
| FLINT standard ndpm=3000 | 481s | 1.1× |
| + ndpm=1000 | 195s | 2.7× |
| + predicted_swap | 176s | 3.0× |
| + IP=300 | 118s | 4.5× |
| + 4 threads | 38s | **14×** |

**Poisson φ₂ s=29 (n=197), 4 threads — no QP (izd=2 prevents it):**

| Config | Wall | vs standard baseline |
|--------|------|--------------------|
| standard ndpm=1000 | 2511s (42 min) | baseline |
| predicted_swap ndpm=1000 | 1183s (20 min) | 2.1× |
| predicted_swap ndpm=800 | 1088s (18 min) | **2.3×** |
| predicted_swap ndpm=600 | 1088s (18 min) | 2.3× |

**phi_s29 cumulative speedup (n=197, no QP available):**

| Config | Time | vs standard 1-thread |
|--------|------|---------------------|
| standard ndpm=1000, 1 thread | 7387s (2.1 hr) | baseline |
| + predicted_swap | 3596s (60 min) | 2.1× |
| + 4 threads | 1183s (20 min) | 6.2× |
| + ndpm=800 | 1088s (18 min) | **6.8×** |

phi_s29 gains are smaller than psi_s24 because QP is unavailable (izd=2 > 0) and fullMP dominates at n=197 (~34s per fullMP at 1 thread). The predicted_swap iteration reduction (4.8× fewer) still translates to 2.1× wall time through fewer mxmdm calls, and threading adds another 3.0×.

## When Each Optimization Helps

| Optimization | Best when | Doesn't help when |
|-------------|-----------|-------------------|
| Lower ndpm | Small n (mxmdm dominates) | Large n (fullMP dominates) |
| Predicted_swap | Always reduces iterations | Smallest wall-time gain on easy problems |
| IP/QP | H_spread > 370 AND izd=2 = 0 | Narrow H_spread or izd=2 > 0 |
| Threading | Always helps to 8-16 threads | Diminishing returns past 16 for n < 200 |

## Problem Families Tested

| Problem | Type | n | Degree | Description |
|---------|------|---|--------|-------------|
| BBP40 | BBP | 40 | — | π via Σ 1/16^k terms |
| H2 | Algebraic | 111 | 110 | 3^(1/10) − 2^(1/11) minimal polynomial |
| psi_s24 | Poisson ψ₂ | 65 | 64 | ψ₂(2/24, 5/24) even-index lattice sum |
| phi_s29 | Poisson φ₂ | 197 | 196 | φ₂(1/29, 1/29) standard lattice sum |

The two Poisson families differ in their lattice sum structure:
- **ψ₂ (psi)**: Even-index lattice sum. α = exp(−8πs·ψ₂(p/s, q/s)). Tends to have wide H_spread → IP viable.
- **φ₂ (phi)**: Standard lattice sum. α = exp(8π·φ₂(1/s, 1/s)). Tends to have izd=2 > 0 → IP not viable.

---

## Running It

### Prerequisites

- FLINT 3.5+ (`brew install flint` on macOS) — includes GMP, MPFR
- Python 3.11+ with mpmath (`pip install mpmath`)

### Build

```bash
cd pslq/production
make
```

### Quick Test

```bash
cd test
python3 test.py --quick       # BBP40 smoke test (~1s)
python3 test.py               # Full H2 correctness matrix (~3 min)
```

### Performance Benchmarks

```bash
cd test
python3 performance.py psi_s24                  # ψ₂ s=24, 1 thread (~20 min)
python3 performance.py psi_s24 --threads 4      # ψ₂ s=24, 4 threads (~7 min)
python3 performance.py phi_s29 --threads 4      # φ₂ s=29, 4 threads (~hours)
python3 performance.py all --threads 4          # Both
```

Results are saved to `test/performance_results.json` and summarized in `test/results.md`.

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `STRATEGY` | `predicted_swap` | `standard` (no nudge) or `predicted_swap` (Givens-aware sort) |
| `THREADS` | `1` | FLINT thread count for CRT parallelism |
| `IP_BITS` | `0` (disabled) | Intermediate precision layer in bits. Set to H_spread/6 when viable |
| `IPM_OVERRIDE` | `10` | DP iterations per MPM checkpoint batch |

### Source Files

| File | Lines | Description |
|------|-------|-------------|
| `pslqm3.c` | 1769 | Main PSLQ engine — DP/MPM/fullMP cycle, strategy dispatch, entry point |
| `pslq_matmul.c` | 359 | CRT matrix multiply (mxmdm, mxm) + ColPar column-parallel threading |
| `pslq_4level.c` | 249 | Intermediate precision layer — IP init, update, flush to MPM |
| `pslq_arb.h` | 212 | MPFUN-compatible arb wrappers (round-to-nearest matching Fortran semantics) |
| `pslq_sort.c` | 171 | Quicksort for pair selection (exact Fortran translation) |

### Using from Python (ctypes)

```python
import ctypes

lib = ctypes.CDLL("pslqm3.dylib")
lib.pslqm3_c.restype = ctypes.c_int
lib.pslqm3_c.argtypes = [
    ctypes.c_int,                        # n
    ctypes.POINTER(ctypes.c_char_p),     # x_strs (decimal strings)
    ctypes.c_int,                        # digits (full precision)
    ctypes.c_int,                        # ndpm (medium precision digits)
    ctypes.c_int,                        # ndr (dynamic range requirement)
    ctypes.c_int,                        # nrb (norm bound limit)
    ctypes.c_int,                        # itm (iteration limit)
    ctypes.c_double,                     # unused
    ctypes.POINTER(ctypes.c_int64),      # result (output coefficients)
    ctypes.c_int,                        # strat (0=standard, 1=predicted_swap)
]

# Returns iq: 1 = relation found, 0 = not found
# Coefficients written to result array and /tmp/v29_rel.txt
```
