# PSLQ Performance Results

## Correctness Tests

All passed: BBP40 (both strategies), H2 (8 configs: 2 strategies × 4 ndpm × QP on/off).

## Poisson ψ₂ s=24 — 1 thread

Problem: ψ₂(2/24, 5/24), even-index lattice sum, degree=64, n=65.

| Config | Time | Iterations | FullMP | vs baseline |
|--------|------|-----------|--------|-------------|
| standard ndpm=3000 (Bailey default) | 532s | 94,272 | 11 | baseline |
| predicted_swap ndpm=3000 | 481s | 27,594 | 7 | 1.1× |
| standard ndpm=1000 | 195s | 94,083 | 20 | 2.7× |
| predicted_swap ndpm=1000 | 176s | 27,603 | 18 | 3.0× |
| predicted_swap ndpm=700 | 178s | 27,603 | 26 | 3.0× |
| predicted_swap ndpm=1500 | 222s | 27,603 | 13 | 2.4× |
| predicted_swap ndpm=1000 +ip200 | 122s | 27,603 | 18 | 4.4× |
| standard ndpm=1000 +ip200 | 129s | 94,263 | 18 | 4.1× |
| predicted_swap ndpm=1000 +ip130 | 136s | 27,603 | 18 | 3.9× |
| **predicted_swap ndpm=1000 +ip300** | **118s** | 27,806 | 18 | **4.5×** |

## Poisson ψ₂ s=24 — 4 threads

| Config | Time | Iterations | FullMP | vs baseline |
|--------|------|-----------|--------|-------------|
| standard ndpm=3000 (Bailey default) | 93.3s | 94,272 | 11 | baseline |
| predicted_swap ndpm=3000 | 82.8s | 27,594 | 7 | 1.1× |
| standard ndpm=1000 | 51.2s | 94,083 | 20 | 1.8× |
| predicted_swap ndpm=1000 | 46.1s | 27,603 | 18 | 2.0× |
| predicted_swap ndpm=700 | 45.5s | 27,603 | 26 | 2.1× |
| predicted_swap ndpm=1500 | 52.0s | 27,603 | 13 | 1.8× |
| predicted_swap ndpm=1000 +ip200 | 38.2s | 27,603 | 18 | 2.4× |
| standard ndpm=1000 +ip200 | 43.8s | 94,263 | 18 | 2.1× |
| predicted_swap ndpm=1000 +ip130 | 40.2s | 27,603 | 18 | 2.3× |
| **predicted_swap ndpm=1000 +ip300** | **38.1s** | 27,806 | 18 | **2.4×** |

## Poisson φ₂ s=29 — 4 threads

Problem: φ₂(1/29, 1/29), standard lattice sum, degree=196, n=197.
QP not viable (izd=2 > 0 at all ndpm levels).

| Config | Time | Iterations | FullMP | vs standard ndpm=1000 |
|--------|------|-----------|--------|----------------------|
| standard ndpm=1000 | 2511s | 166,223 | 77 | baseline |
| **predicted_swap ndpm=1000** | **1183s** | **34,434** | **63** | **2.1×** |
| standard ndpm=600 | 1985s | 166,223 | 88 | 1.3× |
| predicted_swap ndpm=600 | 1088s | 34,434 | 75 | 2.3× |
| **predicted_swap ndpm=800** | **1088s** | **34,434** | **68** | **2.3×** |
| predicted_swap ndpm=1500 | 1474s | 34,434 | 58 | 1.7× |

## Cross-thread scaling

**psi_s24 best config (predicted_swap + ndpm=1000 + qp300):**

| Threads | Time | vs MPFUN2020 1-thread |
|---------|------|--------------------|
| 1 | 118s | 4.5× |
| 4 | 38s | **14×** |

**phi_s29:**

| Config | Time | vs standard 1-thread |
|--------|------|---------------------|
| standard ndpm=1000, 1 thread | 7387s (2.1 hr) | baseline |
| predicted_swap ndpm=1000, 1 thread | 3596s (60 min) | 2.1× |
| predicted_swap ndpm=1000, 4 threads | 1183s (20 min) | 6.2× |
| predicted_swap ndpm=800, 4 threads | 1088s (18 min) | **6.8×** |

## Key findings

- **predicted_swap** does 3.4× fewer iterations on psi_s24 (27K vs 94K) and 4.8× fewer on phi_s29 (34K vs 166K)
- **ndpm tuning** (3000→1000): 2.7× on psi_s24 at 1 thread. On phi_s29, ndpm=600–800 is optimal (more fullMPs but cheaper mxmdm)
- **IP=300**: additional 1.5× on psi_s24. Not viable on phi_s29 (izd=2 > 0)
- **4 threads**: 3.1× on psi_s24 (118s→38s), similar scaling on phi_s29
- All optimizations stack: **14× on psi_s24**, **2.3× on phi_s29** vs FLINT standard baseline
- phi_s29 gains are smaller because QP is unavailable and fullMP cost dominates at n=197
