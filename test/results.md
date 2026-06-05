# PSLQ Performance Results

All benchmarks: Apple M1 Max, 64 GB RAM, FLINT 3.5, `cc -O2 -march=native`.

## Correctness Tests

All passed: BBP40 (both strategies), H2 (8 configs: 2 strategies × 4 ndpm × IP on/off).

## Poisson ψ₂ s=24 — 1 thread

Problem: ψ₂(2/24, 5/24), even-index lattice sum, degree=64, n=65.

| Config | Time | Iterations | FullMP | vs baseline |
|--------|------|-----------|--------|-------------|
| standard ndpm=3000 (Bailey default) | 529s | 94,272 | 11 | baseline |
| predicted_swap ndpm=3000 | 476s | 27,594 | 7 | 1.1× |
| standard ndpm=1000 | 194s | 94,083 | 20 | 2.7× |
| predicted_swap ndpm=1000 | 175s | 27,603 | 18 | 3.0× |
| predicted_swap ndpm=700 | 178s | 27,603 | 26 | 3.0× |
| predicted_swap ndpm=1500 | 220s | 27,603 | 13 | 2.4× |
| predicted_swap ndpm=1000 +ip200 | 121s | 27,603 | 18 | 4.4× |
| standard ndpm=1000 +ip200 | 128s | 94,263 | 18 | 4.1× |
| predicted_swap ndpm=1000 +ip130 | 134s | 27,603 | 18 | 4.0× |
| **predicted_swap ndpm=1000 +ip300** | **116s** | 27,806 | 18 | **4.6×** |

## Poisson ψ₂ s=24 — 4 threads

| Config | Time | Iterations | FullMP | vs baseline |
|--------|------|-----------|--------|-------------|
| standard ndpm=3000 (Bailey default) | 91.7s | 94,272 | 11 | baseline |
| predicted_swap ndpm=3000 | 82.0s | 27,594 | 7 | 1.1× |
| standard ndpm=1000 | 50.6s | 94,083 | 20 | 1.8× |
| predicted_swap ndpm=1000 | 44.7s | 27,603 | 18 | 2.1× |
| predicted_swap ndpm=700 | 44.4s | 27,603 | 26 | 2.1× |
| predicted_swap ndpm=1500 | 51.4s | 27,603 | 13 | 1.8× |
| **predicted_swap ndpm=1000 +ip200** | **37.2s** | 27,603 | 18 | **2.5×** |
| standard ndpm=1000 +ip200 | 42.6s | 94,263 | 18 | 2.2× |
| predicted_swap ndpm=1000 +ip130 | 38.8s | 27,603 | 18 | 2.4× |
| predicted_swap ndpm=1000 +ip300 | 37.6s | 27,806 | 18 | 2.4× |

## Poisson φ₂ s=29 — 4 threads

Problem: φ₂(1/29, 1/29), standard lattice sum, degree=196, n=197.
IP not viable (izd=2 > 0 at all ndpm levels).

| Config | Time | Iterations | FullMP | vs standard ndpm=1000 |
|--------|------|-----------|--------|----------------------|
| standard ndpm=1000 | 2511s | 166,223 | 77 | baseline |
| **predicted_swap ndpm=1000** | **1183s** | **34,434** | **63** | **2.1×** |
| standard ndpm=600 | 1985s | 166,223 | 88 | 1.3× |
| predicted_swap ndpm=600 | 1088s | 34,434 | 75 | 2.3× |
| **predicted_swap ndpm=800** | **1088s** | **34,434** | **68** | **2.3×** |
| predicted_swap ndpm=1500 | 1474s | 34,434 | 58 | 1.7× |

## Cross-thread scaling

**psi_s24 best config (predicted_swap + ndpm=1000 + ip300):**

| Threads | Time | vs 1-thread baseline (529s) |
|---------|------|---------------------------|
| 1 | 116s | 4.6× |
| 4 | 37.6s | **14×** |

**phi_s29 (predicted_swap + ndpm=1000):**

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
- **4 threads**: 3.1× on psi_s24 (116s→37.6s), similar scaling on phi_s29
- All optimizations stack: **14× on psi_s24**, **6.8× on phi_s29** vs single-thread standard baseline
- phi_s29 gains are smaller because IP is unavailable and fullMP cost dominates at n=197
