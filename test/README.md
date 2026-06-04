# Test Suite

## Files

| File | Purpose |
|------|---------|
| `test.py` | **Correctness tests.** BBP40 smoke test + H2 full matrix (strategies × ndpm × QP). Validates that the production build finds correct relations. |
| `performance.py` | **Performance benchmarks.** Runs psi_s24 and phi_s29 across all strategy/ndpm/QP/thread combinations. Outputs timing tables and saves JSON results. |
| `results.md` | Latest benchmark results. |

## Running

```bash
# Correctness
python3 test.py --quick              # BBP40 only (~1s)
python3 test.py                      # BBP40 + H2 full matrix (~3 min)
python3 test.py BBP40 H2             # Specific problems

# Performance — Poisson ψ₂ s=24 (medium, minutes)
python3 performance.py psi_s24                  # 1 thread (~20 min)
python3 performance.py psi_s24 --threads 4      # 4 threads (~7 min)

# Performance — Poisson φ₂ s=29 (hard, hours)
python3 performance.py phi_s29 --threads 4      # 4 threads (~2-3 hours)

# Both
python3 performance.py all --threads 4
```

## What the tests cover

### test.py

**Quick smoke (BBP40):** Both strategies on the easiest problem. Validates basic correctness.

**H2 matrix:** 8 configs testing all dimensions:
- Strategies: standard, predicted_swap
- ndpm: 100, 160, 200, 300
- QP: off, 130 bits

All must find the correct degree-110 minimal polynomial of 3^(1/10) − 2^(1/11).

### performance.py

**psi_s24** — Poisson ψ₂(2/24, 5/24), n=65. Tests 10 configs:
- Strategies: standard, predicted_swap
- ndpm: 700, 1000, 1500, 3000
- QP: 0 (off), 130, 200, 300 bits
- IP/QP is viable here (H_spread > 370, izd=2 = 0)

**phi_s29** — Poisson φ₂(1/29, 1/29), n=197. Tests 6 configs:
- Strategies: standard, predicted_swap
- ndpm: 600, 800, 1000, 1500
- No QP (izd=2 > 0 prevents it)

## Cached input vectors

The first run of each problem precomputes the input vector (alpha powers at high precision) and caches it to `digits/`. Subsequent runs load instantly. Cache files:

- `digits/poisson_psi_s24_alpha.dec` — 22000-digit α constant
- `digits/poisson_psi_s24_x_n65_d22000.txt` — Full x vector (65 entries)
- `digits/poisson_phi_s29_x_n197_d20000.txt` — Full x vector (197 entries)

Delete these to force regeneration.

## Output

- **Console:** Timing table with all configs
- **`performance_results.json`:** Machine-readable results for further analysis
- **`/tmp/pslq_perf_*.log`:** Per-config stderr logs with full PSLQ diagnostics
- **`/tmp/v2_stats.json`:** Structured stats from the last run (iterations, fullMP counts, timing breakdowns)
