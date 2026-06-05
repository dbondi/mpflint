#!/usr/bin/env python3
"""Production PSLQ performance benchmarks.

Tests two Poisson problems with all strategy/ndpm/QP combinations:

  poisson_psi_s24:  ψ₂(2/24, 5/24)  n=65   degree=64   — medium (minutes)
  poisson_phi_s29:  φ₂(1/29, 1/29)  n=197  degree=196  — hard (hours)

These are different Poisson summation families:
  - psi (ψ₂): even-index lattice sum, α = exp(-8πs·ψ₂(p/s, q/s))
  - phi (φ₂): standard lattice sum,   α = exp(8π·φ₂(1/s, 1/s))

Usage:
  python3 performance.py                          # Run psi_s24 (minutes)
  python3 performance.py phi_s29                   # Run phi_s29 (hours)
  python3 performance.py all                       # Run both
  python3 performance.py --threads 8               # Override thread count
  python3 performance.py --threads 8 all           # Both with 8 threads
"""

import ctypes, os, sys, time, json

WORKDIR = os.path.dirname(os.path.abspath(__file__))
PROD_DIR = os.path.join(WORKDIR, "..")
LIB = os.path.join(PROD_DIR, "pslqm3.dylib")

sys.path.insert(0, os.path.join(PROD_DIR, ".."))
from problem_catalog import get_problem, poisson_psi_alpha, gen_alpha_powers


import mpmath

DIGITS_DIR = os.path.join(PROD_DIR, "..", "digits")


def gen_poisson_psi_s24(digits=22000, ndpm=3000):
    """Poisson ψ₂(2/24, 5/24): even-index lattice sum.
    α = exp(-8π·24·ψ₂(2/24, 5/24)), minimal polynomial degree 64, n=65.
    Bailey's Fortran default: digits=22000, ndpm=3000.
    Caches full x_strs vector to avoid expensive recomputation."""
    n = 65
    vec_cache = os.path.join(DIGITS_DIR, f"poisson_psi_s24_x_n{n}_d{digits}.txt")
    alpha_cache = os.path.join(DIGITS_DIR, "poisson_psi_s24_alpha.dec")

    if os.path.exists(vec_cache):
        with open(vec_cache) as f:
            x_strs = [line.strip() for line in f if line.strip()]
        if len(x_strs) == n:
            print(f"  (loaded cached x vector from {vec_cache})")
            return x_strs, n, digits, ndpm

    if os.path.exists(alpha_cache):
        with open(alpha_cache) as f:
            alpha_str = f.read().strip()
        mpmath.mp.dps = digits + 200
        alpha = mpmath.mpf(alpha_str[:digits + 100])
        print(f"  (loaded cached alpha from {alpha_cache})")
    else:
        alpha = poisson_psi_alpha(24, digits, p=2, q=5)
        os.makedirs(DIGITS_DIR, exist_ok=True)
        with open(alpha_cache, "w") as f:
            f.write(mpmath.nstr(alpha, digits + 50))
        print(f"  (saved alpha to {alpha_cache})")

    x_strs = gen_alpha_powers(alpha, n, digits)

    os.makedirs(DIGITS_DIR, exist_ok=True)
    with open(vec_cache, "w") as f:
        for s in x_strs:
            f.write(s + "\n")
    print(f"  (saved x vector to {vec_cache})")

    return x_strs, n, digits, ndpm


def run_pslq(x_strs, n, digits, ndpm, ndr, nrb, itm,
             strategy="predicted_swap", threads=1, ip_bits=0,
             stderr_path=None, label=""):
    """Run PSLQ with explicit parameters."""
    lib = ctypes.CDLL(LIB)
    lib.pslqm3_c.restype = ctypes.c_int
    lib.pslqm3_c.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_char_p),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_int64),
        ctypes.c_int,
    ]

    x_arr = (ctypes.c_char_p * n)(*[s.encode() for s in x_strs])
    result = (ctypes.c_int64 * n)()

    os.environ["STRATEGY"] = strategy
    os.environ["THREADS"] = str(threads)
    if ip_bits > 0:
        os.environ["IP_BITS"] = str(ip_bits)
    elif "IP_BITS" in os.environ:
        del os.environ["IP_BITS"]
    if "QP_PREC" in os.environ:
        del os.environ["QP_PREC"]

    strat_int = 1 if strategy == "predicted_swap" else 0

    if stderr_path:
        old_stderr = os.dup(2)
        f = open(stderr_path, "w")
        os.dup2(f.fileno(), 2)

    t0 = time.time()
    iq = lib.pslqm3_c(n, x_arr, digits, ndpm, ndr, nrb, itm, 0.0, result, strat_int)
    elapsed = time.time() - t0

    if stderr_path:
        os.dup2(old_stderr, 2)
        os.close(old_stderr)
        f.close()

    stats = {}
    try:
        with open("/tmp/v2_stats.json") as f2:
            stats = json.load(f2)
    except:
        pass

    return {
        "label": label,
        "strategy": strategy,
        "ndpm": ndpm,
        "ip_bits": ip_bits,
        "threads": threads,
        "iq": iq,
        "elapsed": elapsed,
        "n": n,
        "digits": digits,
        "iterations": stats.get("it", 0),
        "fullmp": stats.get("fullmp", 0),
        "predicted_swaps": stats.get("predicted_swaps", 0),
        "mxmdm_sec": stats.get("mxmdm_sec", 0),
        "mxm_sec": stats.get("mxm_sec", 0),
        "updtmpm_sec": stats.get("updtmpm_sec", 0),
        "iterdp_sec": stats.get("iterdp_sec", 0),
        "izd2": stats.get("izd2", 0),
        "stats": stats,
    }


def print_header():
    print(f"{'Config':<45s} {'Time':>8s} {'It':>7s} {'FMP':>4s} {'Swaps':>8s} {'mxmdm':>7s} {'mxm':>7s} {'Status':>6s}")
    print("-" * 100)


def print_row(r):
    status = "OK" if r["iq"] == 1 else "FAIL"
    qp = f"+ip{r['ip_bits']}" if r["ip_bits"] > 0 else ""
    thr = f" t{r['threads']}" if r["threads"] > 1 else ""
    config = f"{r['strategy']:<16s} ndpm={r['ndpm']:<5d}{qp}{thr}"
    print(f"{config:<45s} {r['elapsed']:7.1f}s {r['iterations']:>7d} {r['fullmp']:>4d} "
          f"{r['predicted_swaps']:>8d} {r['mxmdm_sec']:>6.1f}s {r['mxm_sec']:>6.1f}s {status:>6s}")


def run_poisson_psi_s24_benchmarks(threads=1):
    """Poisson ψ₂ s=24 benchmark: strategies × ndpm × QP.

    ψ₂ (psi) is the even-index lattice sum variant.
    α = exp(-8π·24·ψ₂(2/24, 5/24)), degree 64, n=65.
    QP is viable here (H_spread > 370, baseline izd=2 = 0).
    """
    print("=" * 100)
    print("POISSON ψ₂ s=24 PERFORMANCE BENCHMARKS")
    print(f"  Function: ψ₂(2/24, 5/24) — even-index lattice sum")
    print(f"  α = exp(-8π·24·ψ₂), minimal polynomial degree=64, n=65")
    print(f"  Threads: {threads}")
    print()

    print("Generating poisson_psi_s24 input vector (digits=22000)...", flush=True)
    t0 = time.time()
    x_strs, n, digits, _ = gen_poisson_psi_s24(digits=22000, ndpm=3000)
    print(f"  Done in {time.time()-t0:.1f}s\n")

    configs = [
        # (strategy, ndpm, ip_bits, description)
        ("standard",       3000, 0,   "Bailey default ndpm"),
        ("predicted_swap", 3000, 0,   "Bailey default ndpm"),
        ("standard",       1000, 0,   "optimal ndpm"),
        ("predicted_swap", 1000, 0,   "optimal ndpm"),
        ("predicted_swap",  700, 0,   "low ndpm"),
        ("predicted_swap", 1500, 0,   "mid ndpm"),
        ("predicted_swap", 1000, 200, "optimal ndpm + QP 200"),
        ("standard",       1000, 200, "optimal ndpm + QP 200"),
        ("predicted_swap", 1000, 130, "optimal ndpm + QP 130"),
        ("predicted_swap", 1000, 300, "optimal ndpm + QP 300"),
    ]

    print_header()
    results = []
    all_pass = True

    for strat, ndpm, qp, desc in configs:
        tag = f"psi_s24_{strat}_{ndpm}_{qp}"
        r = run_pslq(x_strs, n, digits, ndpm,
                     ndr=100, nrb=600, itm=1000000,
                     strategy=strat, threads=threads, ip_bits=qp,
                     stderr_path=f"/tmp/pslq_perf_{tag}.log",
                     label=f"psi_s24 {desc}")
        results.append(r)
        print_row(r)
        if r["iq"] != 1:
            all_pass = False

    print()
    return all_pass, results


def run_poisson_phi_s29_benchmarks(threads=1):
    """Poisson φ₂ s=29 benchmark: strategies × ndpm.

    φ₂ (phi) is the standard lattice sum variant.
    α = exp(8π·φ₂(1/29, 1/29)), degree 196, n=197.
    QP is NOT viable (izd=2 > 0 at all ndpm levels).
    """
    print("=" * 100)
    print("POISSON φ₂ s=29 PERFORMANCE BENCHMARKS")
    print(f"  Function: φ₂(1/29, 1/29) — standard lattice sum")
    print(f"  α = exp(8π·φ₂), minimal polynomial degree=196, n=197")
    print(f"  Threads: {threads}")
    print(f"  NOTE: QP not viable (izd=2 > 0 at all ndpm levels)")
    print()

    n = 197
    digits = 20000
    default_ndpm = 1000
    nrb = 200
    itm = 1000000
    vec_cache = os.path.join(DIGITS_DIR, f"poisson_phi_s29_x_n{n}_d{digits}.txt")

    print("Generating poisson_phi_s29 input vector (digits=20000)...", flush=True)
    t0 = time.time()
    if os.path.exists(vec_cache):
        with open(vec_cache) as f:
            x_strs = [line.strip() for line in f if line.strip()]
        if len(x_strs) == n:
            print(f"  (loaded cached x vector from {vec_cache})")
        else:
            x_strs = None
    else:
        x_strs = None

    if x_strs is None:
        tup = get_problem("S29")
        x_strs = tup[0]
        os.makedirs(DIGITS_DIR, exist_ok=True)
        with open(vec_cache, "w") as f:
            for s in x_strs:
                f.write(s + "\n")
        print(f"  (saved x vector to {vec_cache})")
    print(f"  Done in {time.time()-t0:.1f}s  n={n} digits={digits} default_ndpm={default_ndpm}\n")

    configs = [
        # (strategy, ndpm, description)
        ("standard",       1000, "default ndpm"),
        ("predicted_swap", 1000, "default ndpm"),
        ("standard",        600, "low ndpm"),
        ("predicted_swap",  600, "low ndpm"),
        ("predicted_swap",  800, "mid ndpm"),
        ("predicted_swap", 1500, "high ndpm"),
    ]

    print_header()
    results = []
    all_pass = True

    for strat, ndpm, desc in configs:
        tag = f"phi_s29_{strat}_{ndpm}"
        r = run_pslq(x_strs, n, digits, ndpm,
                     ndr=100, nrb=200, itm=itm,
                     strategy=strat, threads=threads, ip_bits=0,
                     stderr_path=f"/tmp/pslq_perf_{tag}.log",
                     label=f"phi_s29 {desc}")
        results.append(r)
        print_row(r)
        if r["iq"] != 1:
            all_pass = False

    print()
    return all_pass, results


def save_results(results, path):
    """Save results to JSON for later analysis."""
    out = []
    for r in results:
        d = dict(r)
        d.pop("stats", None)
        out.append(d)
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Results saved to {path}")


if __name__ == "__main__":
    args = sys.argv[1:]

    if not os.path.exists(LIB):
        print(f"ERROR: {LIB} not found. Run 'make' in production/ first.")
        sys.exit(1)

    threads = 1
    run_psi_s24 = True
    run_phi_s29 = False

    i = 0
    while i < len(args):
        if args[i] == "--threads":
            threads = int(args[i+1])
            i += 2
        elif args[i].lower() in ("phi_s29", "s29"):
            run_psi_s24 = False
            run_phi_s29 = True
            i += 1
        elif args[i].lower() in ("psi_s24", "s24"):
            run_psi_s24 = True
            run_phi_s29 = False
            i += 1
        elif args[i].lower() == "all":
            run_psi_s24 = True
            run_phi_s29 = True
            i += 1
        else:
            i += 1

    all_results = []
    all_pass = True

    if run_psi_s24:
        ok, results = run_poisson_psi_s24_benchmarks(threads=threads)
        all_results.extend(results)
        if not ok: all_pass = False

    if run_phi_s29:
        ok, results = run_poisson_phi_s29_benchmarks(threads=threads)
        all_results.extend(results)
        if not ok: all_pass = False

    if all_results:
        save_results(all_results, os.path.join(WORKDIR, "performance_results.json"))

    print("=" * 100)
    print("ALL PASSED" if all_pass else "SOME FAILED")
    sys.exit(0 if all_pass else 1)
