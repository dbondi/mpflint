#!/usr/bin/env python3
"""Production PSLQ correctness test.

Runs all strategies (standard + predicted_swap), QP on/off, and different
ndpm levels on H2 to validate the production build.

Usage:
  python3 test.py              # Run full H2 matrix
  python3 test.py --quick      # Quick smoke test (BBP40 only)
  python3 test.py BBP40 H2     # Run specific problems with default settings
"""

import ctypes, os, sys, time, json, subprocess

WORKDIR = os.path.dirname(os.path.abspath(__file__))
PROD_DIR = os.path.join(WORKDIR, "..")
LIB = os.path.join(PROD_DIR, "pslqm3.dylib")

sys.path.insert(0, os.path.join(PROD_DIR, ".."))
from problem_catalog import get_problem


def run_pslq(problem_name, strategy="predicted_swap", threads=1, ip_bits=0,
             ndpm_override=None, stderr_path=None):
    """Run PSLQ on a catalog problem. Returns dict with results."""
    tup = get_problem(problem_name)
    x_strs = tup[0]
    n = tup[1]
    digits = tup[2]
    ndpm = ndpm_override if ndpm_override else tup[3]
    nrb = tup[4]
    itm = tup[5]
    ndr = 100

    return _run_pslq_raw(x_strs, n, digits, ndpm, ndr, nrb, itm,
                         strategy, threads, ip_bits, stderr_path,
                         label=problem_name)


def _run_pslq_raw(x_strs, n, digits, ndpm, ndr, nrb, itm,
                  strategy="predicted_swap", threads=1, ip_bits=0,
                  stderr_path=None, label=""):
    """Run PSLQ with explicit parameters. Returns dict with results."""
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

    coeffs = [result[i] for i in range(n)]

    stats = {}
    try:
        with open("/tmp/v2_stats.json") as f:
            stats = json.load(f)
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
        "coeffs": coeffs,
        "n": n,
        "digits": digits,
        "iterations": stats.get("it", 0),
        "fullmp": stats.get("fullmp", 0),
        "nudge_swaps": stats.get("nudge_swaps", 0),
        "mxmdm_sec": stats.get("mxmdm_sec", 0),
        "updtmp_sec": stats.get("updtmpm_sec", 0) + stats.get("mxm_sec", 0),
        "stats": stats,
    }


def print_result_row(r):
    status = "OK" if r["iq"] == 1 else "FAIL"
    qp_str = f"ip={r['ip_bits']}" if r["ip_bits"] > 0 else "no-ip"
    print(f"  [{status}] {r['strategy']:16s} ndpm={r['ndpm']:<5d} {qp_str:<10s} "
          f"{r['elapsed']:7.1f}s  it={r['iterations']:<7d} fmp={r['fullmp']:<3d} "
          f"nudge={r['nudge_swaps']}")


def run_quick_smoke():
    """Quick smoke test: BBP40 with both strategies."""
    print("=== Quick Smoke Test (BBP40) ===\n")
    all_pass = True
    for strat in ["standard", "predicted_swap"]:
        r = run_pslq("BBP40", strategy=strat,
                     stderr_path=f"/tmp/pslq_test_bbp40_{strat}.log")
        print_result_row(r)
        if r["iq"] != 1:
            all_pass = False
    return all_pass


def run_h2_matrix():
    """Full H2 test matrix: all strategies × QP × ndpm levels."""
    print("=== H2 Correctness Matrix ===")
    print(f"  n=111, digits=2500\n")

    configs = [
        # (strategy, ndpm, ip_bits, description)
        ("standard",       200, 0,   "baseline"),
        ("predicted_swap", 200, 0,   "baseline"),
        ("standard",       160, 0,   "optimal ndpm"),
        ("predicted_swap", 160, 0,   "optimal ndpm"),
        ("predicted_swap", 100, 0,   "low ndpm"),
        ("predicted_swap", 300, 0,   "high ndpm"),
        ("predicted_swap", 160, 130, "optimal ndpm + QP"),
        ("standard",       160, 130, "optimal ndpm + QP"),
    ]

    all_pass = True
    results = []
    for strat, ndpm, qp, desc in configs:
        tag = f"h2_{strat}_{ndpm}_{qp}"
        r = run_pslq("H2", strategy=strat, ndpm_override=ndpm, ip_bits=qp,
                     stderr_path=f"/tmp/pslq_test_{tag}.log")
        r["desc"] = desc
        results.append(r)
        print_result_row(r)
        if r["iq"] != 1:
            all_pass = False

    return all_pass, results


if __name__ == "__main__":
    args = sys.argv[1:]

    if not os.path.exists(LIB):
        print(f"ERROR: {LIB} not found. Run 'make' in production/ first.")
        sys.exit(1)

    if "--quick" in args:
        ok = run_quick_smoke()
        print("\n" + ("ALL PASSED" if ok else "SOME FAILED"))
        sys.exit(0 if ok else 1)

    # If specific problems given, just run them with default settings
    explicit = [a for a in args if not a.startswith("--")]
    if explicit:
        print(f"=== Running: {', '.join(explicit)} ===\n")
        all_pass = True
        for pname in explicit:
            for strat in ["standard", "predicted_swap"]:
                r = run_pslq(pname, strategy=strat,
                             stderr_path=f"/tmp/pslq_test_{pname}_{strat}.log")
                print_result_row(r)
                if r["iq"] != 1:
                    all_pass = False
        print("\n" + ("ALL PASSED" if all_pass else "SOME FAILED"))
        sys.exit(0 if all_pass else 1)

    # Default: run full H2 matrix
    ok_smoke = run_quick_smoke()
    print()
    ok_h2, _ = run_h2_matrix()

    print("\n" + "=" * 70)
    all_ok = ok_smoke and ok_h2
    print("ALL PASSED" if all_ok else "SOME FAILED")
    sys.exit(0 if all_ok else 1)
