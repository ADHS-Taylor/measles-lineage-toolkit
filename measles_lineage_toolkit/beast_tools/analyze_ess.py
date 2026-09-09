#!/usr/bin/env python3
"""
analyze_ess.py
==============
Bulk convergence analysis for BEAST2 runs — quantitative, no GUI.

Wraps BEAST2's own loganalyser to get Tracer-identical ESS values, then
produces a consolidated PASS/FAIL table across all runs. Also compares R0
estimates across baseline and sensitivity runs.

Usage:
    python analyze_ess.py                       # Analyze all runs in beast_runs/
    python analyze_ess.py --burnin 10           # Burnin percent (default 10)
    python analyze_ess.py --csv results.csv     # Export CSV
    python analyze_ess.py --sensitivity         # Group + compare sensitivity runs

Requires: beast/bin/loganalyser (BEAST2 install)
"""

import argparse
import glob
import os
import re
import subprocess
import sys

LOGANALYSER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "beast", "bin", "loganalyser")
KEY_PARAMS = ["posterior", "likelihood", "R0", "gamma", "eta", "initS0", "initI0", "srcSize"]
ESS_THRESHOLD = 200


def run_loganalyser(logfile, burnin_pct):
    """Run BEAST loganalyser, return dict {param: {mean, hpd_lo, hpd_up, ess}}."""
    try:
        out = subprocess.run(
            [LOGANALYSER, "-b", str(burnin_pct), logfile],
            capture_output=True, text=True, timeout=300
        ).stdout
    except Exception as e:
        return None

    results = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) < 9:
            continue
        if parts[0] == "item":
            continue
        name = parts[0]
        try:
            results[name] = {
                "mean": float(parts[1]),
                "hpd_lo": float(parts[5]),
                "hpd_up": float(parts[6]),
                "ess": float(parts[8]) if parts[8] != "NaN" else float("nan"),
            }
        except (ValueError, IndexError):
            continue
    return results


def main():
    ap = argparse.ArgumentParser(description="Bulk BEAST2 convergence analysis")
    ap.add_argument("--dir", default="beast_runs")
    ap.add_argument("--burnin", type=int, default=10, help="Burnin percent (default 10)")
    ap.add_argument("--csv", help="Export CSV")
    ap.add_argument("--sensitivity", action="store_true", help="Group sensitivity runs for comparison")
    ap.add_argument("--pattern", default="*_m1.log")
    args = ap.parse_args()

    logs = sorted(glob.glob(os.path.join(args.dir, args.pattern)))
    if not logs:
        print(f"No logs found in {args.dir}/{args.pattern}")
        sys.exit(1)

    all_results = {}
    for log in logs:
        run = os.path.basename(log).replace("_m1.log", "")
        res = run_loganalyser(log, args.burnin)
        if res:
            all_results[run] = res

    # === Main ESS table ===
    print("=" * 130)
    print(f"CONVERGENCE ANALYSIS — ESS threshold {ESS_THRESHOLD}, burnin {args.burnin}%  (via BEAST loganalyser)")
    print("=" * 130)
    hdr = f"{'Run':<34}"
    for p in KEY_PARAMS:
        hdr += f"{p[:7]:>8}"
    hdr += f"{'R0 mean':>9}{'R0 95% HPD':>18}{'VERDICT':>9}"
    print(hdr)
    print("-" * 130)

    baseline_pass, baseline_fail = [], []
    sens_runs = {}

    for run in sorted(all_results.keys()):
        res = all_results[run]
        row = f"{run:<34}"
        min_ess = float("inf")
        for p in KEY_PARAMS:
            if p in res:
                ess = res[p]["ess"]
                row += f"{ess:>8.0f}" if ess == ess else f"{'NaN':>8}"  # NaN check
                if ess == ess and ess < min_ess:
                    min_ess = ess
            else:
                row += f"{'--':>8}"
        # R0
        if "R0" in res:
            r0 = res["R0"]
            row += f"{r0['mean']:>9.2f}"
            row += f" [{r0['hpd_lo']:>6.2f},{r0['hpd_up']:>7.2f}]"
        else:
            row += f"{'--':>9}{'--':>18}"

        verdict = "PASS" if min_ess >= ESS_THRESHOLD else "FAIL"
        row += f"{verdict:>9}"
        print(row)

        if run.startswith("sens_"):
            sens_runs[run] = res
        elif verdict == "PASS":
            baseline_pass.append(run)
        else:
            baseline_fail.append((run, min_ess))

    print("-" * 130)

    # === Baseline summary ===
    n_baseline = len([r for r in all_results if not r.startswith("sens_")])
    print(f"\nBASELINE: {len(baseline_pass)}/{n_baseline} runs PASS (all key params ESS >= {ESS_THRESHOLD})")
    if baseline_fail:
        print(f"\n{len(baseline_fail)} baseline runs need attention (sorted by lowest ESS):")
        for run, ess in sorted(baseline_fail, key=lambda x: x[1]):
            print(f"   {run:<40} lowest key-param ESS = {ess:.0f}")

    # === Sensitivity comparison ===
    if sens_runs:
        print("\n" + "=" * 90)
        print("SENSITIVITY ANALYSIS — R0 stability across parameter perturbations")
        print("=" * 90)
        # Group by alignment (Oct31 / Dec31)
        groups = {}
        for run, res in sens_runs.items():
            # sens_Oct31_initS0_2000 -> alignment=Oct31, param=initS0, val=2000
            m = re.match(r"sens_(\w+?)_(\w+?)_([\d.]+)", run)
            if m:
                aln, param, val = m.groups()
                groups.setdefault(aln, {}).setdefault(param, []).append((val, res.get("R0", {})))

        for aln in sorted(groups):
            print(f"\n  Alignment: {aln}")
            print(f"  {'Parameter':<14}{'Value':>10}{'R0 mean':>10}{'R0 95% HPD':>20}{'ESS':>8}")
            print("  " + "-" * 62)
            for param in sorted(groups[aln]):
                for val, r0 in sorted(groups[aln][param], key=lambda x: float(x[0])):
                    if r0:
                        ess = r0.get("ess", float("nan"))
                        print(f"  {param:<14}{val:>10}{r0['mean']:>10.2f}"
                              f" [{r0['hpd_lo']:>6.2f},{r0['hpd_up']:>7.2f}]{ess:>8.0f}")

        print("\n  INTERPRETATION:")
        print("  - If R0 mean stays similar across initS0 values (2000 vs 15000), R0 is data-driven (robust).")
        print("  - If R0 shifts substantially with initS0, the estimate is confounded with population size.")

    # === CSV ===
    if args.csv:
        import csv
        with open(args.csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["run"] + [f"ESS_{p}" for p in KEY_PARAMS] + ["R0_mean", "R0_hpd_lo", "R0_hpd_up", "verdict"])
            for run in sorted(all_results):
                res = all_results[run]
                min_ess = min([res[p]["ess"] for p in KEY_PARAMS if p in res and res[p]["ess"] == res[p]["ess"]] or [0])
                verdict = "PASS" if min_ess >= ESS_THRESHOLD else "FAIL"
                row = [run] + [f"{res[p]['ess']:.1f}" if p in res else "" for p in KEY_PARAMS]
                r0 = res.get("R0", {})
                row += [f"{r0.get('mean',''):.3f}" if r0 else "", f"{r0.get('hpd_lo',''):.3f}" if r0 else "",
                        f"{r0.get('hpd_up',''):.3f}" if r0 else "", verdict]
                w.writerow(row)
        print(f"\nCSV exported to {args.csv}")


if __name__ == "__main__":
    main()
