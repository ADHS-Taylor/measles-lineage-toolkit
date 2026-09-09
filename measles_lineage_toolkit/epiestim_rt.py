#!/usr/bin/env python3
"""
epiestim_rt.py
==============
Estimate time-varying reproduction number Rt from an incidence time series using
the Cori et al. (2013) method (the algorithm behind the EpiEstim R package).

Deidentified: reads a generic line list (id, onset_date, region). Optionally
filter to one region. No case identifiers are used beyond counting.

Method (Cori et al. 2013, Am J Epidemiol):
  Rt over window [t-tau+1, t] with a Gamma prior; incidence I_t ~ Poisson(Rt*Lambda_t),
  Lambda_t = sum_s I_{t-s} * w_s, w_s = discretized serial interval.
  Posterior: Gamma(a + sum I, 1/(1/b + sum Lambda)).

Example:
    python epiestim_rt.py --linelist ../data_synthetic/linelist.csv --region R1 \
        --si-mean 11.7 --si-sd 2.0 --csv rt.csv
"""
import argparse
import csv
from datetime import datetime, timedelta

import numpy as np
from scipy import stats

import io_utils as io


def build_incidence(dates):
    parsed = []
    for d in dates:
        for fmt in ("%Y-%m-%d", "%Y-%m", "%m/%d/%Y", "%Y/%m/%d"):
            try:
                parsed.append(datetime.strptime(d, fmt)); break
            except ValueError:
                continue
    if not parsed:
        raise ValueError("No parseable dates")
    start, end = min(parsed), max(parsed)
    n = (end - start).days + 1
    counts = np.zeros(n)
    for p in parsed:
        counts[(p - start).days] += 1
    days = [start + timedelta(days=i) for i in range(n)]
    return days, counts


def discretize_si(mean, sd, max_days=40):
    shape = (mean / sd) ** 2
    scale = sd ** 2 / mean
    dist = stats.gamma(a=shape, scale=scale)
    w = np.zeros(max_days + 1)
    for s in range(1, max_days + 1):
        w[s] = dist.cdf(s + 0.5) - dist.cdf(s - 0.5)
    return w / w.sum()


def estimate_rt(incidence, w, window, prior_mean=5.0, prior_sd=5.0):
    a = (prior_mean / prior_sd) ** 2
    b = prior_sd ** 2 / prior_mean
    n = len(incidence); max_s = len(w) - 1
    lam = np.zeros(n)
    for t in range(n):
        lam[t] = sum(incidence[t - s] * w[s] for s in range(1, min(t, max_s) + 1))
    out = []
    for t in range(n):
        s0 = t - window + 1
        if s0 < 1 or lam[s0:t + 1].sum() <= 0:
            out.append(None); continue
        a_post = a + incidence[s0:t + 1].sum()
        b_post = 1.0 / (1.0 / b + lam[s0:t + 1].sum())
        post = stats.gamma(a=a_post, scale=b_post)
        out.append({"median": post.median(), "lo": post.ppf(0.025), "hi": post.ppf(0.975)})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--linelist", required=True)
    ap.add_argument("--region", help="filter to a single region (optional)")
    ap.add_argument("--window", type=int, default=7)
    ap.add_argument("--si-mean", type=float, default=11.7, help="serial interval mean (measles ~11.7)")
    ap.add_argument("--si-sd", type=float, default=2.0)
    ap.add_argument("--csv", help="write daily Rt to CSV")
    args = ap.parse_args()

    rows = io.read_linelist(args.linelist)
    if args.region:
        rows = [r for r in rows if r["region"] == args.region]
    dates = [r["onset_date"] for r in rows if r["onset_date"]]
    days, incidence = build_incidence(dates)
    w = discretize_si(args.si_mean, args.si_sd)
    results = estimate_rt(incidence, w, args.window)

    scope = f"region={args.region}" if args.region else "all regions"
    print(f"Incidence: {int(incidence.sum())} cases over {len(incidence)} days ({scope})")
    print(f"  {days[0].date()} .. {days[-1].date()}")
    early = [r["median"] for r in results[:35] if r]
    if early:
        print(f"Early-phase R (first ~4 weeks): {np.mean(early):.2f}")

    # monthly summary
    from collections import defaultdict
    md = defaultdict(list)
    for i, r in enumerate(results):
        if r:
            md[days[i].strftime("%Y-%m")].append(r["median"])
    print("\n  Month      mean Rt")
    for mo in sorted(md):
        print(f"  {mo}   {np.mean(md[mo]):.2f}")

    if args.csv:
        with open(args.csv, "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["date", "incidence", "Rt_median", "Rt_lo95", "Rt_hi95"])
            for i, r in enumerate(results):
                if r:
                    wr.writerow([days[i].date(), int(incidence[i]),
                                 f"{r['median']:.3f}", f"{r['lo']:.3f}", f"{r['hi']:.3f}"])
                else:
                    wr.writerow([days[i].date(), int(incidence[i]), "", "", ""])
        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
