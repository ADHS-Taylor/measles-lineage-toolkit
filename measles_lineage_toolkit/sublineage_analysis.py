#!/usr/bin/env python3
"""
sublineage_analysis.py
======================
Tests whether genotype sub-lineages represent distinct introductions vs. one
diversifying chain, using:
  1. Within- vs between-sub-lineage pairwise SNP distances (bimodality test)
  2. Clade-based monophyly / nesting on an ML tree
  3. Region composition and timing per sub-lineage

Deidentified: works on generic id / region / sublineage schema. No case
identifiers. All inputs are explicit paths (see io_utils schema).

Example:
    python sublineage_analysis.py \
        --typing ../data_synthetic/typing.csv \
        --metadata ../data_synthetic/metadata.csv \
        --snp-matrix ../data_synthetic/snp_matrix.tsv \
        --tree ../data_synthetic/tree.nwk
"""
import argparse
import itertools
import statistics
from collections import Counter, defaultdict

import io_utils as io


def hist(values, maxbin=30):
    if not values:
        return "  (none)"
    h = Counter(values)
    mx = max(h.values())
    lines = []
    for d in range(0, min(max(values) + 1, maxbin)):
        if h[d] > 0:
            bar = "#" * int(50 * h[d] / mx)
            lines.append(f"    {d:>3} SNPs: {bar} {h[d]}")
    return "\n".join(lines)


def within_between(typing, ids, dmat):
    within, between = [], []
    within_by = defaultdict(list)
    for a, b in itertools.combinations(ids, 2):
        d = io.pair_dist(dmat, a, b)
        if d is None:
            continue
        if typing.get(a) == typing.get(b):
            within.append(d); within_by[typing[a]].append(d)
        else:
            between.append(d)
    return within, between, within_by


def monophyly(typing, tree, focus_subs):
    node_tips = {id(n): set(t for t in io.tip_names(n)) for n in io.all_nodes(tree)}

    def mrca(members):
        best, bs = None, None
        for n in io.all_nodes(tree):
            ts = node_tips[id(n)]
            if members.issubset(ts) and (bs is None or len(ts) < bs):
                best, bs = n, len(ts)
        return best

    results = {}
    all_tips = node_tips[id(tree)]
    for sub in focus_subs:
        members = set(t for t in all_tips if typing.get(t) == sub)
        if not members:
            continue
        clade = mrca(members)
        ct = node_tips[id(clade)]
        intruders = [t for t in ct if typing.get(t) != sub and t in typing]
        results[sub] = {
            "n": len(members),
            "clade_size": len(ct),
            "intruders": len(intruders),
            "intruder_subs": dict(Counter(typing.get(t) for t in intruders)),
        }
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--typing", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--snp-matrix", required=True)
    ap.add_argument("--tree", required=False)
    ap.add_argument("--min-count", type=int, default=3,
                    help="min sub-lineage size for distance/monophyly tests")
    args = ap.parse_args()

    typing = io.read_typing(args.typing)
    meta = io.read_metadata(args.metadata)
    ids, dmat = io.read_snp_matrix(args.snp_matrix)
    typed_ids = [i for i in ids if i in typing]

    print("=" * 70)
    print("1. SUB-LINEAGE COMPOSITION (in SNP matrix)")
    print("=" * 70)
    counts = Counter(typing[i] for i in typed_ids)
    for sub, n in counts.most_common():
        print(f"  {sub:<12} {n}")
    focus = [s for s, n in counts.items() if n >= args.min_count]

    print("\n" + "=" * 70)
    print("2. WITHIN vs BETWEEN SUB-LINEAGE SNP DISTANCES")
    print("=" * 70)
    within, between, within_by = within_between(typing, typed_ids, dmat)

    def summ(name, v):
        if not v:
            print(f"  {name}: no pairs"); return
        print(f"  {name}: n={len(v)} mean={statistics.mean(v):.1f} "
              f"median={statistics.median(v)} range={min(v)}-{max(v)}")
    summ("WITHIN same sub-lineage", within)
    summ("BETWEEN sub-lineages", between)
    for sub in sorted(within_by):
        summ(f"  within {sub}", within_by[sub])
    print("\n  WITHIN distribution:"); print(hist(within))
    print("\n  BETWEEN distribution:"); print(hist(between))
    print("\n  Interpretation: strong bimodality (a gap between within and between)")
    print("  indicates deep, separate introductions. Heavy overlap indicates")
    print("  shallow sub-lineages of one diversifying chain.")

    if args.tree:
        print("\n" + "=" * 70)
        print("3. MONOPHYLY / NESTING ON ML TREE")
        print("=" * 70)
        tree = io.parse_newick(open(args.tree).read())
        res = monophyly(typing, tree, focus)
        for sub, r in res.items():
            verdict = ("MONOPHYLETIC" if r["intruders"] == 0
                       else f"NOT monophyletic ({r['intruders']} intruders: {r['intruder_subs']})")
            print(f"  {sub}: {r['n']} members, smallest clade={r['clade_size']} tips -> {verdict}")

    print("\n" + "=" * 70)
    print("4. REGION COMPOSITION & TIMING PER SUB-LINEAGE")
    print("=" * 70)
    by_sub_region = defaultdict(Counter)
    by_sub_month = defaultdict(Counter)
    for i, m in meta.items():
        sub = m["sublineage"]
        by_sub_region[sub][m["region"] or "unknown"] += 1
        if m["date"] and len(m["date"]) >= 7:
            by_sub_month[sub][m["date"][:7]] += 1
    for sub in sorted(by_sub_region, key=lambda x: -sum(by_sub_region[x].values())):
        regions = ", ".join(f"{k}={v}" for k, v in by_sub_region[sub].most_common())
        print(f"\n  {sub}: regions -> {regions}")
        months = by_sub_month[sub]
        if months:
            ds = sorted(months)
            print(f"    dated span {ds[0]} .. {ds[-1]}")
            for mo in ds:
                print(f"      {mo}: {'*'*months[mo]} {months[mo]}")


if __name__ == "__main__":
    main()
