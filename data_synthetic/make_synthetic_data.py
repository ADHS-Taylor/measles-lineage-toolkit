#!/usr/bin/env python3
"""
make_synthetic_data.py
======================
Generate a small, fully SYNTHETIC dataset that demonstrates the measles
lineage toolkit end to end. No real sequences, dates, or case information are
used or reproduced. All identifiers, dates, regions, and sub-lineage labels are
fabricated for demonstration only.

The synthetic scenario mirrors the analysis structure (NOT the real data):
  - Two sub-lineages ("SUBLIN-A" dominant, "SUBLIN-B" smaller) diverged from a
    common ancestor, plus a rare unrelated genotype ("OTHER-1").
  - SUBLIN-B is seeded in region "R1", then exported to regions "R2"/"R3".
  - Outputs use the toolkit's deidentified schema:
        id, date, region, sublineage   (metadata)
        id, sublineage, total_snps      (typing)
        pairwise SNP matrix (TSV)
        a newick tree
        a line list: id, onset_date, region

Usage:
    python make_synthetic_data.py --out ../data_synthetic --n-a 40 --n-b 15 --seed 1
"""
import argparse
import csv
import os
import random
from datetime import date, timedelta


def rand_id(prefix, rng):
    return f"{prefix}{rng.randint(100000, 999999)}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="../data_synthetic")
    ap.add_argument("--n-a", type=int, default=40, help="SUBLIN-A (dominant) count")
    ap.add_argument("--n-b", type=int, default=15, help="SUBLIN-B (exported) count")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    os.makedirs(args.out, exist_ok=True)

    genome_len = 15900
    start = date(2025, 1, 1)  # arbitrary synthetic start

    # ---- build samples ----
    samples = []   # (id, sublineage, region, day_offset, snp_profile:set(positions))

    # A common ancestral profile; sub-lineages add defining mutations.
    def profile(base_positions, extra):
        s = set(base_positions)
        for _ in range(extra):
            s.add(rng.randint(1, genome_len))
        return s

    base_a = set(rng.sample(range(1, genome_len), 3))
    base_b = base_a | set(rng.sample(range(1, genome_len), 4))  # B derived from A backbone

    # SUBLIN-A: dominant, spread across regions but mostly R1, over ~8 months
    for i in range(args.n_a):
        sid = rand_id("S", rng)
        region = rng.choices(["R1", "R2", "R3", "R4"], weights=[6, 1, 1, 1])[0]
        day = rng.randint(20, 240)
        samples.append([sid, "SUBLIN-A", region, day, profile(base_a, rng.randint(0, 5))])

    # SUBLIN-B: smaller, seeded R1 early then exported R2/R3
    for i in range(args.n_b):
        sid = rand_id("S", rng)
        if i == 0:
            region, day = "R1", 15                      # earliest, basal
        elif i < 3:
            region, day = rng.choice(["R2", "R3"]), rng.randint(60, 90)  # early exports
        else:
            region, day = rng.choices(["R1", "R2", "R3"], weights=[3, 1, 1])[0], rng.randint(90, 200)
        samples.append([sid, "SUBLIN-B", region, day, profile(base_b, rng.randint(0, 3))])

    # one unrelated genotype
    other = set(rng.sample(range(1, genome_len), 60))
    samples.append([rand_id("S", rng), "OTHER-1", "R1", rng.randint(30, 120), other])

    rng.shuffle(samples)

    # ---- metadata + typing + line list ----
    meta_path = os.path.join(args.out, "metadata.csv")
    typ_path = os.path.join(args.out, "typing.csv")
    ll_path = os.path.join(args.out, "linelist.csv")
    snp_path = os.path.join(args.out, "snp_matrix.tsv")
    tree_path = os.path.join(args.out, "tree.nwk")

    with open(meta_path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["id", "date", "region", "sublineage"])
        for sid, sub, region, day, prof in samples:
            d = (start + timedelta(days=day)).isoformat()
            # leave ~20% dates blank to mimic real-world missingness
            w.writerow([sid, d if rng.random() > 0.2 else "", region, sub])

    with open(typ_path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["id", "sublineage", "total_snps"])
        for sid, sub, region, day, prof in samples:
            w.writerow([sid, sub, rng.randint(0, 2)])

    with open(ll_path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["id", "onset_date", "region"])
        # a larger synthetic line list (cases, only some sequenced)
        for sid, sub, region, day, prof in samples:
            d = (start + timedelta(days=day)).isoformat()
            w.writerow([sid, d, region])
        for _ in range(120):  # extra unsequenced cases for Rt demo
            region = rng.choices(["R1", "R2", "R3"], weights=[6, 1, 1])[0]
            day = rng.randint(1, 250)
            w.writerow([rand_id("C", rng), (start + timedelta(days=day)).isoformat(), region])

    # ---- pairwise SNP matrix (symmetric Hamming of profiles) ----
    ids = [s[0] for s in samples]
    prof = {s[0]: s[4] for s in samples}
    with open(snp_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([""] + ids)
        for a in ids:
            row = [a]
            for b in ids:
                row.append(len(prof[a] ^ prof[b]))
            w.writerow(row)

    # ---- a simple newick tree (neighbor-ish by sublineage grouping) ----
    # Not a real inference -- a demo topology grouping sublineages.
    def clade(members):
        if len(members) == 1:
            return f"{members[0]}:0.0001"
        mid = len(members) // 2
        return f"({clade(members[:mid])},{clade(members[mid:])}):0.0001"
    a_ids = [s[0] for s in samples if s[1] == "SUBLIN-A"]
    b_ids = [s[0] for s in samples if s[1] == "SUBLIN-B"]
    o_ids = [s[0] for s in samples if s[1] == "OTHER-1"]
    nwk = f"(({clade(b_ids)},{clade(a_ids)}):0.001,{clade(o_ids)}):0.0;"
    with open(tree_path, "w") as f:
        f.write(nwk + "\n")

    print(f"Wrote synthetic demo data to {args.out}/")
    for p in (meta_path, typ_path, ll_path, snp_path, tree_path):
        print(f"  {os.path.basename(p)}")
    print(f"  samples: {len(samples)} (SUBLIN-A={len(a_ids)}, SUBLIN-B={len(b_ids)}, OTHER-1={len(o_ids)})")


if __name__ == "__main__":
    main()
