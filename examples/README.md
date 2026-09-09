# Examples

All commands assume you have generated the synthetic demo data first:

```bash
cd data_synthetic
python make_synthetic_data.py --out . --seed 1
cd ../measles_lineage_toolkit
```

## 1. Sub-lineage structure / import-vs-local test

```bash
python sublineage_analysis.py \
    --typing     ../data_synthetic/typing.csv \
    --metadata   ../data_synthetic/metadata.csv \
    --snp-matrix ../data_synthetic/snp_matrix.tsv \
    --tree       ../data_synthetic/tree.nwk
```

Expected (synthetic, seed 1) — illustrative structure, not real results:
- Sub-lineage composition: `SUBLIN-A` dominant, `SUBLIN-B` smaller, `OTHER-1` rare.
- Within-`SUBLIN-B` distances tighter than within-`SUBLIN-A`.
- Both `SUBLIN-A` and `SUBLIN-B` come out MONOPHYLETIC on the demo tree.
- `SUBLIN-B` concentrated in region `R1` with some `R2`/`R3` (the "export" motif).

## 2. Rt from the line list

```bash
# whole line list
python epiestim_rt.py --linelist ../data_synthetic/linelist.csv --csv rt_all.csv

# filtered to one region
python epiestim_rt.py --linelist ../data_synthetic/linelist.csv --region R1 --csv rt_R1.csv

# custom serial interval / window
python epiestim_rt.py --linelist ../data_synthetic/linelist.csv \
    --si-mean 11.7 --si-sd 2.0 --window 7
```

## 3. BEAST2 convergence check (requires a BEAST2 install + real runs)

```bash
# Not runnable on the synthetic data (no BEAST run files); shown for reference.
python beast_tools/analyze_ess.py --burnin 10 --csv ess_results.csv
```
