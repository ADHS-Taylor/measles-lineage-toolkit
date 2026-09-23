# Measles Lineage Toolkit

Reusable, deidentified code for analyzing within-outbreak measles whole-genome
sequences: genotype sub-lineage structure, transmission-lineage / import-vs-local
questions, and case-based reproduction number (Rt) estimation.

**Author:** Taylor Martins

> **No real data is included.** This repository contains code and a fully
> **synthetic** demonstration dataset only. All identifiers, dates, regions, and
> sub-lineage labels in `data_synthetic/` are fabricated. See
> [Data availability](#data-availability).

---

## Motivation

This toolkit was built during a measles outbreak genomic investigation. It
captures two things worth reusing:

1. **A working case-based Rt estimator** (Cori et al. 2013 / EpiEstim method) and
   **sub-lineage transmission-structure analysis** that produced defensible
   results.
2. **A documented negative result:** Bayesian phylodynamic estimation of R0 from
   within-outbreak measles WGS (BEAST2 + PhyDyn coalescent) **does not converge**,
   because within-outbreak measles sequences carry insufficient phylogenetic
   signal. The BEAST tooling is included so others can reproduce the diagnosis —
   not because the approach is recommended.

---

## What's here

```
measles_lineage_toolkit/
  io_utils.py             Shared deidentified I/O + a dependency-free newick parser
  sublineage_analysis.py  Within/between sub-lineage SNP distances, tree monophyly,
                          region composition & timing (import-vs-local test)
  epiestim_rt.py          Time-varying Rt from a line list (Cori et al. 2013)
  beast_tools/
    generate_beast_xmls.py  BEAST2 PhyDyn XML generator (for the negative-result work)
    analyze_ess.py          Bulk BEAST2 convergence (ESS) analyzer via loganalyser
data_synthetic/
  make_synthetic_data.py  Generates the synthetic demo dataset below
  metadata.csv            id, date, region, sublineage
  typing.csv              id, sublineage, total_snps
  linelist.csv            id, onset_date, region
  snp_matrix.tsv          square pairwise SNP matrix
  tree.nwk                demo topology
examples/                 example command lines / expected output
```

### Deidentified data schema

The toolkit uses a **generic, non-identifying schema** — no case IDs, no
registry identifiers, no PII:

| File | Columns |
|------|---------|
| `metadata.csv` | `id`, `date`, `region`, `sublineage` |
| `typing.csv`   | `id`, `sublineage`, `total_snps` |
| `linelist.csv` | `id`, `onset_date`, `region` |
| `snp_matrix.tsv` | first column `id`, then one column per `id` (integer SNP distances) |

`id` is an arbitrary sequence label. `region` is an arbitrary grouping label.
Substitute your own values behind these column names.

---

## Quick start

```bash
pip install -r requirements.txt

# 1. Generate the synthetic demo data
cd data_synthetic
python make_synthetic_data.py --out . --seed 1

# 2. Sub-lineage structure / import-vs-local test
cd ../measles_lineage_toolkit
python sublineage_analysis.py \
    --typing     ../data_synthetic/typing.csv \
    --metadata   ../data_synthetic/metadata.csv \
    --snp-matrix ../data_synthetic/snp_matrix.tsv \
    --tree       ../data_synthetic/tree.nwk

# 3. Rt from the line list (optionally filter to one region)
python epiestim_rt.py --linelist ../data_synthetic/linelist.csv --region R1 --csv rt.csv
```

The synthetic run reproduces the analysis *structure* (a dominant sub-lineage, a
smaller one concentrated in one region and exported to others) — not the real
findings.

---

## Methods

### Sub-lineage transmission structure (`sublineage_analysis.py`)

Genotype-level typing (e.g., measles D8) is often too coarse to separate
co-circulating transmission chains. Given a finer **sub-lineage** typing layer
(for example, N450-based sub-lineage assignment), this tool tests whether
sub-lineages represent **distinct introductions** or **one diversifying chain**:

1. **Within- vs between-sub-lineage SNP distances.** Deep, independent
   introductions produce a *bimodal* distribution (a gap between within-cluster
   and between-cluster distances). Shallow sub-lineages of one chain produce
   *overlapping* distributions.
2. **Tree monophyly / nesting.** For each sub-lineage it finds the smallest clade
   (MRCA) containing all its members and counts "intruders" from other
   sub-lineages. Zero intruders = monophyletic (a coherent clade); a sub-lineage
   whose MRCA clade engulfs another is *paraphyletic background*.
3. **Region composition & timing.** Summarizes where and when each sub-lineage
   appears — useful for import-vs-local and interstate-export questions.

> **Signal limit (important for measles):** within-outbreak measles WGS differ by
> only a handful of SNPs across the ~15.9 kb genome. Genome-wide distance/tree
> summaries can **under-detect** fine structure that targeted sub-lineage typing
> (e.g., N450) resolves. Run both.

### Rt estimation (`epiestim_rt.py`)

Implements the Cori et al. (2013) sliding-window method (the algorithm behind the
EpiEstim R package) directly in Python (numpy/scipy). Default serial interval is
measles-appropriate (mean 11.7 d, SD 2.0 d). Validated against a synthetic
renewal-process epidemic with a known R.

### BEAST2 / PhyDyn tooling and the negative result (`beast_tools/`)

`generate_beast_xmls.py` builds BEAST2 XMLs for a PhyDyn SIR coalescent across
rolling time windows; `analyze_ess.py` wraps BEAST2's `loganalyser` to produce a
bulk PASS/FAIL convergence (ESS) table without a GUI.

**Finding:** across many runs and extensive runtime, R0 effective sample sizes
stayed in the single digits (target ≥ 200). This is a **data/method mismatch, not
a tuning problem** — measles evolves ~10 substitutions/genome/year, so
within-outbreak sequences (mean ~4–8 pairwise SNPs) don't contain enough
phylogenetic signal for a coalescent model to estimate population dynamics.
**Recommendation: use case-based Rt (EpiEstim) for outbreak R; reserve genomics
for transmission-lineage / import-vs-local questions, where it is informative.**

---

## Data availability

No sequence data, case data, or line lists are included in this repository, by
design. Real measles sequences are available through public repositories
(e.g., NCBI / Pathoplexus / GISAID) under their respective terms. Case-level
epidemiologic data are not public. The `data_synthetic/` dataset is fabricated
solely to demonstrate the code.

## License

Licensed under the Apache License, Version 2.0 — see [`LICENSE`](LICENSE).
Copyright 2026 Taylor Martins.

## Citation / method references

- Cori A, Ferguson NM, Fraser C, Cauchemez S. *A new framework and software to
  estimate time-varying reproduction numbers during epidemics.* Am J Epidemiol.
  2013;178(9):1505–1512.
- Volz EM, Siveroni I. *Bayesian phylodynamic inference with complex models.*
  PLoS Comput Biol. 2018 (PhyDyn).
