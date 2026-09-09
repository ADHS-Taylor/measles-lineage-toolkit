#!/usr/bin/env python3
"""
generate_beast_xmls.py
======================
Generates BEAST2 XML files for PhyDyn SIR 2-deme R0 estimation from:
  - FASTA alignment
  - Tab-delimited date file (id\tDate)
  - Newick starting tree

Appends '_I' deme suffix to all taxon identifiers (required by PhyDyn).
Uses a SIR 2-deme template XML as the base template.

Usage:
    python generate_beast_xmls.py --fasta <alignment.fasta> --dates <dates.txt> --tree <tree.nwk> [options]
    python generate_beast_xmls.py --batch <batch_config.tsv>  # Run multiple at once

Example:
    python generate_beast_xmls.py \\
        --fasta alignment.fasta \\
        --dates dates.txt \\
        --tree tree.nwk \\
        --run-name window_01 \\
        --outdir beast_runs/

Author: Taylor
"""

import argparse
import os
import re
import sys
from pathlib import Path


# ============================================================================
# Default parameter values
# ============================================================================
DEFAULTS = {
    # MCMC
    "chain_length": 1000000,
    # Clock rate
    "clock_rate": 6.61e-4,
    "clock_rate_lower": 0.0,
    "clock_rate_upper": 0.01,
    # Model timing
    "start_date": 2025.597,       # Aug 6, 2025 in decimal year
    "starting_tmrca": 2024.1,     # t0 parameter
    # SIR parameters: initI0
    "initI0": 4.0,
    "initI0_lower": 0.0,
    "initI0_upper": 100.0,
    # initS0
    "initS0": 5000.0,
    "initS0_lower": 0.0,
    # R0
    "R0": 9.9,
    "R0_lower": 0.0,
    "R0_upper": 100.0,
    # gamma (recovery rate = 1/infectious_period * 365)
    "gamma": 40.556,
    "gamma_lower": 0.0,
    "gamma_upper": 500.0,
    # eta (migration rate)
    "eta": 2.0,
    "eta_lower": 0.0,
    "eta_upper": 100.0,
    # srcSize (reservoir effective population size)
    "srcSize": 4000.0,
    "srcSize_lower": 0.0,
    "srcSize_upper": 100000.0,
}


# ============================================================================
# Helper functions
# ============================================================================

def read_fasta(fasta_path):
    """Read FASTA file, return list of (name, sequence) tuples."""
    sequences = []
    current_name = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip().replace('\r', '')
            if line.startswith('>'):
                if current_name is not None:
                    sequences.append((current_name, ''.join(current_seq)))
                current_name = line[1:].split()[0]  # Take first word only
                current_seq = []
            elif line:
                current_seq.append(line)
        if current_name is not None:
            sequences.append((current_name, ''.join(current_seq)))

    return sequences


def read_dates(dates_path):
    """Read tab-delimited date file (id\\tDate), return dict."""
    dates = {}
    with open(dates_path, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            line = line.strip().replace('\r', '')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                dates[parts[0]] = parts[1]
    return dates


def read_newick(tree_path):
    """Read newick tree from file."""
    with open(tree_path, 'r') as f:
        tree = f.read().strip().replace('\r', '')
    return tree


def add_deme_suffix(name, suffix="_I"):
    """Add deme suffix to taxon name if not already present."""
    if not name.endswith(suffix):
        return name + suffix
    return name


def suffix_newick(newick_str, suffix="_I"):
    """
    Add deme suffix to all leaf names in a newick string.
    Leaf names are sequences of non-special characters not followed by '('.
    """
    # Match leaf names: sequences of word characters (digits, letters, underscores, dots, hyphens)
    # that are NOT already suffixed with _I
    # Newick special chars: ( ) , ; : [ ]
    def replace_leaf(match):
        name = match.group(0)
        if name.endswith(suffix):
            return name
        return name + suffix

    # Pattern: match a leaf name (alphanumeric, dots, hyphens, underscores)
    # preceded by ( or , or start-of-string and followed by : or , or ) or ;
    result = re.sub(
        r'(?<=[(,])([A-Za-z0-9._-]+?)(?=[):,;])',
        replace_leaf,
        newick_str
    )
    # Also handle the very first leaf if tree starts without (
    result = re.sub(
        r'^([A-Za-z0-9._-]+?)(?=[):,;])',
        replace_leaf,
        result
    )
    return result


def format_sequence_data(sequences):
    """Format sequences as BEAST2 XML alignment entries."""
    lines = []
    for name, seq in sequences:
        taxon_name = add_deme_suffix(name)
        lines.append(
            f'        <sequence id="seq_{taxon_name}" spec="Sequence" '
            f'taxon="{taxon_name}" totalcount="4" value="{seq}"/>'
        )
    return '\n'.join(lines)


def format_dates(dates_dict, sequences):
    """
    Format dates as BEAST2 trait value string.
    Format: taxon1=date1,taxon2=date2,...
    Only includes taxa present in the alignment.
    """
    seq_names = {name for name, _ in sequences}
    entries = []
    for name in seq_names:
        taxon_name = add_deme_suffix(name)
        # Try to find date by original name
        date = dates_dict.get(name)
        if date is None:
            # Try with suffix already
            date = dates_dict.get(taxon_name)
        if date is None:
            print(f"WARNING: No date found for taxon '{name}', skipping")
            continue
        entries.append(f"{taxon_name}={date}")
    return ',\n                    '.join(entries)


def compute_alignment_length(sequences):
    """Return the alignment length from first sequence."""
    if sequences:
        return len(sequences[0][1])
    return 15897  # fallback


def generate_xml(
    fasta_path,
    dates_path,
    tree_path,
    run_name,
    outdir=".",
    params=None,
    template_path=None,
):
    """
    Generate a BEAST2 XML file from inputs.

    Parameters
    ----------
    fasta_path : str
        Path to FASTA alignment
    dates_path : str
        Path to tab-delimited date file
    tree_path : str
        Path to newick tree file
    run_name : str
        Name for this run (used in filenames and XML IDs)
    outdir : str
        Output directory for XML file
    params : dict
        Override default parameters (any key from DEFAULTS)
    template_path : str
        Path to XML template (default: auto-detect in script directory)

    Returns
    -------
    str : Path to generated XML file
    """
    # Merge params with defaults
    p = dict(DEFAULTS)
    if params:
        p.update(params)

    # Find template
    if template_path is None:
        script_dir = Path(__file__).parent
        template_path = script_dir / "SIR_2deme_template.xml"
    if not os.path.exists(template_path):
        raise FileNotFoundError(f"Template not found: {template_path}")

    # Read inputs
    print(f"Reading alignment: {fasta_path}")
    sequences = read_fasta(fasta_path)
    print(f"  {len(sequences)} sequences, {compute_alignment_length(sequences)} bp")

    print(f"Reading dates: {dates_path}")
    dates = read_dates(dates_path)
    print(f"  {len(dates)} date entries")

    print(f"Reading tree: {tree_path}")
    newick = read_newick(tree_path)

    # Add _I suffix to newick leaf names
    newick_suffixed = suffix_newick(newick)

    # Reconcile: only keep sequences that are in the tree
    tree_taxa = set(re.findall(r'(?<=[(,])([A-Za-z0-9._-]+?)(?=[):,;])', newick_suffixed))
    # Also check start of string
    start_match = re.match(r'^([A-Za-z0-9._-]+?)(?=[):,;])', newick_suffixed)
    if start_match:
        tree_taxa.add(start_match.group(1))

    original_count = len(sequences)
    sequences = [(name, seq) for name, seq in sequences if add_deme_suffix(name) in tree_taxa]
    if len(sequences) < original_count:
        dropped = original_count - len(sequences)
        print(f"  NOTE: Dropped {dropped} taxa not found in tree ({original_count} -> {len(sequences)})")

    # Compute alignment length for weightvector
    aln_length = compute_alignment_length(sequences)

    # Format data for XML
    seq_data = format_sequence_data(sequences)
    date_str = format_dates(dates, sequences)

    # Model name derived from run name
    model_name = f"SIR_{run_name}"
    alignment_name = run_name

    # Read template
    with open(template_path, 'r') as f:
        xml = f.read()

    # Perform replacements
    replacements = {
        "[INSERT ALIGNMENT NAME]": alignment_name,
        "[INSERT SEQUENCE DATA]": seq_data,
        "[INSERT MODEL NAME]": model_name,
        "[INSERT START DATE]": str(p["start_date"]),
        "[INSERT STARTING TMRCA]": str(p["starting_tmrca"]),
        "[INSERT MCMC CHAIN LENGTH]": str(p["chain_length"]),
        "[INSERT CLOCK RATE PRIOR LOWER BOUND]": str(p["clock_rate_lower"]),
        "[INSERT CLOCK RATE PRIOR UPPER BOUND]": str(p["clock_rate_upper"]),
        "[INSERT CLOCK RATE PRIOR]": str(p["clock_rate"]),
        "[INSERT DATES]": date_str,
        "[INSERT NEWICK]": newick_suffixed,
        # Parameter priors
        "[INSERT PRIOR INITI0 LOWER BOUND]": str(p["initI0_lower"]),
        "[INSERT PRIOR INITI0 UPPER BOUND]": str(p["initI0_upper"]),
        "[INSERT INITI0 PRIOR]": str(p["initI0"]),
        "[INSERT INITS0 PRIOR LOWER BOUND]": str(p["initS0_lower"]),
        "[INSERT INITS0 PRIOR]": str(p["initS0"]),
        "[INSERT R0 PRIOR LOWER BOUND]": str(p["R0_lower"]),
        "[INSERT R0 PRIOR UPPER BOUND]": str(p["R0_upper"]),
        "[INSERT R0 PRIOR]": str(p["R0"]),
        "[INSERT GAMMA PRIOR LOWER BOUND]": str(p["gamma_lower"]),
        "[INSERT GAMMA PRIOR UPPER BOUND]": str(p["gamma_upper"]),
        "[INSERT GAMMA PRIOR]": str(p["gamma"]),
        "[INSERT ETA PRIOR LOWER BOUND]": str(p["eta_lower"]),
        "[INSERT ETA PRIOR UPPER BOUND]": str(p["eta_upper"]),
        "[INSERT ETA PRIOR]": str(p["eta"]),
        "[INSERT NE PRIOR LOWER BOUND]": str(p["srcSize_lower"]),
        "[INSERT NE PRIOR UPPER BOUND]": str(p["srcSize_upper"]),
        "[INSERT NE PRIOR]": str(p["srcSize"]),
    }

    for placeholder, value in replacements.items():
        xml = xml.replace(placeholder, value)

    # Fix weightvector to match alignment length
    xml = re.sub(
        r'(<weightvector[^>]*>)\d+(</weightvector>)',
        f'\\g<1>{aln_length}\\2',
        xml
    )

    # Fix log filenames: replace hardcoded lfs2_m1 with run name
    xml = xml.replace('lfs2_m1.trees', f'{alignment_name}_m1.trees')
    xml = xml.replace('lfs2_m1.traj', f'{alignment_name}_m1.traj')

    # Remove XML comments (template comments contain "--" which is invalid XML)
    xml = re.sub(r'<!--.*?-->', '', xml, flags=re.DOTALL)

    # Write output
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, f"{alignment_name}_SIR.xml")
    with open(out_path, 'w') as f:
        f.write(xml)

    print(f"Generated: {out_path}")

    # Verify no remaining placeholders
    remaining = re.findall(r'\[INSERT [^\]]+\]', xml)
    if remaining:
        unique = set(remaining)
        print(f"WARNING: {len(unique)} unfilled placeholders remaining: {unique}")

    return out_path


# ============================================================================
# Batch mode
# ============================================================================

def run_batch(batch_file, outdir="beast_runs", params=None):
    """
    Run batch generation from a TSV config file.

    TSV format (tab-delimited, header row):
    run_name\tfasta\tdates\ttree

    Lines starting with # are comments.
    """
    generated = []
    with open(batch_file, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                print(f"Skipping malformed line: {line}")
                continue
            run_name, fasta, dates, tree = parts[0], parts[1], parts[2], parts[3]

            # Allow per-run parameter overrides in columns 5+
            run_params = dict(params) if params else {}
            if len(parts) > 4:
                # Format: key=value pairs separated by semicolons
                for kv in parts[4].split(';'):
                    if '=' in kv:
                        k, v = kv.split('=', 1)
                        try:
                            run_params[k.strip()] = float(v.strip())
                        except ValueError:
                            run_params[k.strip()] = v.strip()

            try:
                path = generate_xml(fasta, dates, tree, run_name, outdir, run_params)
                generated.append(path)
            except Exception as e:
                print(f"ERROR generating {run_name}: {e}")

    print(f"\nBatch complete: {len(generated)} XMLs generated in {outdir}/")
    return generated


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate BEAST2 XML files for PhyDyn SIR 2-deme R0 estimation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single run:
  python generate_beast_xmls.py \\
      --fasta alignment.fasta \\
      --dates dates.txt \\
      --tree tree.nwk \\
      --run-name window_01

  # Batch mode:
  python generate_beast_xmls.py --batch batch_config.tsv --outdir beast_runs/

  # Override parameters:
  python generate_beast_xmls.py \\
      --fasta alignment.fasta --dates dates.txt --tree tree.nwk \\
      --run-name test_run \\
      --chain-length 10000000 --R0 12.0 --initS0 8000
        """
    )

    # Input mode
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--batch', help='Batch config TSV file')
    mode.add_argument('--fasta', help='FASTA alignment file')

    # Single-run inputs
    parser.add_argument('--dates', help='Tab-delimited date file (id\\tDate)')
    parser.add_argument('--tree', help='Newick tree file')
    parser.add_argument('--run-name', help='Name for this run')

    # Output
    parser.add_argument('--outdir', default='beast_runs',
                        help='Output directory (default: beast_runs/)')
    parser.add_argument('--template', default=None,
                        help='Path to XML template (auto-detected if not specified)')

    # Parameter overrides
    param_group = parser.add_argument_group('Parameter overrides (defaults from Run 1)')
    param_group.add_argument('--chain-length', type=int, help=f'MCMC chain length (default: {DEFAULTS["chain_length"]})')
    param_group.add_argument('--clock-rate', type=float, help=f'Clock rate prior (default: {DEFAULTS["clock_rate"]})')
    param_group.add_argument('--start-date', type=float, help=f'Model start date in decimal year (default: {DEFAULTS["start_date"]})')
    param_group.add_argument('--starting-tmrca', type=float, help=f'Starting TMRCA / t0 (default: {DEFAULTS["starting_tmrca"]})')
    param_group.add_argument('--R0', type=float, help=f'R0 prior (default: {DEFAULTS["R0"]})')
    param_group.add_argument('--gamma', type=float, help=f'Gamma/recovery rate prior (default: {DEFAULTS["gamma"]})')
    param_group.add_argument('--eta', type=float, help=f'Eta/migration rate prior (default: {DEFAULTS["eta"]})')
    param_group.add_argument('--initI0', type=float, help=f'Initial infected count (default: {DEFAULTS["initI0"]})')
    param_group.add_argument('--initS0', type=float, help=f'Initial susceptible count (default: {DEFAULTS["initS0"]})')
    param_group.add_argument('--srcSize', type=float, help=f'Source reservoir size (default: {DEFAULTS["srcSize"]})')

    args = parser.parse_args()

    # Build parameter overrides from CLI args
    params = {}
    param_map = {
        'chain_length': args.chain_length,
        'clock_rate': args.clock_rate,
        'start_date': args.start_date,
        'starting_tmrca': args.starting_tmrca,
        'R0': args.R0,
        'gamma': args.gamma,
        'eta': args.eta,
        'initI0': args.initI0,
        'initS0': args.initS0,
        'srcSize': args.srcSize,
    }
    for k, v in param_map.items():
        if v is not None:
            params[k] = v

    if args.batch:
        run_batch(args.batch, args.outdir, params if params else None)
    else:
        # Validate single-run inputs
        if not args.dates or not args.tree:
            parser.error("--dates and --tree are required for single-run mode")
        if not args.run_name:
            # Auto-generate from fasta filename
            args.run_name = Path(args.fasta).stem.replace('.', '_').replace('-', '_')

        generate_xml(
            args.fasta,
            args.dates,
            args.tree,
            args.run_name,
            args.outdir,
            params if params else None,
            args.template,
        )


if __name__ == '__main__':
    main()
