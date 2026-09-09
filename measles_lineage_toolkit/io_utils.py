"""
io_utils.py
===========
Shared, deidentified I/O for the measles lineage toolkit.

Deidentified schema (no case identifiers, no PII):
  metadata.csv : id, date, region, sublineage
  typing.csv   : id, sublineage, total_snps
  linelist.csv : id, onset_date, region
  snp_matrix.tsv : square pairwise SNP matrix, first column = id
  tree.nwk     : newick

All functions take explicit file paths (no hardcoded locations).
"""
import csv
import re


def read_metadata(path):
    """Return {id: {date, region, sublineage}}."""
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            out[r["id"]] = {
                "date": r.get("date", "") or "",
                "region": r.get("region", "") or "",
                "sublineage": r.get("sublineage", "") or "",
            }
    return out


def read_typing(path):
    """Return {id: sublineage}."""
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            out[r["id"]] = r["sublineage"]
    return out


def read_linelist(path):
    """Return list of {id, onset_date, region}."""
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append({"id": r["id"], "onset_date": r.get("onset_date", ""),
                         "region": r.get("region", "")})
    return rows


def read_snp_matrix(path):
    """Return (ids, {(a,b): dist})."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        ids = header[1:]
        dist = {}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            a = parts[0]
            for b, v in zip(ids, parts[1:]):
                dist[(a, b)] = int(v)
    return ids, dist


def pair_dist(dmat, a, b):
    if (a, b) in dmat:
        return dmat[(a, b)]
    if (b, a) in dmat:
        return dmat[(b, a)]
    return None


# --- minimal newick parser (topology + branch lengths, no deps) ---
class Node:
    __slots__ = ("children", "name", "parent")
    def __init__(self):
        self.children = []
        self.name = None
        self.parent = None


def parse_newick(text):
    s = text.strip().rstrip(";")
    pos = 0

    def rec():
        nonlocal pos
        n = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                c = rec(); c.parent = n; n.children.append(c)
                if s[pos] == ",":
                    pos += 1; continue
                if s[pos] == ")":
                    pos += 1; break
        m = re.match(r"[^,()\:]+", s[pos:])
        if m and s[pos] not in ":":
            n.name = m.group(0); pos += len(m.group(0))
        if pos < len(s) and s[pos] == ":":
            m2 = re.match(r":[0-9.eE\-]+", s[pos:])
            pos += len(m2.group(0))
        return n
    return rec()


def tip_names(node):
    if not node.children:
        return [node.name] if node.name else []
    out = []
    for c in node.children:
        out.extend(tip_names(c))
    return out


def all_nodes(node):
    yield node
    for c in node.children:
        yield from all_nodes(c)
