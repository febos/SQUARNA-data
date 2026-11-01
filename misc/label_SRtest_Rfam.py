#!/usr/bin/env python3
"""Label SRtest sequences by their Rfam family hits relative to SRtrain.

Labels:
  in_train      - matches an Rfam family that also hits some SRtrain sequence
  not_in_train  - matches Rfam family/families, none of which hit SRtrain
  no_rfam       - no Rfam hit

Two modes:
  inclusion - only hits above the inclusion threshold ("!")
  any       - any reported hit ("!" or "?")

Usage:
  python label_SRtest_Rfam.py SRtestSEQ.fas SRtest.fams SRtrain.fams [-m inclusion|any|both] [-o out.tsv]
"""
import argparse
import re
from collections import Counter

LABELS = ("in_train", "not_in_train", "no_rfam")

HIT_RE = re.compile(r"^\s*\(\d+\)\s+([!?])\s+\S+\s+\S+\s+\S+\s+(\S+)")


def read_fasta_headers(path):
    """Return FASTA headers (without '>') in file order."""
    with open(path) as f:
        return [line[1:].strip() for line in f if line.startswith(">")]


def parse_cmscan(path):
    """Parse cmscan standard output.

    Returns a list of (query_header, hits) in file order, where hits is a list
    of (family, included) and query_header is "name description".
    Order is kept because query names are not unique.
    """
    queries = []
    name = desc = None
    hits = None
    in_scores = False
    with open(path) as f:
        for line in f:
            if line.startswith("Query:"):
                name = line.split()[1]
                desc = ""
                hits = []
                in_scores = False
            elif line.startswith("Description:"):
                desc = line[len("Description:"):].strip()
            elif line.startswith("Hit scores:"):
                in_scores = True
            elif line.startswith("Hit alignments:"):
                in_scores = False
            elif in_scores:
                m = HIT_RE.match(line)
                if m:
                    hits.append((m.group(2), m.group(1) == "!"))
            elif line.startswith("//"):
                queries.append((f"{name} {desc}".strip(), hits))
                name = desc = hits = None
    return queries


def families(hits, inclusion_only):
    return {fam for fam, inc in hits if inc or not inclusion_only}


def label_test(test_headers, test_res, train_res, inclusion_only):
    """Return a list of labels, one per SRtest FASTA entry (in FASTA order)."""
    if len(test_headers) != len(test_res):
        raise ValueError(f"{len(test_headers)} FASTA entries vs {len(test_res)} cmscan queries")
    for h, (q, _) in zip(test_headers, test_res):
        if " ".join(h.split()) != " ".join(q.split()):
            raise ValueError(f"FASTA/cmscan order mismatch: '{h}' vs '{q}'")

    train_fams = set()
    for _, hits in train_res:
        train_fams |= families(hits, inclusion_only)

    labels = []
    for _, hits in test_res:
        fams = families(hits, inclusion_only)
        if not fams:
            labels.append("no_rfam")
        elif fams & train_fams:
            labels.append("in_train")
        else:
            labels.append("not_in_train")
    return labels


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("test_fasta")
    ap.add_argument("test_fams")
    ap.add_argument("train_fams")
    ap.add_argument("-m", "--mode", choices=["inclusion", "any", "both"], default="both")
    ap.add_argument("-o", "--out", help="write per-sequence labels to this TSV")
    args = ap.parse_args()

    test_headers = read_fasta_headers(args.test_fasta)
    test_res = parse_cmscan(args.test_fams)
    train_res = parse_cmscan(args.train_fams)

    modes = ["inclusion", "any"] if args.mode == "both" else [args.mode]
    results = {}
    for mode in modes:
        labels = label_test(test_headers, test_res, train_res, mode == "inclusion")
        results[mode] = labels
        counts = Counter(labels)
        print(f"== {mode} hits (n={len(labels)})")
        for lab in LABELS:
            print(f"  {lab:<14}{counts[lab]}")
        print(f"labels_{mode} = {labels!r}\n")

    if args.out:
        with open(args.out, "w") as f:
            f.write("header\t" + "\t".join(f"label_{m}" for m in modes) + "\n")
            for i, h in enumerate(test_headers):
                f.write(h + "\t" + "\t".join(results[m][i] for m in modes) + "\n")


if __name__ == "__main__":
    main()
