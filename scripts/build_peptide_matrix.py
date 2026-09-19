#!/usr/bin/env python3
"""Reduce a Sage lfq.tsv or a MaxQuant peptides.txt to one tidy shape.

Output columns: peptide, proteins, then one column per run, named with the
`run` values from config/samples.tsv. differential.R consumes exactly this
and nothing else, which is what lets the identical statistical workflow run
on both search engines. If the two pipelines differed in their statistics as
well as their search, stage 09 would be comparing two things at once and
could not attribute a difference to either.

Sage lfq.tsv columns are, from the writer in crates/sage-cli/src/output.rs:
    peptide, charge, proteins, q_value, score, spectral_angle, <one per file>
Decoys are already excluded by Sage from this file. The per-file columns are
named with the mzML basename.

MaxQuant peptides.txt carries its quantification in "Intensity <experiment>"
columns, with Reverse and "Potential contaminant" as separate flag columns.
Those flags are folded into the accession string as rev_ and Cont_ prefixes
so that the single filter in differential.R catches both sources.
"""
import argparse
import csv
import os
import re
import sys

csv.field_size_limit(10 ** 8)


def read_samples(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build_sage(lfq_path, samples):
    with open(lfq_path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames or []
        fixed = ["peptide", "charge", "proteins", "q_value", "score", "spectral_angle"]
        missing = [c for c in fixed if c not in cols]
        if missing:
            sys.exit(f"build_peptide_matrix: lfq.tsv is missing {missing}")
        quant_cols = [c for c in cols if c not in fixed]
        # Sage names the quantification columns after the mzML files. Map
        # them back to run names by stripping the extensions.
        def to_run(c):
            b = os.path.basename(c)
            for ext in (".mzML.gz", ".mzml.gz", ".mzML", ".mzml"):
                if b.endswith(ext):
                    return b[: -len(ext)]
            return b
        colmap = {c: to_run(c) for c in quant_cols}
        known = {s["run"] for s in samples}
        unknown = sorted(set(colmap.values()) - known)
        if unknown:
            print(f"build_peptide_matrix: warning, {len(unknown)} quant column(s) "
                  f"not in the sample sheet, e.g. {unknown[:3]}", file=sys.stderr)
        rows = []
        for r in rd:
            out = {"peptide": r["peptide"], "proteins": r["proteins"]}
            for c in quant_cols:
                out[colmap[c]] = r[c]
            rows.append(out)
    runs = [s["run"] for s in samples if s["run"] in set(colmap.values())]
    return rows, runs


def annotation_from_maxquant(pep_path):
    """MaxQuant carries a "Protein names" column parallel to "Proteins"."""
    ann = {}
    with open(pep_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            accs = r.get("Proteins", "").split(";")
            names = r.get("Protein names", "").split(";")
            for i, a in enumerate(accs):
                if a and a not in ann:
                    ann[a] = names[i] if i < len(names) else (names[0] if names else "")
    return ann


def annotation_from_fasta(path):
    """UniProt headers: >db|ACC|ENTRY Description OS=... so cut at OS=."""
    ann = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            hdr = line[1:].rstrip()
            acc = hdr.split()[0]
            desc = hdr[len(acc):].strip()
            desc = re.split(r"\s+(?:OS|OX|GN|PE|SV)=", desc)[0].strip()
            ann[acc] = desc
            # Sage reports the whole header token, so also key on the bare
            # UniProt accession between the pipes.
            parts = acc.split("|")
            if len(parts) >= 2 and parts[1] not in ann:
                ann[parts[1]] = desc
    return ann


def build_maxquant(pep_path, samples):
    with open(pep_path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames or []
        seq_col = "Sequence" if "Sequence" in cols else cols[0]
        prot_col = next((c for c in ("Proteins", "Leading razor protein", "Protein group IDs")
                         if c in cols), None)
        if prot_col is None:
            sys.exit("build_peptide_matrix: no protein column in peptides.txt")
        rev_col = "Reverse" if "Reverse" in cols else None
        con_col = next((c for c in ("Potential contaminant", "Contaminant") if c in cols), None)

        # "Intensity <experiment>" maps to the sample sheet's `sample` field.
        by_sample = {s["sample"]: s["run"] for s in samples}
        colmap = {}
        for c in cols:
            if c.startswith("Intensity ") and c != "Intensity":
                exp = c[len("Intensity "):]
                if exp in by_sample:
                    colmap[c] = by_sample[exp]
        if not colmap:
            sys.exit("build_peptide_matrix: no Intensity columns matched the sample sheet")

        rows = []
        for r in rd:
            acc = r.get(prot_col, "") or "unknown"
            # Fold MaxQuant's flags into the accession so one filter handles
            # both engines. A peptide flagged Reverse is a decoy hit; a
            # peptide flagged as a potential contaminant is a contaminant.
            if rev_col and r.get(rev_col, "").strip() == "+":
                acc = ";".join("rev_" + a for a in acc.split(";") if a) or "rev_unknown"
            elif con_col and r.get(con_col, "").strip() == "+":
                acc = ";".join("Cont_" + a for a in acc.split(";") if a) or "Cont_unknown"
            out = {"peptide": r[seq_col], "proteins": acc}
            for c, run in colmap.items():
                out[run] = r[c]
            rows.append(out)
    runs = [s["run"] for s in samples if s["run"] in set(colmap.values())]
    return rows, runs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True, choices=["sage", "maxquant"])
    ap.add_argument("--input", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--annotation", help="optional protein -> description TSV")
    ap.add_argument("--fasta", help="FASTA to take descriptions from (sage source)")
    args = ap.parse_args()

    samples = read_samples(args.samples)
    if args.source == "sage":
        rows, runs = build_sage(args.input, samples)
    else:
        rows, runs = build_maxquant(args.input, samples)

    if not runs:
        sys.exit("build_peptide_matrix: no runs matched between input and sample sheet")

    header = ["peptide", "proteins"] + runs
    n_written = 0
    with open(args.output, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r.get("peptide", ""), r.get("proteins", "")] +
                       [r.get(run, "") or "0" for run in runs])
            n_written += 1
    print(f"build_peptide_matrix: {n_written} peptides x {len(runs)} runs -> {args.output}",
          file=sys.stderr)

    if args.annotation:
        if args.source == "maxquant":
            ann = annotation_from_maxquant(args.input)
        elif args.fasta:
            ann = annotation_from_fasta(args.fasta)
        else:
            ann = {}
            print("build_peptide_matrix: no --fasta given, annotation will be empty",
                  file=sys.stderr)
        with open(args.annotation, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["protein", "description"])
            for a in sorted(ann):
                w.writerow([a, ann[a]])
        print(f"build_peptide_matrix: {len(ann)} annotations -> {args.annotation}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
