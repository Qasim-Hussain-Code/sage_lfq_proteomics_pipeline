#!/usr/bin/env python3
"""Estimate the true FDP from an entrapment search and compare it to nominal.

All four estimators discussed by Wen B, Freestone J, Riffle M, MacCoss MJ,
Noble WS, Keich U. "Assessment of false discovery rate control in tandem mass
spectrometry analysis using entrapment." Nat Methods 2025;22(7):1454-1463,
PMID 40524023, DOI 10.1038/s41592-025-02719-x. Equation numbers below are
theirs, taken from the preprint text (bioRxiv 2024.06.01.596967v2).

  (1) combined     FDP = N_E (1 + 1/r) / (N_T + N_E)
                   Valid upper bound. Tends to overestimate, so it is
                   conservative and, in the paper's words, underpowered.
                   At r = 1 it reduces to Elias and Gygi's concatenated
                   target-decoy estimate.

  (2) lower bound  FDP = N_E / (N_T + N_E)
                   The combined formula with the 1/r term dropped. This is a
                   lower bound, so it can demonstrate that a tool FAILS to
                   control the FDR but can never demonstrate that it does.
                   Reported here only because it is the number most often
                   quoted in the literature, and seeing the gap between it
                   and (1) is the point.

  (3) sample       FDP = N_E (1/r) / N_T
                   Estimates the FDP among original target discoveries only.
                   Wen et al. call it inherently flawed: it usually
                   underestimates and can occasionally overestimate, so it
                   supports no conclusion in either direction. Reported for
                   completeness and labelled as unusable.

  (4) paired       FDP = (N_E + N_{E>=s>T} + 2 N_{E>T>=s}) / (N_T + N_E)
                   The estimator the paper proposes. A tighter upper bound
                   than (1). Requires r = 1 with each target peptide paired
                   to a unique entrapment peptide, which is exactly what
                   build_entrapment_db.py constructs.
                     N_{E>=s>T} : discovered entrapment peptides whose paired
                                  target scored below the cutoff s
                     N_{E>T>=s} : discovered entrapment peptides whose paired
                                  target scored lower but was also discovered

The headline number is (4) at the peptide level, with (1) alongside it. If
the assessed FDP exceeds nominal, that is the result and it leads the report.
"""
import argparse
import csv
import json
import re
import sys

MOD = re.compile(r"\[[^\]]*\]|\([^)]*\)|[+-]?\d+\.\d+")


def strip_mods(pep):
    """Sage writes modifications inline, e.g. M[+15.9949]PEPTIDER."""
    return MOD.sub("", pep).replace(".", "").strip().upper()


def classify(proteins, ent_prefix, decoy_tag):
    """entrapment / target / decoy for a semicolon-delimited protein list."""
    accs = [a for a in re.split(r"[;,]", proteins) if a]
    if not accs:
        return "target"
    if all(a.startswith(decoy_tag) for a in accs):
        return "decoy"
    stripped = [a[len(decoy_tag):] if a.startswith(decoy_tag) else a for a in accs]
    # A peptide is entrapment only if every protein it maps to is entrapment.
    # Anything shared with a real sequence is treated as a target, which is
    # the conservative direction: it can only lower the entrapment count.
    if all(a.startswith(ent_prefix) for a in stripped):
        return "entrapment"
    return "target"


def estimators(n_t, n_e, r, paired=None):
    out = {}
    denom = n_t + n_e
    out["combined_eq1"] = (n_e * (1 + 1 / r) / denom) if denom else float("nan")
    out["lower_bound_eq2"] = (n_e / denom) if denom else float("nan")
    out["sample_eq3"] = (n_e * (1 / r) / n_t) if n_t else float("nan")
    if paired is not None and denom:
        n_e_ge_s_gt_t, n_e_gt_t_ge_s = paired
        out["paired_eq4"] = (n_e + n_e_ge_s_gt_t + 2 * n_e_gt_t_ge_s) / denom
    else:
        out["paired_eq4"] = None
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="results.sage.tsv from the entrapment search")
    ap.add_argument("--pairing", required=True, help="target/entrapment peptide pairing tsv")
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--ratio", type=float, default=1.0, help="entrapment to target ratio r")
    ap.add_argument("--entrapment-prefix", default="ENT_")
    ap.add_argument("--decoy-tag", default="rev_")
    ap.add_argument("--thresholds", default="0.01,0.05")
    args = ap.parse_args()

    pairing = {}
    with open(args.pairing) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            pairing[row["target_peptide"]] = row["entrapment_peptide"]
    rev_pairing = {v: k for k, v in pairing.items()}

    rows = []
    with open(args.results) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames or []
        prot_col = "proteins" if "proteins" in cols else ("protein" if "protein" in cols else None)
        if prot_col is None:
            sys.exit("assess_entrapment: no protein column in results file")
        score_col = "sage_discriminant_score" if "sage_discriminant_score" in cols else "hyperscore"
        for row in rd:
            rows.append(row)
    if not rows:
        sys.exit("assess_entrapment: results file has no rows")

    print(f"assess_entrapment: {len(rows)} PSM rows, scoring on {score_col}", file=sys.stderr)

    for row in rows:
        row["_class"] = classify(row[prot_col], args.entrapment_prefix, args.decoy_tag)
        row["_pep"] = strip_mods(row.get("peptide", ""))
        try:
            row["_score"] = float(row.get(score_col, "nan"))
        except ValueError:
            row["_score"] = float("nan")

    # Best score per stripped peptide sequence, for the paired estimator and
    # for peptide level counting.
    best = {}
    for row in rows:
        if row["_class"] == "decoy":
            continue
        p = row["_pep"]
        if p not in best or row["_score"] > best[p]["_score"]:
            best[p] = row

    results = {}
    tsv_rows = []
    for thr in [float(t) for t in args.thresholds.split(",")]:
        level_out = {}

        # ---- PSM level, on spectrum_q ----
        disc = [r for r in rows if r["_class"] != "decoy" and _q(r, "spectrum_q") <= thr]
        n_t = sum(1 for r in disc if r["_class"] == "target")
        n_e = sum(1 for r in disc if r["_class"] == "entrapment")
        level_out["psm"] = dict(n_target=n_t, n_entrapment=n_e,
                                **estimators(n_t, n_e, args.ratio))

        # ---- peptide level, on peptide_q, with the paired estimator ----
        pep_disc = {p: r for p, r in best.items() if _q(r, "peptide_q") <= thr}
        n_t = sum(1 for r in pep_disc.values() if r["_class"] == "target")
        n_e = sum(1 for r in pep_disc.values() if r["_class"] == "entrapment")
        # s is the discovery cutoff: the lowest score still called at thr.
        scores_at_thr = [r["_score"] for r in pep_disc.values()]
        s = min(scores_at_thr) if scores_at_thr else float("inf")
        n_e_ge_s_gt_t = n_e_gt_t_ge_s = 0
        n_unpairable = 0
        for p, r in pep_disc.items():
            if r["_class"] != "entrapment":
                continue
            tgt = rev_pairing.get(p)
            if tgt is None:
                n_unpairable += 1
                continue
            t_row = best.get(tgt)
            t_score = t_row["_score"] if t_row else float("-inf")
            if t_score < s:
                n_e_ge_s_gt_t += 1
            elif t_score < r["_score"]:
                n_e_gt_t_ge_s += 1
        paired = (n_e_ge_s_gt_t, n_e_gt_t_ge_s) if args.ratio == 1.0 else None
        level_out["peptide"] = dict(
            n_target=n_t, n_entrapment=n_e, cutoff_score=s,
            n_E_ge_s_gt_T=n_e_ge_s_gt_t, n_E_gt_T_ge_s=n_e_gt_t_ge_s,
            n_entrapment_unpairable=n_unpairable,
            **estimators(n_t, n_e, args.ratio, paired))

        # ---- protein level, on protein_q ----
        # The paired estimator is defined for peptides. A shuffled protein is
        # not "the same protein" as its target in any sense that makes the
        # pairing argument work, so only (1) to (3) are reported here.
        prots = {}
        for r in rows:
            if r["_class"] == "decoy" or _q(r, "protein_q") > thr:
                continue
            for a in re.split(r"[;,]", r[prot_col]):
                if a:
                    prots[a] = r["_class"]
        n_t = sum(1 for c in prots.values() if c == "target")
        n_e = sum(1 for c in prots.values() if c == "entrapment")
        level_out["protein"] = dict(n_target=n_t, n_entrapment=n_e,
                                    **estimators(n_t, n_e, args.ratio))

        results[f"{thr}"] = level_out
        for lvl, d in level_out.items():
            tsv_rows.append(dict(
                nominal_fdr=thr, level=lvl,
                n_target=d["n_target"], n_entrapment=d["n_entrapment"],
                combined_eq1=_f(d["combined_eq1"]), lower_bound_eq2=_f(d["lower_bound_eq2"]),
                sample_eq3=_f(d["sample_eq3"]),
                paired_eq4=_f(d["paired_eq4"]) if d.get("paired_eq4") is not None else "NA",
                controlled=_verdict(d, thr)))

    payload = dict(entrapment_ratio_r=args.ratio, score_column=score_col,
                   n_psm_rows=len(rows), results=results)
    with open(args.out_json, "w") as fh:
        json.dump(payload, fh, indent=2)
    with open(args.out_tsv, "w") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(tsv_rows[0].keys()))
        w.writeheader()
        w.writerows(tsv_rows)

    for r in tsv_rows:
        print("assess_entrapment: " + "  ".join(f"{k}={v}" for k, v in r.items()), file=sys.stderr)


def _q(row, col):
    try:
        return float(row.get(col, "nan"))
    except ValueError:
        return float("nan")


def _f(x):
    return "NA" if x is None or x != x else f"{x:.5f}"


def _verdict(d, thr):
    """The only estimator that can establish control is a valid upper bound."""
    ub = d.get("paired_eq4")
    if ub is None:
        ub = d.get("combined_eq1")
    if ub is None or ub != ub:
        return "undetermined"
    return "consistent_with_control" if ub <= thr else "NOT_CONTROLLED"


if __name__ == "__main__":
    main()
