#!/usr/bin/env python3
"""Build an entrapment database with exact one-to-one peptide pairing.

Called by 06_entrapment.sh. Python rather than awk because a constrained
shuffle over an in-silico tryptic digest is the kind of thing awk makes
unreadable. Standard library only, so it runs in any of the project's
environments.

Why shuffled and not a foreign proteome
---------------------------------------
The brief this repository was built from asked for a foreign entrapment
proteome. Wen et al. 2025, the paper the brief also asks me to follow, do the
opposite: they "primarily employ the shuffled entrapment approach,
highlighting potential pitfalls of using foreign entrapment sequences in
Supplementary Note S3". More decisively, their proposed paired estimator
requires r = 1 with each target peptide paired to exactly one entrapment
peptide, "which in practice means a shuffling or reversal". A foreign
proteome cannot supply that pairing, so choosing foreign entrapment would
mean giving up the estimator the paper recommends and falling back to the
conservative combined one. Shuffled is therefore the default. Foreign is
still available through 02_fetch_fasta.sh --foreign for anyone who wants the
comparison.

How the pairing is kept exact
-----------------------------
An entrapment "protein" is the concatenation of the shuffled partners of one
target protein's tryptic peptides. For the digest of that entrapment protein
to reproduce those same peptides, the shuffle must not move any residue that
defines a cleavage site. So:

  - every K and R is held fixed, because those are the cleavage sites,
  - every P is held fixed too. This one is easy to miss: trypsin does not cut
    K or R when a proline follows, so moving a proline next to a fixed K
    destroys a cleavage site that exists in the target, and the entrapment
    protein then digests into different peptides. Holding K, R and P all in
    place makes the two digests identical by construction,
  - everything else is permuted.

Every sequence in the search database is shuffled, contaminants included.
That is deliberate: the ratio r in the estimators is entrapment database size
over original target database size, and the paired estimator needs r = 1
exactly. Leaving the 381 contaminants unpaired would make r = 1719/2100 and
quietly invalidate equation (4). A shuffled keratin peptide is not in the
tube either, so it is a legitimate entrapment sequence.

The result digests back to peptides of identical length, identical mass and
identical cleavage boundaries, differing only in the order of the interior
residues. That is what makes an entrapment hit interpretable as a false
positive rather than as an artefact of a different digest.
"""
import argparse
import json
import random
import sys
from collections import OrderedDict


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def tryptic_peptides(seq, missed_cleavages, min_len, max_len):
    """Cleave after K/R not followed by P, then expand missed cleavages.

    Returns (start, end) index pairs into seq so the caller can reconstruct
    the protein from its peptides without losing anything.
    """
    sites = [0]
    for i, aa in enumerate(seq):
        if aa in "KR" and not (i + 1 < len(seq) and seq[i + 1] == "P"):
            sites.append(i + 1)
    if sites[-1] != len(seq):
        sites.append(len(seq))
    out = []
    for i in range(len(sites) - 1):
        for j in range(i + 1, min(i + 2 + missed_cleavages, len(sites))):
            s, e = sites[i], sites[j]
            if min_len <= e - s <= max_len:
                out.append((s, e))
    return sites, out


def shuffle_peptide(pep, rng, forbidden, attempts=40):
    """Permute the sequence, holding every cleavage-defining residue in place."""
    if len(pep) <= 2:
        return None
    fixed = {len(pep) - 1}
    for i, aa in enumerate(pep):
        if aa in "KRP":
            fixed.add(i)
    movable = [i for i in range(len(pep)) if i not in fixed]
    if len(movable) < 2:
        return None
    chars = list(pep)
    pool = [chars[i] for i in movable]
    for _ in range(attempts):
        rng.shuffle(pool)
        cand = list(chars)
        for idx, i in enumerate(movable):
            cand[i] = pool[idx]
        cand = "".join(cand)
        # A shuffle that lands on a real target peptide is not an entrapment
        # sequence at all, it is a duplicate target, and counting a hit to it
        # as false would be wrong.
        if cand != pep and cand not in forbidden:
            return cand
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target-fasta", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-pairing", required=True)
    ap.add_argument("--out-stats", required=True)
    ap.add_argument("--entrapment-prefix", default="ENT_")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--missed-cleavages", type=int, default=2)
    ap.add_argument("--min-len", type=int, default=7)
    ap.add_argument("--max-len", type=int, default=50)
    ap.add_argument("--skip-prefix", default="",
                    help="accession substring to leave unpaired. Empty by "
                         "default: everything is shuffled so that r == 1")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    proteins = list(read_fasta(args.target_fasta))

    # Every fully- and partially-cleaved target peptide, so a shuffle can be
    # rejected if it collides with any of them, not just with its own partner.
    def skipped(acc):
        return bool(args.skip_prefix) and args.skip_prefix in acc

    target_peptides = set()
    for hdr, seq in proteins:
        if skipped(hdr.split()[0]):
            continue
        _, peps = tryptic_peptides(seq, args.missed_cleavages, args.min_len, args.max_len)
        for s, e in peps:
            target_peptides.add(seq[s:e])

    pairing = OrderedDict()
    out_records = []
    n_fail = 0
    n_skipped_contaminant = 0

    n_site_mismatch = 0
    for hdr, seq in proteins:
        if skipped(hdr.split()[0]):
            n_skipped_contaminant += 1
            out_records.append((hdr, seq))
            continue
        # Fully cleaved peptides tile the protein exactly, so shuffling each
        # one and concatenating rebuilds a protein of the same length.
        sites, _ = tryptic_peptides(seq, 0, 1, 10 ** 9)
        sites_full, _ = tryptic_peptides(seq, args.missed_cleavages,
                                         args.min_len, args.max_len)
        pieces = []
        for i in range(len(sites) - 1):
            frag = seq[sites[i]:sites[i + 1]]
            shuffled = shuffle_peptide(frag, rng, target_peptides)
            if shuffled is None:
                # Too short to permute, or every permutation collided. Keep
                # the original residues so protein length and mass are
                # preserved; these fragments are mostly below min_len and so
                # are not searched anyway.
                shuffled = frag
                if len(frag) >= args.min_len:
                    n_fail += 1
            pieces.append(shuffled)
        ent_seq = "".join(pieces)
        assert len(ent_seq) == len(seq), "entrapment protein changed length"

        # Verify rather than assume: the entrapment protein must expose the
        # same cleavage sites as the target, or the coordinate pairing below
        # is meaningless.
        ent_sites, ent_peps = tryptic_peptides(ent_seq, args.missed_cleavages,
                                               args.min_len, args.max_len)
        if ent_sites != sites_full:
            n_site_mismatch += 1
            continue

        for s_i, e_i in ent_peps:
            tp, ep = seq[s_i:e_i], ent_seq[s_i:e_i]
            if tp != ep:
                pairing[tp] = ep
        out_records.append((args.entrapment_prefix + hdr, ent_seq))

    with open(args.out_fasta, "w") as fh:
        for hdr, seq in out_records:
            fh.write(">" + hdr + "\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    with open(args.out_pairing, "w") as fh:
        fh.write("target_peptide\tentrapment_peptide\n")
        for t, e in pairing.items():
            fh.write(f"{t}\t{e}\n")

    stats = {
        "seed": args.seed,
        "entrapment_mode": "shuffled",
        "entrapment_ratio_r": 1,
        "entrapment_prefix": args.entrapment_prefix,
        "target_proteins_shuffled": len(out_records),
        "contaminant_proteins_skipped": n_skipped_contaminant,
        "paired_peptides": len(pairing),
        "distinct_target_peptides": len(target_peptides),
        "fully_cleaved_fragments_unshufflable": n_fail,
        "proteins_dropped_for_site_mismatch": n_site_mismatch,
        "missed_cleavages": args.missed_cleavages,
        "min_len": args.min_len,
        "max_len": args.max_len,
    }
    with open(args.out_stats, "w") as fh:
        json.dump(stats, fh, indent=2)
    for k, v in stats.items():
        print(f"build_entrapment_db: {k} = {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
