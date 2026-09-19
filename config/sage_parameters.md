# Why each Sage parameter is set the way it is

The authors' own MaxQuant settings are recoverable: `parameters.txt` inside
`MaxQuantOutput.tar.gz` in the PRIDE submission records what they actually
ran. That file, not the paper's prose, is the source for everything below
marked "from the authors". Where it is silent I chose a value and said so.
Nothing here is presented as coming from the paper when it does not.

| Parameter | Value | Source |
|---|---|---|
| enzyme | trypsin, `cleave_at KR`, `restrict P` | from the authors (Special AAs KR) |
| missed cleavages | 2 | my choice, MaxQuant's default, not stated in parameters.txt |
| min peptide length | 7 | from the authors (Min. peptide Length 7) |
| fixed modification | Carbamidomethyl C, +57.02146 | from the authors |
| variable modifications | Oxidation M +15.99491, Acetyl protein N-term +42.01057 | from the authors |
| max variable mods | 2 | my choice, Sage default |
| fragment tolerance | 20 ppm | from the authors (MS/MS tol. FTMS 20 ppm) |
| precursor tolerance | 20 ppm | **my choice.** See below |
| isotope errors | -1 to +3 | my choice, standard for Orbitrap precursors |
| PSM/peptide/protein FDR | 1 percent | from the authors (PSM FDR 0.01, Protein FDR 0.01) |
| decoys | generated internally by Sage, `decoy_tag rev_` | see below |

## The precursor tolerance is the parameter I could not recover

`parameters.txt` records the MS/MS tolerance but not the precursor tolerance,
because MaxQuant 1.4 does precursor matching in two passes (a wide first
search, then a recalibrated narrow main search) and writes neither number to
that file. So the precursor window is genuinely undetermined by the deposited
metadata.

I used 20 ppm symmetric. That is wider than MaxQuant's 6 ppm main search and
narrower than its 20 ppm first search, and Sage's linear discriminant
rescoring uses the observed precursor error as a feature, so a generous
window costs specificity less than it would in a tool that hard-filters.
It is still a difference I did not control, and it belongs in the limitations
section rather than being quietly absorbed.

## Decoys

`generate_decoys: true` and `decoy_tag: "rev_"`. Sage reverses tryptic
peptides rather than whole proteins, which is what makes the picked-peptide
approach to FDR available (Savitski et al., PMID 36166314 is the reference
Sage's own documentation gives). The distinction matters: reversing a whole
protein scrambles the peptide boundaries, so a decoy peptide no longer has a
one-to-one target partner and the picked approach, which competes each target
peptide against its own reverse, has nothing to pair. Peptide-level reversal
keeps that pairing.

When `generate_decoys` is true Sage *ignores* any FASTA entry whose accession
matches `decoy_tag`. If the database already contained sequences beginning
"rev_" they would be silently dropped. 05_run_sage.sh checks for that before
searching rather than trusting that it cannot happen.

## predict_rt and LFQ

The brief I worked from said `predict_rt` is incompatible with LFQ. The
binary says the opposite, and it is worth being precise because the error
runs the other way: Sage 0.14.6 contains the message

    `predict_rt: false` and `lfq: true` are incompatible. Setting `predict_rt: true`

So turning LFQ on requires retention time prediction, and Sage will switch it
back on for you if you set it false. The config therefore sets
`predict_rt: true` explicitly, so the file matches what actually runs instead
of recording a preference the tool overrode.

## parallel and batch size

`parallel: false` in the config and `--batch-size 1` on the command line.
Sage's documentation says parallel false "can reduce memory usage at the cost
of running slower", which is the trade this pipeline wants: the target
machine is a laptop, and an out-of-memory kill eighteen files into a search
costs more than a slow search does. The Francisella database is 2,100
sequences, so the fragment index is small and the penalty is modest. Both are
exposed as flags on 05_run_sage.sh for anyone running on a real server.
