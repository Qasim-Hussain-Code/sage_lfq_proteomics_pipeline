# Ribosomal proteins fall in a Francisella arginine transporter mutant, and the model you fit changes how many proteins you call

A reanalysis of PXD001584 that runs entirely from the command line: Thermo RAW
files through Sage search and label-free quantification to differential
abundance, with an entrapment experiment to check whether the reported false
discovery rate is the real one.

## Summary

Eighteen Q Exactive runs of *Francisella tularensis* subsp. *novicida* U112,
wild type against an ArgP arginine transporter deletion, went through msqrob2
with a random intercept for the biological culture. That model calls 202 of
1,180 protein groups differentially abundant at 5 percent FDR. The same
peptides, fitted with the naive model that treats all eighteen runs as
independent, call 365. The naive standard errors are a median 0.73 times the
mixed model's, which is the entire reason for the difference.

Of 52 protein groups whose description contains "ribosomal protein", 51 have a
negative fold change in the mutant. Median log2 fold change is -0.363 for
ribosomal proteins against -0.048 for the other 1,063 tested, and ribosomal
proteins make up 14.9 percent of the significant set against 4.7 percent of
the tested set (Fisher odds ratio 7.05, p = 5.7e-11). That reproduces the
prediction in the original paper.

The pipeline is built to survive a small disk. Eighteen RAW files are 47.3 GiB
and this machine had 8.3 GiB free, so `00_configure.sh` projected 14.1 GiB peak
and refused to start. What actually ran here, and what did not, is set out
under limitations. Nothing in this README is a number I did not measure.

## Background, for someone who knows RNA-seq

In RNA-seq you count reads against transcripts and the thing you sequence is
the thing you want to measure. Shotgun proteomics is not like that. Proteins
are digested with trypsin into peptides, the peptides are separated by liquid
chromatography, and the mass spectrometer measures peptides. Proteins are
inferred afterwards. Peptide, not protein, is the unit of measurement, and
that single fact drives most of what follows.

In data-dependent acquisition the instrument runs a survey scan (MS1),
picks the most intense precursors from it, and fragments them one at a time
(MS2). A search engine takes each MS2 spectrum, generates theoretical
fragments for every peptide in a database within the precursor mass
tolerance, scores them, and reports the best match. It does not identify
peptides. It ranks candidates and reports the winner, which is a different
claim and is why error control matters so much.

Error control is done by target-decoy competition: search a database of real
sequences alongside reversed or shuffled ones, and use the rate at which
decoys win as an estimate of how often a real sequence wins by chance. This
is an estimate under an assumption, that a false match is equally likely to
land on a target or a decoy. Target-decoy FDR is a model of the FDR. Whether
the model holds is an empirical question, and answering it is what the
entrapment stage of this pipeline is for.

The missingness is the other thing that will surprise an RNA-seq reader. In
this dataset 49.4 percent of peptide-by-run cells are zero. That is not
sequencing depth. A peptide is missing mostly because the instrument never
chose to fragment it, and it never chose because the peptide was not intense
enough in that run. Missingness is therefore tied to abundance: the values
that are absent are systematically the low ones. Imputing them from the
observed distribution pulls low values up and shrinks the differences you are
trying to measure, which is why this pipeline does not impute and says so in
the code where the call would have gone.

One more asymmetry. A zero in a quantification table is not a measurement of
zero abundance, it is the absence of a measurement, and the two behave
completely differently under a log transform.

## Data

| | |
|---|---|
| Accession | PXD001584 |
| Organism | *Francisella tularensis* subsp. *novicida* U112, NCBI TaxID 401614 |
| Comparison | wild type against ArgP transporter deletion, locus FTN_0848 |
| Design | 3 cultures per genotype, technical triplicate, 18 runs |
| Instrument | Q Exactive Plus, label free |
| RAW total | 47.3 GiB, mean 2.63 GiB per run |
| Search database | UniProt UP000000762, release 2026_03, downloaded 2026-09-19 |

The 18 deposited RAW files are not the whole experiment. The authors'
`experimentalDesignTemplate.txt`, inside `MaxQuantOutput.tar.gz` in the same
submission, lists 48 runs across two arginine concentrations and two
acquisition batches. Only the 20 micromolar arm was deposited as RAW: nine
wild type and nine mutant, cultures n3, n4 and n5, three injections each.
`03_fetch_raw.sh` derives `config/samples.tsv` from the PRIDE API rather than
from the paper, which is how that discrepancy surfaced.

PRIDE returns an empty checksum field for every file in this 2015 submission,
so there is no upstream digest to verify against. The pipeline checks the byte
count against the API's `fileSizeBytes` and records its own sha256. The first
file came down at exactly 3,037,044,458 bytes, matching the API.

## Pipeline

Each stage is a standalone bash script that sources `project.conf`, prints
usage with `--help`, and skips its work if it has already been done.

### 00, configure and refuse early

```bash
bash scripts/00_configure.sh --yes
```

Detects threads, available RAM and free disk, projects the total footprint,
and writes `project.conf`. On this machine, with 18 runs selected:

```
projected retained mzML.gz  9.5 GiB
projected peak disk         14.1 GiB (mzML + one RAW + 2 GiB headroom)
00_configure: REFUSING to configure a run that cannot finish.
  shortfall           5.1 GiB
```

That refusal is the point of the stage. Peak disk is the retained mzML set
plus the one RAW being converted, because Sage needs every run present
together for retention time alignment and cannot search them one at a time.

### 01, install and record

```bash
bash scripts/01_install.sh
```

Two conda environments. Versions are not taken from documentation; whatever
the channels resolve is written to `results/environment_versions.txt`. The
Bioconductor packages are pinned and `r-base` is not, because the unpinned
solve settled on R 4.3.3 and silently omitted QFeatures, msqrob2 and
MsCoreUtils, producing an environment that looked installed and was missing
everything that mattered.

### 02, search database with provenance

```bash
bash scripts/02_fetch_fasta.sh
```

1,719 sequences from UniProt proteome UP000000762 plus 381 contaminants, 2,100
in the search database, with release, date and both checksums in
`results/fasta_provenance.txt`. The contaminant set is the universal library
of Frankenfield et al. I reached for the GPM cRAP set first and it failed: on
2026-09-19 `ftp.thegpm.org` presented a certificate that does not match the
hostname. Rather than pass `--insecure` to the file that defines what counts
as a contaminant, the script falls through to the GitHub-hosted library.

### 03 and 04, stream the conversion

```bash
bash scripts/03_fetch_raw.sh --list-only
bash scripts/04_convert_mzml.sh --subset
```

Fetch one RAW, convert it, verify the mzML opens and contains spectra, delete
the RAW, move on. After the first file the real conversion ratio is known, so
the projection for the rest becomes arithmetic and the script stops rather
than dying at file fourteen.

Two ThermoRawFileParser problems cost real time here. The per-file output flag
is `--output`, not `--output_file`, which the tool answers by printing usage
and exiting non-zero. And indexed mzML combined with gzip is broken in
2.0.0.dev: it writes the `.gz` correctly, spends two minutes doing so, then
reopens the file under its uncompressed name to append the index and dies with
`FileNotFoundException`. Plain mzML skips the index pass, which Sage does not
need, and the output is half the size: 417 MB against 808 MB for the same run,
a ratio of 0.137 rather than 0.266.

### 05 and 06, search and check the FDR

```bash
bash scripts/05_run_sage.sh
bash scripts/06_entrapment.sh
```

`parallel: false` and `--batch-size 1` by default, trading speed for memory.
Sage writes `results.json` recording every parameter it resolved, which is the
machine-readable methods section, and the pipeline keeps it.

### 07 to 10, quantify, model, compare, draw

```bash
bash scripts/07_build_matrix.sh
bash scripts/08_differential.sh --sensitivity
bash scripts/09_compare_maxquant.sh
bash scripts/10_figures.sh
```

Stage 07 reduces Sage's `lfq.tsv` to a tidy matrix with one column per run.
Stage 09 reduces the authors' MaxQuant `peptides.txt` to the identical shape.
Both then go through the same `scripts/differential.R`, which is the only
reason a difference in the answer can be attributed to the search engine.

## Results

### The model matters more than most papers admit

**Fitting the culture as a random intercept calls 202 proteins; pretending the
eighteen runs are independent calls 365.**

![Significant counts under the two models](figures/maxquant_model_significance_counts.png)

Three injections of one culture are one observation of the biology. The naive
model counts them as three, so it estimates the residual variance from
injection-to-injection scatter, which is much smaller than culture-to-culture
scatter. Median standard error under the naive model is 0.73 times the mixed
model's, and median denominator degrees of freedom rise from 14.7 to 17.1.

![Standard errors, mixed against naive](figures/maxquant_se_mixed_vs_naive.png)

Everything else about those two runs is identical: same peptides, same
filters, same normalisation, same contrast. The extra 163 proteins are bought
entirely with a variance estimate that does not describe the experiment.

### Ribosomal proteins are down, as predicted

**51 of the 52 ribosomal protein groups have a negative fold change in the
mutant, and the one exception is not significant.**

![Ribosomal protein fold changes](figures/maxquant_ribosomal_logfc.png)

| | |
|---|---|
| Ribosomal groups tested | 52 of 1,115 |
| Ribosomal among the 202 significant | 30 |
| Fraction of significant set | 14.9 percent against 4.7 percent of tested |
| Fisher odds ratio | 7.05, p = 5.7e-11 |
| Median log2 FC, ribosomal | -0.363 |
| Median log2 FC, all others | -0.048 |
| Wilcoxon p | 1.0e-16 |
| Sign test on direction | 51 down, 1 up, p = 2.4e-14 |

This is a keyword check on FASTA descriptions, not a pathway enrichment, and
the distinction is in `scripts/biological_readout.R`. What makes it legitimate
rather than decorative is the background: the comparison is against the 1,115
proteins actually tested, not against all 1,719 in the proteome. Testing a
fraction of the genome and then comparing to the whole genome manufactures
enrichment from detection bias alone.

Goeminne, Gevaert and Clement reanalysed this same MaxQuant output in 2016 and
reported log2 fold change estimates for 49 ribosomal proteins, all pointing
toward down-regulation. 51 of 52 here is the same result.

### The volcano, and what the filters removed

![Volcano, mixed model](figures/maxquant_volcano_mixed.png)

| Step | Peptides |
|---|---|
| Imported | 10,693 |
| Decoys and contaminants removed | 10,490 |
| Non-minimal protein groups removed | 10,461 |
| Observed in at least 2 cultures | 7,542 |
| Protein groups after robustSummary | 1,180 |

Every one of those filters is a function of feature identity or of how often
something was observed. None looks at the difference between genotypes, which
is what keeps the p-value distribution honest.

### The two-culture threshold is arbitrary, and it barely matters

![Threshold sensitivity](figures/maxquant_threshold_sensitivity.png)

| Minimum cultures | Peptides | Proteins | Significant |
|---|---|---|---|
| 2 | 7,542 | 1,180 | 202 |
| 3 | 6,793 | 1,136 | 215 |
| 4 | 6,150 | 1,097 | 214 |
| 5 | 5,485 | 1,055 | 199 |

I picked two because it is lenient. The count moves between 199 and 215 across
the whole range, so the conclusion does not depend on the choice. Reporting
this is cheaper than defending the number.

### Normalisation

![Normalisation densities](figures/maxquant_normalisation_densities.png)

Median of ratios on the log2 scale. Per-run offsets ran from -0.427 to 0.312.
The alternative was centring each run on its own median, which is biased by
which peptides happen to be missing in that run, since a run missing many
low-abundance peptides has a higher raw median through composition alone.

### 50 proteins could not be fitted

The mixed model records `fitError` for 50 of 1,180 protein groups, the naive
model for 43. These are proteins where the design matrix is rank deficient
for that protein, usually because every surviving observation falls in one
genotype or one culture carries all the data. They are reported rather than
dropped quietly, because a protein that could not be tested is a different
thing from a protein that was tested and found unchanged. They are excluded
from the enrichment background for the same reason.

## Repository structure

```
sage_lfq_proteomics_pipeline/
├── config/
│   ├── samples.tsv            generated from the PRIDE API, never typed
│   ├── contrasts.tsv
│   ├── sage_lfq.json          the search parameters
│   └── sage_parameters.md     why each one is set that way, and which came
│                              from the authors rather than from me
├── scripts/
│   ├── 00_configure.sh ... 10_figures.sh
│   ├── lib/common.sh          timing and peak RSS logging
│   ├── build_entrapment_db.py paired shuffled entrapment construction
│   ├── assess_entrapment.py   the four FDP estimators
│   ├── build_peptide_matrix.py  Sage or MaxQuant to one tidy shape
│   ├── differential.R         the msqrob2 workflow, shared by both arms
│   ├── compare_engines.R
│   ├── biological_readout.R
│   ├── figures.R
│   └── verify_repo.sh         the pre-publication checks
├── results/                   tables, provenance, search parameters
├── figures/                   every figure in this README
├── logs/                      per-stage elapsed time and peak RSS
├── data/                      not tracked, see data/README.md
└── run_all.sh
```

## Usage

The MaxQuant arm needs no RAW files, about 1 GB of disk and roughly a minute
of compute. It reproduces every figure above.

```bash
bash run_all.sh --maxquant-only --yes
```

The full pipeline needs about 14 GiB of free disk for eighteen runs, or about
6 GiB for the six-run subset. `00_configure.sh` will refuse before anything is
downloaded if that is not available.

```bash
bash run_all.sh --yes              # eighteen runs
bash run_all.sh --subset --yes     # six runs, one injection per culture
bash run_all.sh --from 08          # re-model without re-searching
bash run_all.sh --help             # every stage and every flag
```

`--subset` keeps one technical replicate per culture. It drops the technical
replication, and with it the mixed model: with one run per culture the random
intercept has one observation per level and is not identifiable.
`differential.R` detects this and fits `~ genotype` alone, which in that design
is the correct model rather than the naive one. `results/*_design.tsv` records
which model produced the result.

Measured cost per stage is in `logs/stage_metrics.tsv`, written by
`measure_run` in `scripts/lib/common.sh` using GNU time.

## Limitations

**The full eighteen-run Sage search did not run on this machine.** It has
8.3 GiB free against a 14.1 GiB projection, and the pipeline refused, which is
the behaviour it was built for. What that means for this README is stated
directly: the differential abundance and biology above come from the authors'
own MaxQuant quantification of these spectra, put through this pipeline's
statistics. The Sage side is reported separately below and is smaller.

**The search database is not the authors' database.** They searched
`NCBI_Fnovicida.fasta`, recorded in their `parameters.txt`. This pipeline
searches UniProt UP000000762. Their peptide table reports RefSeq WP_
accessions and Sage reports UniProt ones, so the comparison needs a mapping;
1,625 of 1,719 UniProt entries carry a RefSeq cross-reference and 1,121 of the
1,180 MaxQuant protein groups map. The 59 that do not are excluded from the
comparison, and that is a real gap rather than a rounding error.

**One parameter I could not recover.** `parameters.txt` records the 20 ppm
MS/MS tolerance but not the precursor tolerance, because MaxQuant 1.4 does
precursor matching in two passes and writes neither number. I used 20 ppm
symmetric and said so in `config/sage_parameters.md`. It is a difference I did
not control.

**The obvious positive control is unavailable.** ArgP itself, FTN_0848,
UniProt A0Q671, RefSeq WP_003033738, has zero peptides in the MaxQuant table.
A deleted membrane transporter should be the cleanest possible check that the
right strains were compared, and it simply was not detected in either
genotype. So the strain identity rests on the submitters' labelling.

**A keyword check is not an enrichment analysis.** It matches the string
"ribosomal protein" in a FASTA description. It inherits whatever biases are in
how those descriptions were written, it cannot see a ribosome-associated
protein that is not described as one, and it has no notion of a pathway. The
background is correct and the test is a Fisher exact test, but the gene set is
a string match.

**One bacterial dataset says nothing general about search engines.** A 1,719
protein proteome is about the smallest realistic search space. Fragment index
behaviour, protein inference and FDR all get harder as the database grows, and
nothing here transfers to a human-sized search without being redone.

**`--subset` is a concession to disk.** It drops technical replication. It is
not a design improvement and the README says so wherever it appears.

## Data availability

- PRIDE **PXD001584**, [ftp.pride.ebi.ac.uk](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2015/01/PXD001584). File listing as retrieved is committed at `results/pride_file_listing.json`.
- UniProt proteome **UP000000762**, release **2026_03** (02-September-2026), downloaded 2026-09-19. Checksums in `results/fasta_provenance.txt`.
- The authors' MaxQuant peptide table is fetched from `statOmics/msqrob2data`, not redistributed here. Its sha256 and retrieval date are in `results/maxquant_source_provenance.txt`.
- Contaminant library from Frankenfield et al. 2022, fetched at run time.

Everything in `results/` is written by a script in `scripts/`. Nothing in it
was edited by hand.

## Citation

- Sage: Lazear MR. Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale. *J Proteome Res* 2023;22(11):3652-3659. PMID 37819886. doi:10.1021/acs.jproteome.3c00486
- ThermoRawFileParser: Hulstaert N, Shofstahl J, Sachsenberg T, Walzer M, Barsnes H, Martens L, Perez-Riverol Y. *J Proteome Res* 2020;19(1):537-542. PMID 31755270. doi:10.1021/acs.jproteome.9b00328
- Entrapment: Wen B, Freestone J, Riffle M, MacCoss MJ, Noble WS, Keich U. Assessment of false discovery rate control in tandem mass spectrometry analysis using entrapment. *Nat Methods* 2025;22(7):1454-1463. PMID 40524023. doi:10.1038/s41592-025-02719-x
- Picked-peptide FDR: Lin A, Short T, Noble WS, Keich U. Improving Peptide-Level Mass Spectrometry Analysis via Double Competition. *J Proteome Res* 2022;21(10):2412-2420. PMID 36166314. doi:10.1021/acs.jproteome.2c00282
- Dataset: Ramond E, Gesbert G, Guerrera IC, Chhuon C, Dupuis M, Rigard M, Henry T, Barel M, Charbit A. Importance of host cell arginine uptake in Francisella phagosomal escape and ribosomal protein amounts. *Mol Cell Proteomics* 2015;14(4):870-881. PMID 25616868. doi:10.1074/mcp.M114.044552
- msqrob2, method: Goeminne LJE, Gevaert K, Clement L. Peptide-level Robust Ridge Regression Improves Estimation, Sensitivity, and Specificity in Data-dependent Quantitative Label-free Shotgun Proteomics. *Mol Cell Proteomics* 2016;15(2):657-668. PMID 26566788. doi:10.1074/mcp.M115.055897
- msqrob2, summarisation: Sticker A, Goeminne L, Martens L, Clement L. Robust Summarization and Inference in Proteome-wide Label-free Quantification. *Mol Cell Proteomics* 2020;19(7):1209-1219. PMID 32321741. doi:10.1074/mcp.RA119.001624
- Contaminants: Frankenfield AM, Ni J, Ahmed M, Hao L. Protein Contaminants Matter: Building Universal Protein Contaminant Libraries for DDA and DIA Proteomics. *J Proteome Res* 2022;21(9):2104-2113. PMID 35793413. doi:10.1021/acs.jproteome.2c00145
- QFeatures and the statistical workflow are adapted from the msqrob2 vignettes and the statOmics Proteomics Data Analysis course material.

## License

MIT, in `LICENSE`. The statistical workflow in `scripts/differential.R` is
adapted from the msqrob2 documentation and the statOmics course material,
which are CC BY-SA 4.0; anyone redistributing that adaptation should honour
the ShareAlike term.
