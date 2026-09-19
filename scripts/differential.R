#!/usr/bin/env Rscript
# Peptide level differential abundance with msqrob2.
#
# One script, two inputs. It runs on the Sage LFQ table and on the authors'
# MaxQuant peptides.txt without changing a line, because 07_build_matrix.sh
# and 09_compare_maxquant.sh both reduce their source to the same tidy shape
# first. That is the only way the comparison in stage 09 means anything: if
# the two pipelines differed in their statistics as well as their search, a
# difference in the answer would tell you nothing about the search engines.
#
# Adapted from the msqrob2 vignettes and the statOmics Proteomics Data
# Analysis course (CC BY-SA 4.0). What changed, and why:
#
#   - The book reads a MaxQuant peptides.txt and locates its quantification
#     columns by the "Intensity " prefix. Sage writes a different table:
#     columns peptide, charge, proteins, q_value, score, spectral_angle, then
#     one column per mzML file. So the import is driven by an explicit sample
#     sheet instead of a column name prefix, and both sources are mapped onto
#     it upstream.
#   - The book's contaminant filter looks for MaxQuant's "REV__" and "CON__"
#     flags. Sage has neither concept, so filtering is by accession prefix
#     against the tags this pipeline controls: rev_ for Sage decoys, Cont_
#     for the contaminant library, ENT_ for entrapment sequences.
#   - msqrobCollect does not exist in msqrob2 1.18.0. Results are collected
#     with topFeatures, which is the current API.
#   - Both the mixed model and the naive model are fitted deliberately, and
#     both are reported. The book fits the correct one.

suppressPackageStartupMessages({
  library(QFeatures); library(msqrob2); library(MsCoreUtils)
  library(SummarizedExperiment); library(BiocParallel); library(limma)
})

# Everything downstream of here is deterministic, but set the seed anyway so
# that any future change which introduces sampling does not silently become
# irreproducible.
set.seed(20260919)
register(SerialParam())

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  args[i + 1]
}
if ("--help" %in% args || length(args) == 0) {
  cat("Usage: differential.R --matrix <tsv> --samples <tsv> --outdir <dir> --prefix <name>\n")
  cat("                      [--min-cultures N] [--contrast <coef>] [--sensitivity]\n")
  quit(status = 0)
}

matrix_file <- parse_arg("--matrix")
samples_file <- parse_arg("--samples")
outdir <- parse_arg("--outdir", "results")
prefix <- parse_arg("--prefix", "sage")
min_cultures <- as.integer(parse_arg("--min-cultures", "2"))
run_sensitivity <- "--sensitivity" %in% args
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(sprintf("differential.R: %s\n", paste0(...)))

# ---------------------------------------------------------------- import ----
samples <- read.delim(samples_file, stringsAsFactors = FALSE)
pep <- read.delim(matrix_file, stringsAsFactors = FALSE, check.names = FALSE)

quant_cols <- intersect(samples$run, colnames(pep))
if (length(quant_cols) == 0) stop("no sample columns from the sample sheet are present in the matrix")
samples <- samples[match(quant_cols, samples$run), , drop = FALSE]
msg(nrow(pep), " peptide rows, ", length(quant_cols), " runs")

coldata <- data.frame(
  row.names = quant_cols,
  run       = samples$run,
  sample    = samples$sample,
  genotype  = factor(samples$genotype, levels = c("WT", "argP_KO")),
  biorep    = factor(samples$biorep),
  techrep   = factor(samples$techrep),
  stringsAsFactors = FALSE
)
# Culture identity has to be unique across genotypes. The sheet labels
# cultures n3/n4/n5 within each genotype, so WT n3 and mutant n3 would
# otherwise collapse into one random effect level and the model would be
# pooling two unrelated cultures.
coldata$culture <- factor(paste(coldata$genotype, coldata$biorep, sep = "_"))

# The SummarizedExperiment is assembled by hand rather than through
# readQFeatures. QFeatures 1.20.0 changed that function's contract: passing a
# colData now requires a "quantCols" column inside the colData itself, and
# the older positional form fails with a message that does not say so.
# Constructing the object directly is a few more lines and does not move
# between releases.
int_mat <- as.matrix(pep[, quant_cols, drop = FALSE])
mode(int_mat) <- "numeric"
rownames(int_mat) <- make.unique(as.character(pep$peptide))
se <- SummarizedExperiment(
  assays  = list(assay = int_mat),
  rowData = DataFrame(peptide_id = pep$peptide, proteins = pep$proteins),
  colData = DataFrame(coldata))
qf <- QFeatures(list(peptide = se), colData = DataFrame(coldata))

waterfall <- data.frame(step = "imported", features = nrow(qf[["peptide"]]))
add_step <- function(w, label, n) rbind(w, data.frame(step = label, features = n))

# ------------------------------------------------------- zeros are not 0 ----
# A zero in an LFQ table is not a measurement of zero abundance. It means no
# MS1 peak was integrated for that precursor in that run, which happens when
# the peptide was not selected for fragmentation, when its peak failed the
# spectral angle cutoff, or when it genuinely was not there. Those are
# different events and none of them is "the abundance was zero". Left as 0 and
# log transformed they become -Inf; replaced by a small number they become a
# confident measurement of almost nothing, which is worse. NA says what is
# actually known, and msqrob2's models drop missing observations per protein
# rather than pretending.
zero_n <- sum(assay(qf[["peptide"]]) == 0, na.rm = TRUE)
qf <- zeroIsNA(qf, "peptide")
msg(zero_n, " zero intensities recoded as NA (",
    round(100 * zero_n / length(assay(qf[["peptide"]])), 1), " percent of cells)")

qf <- logTransform(qf, base = 2, i = "peptide", name = "log2")
waterfall <- add_step(waterfall, "log2 transformed", nrow(qf[["log2"]]))

# --------------------------------------------------------------- filters ----
# Every filter below is a function of the identity of a feature or of how
# often it was observed. None of them looks at the difference between
# genotypes. That independence is the whole point: a filter that used the
# observed effect size, or that kept peptides because they looked different
# between groups, would enrich the surviving set for apparent differences and
# the resulting p values would no longer be uniform under the null. Filtering
# on annotation and on total observation count is safe because it is blind to
# the contrast being tested.
rd <- rowData(qf[["log2"]])
acc <- rd$proteins

is_decoy  <- grepl("(^|;)rev_",  acc)
is_contam <- grepl("(^|;)(Cont_|CON__|sp\\|Cont_)", acc) | grepl("Cont_", acc)
is_entrap <- grepl("(^|;)ENT_",  acc)
keep <- !(is_decoy | is_contam | is_entrap)
msg(sum(is_decoy), " decoy, ", sum(is_contam), " contaminant, ", sum(is_entrap),
    " entrapment peptides removed")
qf <- qf[keep, , ]
waterfall <- add_step(waterfall, "decoys/contaminants removed", nrow(qf[["log2"]]))

# Razor and shared peptides. A peptide matching more than one protein group
# cannot be attributed to one of them, and summarising it into every group it
# touches propagates one measurement into several supposedly independent
# protein estimates. smallestUniqueGroups keeps the smallest set of protein
# groups that explains the peptides.
prot <- rowData(qf[["log2"]])$proteins
keep_unique <- prot %in% smallestUniqueGroups(prot)
qf <- qf[keep_unique, , ]
msg(sum(!keep_unique), " peptides from non-minimal protein groups removed")
waterfall <- add_step(waterfall, "shared peptides removed", nrow(qf[["log2"]]))

# Observed in at least min_cultures distinct biological cultures. Counting
# cultures rather than runs matters here: three technical injections of one
# culture are one observation of the biology, so a peptide seen in three runs
# of a single culture has been seen once, not three times. The threshold is
# arbitrary. Two of six cultures is lenient; --sensitivity re-runs the whole
# analysis across a range of it so the reader can see how much the answer
# moves.
count_cultures <- function(assay_mat, cultures) {
  obs <- !is.na(assay_mat)
  vapply(seq_len(nrow(obs)), function(i)
    length(unique(cultures[obs[i, ]])), integer(1))
}
cult <- colData(qf)$culture
ncult <- count_cultures(assay(qf[["log2"]]), cult)
qf <- qf[ncult >= min_cultures, , ]
msg(sum(ncult < min_cultures), " peptides seen in fewer than ", min_cultures,
    " cultures removed; ", nrow(qf[["log2"]]), " remain")
waterfall <- add_step(waterfall, paste0("observed in >= ", min_cultures, " cultures"),
                      nrow(qf[["log2"]]))

# --------------------------------------------------------- normalisation ----
# Median of ratios, computed on the log2 scale. For each run, the offset is
# the median across peptides of (log2 intensity in that run minus the median
# log2 intensity of that peptide across runs). Subtracting it puts every run
# on a common scale without assuming the runs share a total.
#
# The alternative was QFeatures' center.median, which subtracts each run's own
# median intensity. That is simpler but it is biased by which peptides happen
# to be missing in a run: a run with many missing low abundance peptides has a
# higher raw median through composition alone. Median of ratios compares each
# run against a per-peptide reference, so a peptide only contributes where it
# was actually measured. With no missing values the two coincide.
m <- assay(qf[["log2"]])
ref <- apply(m, 1, median, na.rm = TRUE)
offsets <- apply(m - ref, 2, median, na.rm = TRUE)
norm_assay <- sweep(m, 2, offsets, "-")
qf <- addAssay(qf, SummarizedExperiment(
  assays = list(assay = norm_assay),
  rowData = rowData(qf[["log2"]]), colData = colData(qf)), name = "norm")
qf <- addAssayLinkOneToOne(qf, from = "log2", to = "norm")
msg("median-of-ratios offsets range ", round(min(offsets), 3), " to ", round(max(offsets), 3))
write.table(data.frame(run = names(offsets), offset_log2 = as.numeric(offsets)),
            file.path(outdir, paste0(prefix, "_normalisation_offsets.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Long format intensities before and after, for the density figure. Written
# rather than plotted here so that figures.R owns all the drawing.
dens <- rbind(
  data.frame(run = rep(colnames(m), each = nrow(m)), stage = "before",
             log2_intensity = as.vector(m)),
  data.frame(run = rep(colnames(norm_assay), each = nrow(norm_assay)), stage = "after",
             log2_intensity = as.vector(norm_assay)))
dens <- dens[!is.na(dens$log2_intensity), ]
write.table(dens, file.path(outdir, paste0(prefix, "_normalisation_densities.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------- peptides to protein ----
# robustSummary fits an M-estimator across the peptides of a protein, so a
# single badly behaved peptide moves the protein estimate less than it would
# under a mean or a median-polish.
#
# No imputation anywhere in this script. The call would go here:
#   qf <- impute(qf, i = "norm", method = "MinDet")
# Leaving it out is the default because missingness in a label free DDA
# experiment is informative: a peptide is missing largely because it was too
# low in abundance to be picked for fragmentation, so imputing it from the
# observed distribution of that protein pulls low values upward and shrinks
# exactly the differences being tested. Someone might want it if they need a
# complete matrix for a method that cannot tolerate NA, such as PCA or most
# clustering. The assumption they take on is that values are missing at
# random given the observed data, which for this assay is close to known to
# be false.
# Snapshot before aggregation. Once a protein assay is attached, indexing
# the whole object with a peptide-length logical vector is out of bounds on
# that assay, which is what broke the first version of the sweep below.
qf_pep <- qf
qf <- aggregateFeatures(qf, i = "norm", fcol = "proteins",
                        name = "protein", fun = robustSummary, na.rm = TRUE)
n_prot <- nrow(qf[["protein"]])
msg(n_prot, " protein groups after robustSummary")
waterfall <- add_step(waterfall, "protein groups", n_prot)
write.table(waterfall, file.path(outdir, paste0(prefix, "_filtering_waterfall.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------- models ----
# Two models on the same data.
#
# The mixed model is the honest one. Eighteen runs are three cultures per
# genotype injected three times each, so the runs are not eighteen
# independent samples. (1 | culture) says that injections of one culture
# share a random offset, and the genotype effect is then judged against
# variation between cultures rather than between injections.
#
# The naive model is fitted on purpose, to be wrong in the specific way that
# is easy to be wrong. It treats all eighteen runs as independent replicates.
# Technical injections are far less variable than separate cultures, so the
# residual variance it estimates is too small, the standard errors are too
# narrow, and it calls far too much significant. This is pseudo-replication,
# and showing what it costs on real data is more use than a warning.
fit_one <- function(qf, formula, label) {
  ok <- TRUE
  out <- tryCatch(msqrob(object = qf, i = "protein", formula = formula,
                         overwrite = TRUE, modelColumnName = label),
                  error = function(e) { ok <<- FALSE; message("  model failed: ", conditionMessage(e)); NULL })
  if (!ok) return(NULL)
  out
}

contrast_name <- parse_arg("--contrast", "genotypeargP_KO")
msg("contrast: ", contrast_name, " = 0")

results_for <- function(qf, formula, label) {
  qf2 <- fit_one(qf, formula, label)
  if (is.null(qf2)) return(NULL)
  L <- makeContrast(paste0(contrast_name, " = 0"), parameterNames = contrast_name)
  qf2 <- hypothesisTest(object = qf2, i = "protein", contrast = L,
                        modelColumn = label, resultsColumnNamePrefix = paste0(label, "_"),
                        overwrite = TRUE)
  rd <- rowData(qf2[["protein"]])
  res <- rd[[paste0(label, "_", contrast_name)]]
  res$protein <- rownames(rd)
  models <- rd[[label]]
  res$fit_status <- vapply(models, getFitMethod, character(1))
  res$model <- label
  list(qf = qf2, res = as.data.frame(res))
}

# Is the random intercept identifiable at all? It needs at least one culture
# measured more than once. In --subset mode there is a single injection per
# culture, so (1 | culture) has one observation per level, the variance
# component has nothing to estimate from, and every fit is singular.
#
# The distinction matters for how the second model should be read. With
# eighteen runs, three injections of each culture, ~ genotype is the naive
# model and it is wrong: it counts injections as replicates. With six runs,
# one injection per culture, ~ genotype is the correct model, because then
# the runs really are six independent cultures. Same formula, opposite
# status, decided by the design and not by preference.
runs_per_culture <- table(coldata$culture)
mixed_identifiable <- any(runs_per_culture > 1)
msg("runs per culture: ", paste(sprintf("%s=%d", names(runs_per_culture),
    as.integer(runs_per_culture)), collapse = ", "))

if (mixed_identifiable) {
  msg("fitting ~ genotype + (1 | culture) and, for contrast, the naive ~ genotype")
  mixed <- results_for(qf, ~ genotype + (1 | culture), "mixed")
  naive <- results_for(qf, ~ genotype, "naive")
  second_model_role <- "naive, pseudo-replicated"
} else {
  msg("every culture contributes one run, so (1 | culture) is not identifiable.")
  msg("fitting ~ genotype only. Here that is the correct model, not the naive one.")
  mixed <- NULL
  naive <- results_for(qf, ~ genotype, "naive")
  second_model_role <- "correct: one run per culture, no technical replication"
}
writeLines(c(paste0("mixed_model_identifiable\t", mixed_identifiable),
             paste0("second_model_role\t", second_model_role),
             paste0("runs_per_culture\t",
                    paste(sprintf("%s=%d", names(runs_per_culture),
                                  as.integer(runs_per_culture)), collapse = ";"))),
           file.path(outdir, paste0(prefix, "_design.tsv")))

summarise_fit <- function(r, label) {
  if (is.null(r)) return(NULL)
  d <- r$res
  d$adjPval_BH <- p.adjust(d$pval, method = "BH")
  n_err <- sum(d$fit_status == "fitError", na.rm = TRUE)
  msg(label, ": ", nrow(d), " proteins, ", n_err, " fitError, ",
      sum(d$adjPval_BH < 0.05, na.rm = TRUE), " significant at 5 percent FDR")
  d
}
mixed_d <- summarise_fit(mixed, "mixed")
naive_d <- summarise_fit(naive, "naive")

# fitError proteins. msqrob2 records fitError when the model could not be
# estimated for that protein: usually every observation falls in one genotype,
# or one culture contributes all of the data, so the design matrix is rank
# deficient or the random effect has nothing to estimate its variance from.
# They are reported rather than dropped silently, because a protein that
# cannot be tested is a different thing from a protein that was tested and
# found unchanged.
for (nm in c("mixed", "naive")) {
  d <- get(paste0(nm, "_d"))
  if (is.null(d)) next
  write.table(d[order(d$pval), ],
              file.path(outdir, paste0(prefix, "_protein_results_", nm, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
# Stages 09 and 10 read the "_mixed" file as the primary result. When the
# mixed model was not identifiable, the single fitted model is the primary
# result, so it is written under that name as well. _design.tsv records
# which model actually produced it, so the file name cannot mislead on its
# own.
if (is.null(mixed_d) && !is.null(naive_d)) {
  write.table(naive_d[order(naive_d$pval), ],
              file.path(outdir, paste0(prefix, "_protein_results_mixed.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

if (is.null(mixed_d) && !is.null(naive_d)) {
  summ <- data.frame(
    metric = c("proteins_tested", "sig_single_model_5pct", "fitError_single_model",
               "median_df_single_model", "mixed_model_identifiable"),
    value = c(nrow(naive_d),
              sum(naive_d$adjPval_BH < 0.05, na.rm = TRUE),
              sum(naive_d$fit_status == "fitError", na.rm = TRUE),
              round(median(naive_d$df, na.rm = TRUE), 2),
              FALSE))
  write.table(summ, file.path(outdir, paste0(prefix, "_model_summary.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  print(summ)
}

if (!is.null(mixed_d) && !is.null(naive_d)) {
  common <- intersect(mixed_d$protein, naive_d$protein)
  mm <- mixed_d[match(common, mixed_d$protein), ]
  nn <- naive_d[match(common, naive_d$protein), ]
  cmp <- data.frame(
    protein = common,
    logFC_mixed = mm$logFC, logFC_naive = nn$logFC,
    se_mixed = mm$se, se_naive = nn$se,
    pval_mixed = mm$pval, pval_naive = nn$pval,
    adjP_mixed = mm$adjPval_BH, adjP_naive = nn$adjPval_BH,
    df_mixed = mm$df, df_naive = nn$df)
  write.table(cmp, file.path(outdir, paste0(prefix, "_model_comparison.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  se_ratio <- median(cmp$se_naive / cmp$se_mixed, na.rm = TRUE)
  summ <- data.frame(
    metric = c("proteins_tested", "sig_mixed_5pct", "sig_naive_5pct",
               "fitError_mixed", "fitError_naive", "median_se_ratio_naive_over_mixed",
               "median_df_mixed", "median_df_naive"),
    value = c(length(common),
              sum(cmp$adjP_mixed < 0.05, na.rm = TRUE),
              sum(cmp$adjP_naive < 0.05, na.rm = TRUE),
              sum(mixed_d$fit_status == "fitError", na.rm = TRUE),
              sum(naive_d$fit_status == "fitError", na.rm = TRUE),
              round(se_ratio, 4),
              round(median(cmp$df_mixed, na.rm = TRUE), 2),
              round(median(cmp$df_naive, na.rm = TRUE), 2)))
  write.table(summ, file.path(outdir, paste0(prefix, "_model_summary.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  print(summ)
}

# --------------------------------------------------- threshold sensitivity ----
# The culture threshold above is arbitrary, so report how much it matters
# instead of defending the particular number.
if (run_sensitivity) {
  msg("sensitivity sweep over --min-cultures")
  nc_all <- count_cultures(assay(qf_pep[["norm"]]), cult)
  sens <- do.call(rbind, lapply(2:5, function(k) {
    sub <- tryCatch(qf_pep[nc_all >= k, , ], error = function(e) NULL)
    if (is.null(sub)) return(NULL)
    sub <- tryCatch(aggregateFeatures(sub, i = "norm", fcol = "proteins",
                                      name = "protein_s", fun = robustSummary, na.rm = TRUE),
                    error = function(e) NULL)
    if (is.null(sub)) return(NULL)
    fml <- if (mixed_identifiable) ~ genotype + (1 | culture) else ~ genotype
    sub <- tryCatch(msqrob(sub, i = "protein_s", formula = fml,
                           overwrite = TRUE), error = function(e) NULL)
    if (is.null(sub)) return(NULL)
    L <- makeContrast(paste0(contrast_name, " = 0"), parameterNames = contrast_name)
    sub <- tryCatch(hypothesisTest(sub, i = "protein_s", contrast = L, overwrite = TRUE),
                    error = function(e) NULL)
    if (is.null(sub)) return(NULL)
    d <- as.data.frame(rowData(sub[["protein_s"]])[[contrast_name]])
    data.frame(min_cultures = k,
               peptides = sum(nc_all >= k),
               proteins = nrow(d),
               significant_5pct = sum(p.adjust(d$pval, "BH") < 0.05, na.rm = TRUE))
  }))
  write.table(sens, file.path(outdir, paste0(prefix, "_threshold_sensitivity.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  print(sens)
}

saveRDS(qf, file.path(outdir, paste0(prefix, "_qfeatures.rds")))
msg("done")
