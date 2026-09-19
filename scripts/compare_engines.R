#!/usr/bin/env Rscript
# Sage against MaxQuant on the same spectra.
#
# Same eighteen raw files, same statistical workflow, two different search
# engines and two different quantifications. Anything that differs in the
# answer is attributable to the search and the quantification, because
# nothing else was allowed to vary: both sides go through the identical
# differential.R.
#
# The two engines do not speak the same accession language. Sage searched a
# UniProt proteome and reports UniProt accessions; the authors searched an
# NCBI FASTA and the table reports RefSeq WP_ accessions. The bridge is
# UniProt's own RefSeq cross-reference, fetched in 09_compare_maxquant.sh.
# Proteins that fail to map are counted and reported rather than dropped
# quietly, because an unmappable protein is a real limit on the comparison
# and not a rounding error.

suppressPackageStartupMessages({ library(stats) })

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i)) return(default); args[i + 1]
}
if ("--help" %in% args || length(args) == 0) {
  cat("Usage: compare_engines.R --sage <tsv> --maxquant <tsv> --mapping <tsv> --outdir <dir>\n")
  cat("                         [--sage-matrix <tsv>] [--maxquant-matrix <tsv>] [--fdr 0.05]\n")
  quit(status = 0)
}
sage_file <- parse_arg("--sage")
mq_file   <- parse_arg("--maxquant")
map_file  <- parse_arg("--mapping")
outdir    <- parse_arg("--outdir", "results")
sage_mat  <- parse_arg("--sage-matrix")
mq_mat    <- parse_arg("--maxquant-matrix")
fdr_cut   <- as.numeric(parse_arg("--fdr", "0.05"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(sprintf("compare_engines.R: %s\n", paste0(...)))

sage <- read.delim(sage_file, stringsAsFactors = FALSE)
mq   <- read.delim(mq_file, stringsAsFactors = FALSE)
map  <- read.delim(map_file, stringsAsFactors = FALSE)

# map has columns uniprot, refseq (one row per pair, already exploded).
# Sage protein groups are semicolon-delimited header tokens like
# sp|A0Q671|A0Q671_FRATN; reduce each to a set of UniProt accessions and then
# to the set of RefSeq ids they correspond to.
uniprot_of <- function(x) {
  toks <- unlist(strsplit(x, ";"))
  acc <- sub("^(sp|tr)\\|([^|]+)\\|.*$", "\\2", toks)
  acc <- sub("^(sp|tr)\\|", "", acc)
  unique(acc[nzchar(acc)])
}
u2r <- split(map$refseq, map$uniprot)

sage_key <- vapply(sage$protein, function(p) {
  us <- uniprot_of(p)
  rs <- unique(unlist(u2r[us]))
  if (length(rs) == 0) NA_character_ else paste(sort(rs), collapse = ";")
}, character(1))

# MaxQuant accessions are WP_ with a version suffix in UniProt's xref but
# usually without one in the table, so compare on the unversioned stem.
strip_ver <- function(x) gsub("\\.[0-9]+$", "", x)
mq_key <- vapply(mq$protein, function(p) {
  rs <- strip_ver(unlist(strsplit(p, ";")))
  rs <- unique(rs[nzchar(rs)])
  if (length(rs) == 0) NA_character_ else paste(sort(rs), collapse = ";")
}, character(1))

sage$key <- sage_key
mq$key   <- mq_key

n_sage_unmapped <- sum(is.na(sage$key))
msg(n_sage_unmapped, " of ", nrow(sage), " Sage protein groups had no RefSeq cross-reference")

# Compare on single-accession keys. A group that merges several proteins in
# one engine and splits them in the other cannot be matched one to one, and
# forcing it would invent agreement or disagreement that is an artefact of
# grouping rather than of the search.
simple <- function(k) !is.na(k) & !grepl(";", k)
s1 <- sage[simple(sage$key), ]
m1 <- mq[simple(mq$key), ]
s1 <- s1[!duplicated(s1$key), ]
m1 <- m1[!duplicated(m1$key), ]

tested_s <- s1$key[!is.na(s1$pval) & s1$fit_status != "fitError"]
tested_m <- m1$key[!is.na(m1$pval) & m1$fit_status != "fitError"]
shared   <- intersect(tested_s, tested_m)
msg(length(tested_s), " Sage / ", length(tested_m), " MaxQuant testable proteins, ",
    length(shared), " shared")

sig_s <- s1$key[!is.na(s1$adjPval_BH) & s1$adjPval_BH < fdr_cut]
sig_m <- m1$key[!is.na(m1$adjPval_BH) & m1$adjPval_BH < fdr_cut]
sig_shared <- intersect(sig_s, sig_m)

ss <- s1[match(shared, s1$key), ]
mm <- m1[match(shared, m1$key), ]
rho_t  <- suppressWarnings(cor(ss$t, mm$t, method = "spearman", use = "complete.obs"))
rho_fc <- suppressWarnings(cor(ss$logFC, mm$logFC, method = "spearman", use = "complete.obs"))
pear_fc <- suppressWarnings(cor(ss$logFC, mm$logFC, method = "pearson", use = "complete.obs"))

top_n <- function(d, n = 100) d$key[order(d$pval)][seq_len(min(n, nrow(d)))]
top_s <- top_n(s1[!is.na(s1$pval), ])
top_m <- top_n(m1[!is.na(m1$pval), ])
top_overlap <- length(intersect(top_s, top_m))

# Peptide and protein counts from the two input matrices, if given.
pep_counts <- c(sage = NA_integer_, maxquant = NA_integer_)
if (!is.null(sage_mat) && file.exists(sage_mat))
  pep_counts["sage"] <- nrow(read.delim(sage_mat, stringsAsFactors = FALSE, check.names = FALSE))
if (!is.null(mq_mat) && file.exists(mq_mat))
  pep_counts["maxquant"] <- nrow(read.delim(mq_mat, stringsAsFactors = FALSE, check.names = FALSE))

jacc <- function(a, b) if (length(union(a, b)) == 0) NA else length(intersect(a, b)) / length(union(a, b))

out <- data.frame(
  metric = c("peptides_sage", "peptides_maxquant",
             "protein_groups_sage", "protein_groups_maxquant",
             "sage_groups_unmapped_to_refseq",
             "testable_sage", "testable_maxquant", "testable_shared",
             "jaccard_testable",
             "significant_sage", "significant_maxquant", "significant_shared",
             "jaccard_significant",
             "spearman_t_shared", "spearman_logFC_shared", "pearson_logFC_shared",
             "top100_overlap"),
  value = c(pep_counts["sage"], pep_counts["maxquant"],
            nrow(sage), nrow(mq), n_sage_unmapped,
            length(tested_s), length(tested_m), length(shared),
            round(jacc(tested_s, tested_m), 4),
            length(sig_s), length(sig_m), length(sig_shared),
            round(jacc(sig_s, sig_m), 4),
            round(rho_t, 4), round(rho_fc, 4), round(pear_fc, 4),
            top_overlap),
  stringsAsFactors = FALSE)
write.table(out, file.path(outdir, "engine_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
print(out, row.names = FALSE)

paired <- data.frame(key = shared,
                     logFC_sage = ss$logFC, logFC_maxquant = mm$logFC,
                     t_sage = ss$t, t_maxquant = mm$t,
                     adjP_sage = ss$adjPval_BH, adjP_maxquant = mm$adjPval_BH,
                     sig_sage = shared %in% sig_s, sig_maxquant = shared %in% sig_m)
write.table(paired, file.path(outdir, "engine_comparison_per_protein.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
msg("wrote engine_comparison.tsv and engine_comparison_per_protein.tsv")
