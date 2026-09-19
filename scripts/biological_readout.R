#!/usr/bin/env Rscript
# Are ribosomal proteins down in the argP mutant?
#
# This is a keyword check, and the distinction from an enrichment analysis is
# not pedantry. A keyword check asks whether proteins whose description
# contains "ribosomal protein" behave differently from the other proteins in
# the same tested set. It inherits every bias in how those descriptions were
# written, it cannot see a ribosome-associated protein that is not described
# as one, and it has no notion of a pathway. What it does have is a correct
# background: the comparison is against the proteins actually tested here,
# not against every protein in the genome. Testing 1,180 proteins and then
# comparing the significant ones against all 1,719 in the proteome would
# manufacture enrichment out of which proteins happened to be detected.
#
# Ramond et al. 2015 predict ribosomal proteins should be reduced in the
# mutant. Goeminne et al. 2016, reanalysing the same MaxQuant output,
# reported log2 fold change estimates for 49 ribosomal proteins with all of
# them pointing toward down-regulation. That is a specific, falsifiable
# prediction and this script either reproduces it or it does not.

suppressPackageStartupMessages({ library(stats) })

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i)) return(default); args[i + 1]
}
if ("--help" %in% args || length(args) == 0) {
  cat("Usage: biological_readout.R --results <tsv> --annotation <tsv> --out <tsv> [--keyword <regex>] [--fdr 0.05]\n")
  quit(status = 0)
}
res_file <- parse_arg("--results")
ann_file <- parse_arg("--annotation")
out_file <- parse_arg("--out")
keyword  <- parse_arg("--keyword", "ribosomal protein")
fdr_cut  <- as.numeric(parse_arg("--fdr", "0.05"))

res <- read.delim(res_file, stringsAsFactors = FALSE)
ann <- read.delim(ann_file, stringsAsFactors = FALSE)

# A protein group can list several accessions. It counts as matching if any
# of its members does, which is the permissive direction and is stated here
# rather than buried.
desc_for <- function(p) {
  accs <- strsplit(p, ";")[[1]]
  d <- ann$description[match(accs, ann$protein)]
  paste(d[!is.na(d)], collapse = " | ")
}
res$description <- vapply(res$protein, desc_for, character(1))
res$is_keyword <- grepl(keyword, res$description, ignore.case = TRUE)

# Proteins whose model failed were never tested and must not sit in the
# background as though they had been.
tested <- res[!is.na(res$pval) & res$fit_status != "fitError", ]
sig <- tested[!is.na(tested$adjPval_BH) & tested$adjPval_BH < fdr_cut, ]

n_tested    <- nrow(tested)
n_kw_tested <- sum(tested$is_keyword)
n_sig       <- nrow(sig)
n_kw_sig    <- sum(sig$is_keyword)

ct <- matrix(c(n_kw_sig, n_kw_tested - n_kw_sig,
               n_sig - n_kw_sig, n_tested - n_kw_tested - (n_sig - n_kw_sig)),
             nrow = 2, byrow = TRUE,
             dimnames = list(c("keyword", "other"), c("significant", "not")))
ft <- fisher.test(ct)

kw_fc  <- tested$logFC[tested$is_keyword]
oth_fc <- tested$logFC[!tested$is_keyword]
wt <- suppressWarnings(wilcox.test(kw_fc, oth_fc))

n_kw_down <- sum(kw_fc < 0, na.rm = TRUE)
n_kw_up   <- sum(kw_fc > 0, na.rm = TRUE)
# Sign test against the null that direction is a coin flip.
bt <- binom.test(n_kw_down, n_kw_down + n_kw_up, p = 0.5)

out <- data.frame(
  metric = c("keyword", "fdr_threshold",
             "proteins_tested", "keyword_proteins_tested",
             "proteins_significant", "keyword_proteins_significant",
             "keyword_fraction_of_significant", "keyword_fraction_of_tested",
             "fisher_odds_ratio", "fisher_p",
             "median_logFC_keyword", "median_logFC_other", "wilcoxon_p",
             "keyword_down", "keyword_up", "sign_test_p"),
  value = c(keyword, fdr_cut,
            n_tested, n_kw_tested, n_sig, n_kw_sig,
            round(ifelse(n_sig > 0, n_kw_sig / n_sig, NA), 4),
            round(n_kw_tested / n_tested, 4),
            round(unname(ft$estimate), 4), signif(ft$p.value, 4),
            round(median(kw_fc, na.rm = TRUE), 4),
            round(median(oth_fc, na.rm = TRUE), 4), signif(wt$p.value, 4),
            n_kw_down, n_kw_up, signif(bt$p.value, 4)),
  stringsAsFactors = FALSE)

write.table(out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
print(out, row.names = FALSE)

kwtab <- tested[tested$is_keyword, c("protein", "description", "logFC", "se", "pval", "adjPval_BH")]
kwtab <- kwtab[order(kwtab$logFC), ]
write.table(kwtab, sub("\\.tsv$", "_members.tsv", out_file),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("biological_readout.R: %d/%d %s proteins have negative logFC\n",
            n_kw_down, n_kw_down + n_kw_up, keyword))
