#!/usr/bin/env Rscript
# Every figure in the README. Each one is drawn from a file in results/, so
# nothing in the README is a number I typed by hand.
#
# Figures are skipped, with a message, when their input does not exist. That
# matters on a machine which could only run part of the pipeline: the stages
# that did run still produce their plots instead of the whole script failing
# on the first missing file.

suppressPackageStartupMessages({
  library(ggplot2); library(data.table)
})
theme_set(theme_bw(base_size = 11) +
          theme(panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "grey92", colour = NA)))

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i)) return(default); args[i + 1]
}
resdir <- parse_arg("--results", "results")
figdir <- parse_arg("--figures", "figures")
prefix <- parse_arg("--prefix", "maxquant")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

msg  <- function(...) cat(sprintf("figures.R: %s\n", paste0(...)))
rp   <- function(...) file.path(resdir, paste0(...))
have <- function(p) file.exists(p) && file.size(p) > 0

save_fig <- function(plot, name, width = 7, height = 4.5, dpi = 150) {
  out <- file.path(figdir, name)
  ggsave(out, plot, width = width, height = height, dpi = dpi)
  msg("wrote ", out)
}

# --------------------------------------------- 1. normalisation densities ----
f <- rp(prefix, "_normalisation_densities.tsv")
if (have(f)) {
  d <- fread(f)
  d[, stage := factor(stage, levels = c("before", "after"),
                      labels = c("before normalisation", "after normalisation"))]
  p <- ggplot(d, aes(x = log2_intensity, group = run, colour = run)) +
    geom_density(linewidth = 0.35, alpha = 0.8) +
    facet_wrap(~ stage, nrow = 1) +
    guides(colour = "none") +
    labs(x = expression(log[2]~intensity), y = "density",
         title = "Median-of-ratios normalisation",
         subtitle = sprintf("%s arm, %d runs", prefix, length(unique(d$run))))
  save_fig(p, paste0(prefix, "_normalisation_densities.png"), 8, 4)
} else msg("skip densities, no ", f)

# ------------------------------------------------------------ 2. volcano ----
f <- rp(prefix, "_protein_results_mixed.tsv")
if (have(f)) {
  d <- fread(f)
  d <- d[!is.na(pval) & fit_status != "fitError"]
  d[, sig := !is.na(adjPval_BH) & adjPval_BH < 0.05]
  nsig <- sum(d$sig)
  p <- ggplot(d, aes(x = logFC, y = -log10(pval))) +
    geom_point(aes(colour = sig), size = 0.7, alpha = 0.6) +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#b2182b"),
                        labels = c("not significant", "BH < 0.05"), name = NULL) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
    labs(x = expression(log[2]~fold~change~(argP~KO~over~WT)),
         y = expression(-log[10]~p),
         title = "Differential abundance, mixed model",
         subtitle = sprintf("%s arm, %d proteins tested, %d significant at 5 percent FDR",
                            prefix, nrow(d), nsig)) +
    theme(legend.position = "top")
  save_fig(p, paste0(prefix, "_volcano_mixed.png"))
} else msg("skip volcano, no ", f)

# -------------------------------------- 3. mixed against naive, the cost -----
f <- rp(prefix, "_model_comparison.tsv")
if (have(f)) {
  d <- fread(f)
  d <- d[is.finite(se_mixed) & is.finite(se_naive) & se_mixed > 0]
  ratio <- median(d$se_naive / d$se_mixed, na.rm = TRUE)
  p1 <- ggplot(d, aes(x = se_mixed, y = se_naive)) +
    geom_point(size = 0.6, alpha = 0.4, colour = "grey30") +
    geom_abline(slope = 1, intercept = 0, colour = "#b2182b", linewidth = 0.5) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "standard error, mixed model", y = "standard error, naive model",
         title = "Pseudo-replication shrinks the standard error",
         subtitle = sprintf("median naive/mixed ratio %.3f; points below the line are over-confident",
                            ratio))
  save_fig(p1, paste0(prefix, "_se_mixed_vs_naive.png"), 5.5, 5)

  sig <- data.table(
    model = c("mixed\n~ genotype + (1 | culture)", "naive\n~ genotype"),
    n = c(sum(d$adjP_mixed < 0.05, na.rm = TRUE), sum(d$adjP_naive < 0.05, na.rm = TRUE)))
  p2 <- ggplot(sig, aes(x = model, y = n, fill = model)) +
    geom_col(width = 0.55) +
    geom_text(aes(label = n), vjust = -0.35, size = 4) +
    scale_fill_manual(values = c("#2166ac", "#b2182b")) +
    guides(fill = "none") +
    expand_limits(y = max(sig$n) * 1.15) +
    labs(x = NULL, y = "proteins significant at 5 percent FDR",
         title = "What treating injections as replicates costs",
         subtitle = sprintf("%s arm, same peptides, same normalisation, same contrast", prefix))
  save_fig(p2, paste0(prefix, "_model_significance_counts.png"), 5.5, 4.5)
} else msg("skip model comparison, no ", f)

# --------------------------------------------- 4. ribosomal keyword check ----
f  <- rp(prefix, "_protein_results_mixed.tsv")
fm <- rp(prefix, "_ribosomal_check_members.tsv")
if (have(f) && have(fm)) {
  d <- fread(f); d <- d[!is.na(pval) & fit_status != "fitError"]
  mem <- fread(fm)
  d[, ribosomal := protein %in% mem$protein]
  d[, grp := ifelse(ribosomal, "ribosomal protein", "all other tested proteins")]
  med <- d[, .(m = median(logFC, na.rm = TRUE), n = .N), by = grp]
  p <- ggplot(d, aes(x = grp, y = logFC, fill = grp)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_violin(alpha = 0.45, colour = NA, scale = "width") +
    geom_boxplot(width = 0.16, outlier.size = 0.4, fill = "white") +
    scale_fill_manual(values = c("grey60", "#2166ac")) +
    guides(fill = "none") +
    labs(x = NULL, y = expression(log[2]~fold~change~(argP~KO~over~WT)),
         title = "Ribosomal proteins are lower in the argP mutant",
         subtitle = sprintf("median %.3f across %d ribosomal against %.3f across %d others",
                            med[grp == "ribosomal protein", m], med[grp == "ribosomal protein", n],
                            med[grp != "ribosomal protein", m], med[grp != "ribosomal protein", n]))
  save_fig(p, paste0(prefix, "_ribosomal_logfc.png"), 6.5, 4.5)
} else msg("skip ribosomal figure, no ", fm)

# -------------------------------------------- 5. entrapment FDR assessment ---
f <- rp("entrapment_fdr_assessment.tsv")
if (have(f)) {
  d <- fread(f)
  long <- melt(d, id.vars = c("nominal_fdr", "level", "n_target", "n_entrapment"),
               measure.vars = c("combined_eq1", "lower_bound_eq2", "sample_eq3", "paired_eq4"),
               variable.name = "estimator", value.name = "estimate")
  long <- long[!is.na(estimate) & estimate != "NA"]
  long[, estimate := as.numeric(estimate)]
  long[, estimator := factor(estimator,
        levels = c("paired_eq4", "combined_eq1", "lower_bound_eq2", "sample_eq3"),
        labels = c("paired (eq 4, proposed)", "combined (eq 1, valid bound)",
                   "lower bound (eq 2)", "sample (eq 3, unusable)"))]
  p <- ggplot(long, aes(x = level, y = estimate, fill = estimator)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_hline(aes(yintercept = nominal_fdr), linetype = 2, colour = "#b2182b") +
    facet_wrap(~ paste0("nominal FDR ", nominal_fdr), scales = "free_y") +
    scale_fill_brewer(palette = "Blues", direction = -1, name = NULL) +
    labs(x = NULL, y = "estimated false discovery proportion",
         title = "Entrapment assessment of the reported FDR",
         subtitle = "dashed line is the nominal rate; only the upper bounds can establish control") +
    theme(legend.position = "top")
  save_fig(p, "entrapment_fdr_assessment.png", 8.5, 5)
} else msg("skip entrapment figure, no ", f)

# ------------------------------------------------- 6. engine comparison ------
f <- rp("engine_comparison_per_protein.tsv")
if (have(f)) {
  d <- fread(f)
  d[, agree := fifelse(sig_sage & sig_maxquant, "both",
               fifelse(sig_sage, "Sage only",
               fifelse(sig_maxquant, "MaxQuant only", "neither")))]
  rho <- suppressWarnings(cor(d$logFC_sage, d$logFC_maxquant,
                              method = "spearman", use = "complete.obs"))
  p <- ggplot(d, aes(x = logFC_maxquant, y = logFC_sage, colour = agree)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50", linewidth = 0.4) +
    geom_point(size = 0.8, alpha = 0.65) +
    scale_colour_manual(values = c(both = "#b2182b", `Sage only` = "#2166ac",
                                   `MaxQuant only` = "#f4a582", neither = "grey75"),
                        name = NULL) +
    labs(x = expression(log[2]~FC~MaxQuant), y = expression(log[2]~FC~Sage),
         title = "Sage against MaxQuant on the same spectra",
         subtitle = sprintf("%d shared proteins, Spearman rho %.3f", nrow(d), rho)) +
    theme(legend.position = "top")
  save_fig(p, "engine_logfc_scatter.png", 6.5, 5.5)
} else msg("skip engine comparison figure, no ", f)

# ------------------------------------------------ 7. threshold sensitivity ---
f <- rp(prefix, "_threshold_sensitivity.tsv")
if (have(f)) {
  d <- fread(f)
  p <- ggplot(d, aes(x = min_cultures, y = significant_5pct)) +
    geom_line(colour = "#2166ac") +
    geom_point(size = 2, colour = "#2166ac") +
    geom_text(aes(label = paste0(proteins, " prot")), vjust = -0.9, size = 3, colour = "grey30") +
    labs(x = "peptide must be observed in at least N cultures",
         y = "proteins significant at 5 percent FDR",
         title = "How much the arbitrary filter matters",
         subtitle = sprintf("%s arm, mixed model", prefix))
  save_fig(p, paste0(prefix, "_threshold_sensitivity.png"), 6, 4)
} else msg("skip sensitivity figure, no ", f)

msg("done")
