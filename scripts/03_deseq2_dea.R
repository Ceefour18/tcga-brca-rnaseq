################################################################################
# 03_deseq2_dea.R
#
# Project : TCGA-BRCA RNA-seq — Basal-like vs Luminal A
# Purpose : Differential expression between Basal-like and Luminal A tumours.
#             (1) fit the DESeq2 negative-binomial model
#             (2) extract the Basal-vs-LumA contrast (LumA = reference)
#             (3) apply apeglm LFC shrinkage (trustworthy fold-changes)
#             (4) annotate Ensembl IDs -> gene symbols (from rowData)
#             (5) call DEGs at padj < 0.05 & |log2FC| > 1
#             (6) save full + significant tables, volcano and MA plots, and an
#                 annotated results object for the enrichment step (04)
#
# Author  : Clement Akinsola (Ceefor Analytic Hub)
#
# Input   : data/processed/dds_basal_vs_luma.rds          (from 02)
# Outputs : results/tables/deseq2_results_full.csv        (all genes, annotated)
#           results/tables/deseq2_DEGs_significant.csv     (thresholded DEGs)
#           results/figures/volcano_basal_vs_luma.png
#           results/figures/ma_plot_basal_vs_luma.png
#           data/processed/deseq2_results_annotated.rds    (feeds script 04)
#           results/tables/dea_session_info.txt
################################################################################


# ---- 0. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)                # DESeq(), results(), lfcShrink(), plotMA()
  library(SummarizedExperiment)  # rowData() for gene annotation
  library(dplyr)
  library(ggplot2)
  library(apeglm)                # LFC shrinkage estimator (used by lfcShrink)
  library(here)
})
# ggrepel is OPTIONAL — only used to label top genes on the volcano. The script
# runs fine without it (labels are simply skipped).
have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)


# ---- 1. Parameters & paths --------------------------------------------------
PADJ_CUTOFF <- 0.05    # FDR threshold for calling a gene significant
LFC_CUTOFF  <- 1       # |log2FC| > 1  == at least a 2-fold change

dir_proc    <- here("data", "processed")
dir_figures <- here("results", "figures")
dir_tables  <- here("results", "tables")
for (d in c(dir_figures, dir_tables)) dir.create(d, recursive = TRUE, showWarnings = FALSE)


# ---- 2. Load the analysis-ready dds -----------------------------------------
dds <- readRDS(file.path(dir_proc, "dds_basal_vs_luma.rds"))
message("Loaded dds: ", nrow(dds), " genes x ", ncol(dds), " samples")
message("Design: ", deparse(design(dds)),
        " | reference level: ", levels(dds$condition)[1])


# ---- 3. Fit the model -------------------------------------------------------
# DESeq() runs three steps in one call:
#   (a) size-factor estimation  -> corrects for sequencing-depth differences
#   (b) dispersion estimation   -> models each gene's biological variability,
#                                  borrowing strength across genes (shrinkage
#                                  toward a fitted trend) so low-replicate
#                                  variance estimates aren't wild
#   (c) negative-binomial GLM + Wald test for the condition effect
# It operates on RAW counts (never VST/TPM) — the NB model needs the true
# count distribution to weight noisy low-count genes correctly.
dds <- DESeq(dds)
message("Model fitted. resultsNames: ",
        paste(resultsNames(dds), collapse = ", "))


# ---- 4. Extract the Basal-vs-LumA contrast ----------------------------------
# With LumA as reference, the coefficient is 'condition_Basal_vs_LumA', so a
# positive log2FC = higher expression in Basal-like tumours. alpha = PADJ_CUTOFF
# tells results() to tune its independent filtering for that FDR threshold.
coef_name <- grep("^condition", resultsNames(dds), value = TRUE)[1]
res <- results(dds, name = coef_name, alpha = PADJ_CUTOFF)
message("\nResults summary (unshrunken):")
summary(res)


# ---- 5. LFC shrinkage (apeglm) ----------------------------------------------
# Raw fold-changes for low-count genes are noisy and often wildly overstated.
# apeglm shrinks each gene's log2FC toward 0 in proportion to its uncertainty:
# well-powered genes barely move, poorly-estimated ones are pulled in. This
# gives fold-changes that are safe to RANK and PLOT. p-values / padj are the
# original Wald-test values (shrinkage changes effect sizes, not significance).
res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")


# ---- 6. Assemble one annotated results table --------------------------------
# res_shrunk, res and rowData(dds) share the same gene order, so we can combine
# directly. We keep the SHRUNKEN log2FC (for thresholds/plots) AND the Wald
# 'stat' from the unshrunken fit (a good ranking metric for GSEA in script 04).
# Gene symbols come from rowData — the same GENCODE annotation the counts were
# built on — so no external ID-mapping is needed here.
stopifnot(identical(rownames(res_shrunk), rownames(res)),
          identical(rownames(res_shrunk), rownames(dds)))

res_df <- as.data.frame(res_shrunk) %>%
  mutate(
    gene_id   = rownames(res_shrunk),
    stat      = res$stat,                    # Wald statistic (unshrunken)
    symbol    = rowData(dds)$gene_name,
    gene_type = rowData(dds)$gene_type
  ) %>%
  select(gene_id, symbol, gene_type, baseMean,
         log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  arrange(padj)


# ---- 7. Call DEGs & report --------------------------------------------------
# is.na guards matter: padj is NA for genes dropped by independent filtering or
# flagged as count outliers — those must not be counted as significant.
is_sig <- !is.na(res_df$padj) & res_df$padj < PADJ_CUTOFF &
          abs(res_df$log2FoldChange) > LFC_CUTOFF

n_up   <- sum(is_sig & res_df$log2FoldChange > 0)
n_down <- sum(is_sig & res_df$log2FoldChange < 0)

message("\n--- DEG summary (padj < ", PADJ_CUTOFF,
        " & |log2FC| > ", LFC_CUTOFF, ") ---")
message("Total DEGs   : ", sum(is_sig))
message("Up in Basal  : ", n_up)
message("Down in Basal: ", n_down)

deg_df <- res_df[is_sig, ]

write.csv(res_df, file.path(dir_tables, "deseq2_results_full.csv"), row.names = FALSE)
write.csv(deg_df, file.path(dir_tables, "deseq2_DEGs_significant.csv"), row.names = FALSE)


# ---- 8. Volcano plot --------------------------------------------------------
# x = shrunken log2FC (magnitude), y = -log10(pvalue) (evidence). Colour by the
# combined significance call so the reader sees direction and strength at once.
res_df$category <- "Not sig."
res_df$category[is_sig & res_df$log2FoldChange > 0] <- "Up in Basal"
res_df$category[is_sig & res_df$log2FoldChange < 0] <- "Down in Basal"

volcano <- ggplot(res_df, aes(log2FoldChange, -log10(pvalue), colour = category)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_colour_manual(values = c("Up in Basal"   = "#D55E00",
                                 "Down in Basal" = "#0072B2",
                                 "Not sig."      = "grey80")) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  labs(title = "Basal-like vs Luminal A (TCGA-BRCA)",
       subtitle = paste0(sum(is_sig), " DEGs at padj < ", PADJ_CUTOFF,
                         ", |log2FC| > ", LFC_CUTOFF),
       x = "log2 fold-change (Basal vs LumA, apeglm-shrunk)",
       y = "-log10(p-value)", colour = NULL) +
  theme_bw() + theme(legend.position = "top")

# Label the top genes by significance, if ggrepel is available.
if (have_ggrepel) {
  top_genes <- res_df[is_sig, ] %>% arrange(padj) %>% head(20)
  volcano <- volcano +
    ggrepel::geom_text_repel(data = top_genes,
                             aes(label = symbol), size = 3, max.overlaps = 20,
                             show.legend = FALSE)
} else {
  message("(ggrepel not installed — skipping gene labels on volcano. ",
          "Add r-ggrepel to the env to enable.)")
}

ggsave(file.path(dir_figures, "volcano_basal_vs_luma.png"),
       volcano, width = 8, height = 6, dpi = 300)
message("Saved volcano plot.")


# ---- 9. MA plot -------------------------------------------------------------
# Diagnostic: log2FC vs mean expression. On the shrunken results you should see
# high-count genes retain their fold-changes while low-count genes are pulled
# toward 0 — a visual confirmation shrinkage behaved as intended.
png(file.path(dir_figures, "ma_plot_basal_vs_luma.png"),
    width = 8, height = 6, units = "in", res = 300)
plotMA(res_shrunk, ylim = c(-6, 6),
       main = "MA plot (apeglm-shrunk) — Basal vs LumA")
dev.off()
message("Saved MA plot.")


# ---- 10. Save annotated results for downstream ------------------------------
# Script 04 (enrichment) reads this: it has gene_id, symbol, shrunken log2FC,
# and the Wald 'stat' for building the GSEA ranking.
saveRDS(res_df, file.path(dir_proc, "deseq2_results_annotated.rds"))
message("Saved deseq2_results_annotated.rds")


# ---- 11. Reproducibility record ---------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(dir_tables, "dea_session_info.txt"))

message("\nDone. Next: 04_enrichment.R (GO / KEGG over-representation on the ",
        "DEGs, and GSEA on the full ranked list).")
