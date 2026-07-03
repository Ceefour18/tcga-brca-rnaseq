################################################################################
# 02_preprocess.R
#
# Project : TCGA-BRCA RNA-seq — Basal-like vs Luminal A
# Purpose : Turn the full annotated cohort from script 01 into an analysis-ready
#           DESeqDataSet for the two-group contrast:
#             (1) keep primary tumours only
#             (2) keep the two PAM50 groups of interest (Basal, LumA)
#             (3) one sample per patient (independence)
#             (4) define the condition factor with LumA as reference
#             (5) build the DESeqDataSet on RAW counts
#             (6) pre-filter low-count genes
#             (7) QC: VST + PCA to confirm the subtypes separate
#
# Author  : Clement Akinsola (Ceefor Analytic Hub)
#
# Input   : data/processed/brca_se_raw_annotated.rds   (from 01)
# Outputs : data/processed/dds_basal_vs_luma.rds        (analysis-ready dds)
#           results/figures/qc_pca_subtype.png          (QC sanity check)
#           results/tables/preprocess_session_info.txt
################################################################################


# ---- 0. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)                # DESeqDataSet, vst, counts()
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(here)
})


# ---- 1. Paths ---------------------------------------------------------------
dir_proc    <- here("data", "processed")
dir_figures <- here("results", "figures")
dir_tables  <- here("results", "tables")
for (d in c(dir_figures, dir_tables)) dir.create(d, recursive = TRUE, showWarnings = FALSE)


# ---- 2. Load the annotated SE ----------------------------------------------
# This single object carries counts + subtype + PFI, all aligned by 01.
se <- readRDS(file.path(dir_proc, "brca_se_raw_annotated.rds"))
message("Loaded SE: ", nrow(se), " genes x ", ncol(se), " samples")


# ---- 3. Filter to primary tumours + the two PAM50 groups --------------------
# sample_type_code "01" = primary solid tumour. We drop metastatic ("06") and
# normal ("11") so the contrast is tumour-vs-tumour, not confounded by tissue.
# PAM50 %in% {Basal, LumA} both selects our groups AND drops NA labels
# (NA %in% x is FALSE), so unlabelled samples fall away cleanly here.
keep_cols <- colData(se)$sample_type_code == "01" &
             colData(se)$PAM50 %in% c("Basal", "LumA")

se_f <- se[, keep_cols]
message("After tumour + subtype filter: ", ncol(se_f), " samples")
print(table(colData(se_f)$PAM50))


# ---- 4. One sample per patient ----------------------------------------------
# A few patients contributed multiple aliquots (e.g. 01A and 01B). DE assumes
# independent samples, so replicates from the same patient would violate that
# and inflate significance. We keep ONE per patient — the most deeply sequenced
# aliquot (largest library size), which retains the highest-quality data rather
# than choosing arbitrarily.
lib_size <- colSums(assay(se_f, "unstranded"))

keep_barcodes <- as.data.frame(colData(se_f)) %>%
  mutate(lib_size = lib_size) %>%
  group_by(patient_barcode) %>%
  slice_max(lib_size, n = 1, with_ties = FALSE) %>%
  pull(barcode)

se_f <- se_f[, keep_barcodes]
message("After de-duplicating patients: ", ncol(se_f), " samples")


# ---- 5. Define the condition factor (reference = LumA) ----------------------
# Level ORDER sets the DE direction. LumA is the indolent, good-prognosis
# baseline, so we make it the reference level: DESeq2 then reports Basal
# RELATIVE TO LumA, i.e. a positive log2FC = up-regulated in Basal-like tumours.
# That is the biologically natural framing for this comparison.
se_f$condition <- factor(colData(se_f)$PAM50, levels = c("LumA", "Basal"))
print(table(se_f$condition))


# ---- 6. Build the DESeqDataSet on RAW counts --------------------------------
# The SE holds six assays; DESeq2 must model the integer RAW counts, which is
# the "unstranded" layer. We reduce the object to that single assay so DESeq2
# cannot silently pick a normalised layer (TPM/FPKM would be statistically
# wrong here), and coerce to integer as DESeq2 requires. Building from the SE
# (rather than a bare matrix) preserves rowData — the gene symbols/biotypes we
# need for annotation and enrichment downstream.
cts <- assay(se_f, "unstranded")
mode(cts) <- "integer"
assay(se_f, "unstranded") <- cts
assays(se_f) <- assays(se_f)["unstranded"]      # keep only raw counts

dds <- DESeqDataSet(se_f, design = ~ condition)


# ---- 7. Pre-filter low-count genes ------------------------------------------
# Not strictly required (DESeq2 does independent filtering at results time),
# but it removes uninformative near-zero genes, cutting memory/runtime and the
# multiple-testing burden. DESeq2's recommended rule: keep genes with at least
# 10 counts in at least as many samples as the SMALLER group — so a gene must
# be reliably present in enough samples to be one group's real signal.
smallest_group <- min(table(dds$condition))
keep_genes     <- rowSums(counts(dds) >= 10) >= smallest_group

message("Genes before filter: ", nrow(dds),
        " | after filter: ", sum(keep_genes),
        " (smallest group n = ", smallest_group, ")")
dds <- dds[keep_genes, ]


# ---- 8. QC — does the biology separate? -------------------------------------
# Sanity check BEFORE differential expression: variance-stabilising transform
# (blind = TRUE, unsupervised — it does NOT peek at the condition labels, so a
# clean split is honest evidence the groups differ, not an artefact of the
# design), then PCA on the most variable genes. Basal-like and Luminal A are
# the most divergent PAM50 subtypes, so we EXPECT clear separation on PC1. If
# they don't split, something upstream is wrong and it's cheap to catch now.
vsd    <- vst(dds, blind = TRUE)
pca_df <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct    <- round(100 * attr(pca_df, "percentVar"))

p <- ggplot(pca_df, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = paste0("PC1: ", pct[1], "% variance"),
       y = paste0("PC2: ", pct[2], "% variance"),
       title = "QC PCA — Basal-like vs Luminal A (VST, blind)",
       colour = "PAM50") +
  theme_bw()

ggsave(file.path(dir_figures, "qc_pca_subtype.png"),
       plot = p, width = 7, height = 5, dpi = 300)
message("Saved QC PCA to results/figures/qc_pca_subtype.png")


# ---- 9. Save the analysis-ready dds -----------------------------------------
saveRDS(dds, file.path(dir_proc, "dds_basal_vs_luma.rds"))
message("Saved dds_basal_vs_luma.rds  (", nrow(dds), " genes x ",
        ncol(dds), " samples)")


# ---- 10. Reproducibility record ---------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(dir_tables, "preprocess_session_info.txt"))

message("Done. Next: 03_deseq2_dea.R (fit the model, shrink LFCs, extract & ",
        "annotate the Basal-vs-LumA results).")
