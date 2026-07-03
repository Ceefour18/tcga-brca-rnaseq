################################################################################
# 04_enrichment.R
#
# Project : TCGA-BRCA RNA-seq — Basal-like vs Luminal A
# Purpose : Turn the DE gene list into biology, two complementary ways:
#             (A) ORA (over-representation) — GO:BP + KEGG on the SIGNIFICANT
#                 DEGs, run SEPARATELY for up-in-Basal and up-in-LumA so each
#                 pathway list has a clear direction.
#             (B) GSEA — on the FULL ranked list (every gene, ranked by the
#                 Wald statistic). Threshold-free, so it is robust when there
#                 are thousands of DEGs (as here).
#
# Author  : Clement Akinsola (Ceefor Analytic Hub)
#
# Input   : data/processed/deseq2_results_annotated.rds   (from 03)
# Outputs : results/tables/  *_ora_*.csv, gsea_*.csv
#           results/figures/ *_dotplot.png, gsea_*_dotplot.png
#           data/processed/enrichment_results.rds
#           results/tables/enrichment_session_info.txt
#
# NOTE — KEGG (enrichKEGG / gseKEGG) queries the KEGG web API, so it needs an
#   internet connection at run time. GO uses the local org.Hs.eg.db and always
#   works offline. KEGG calls are wrapped in tryCatch so a network failure
#   downgrades gracefully (GO results still produced) instead of crashing.
################################################################################


# ---- 0. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(clusterProfiler)   # enrichGO/KEGG, gseGO/KEGG, bitr, setReadable
  library(org.Hs.eg.db)      # human gene-ID mappings (local)
  library(enrichplot)        # dotplot, ridgeplot, gseaplot2
  library(dplyr)
  library(ggplot2)
  library(here)
})


# ---- 1. Parameters & paths --------------------------------------------------
PADJ_CUTOFF <- 0.05    # must match script 03's DEG definition
LFC_CUTOFF  <- 1
set.seed(42)           # GSEA uses permutations — fix the seed for reproducibility

dir_proc    <- here("data", "processed")
dir_figures <- here("results", "figures")
dir_tables  <- here("results", "tables")
for (d in c(dir_figures, dir_tables)) dir.create(d, recursive = TRUE, showWarnings = FALSE)


# ---- 2. Load results & map IDs (Ensembl -> ENTREZ) --------------------------
res_df <- readRDS(file.path(dir_proc, "deseq2_results_annotated.rds"))

# clusterProfiler / KEGG need ENTREZ IDs. Our IDs are VERSIONED Ensembl
# (ENSG...15); org.Hs.eg.db uses UNVERSIONED IDs, so strip the ".xx" suffix
# first, then map. Not every gene maps (non-coding, retired IDs) — that's
# expected; we report the rate.
res_df$ensembl <- sub("\\..*$", "", res_df$gene_id)

id_map <- bitr(res_df$ensembl, fromType = "ENSEMBL",
               toType = "ENTREZID", OrgDb = org.Hs.eg.db)

res_map <- inner_join(res_df, id_map, by = c("ensembl" = "ENSEMBL"))
message("ID mapping: ", nrow(id_map), " of ", nrow(res_df),
        " genes mapped to ENTREZ (", round(100 * nrow(id_map) / nrow(res_df)), "%).")

# --- Background (universe) for ORA: ALL tested & mapped genes ---
# Using the expressed/tested gene set (not the whole genome) as the background
# is what makes over-representation p-values meaningful.
universe <- unique(res_map$ENTREZID)

# --- Directional DEG sets for ORA ---
sig_up   <- res_map %>%
  filter(!is.na(padj), padj < PADJ_CUTOFF, log2FoldChange >  LFC_CUTOFF) %>%
  pull(ENTREZID) %>% unique()
sig_down <- res_map %>%
  filter(!is.na(padj), padj < PADJ_CUTOFF, log2FoldChange < -LFC_CUTOFF) %>%
  pull(ENTREZID) %>% unique()
message("ORA input: ", length(sig_up), " up-in-Basal, ",
        length(sig_down), " up-in-LumA genes.")

# --- Ranked list for GSEA ---
# Rank by the Wald statistic (blends magnitude + significance). If several
# Ensembl IDs collapse to one ENTREZ, keep the most extreme; drop NA stats;
# sort decreasing (top = most up in Basal, bottom = most up in LumA).
rank_df <- res_map %>%
  filter(!is.na(stat)) %>%
  group_by(ENTREZID) %>%
  slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()
gene_ranks <- rank_df$stat
names(gene_ranks) <- rank_df$ENTREZID
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
message("GSEA ranked list: ", length(gene_ranks), " genes.")


# ---- Small helpers ----------------------------------------------------------
save_table <- function(obj, path) {
  if (is.null(obj)) return(0L)
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return(0L)
  write.csv(df, path, row.names = FALSE)
  nrow(df)
}
save_dotplot <- function(obj, path, title, showCategory = 15) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) {
    message("  (no enriched terms for: ", title, " — skipping plot)")
    return(invisible())
  }
  p <- dotplot(obj, showCategory = showCategory) + ggtitle(title)
  ggsave(path, p, width = 9, height = 7, dpi = 300)
}


# ---- 3. GO:BP over-representation (ORA), by direction -----------------------
# enrichGO tests whether GO Biological Process terms are over-represented among
# the DEGs vs the universe. simplify() then collapses redundant parent/child
# terms (GO is highly nested) for a cleaner, more readable pathway list.
run_go <- function(genes, label) {
  if (length(genes) < 10) { message("  too few genes for GO (", label, ")"); return(NULL) }
  ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                  ont = "BP", universe = universe, pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego <- clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  }
  ego
}

message("\nGO:BP ORA ...")
go_up   <- run_go(sig_up,   "up-in-Basal")
go_down <- run_go(sig_down, "up-in-LumA")

save_table(go_up,   file.path(dir_tables, "go_ora_up_in_basal.csv"))
save_table(go_down, file.path(dir_tables, "go_ora_up_in_luma.csv"))
save_dotplot(go_up,   file.path(dir_figures, "go_ora_up_in_basal_dotplot.png"),
             "GO:BP over-represented — Up in Basal")
save_dotplot(go_down, file.path(dir_figures, "go_ora_up_in_luma_dotplot.png"),
             "GO:BP over-represented — Up in Luminal A")


# ---- 4. KEGG over-representation (ORA), by direction ------------------------
# KEGG pathways = curated signalling/metabolic maps. Needs internet; guarded.
run_kegg <- function(genes, label) {
  if (length(genes) < 10) return(NULL)
  out <- tryCatch({
    ek <- enrichKEGG(gene = genes, organism = "hsa", universe = universe,
                     pAdjustMethod = "BH", pvalueCutoff = 0.05)
    if (!is.null(ek) && nrow(as.data.frame(ek)) > 0)
      ek <- setReadable(ek, org.Hs.eg.db, keyType = "ENTREZID")
    ek
  }, error = function(e) {
    message("  KEGG failed (", label, "): ", conditionMessage(e),
            "  [needs internet — GO results are unaffected]")
    NULL
  })
  out
}

message("KEGG ORA ...")
kegg_up   <- run_kegg(sig_up,   "up-in-Basal")
kegg_down <- run_kegg(sig_down, "up-in-LumA")

save_table(kegg_up,   file.path(dir_tables, "kegg_ora_up_in_basal.csv"))
save_table(kegg_down, file.path(dir_tables, "kegg_ora_up_in_luma.csv"))
save_dotplot(kegg_up,   file.path(dir_figures, "kegg_ora_up_in_basal_dotplot.png"),
             "KEGG over-represented — Up in Basal")
save_dotplot(kegg_down, file.path(dir_figures, "kegg_ora_up_in_luma_dotplot.png"),
             "KEGG over-represented — Up in Luminal A")


# ---- 5. GSEA on GO:BP -------------------------------------------------------
# GSEA walks the FULL ranked list and asks whether each gene set is concentrated
# at the top (positive NES = up in Basal) or bottom (negative NES = up in LumA).
# No DEG threshold — it uses the whole ranking, which is why it stays clean when
# thousands of genes are significant. eps = 0 gives accurate tiny p-values.
message("\nGSEA (GO:BP) ...")
gsea_go <- gseGO(geneList = gene_ranks, OrgDb = org.Hs.eg.db, ont = "BP",
                 keyType = "ENTREZID", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, eps = 0, seed = TRUE, verbose = FALSE)

# Collapse redundant GO terms (as we did for the GO ORA). GO:BP is deeply
# nested, so GSEA returns many near-duplicate cilium/segregation terms that
# crowd the dotplot; simplify() keeps one representative per cluster of
# semantically similar terms, giving a cleaner, more distinct pathway view.
if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
  gsea_go <- tryCatch(
    clusterProfiler::simplify(gsea_go, cutoff = 0.7, by = "p.adjust", select_fun = min),
    error = function(e) { message("  GSEA simplify skipped: ", conditionMessage(e)); gsea_go }
  )
  gsea_go <- setReadable(gsea_go, org.Hs.eg.db, keyType = "ENTREZID")
}
save_table(gsea_go, file.path(dir_tables, "gsea_go_bp.csv"))

# Directional dotplot: activated (up in Basal) vs suppressed (up in LumA).
# Fixes vs the first version: order by significance (orderBy = "p.adjust") so
# the most meaningful terms surface; show fewer terms (8/side) and wrap long
# labels (label_format) so nothing overlaps; taller canvas for breathing room.
if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
  p <- dotplot(gsea_go, showCategory = 8, split = ".sign",
               orderBy = "p.adjust", label_format = 40, font.size = 10) +
       facet_grid(. ~ .sign) +
       ggtitle("GSEA GO:BP — Basal (activated) vs LumA (suppressed)")
  ggsave(file.path(dir_figures, "gsea_go_dotplot.png"), p, width = 12, height = 9, dpi = 300)
}


# ---- 6. GSEA on KEGG --------------------------------------------------------
message("GSEA (KEGG) ...")
gsea_kegg <- tryCatch({
  g <- gseKEGG(geneList = gene_ranks, organism = "hsa", pAdjustMethod = "BH",
               pvalueCutoff = 0.05, eps = 0, seed = TRUE, verbose = FALSE)
  if (!is.null(g) && nrow(as.data.frame(g)) > 0)
    g <- setReadable(g, org.Hs.eg.db, keyType = "ENTREZID")
  g
}, error = function(e) {
  message("  gseKEGG failed: ", conditionMessage(e), "  [needs internet]")
  NULL
})
save_table(gsea_kegg, file.path(dir_tables, "gsea_kegg.csv"))

if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
  p <- dotplot(gsea_kegg, showCategory = 12, split = ".sign") +
       facet_grid(. ~ .sign) +
       ggtitle("GSEA KEGG — Basal (activated) vs LumA (suppressed)")
  ggsave(file.path(dir_figures, "gsea_kegg_dotplot.png"), p, width = 11, height = 7, dpi = 300)
}


# ---- 7. Save all objects for the report -------------------------------------
saveRDS(list(go_up = go_up, go_down = go_down,
             kegg_up = kegg_up, kegg_down = kegg_down,
             gsea_go = gsea_go, gsea_kegg = gsea_kegg,
             gene_ranks = gene_ranks),
        file.path(dir_proc, "enrichment_results.rds"))
message("\nSaved enrichment_results.rds")


# ---- 8. Reproducibility record ----------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(dir_tables, "enrichment_session_info.txt"))

message("Done. Next: 05_survival.R (Kaplan-Meier + age-adjusted Cox on PFI, ",
        "with OS descriptives).")
