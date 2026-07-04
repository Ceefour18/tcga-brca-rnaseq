Transcriptomic and Survival Landscape of Basal-like vs Luminal A Breast Cancer (TCGA-BRCA)

![Analysis](https://img.shields.io/badge/Analysis-RNA--seq-blue)
![Data](https://img.shields.io/badge/Data-TCGA--BRCA-orange)
![DE](https://img.shields.io/badge/DESeq2-1.44-blue)
![Enrichment](https://img.shields.io/badge/clusterProfiler-4.0-blueviolet)
![Survival](https://img.shields.io/badge/Survival-Cox%20%2B%20KM-green)
![Status](https://img.shields.io/badge/Status-Complete-brightgreen)
![License](https://img.shields.io/badge/License-MIT-yellow)

A reproducible transcriptomic and survival analysis comparing PAM50 Basal-like
and Luminal A breast tumours in the TCGA-BRCA cohort. The pipeline runs from
GDC data acquisition through differential expression, functional enrichment, and
survival analysis, with every result traceable to a saved table or figure.


Biological question

Basal-like and Luminal A are the two most transcriptionally divergent PAM50
subtypes. This project asks: (1) which genes and pathways distinguish them, and
(2) does subtype carry prognostic information for progression-free interval (PFI)
and overall survival (OS)?


Key results

All numbers below are the actual output of the pipeline on this cohort. Figures
and tables referenced are in results/.

Cohort


Analysis cohort: 752 primary tumours — 190 Basal-like, 562 Luminal A —
one (deepest-sequenced) sample per patient.
Quality control: blind variance-stabilised PCA separates the two subtypes
along PC1 (27% of variance), confirming they are globally distinct before any
testing.


Show Image

Differential expression (DESeq2)


5,346 differentially expressed genes at padj < 0.05 and |log2FC| > 1:
2,905 up in Basal-like, 2,441 up in Luminal A (LumA as reference;
positive log2FC = higher in Basal). Fold-changes are apeglm-shrunk.
Orientation was validated against canonical markers: luminal genes
(ESR1, GATA3, FOXA1, PGR) are higher in Luminal A, and basal markers
(FOXC1, SOX8) higher in Basal-like — consistent with expected subtype
biology.


Show Image


On these p-values: with 190 vs 562 samples, many genes reach extremely
small p-values. Significance alone is therefore not the whole story — the
results are supported by large effect sizes and by concordance with
established subtype markers, not by p-values in isolation.



Functional enrichment (clusterProfiler: GO/KEGG ORA + GSEA)

Enrichment was run separately by direction (up-in-Basal vs up-in-LumA) and by
two methods (over-representation on the DEGs; GSEA on the full ranked list).
The two methods are concordant:


Up in Basal-like: cell cycle, DNA replication, chromosome segregation,
mitotic division (a proliferation signature); plus adaptive/innate immune and
cytokine terms, and epidermis/keratinization (basal-cytokeratin) programs.
Up in Luminal A: hormone regulation, transport and signalling; cytochrome
P450 drug and xenobiotic metabolism; differentiated epithelial metabolism.


Tables: results/tables/{go_ora_*,kegg_ora_*,gsea_go_bp,gsea_kegg}.csv.
Figures: results/figures/{go_ora_*,kegg_ora_*,gsea_*}_dotplot.png.

Survival (PFI primary; OS secondary — both modelled in parallel)

PFI and OS were analysed identically: Kaplan-Meier + log-rank, unadjusted and
age-adjusted Cox, and a proportional-hazards check.

EndpointEventsMedian follow-upLog-rank pAge-adj. HR (Basal vs LumA)Cox pPH global pPFI96/750 (12.8%)31.3 mo0.0891.46 (95% CI 0.95–2.25)0.0870.0018OS93/750 (12.4%)31.7 mo0.4791.39 (95% CI 0.88–2.20)0.1590.051

What the data supports (and what it does not):


Neither endpoint reached statistical significance. The subtype hazard
ratios point above 1 (higher event hazard in Basal-like) for both endpoints,
but all confidence intervals cross 1 and all p-values exceed 0.05. This is a
non-significant trend, not a demonstrated survival difference.
Age adjustment barely changed the subtype HR (PFI 1.446 → 1.461; OS
1.176 → 1.390), so the subtype estimate is not confounded by age. Age itself
was a statistically significant predictor of OS but with a negligible
per-year effect (HR ≈ 1.0 per year).
The Kaplan-Meier curves cross for both endpoints: Basal-like has more
early events, Luminal A more late events. This crossing corresponds to the
violated proportional-hazards assumption for PFI (Schoenfeld global
p = 0.0018), which means a single Cox HR and the log-rank test — both
time-averaged — understate the early divergence and are only approximate
summaries here (results/figures/km_{pfi,os}_basal_vs_luma.png).
Event counts were near-identical across endpoints (96 vs 93), so PFI is not
favoured on event count; it is retained as primary because progression is the
more subtype-relevant event and OS is more affected by age-related competing
mortality.



On these p-values: a non-significant result does not prove the absence of
a difference — the crossing curves and modest cohort size (few patients at
risk beyond ~14 years) limit power. Late-time separation in the curves rests
on very small numbers at risk and should not be over-interpreted.



Show Image

Data-derived summary: results/tables/survival_summary.txt.


Repository structure
tcga-brca-rnaseq/
├── environment.yml              # conda environment (all dependencies)
├── README.md
├── .gitignore                   # excludes data/ downloads
├── scripts/
│   ├── 01_data_acquisition.R    # GDC STAR counts + PAM50 + PFI/OS + covariates
│   ├── 02_preprocess.R          # subset, dedupe, filter, DESeqDataSet, QC PCA
│   ├── 03_deseq2_dea.R          # differential expression + volcano/MA
│   ├── 04_enrichment.R          # GO/KEGG ORA + GSEA
│   └── 05_survival.R            # KM + Cox (PFI & OS), comparison
├── data/
│   ├── raw/                     # TCGA-CDR xlsx + GDC cache (git-ignored)
│   └── processed/               # saved .rds objects, metadata (git-ignored)
└── results/
    ├── figures/                 # committed figures
    └── tables/                  # committed result tables


Data sources and provenance


Expression: GDC harmonised STAR - Counts (open-access; GENCODE v36,
GRCh38), retrieved with TCGAbiolinks. Raw counts (unstranded assay) are
used for DESeq2. (FASTQ/BAM are dbGaP controlled-access and are not used.)
PAM50 subtypes: TCGAquery_subtype("brca") (Ciriello et al., 2018).
Survival endpoints & covariates (PFI, OS, age, stage, race, gender):
TCGA Pan-Cancer Clinical Data Resource (Liu et al., 2018),
TCGA-CDR-SupplementalTableS1.xlsx — the one file downloaded manually into
data/raw/ (see script 01 header).

# 1. build the environment (installs all R/Bioconductor packages)
conda config --set channel_priority strict
conda env create -f environment.yml
conda activate tcga-brca-rnaseq

# 2. place TCGA-CDR-SupplementalTableS1.xlsx in data/raw/

# 3. run the pipeline in order
Rscript scripts/01_data_acquisition.R
Rscript scripts/02_preprocess.R
Rscript scripts/03_deseq2_dea.R
Rscript scripts/04_enrichment.R
Rscript scripts/05_survival.R


Reproduction

# 1. build the environment (installs all R/Bioconductor packages)
conda config --set channel_priority strict
conda env create -f environment.yml
conda activate tcga-brca-rnaseq

# 2. place TCGA-CDR-SupplementalTableS1.xlsx in data/raw/  (see script 01)

# 3. run the pipeline in order
Rscript scripts/01_data_acquisition.R
Rscript scripts/02_preprocess.R
Rscript scripts/03_deseq2_dea.R
Rscript scripts/04_enrichment.R
Rscript scripts/05_survival.R


Script 01 downloads from the GDC (open-access, cached and idempotent) and
KEGG enrichment in script 04 queries the KEGG web API, so an internet
connection is required for those steps. Package versions are recorded in
results/tables/*_session_info.txt for each stage.

Methods (brief)


Preprocessing: primary tumours only; one sample per patient (largest
library size); genes filtered to ≥10 counts in ≥190 samples (the smaller
group). Luminal A set as the DESeq2 reference level.
Differential expression: DESeq2 negative-binomial model on raw counts;
apeglm log-fold-change shrinkage; DEGs at padj < 0.05 and |log2FC| > 1.
Enrichment: Ensembl→ENTREZ mapping; GO:BP and KEGG over-representation
(background = all expressed, mapped genes) with GO redundancy collapsed via
simplify(); GSEA ranked by the Wald statistic.
Survival: endpoints from TCGA-CDR; Kaplan-Meier + log-rank; Cox
(unadjusted and age-adjusted) fitted on a shared complete-case set;
proportional-hazards checked with Schoenfeld residuals.



Limitations


Survival associations are non-significant trends; the crossing hazards and
the limited number of late events constrain both power and the validity of a
single hazard ratio.
Analysis inherits the GDC's fixed alignment and annotation (a deliberate
choice for a downstream-focused study).
Results describe the TCGA-BRCA cohort and are not independently validated in
an external dataset.



References


Liu et al. (2018) An Integrated TCGA Pan-Cancer Clinical Data Resource…
Cell 173(2):400–416. (TCGA-CDR; PFI/OS definitions.)
Ciriello et al. (2018) Comprehensive Molecular Portraits of Invasive Lobular
Breast Cancer. Cell (PAM50 subtype calls via TCGAbiolinks.)
Love, Huber & Anders (2014) Moderated estimation of fold change and
dispersion for RNA-seq data with DESeq2. Genome Biology.
Zhu, Ibrahim & Love (2019) Heavy-tailed prior distributions for sequence
count data (apeglm). Bioinformatics.
Wu et al. (2021) clusterProfiler 4.0. The Innovation.



Author

Clement Akinsola
