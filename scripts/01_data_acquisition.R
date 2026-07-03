################################################################################
# 01_data_acquisition.R
#
# Project : TCGA-BRCA RNA-seq — Basal-like vs Luminal A
# Purpose : Acquire and assemble the raw inputs for the whole analysis:
#             (1) STAR gene-level RAW COUNTS from the GDC (open-access)
#             (2) PAM50 molecular subtype labels
#             (3) Survival endpoints (TCGA-CDR):
#                   PFI / PFI.time  -> PRIMARY endpoint (modelled in script 05)
#                   OS  / OS.time   -> SECONDARY, descriptive only (see note)
#             (4) Clinical covariates: age, gender, race, tumour stage,
#                 vital_status
#           Merge (2)-(4) into the SummarizedExperiment's colData and save an
#           annotated object that every downstream script reads.
#
# Author  : Clement Akinsola (Ceefor Analytic Hub)
#
# Inputs  : data/raw/TCGA-CDR-SupplementalTableS1.xlsx   (download once — see NOTE)
# Outputs : data/processed/brca_se_raw_annotated.rds     (annotated SE, raw counts)
#           data/processed/sample_metadata.csv           (flat per-sample table)
#           results/tables/acquisition_session_info.txt  (reproducibility record)
#
# WHY STAR counts from the GDC (not self-quantification):
#   TCGA FASTQ/BAM are dbGaP controlled-access. The STAR gene counts are the
#   open-access, reproducible product of GDC's fixed pipeline (GENCODE v36,
#   GRCh38). We inherit their alignment/annotation on purpose — this project is
#   about rigorous DOWNSTREAM analysis; the Nextflow project handles FASTQ->counts.
#
# WHY PFI is primary and OS is descriptive only:
#   BRCA has favourable prognosis, so deaths are sparse within TCGA's follow-up
#   window -> OS is underpowered for modelling. We MODEL PFI, but still carry OS
#   to DESCRIBE that sparsity (event counts, follow-up), which documents and
#   justifies the choice of PFI as the primary endpoint. Acquire generously,
#   model selectively.
#
# NOTE — the survival endpoints and clinical covariates are NOT in the GDC
#   clinical download. Download the TCGA Pan-Cancer Clinical Data Resource
#   supplement once and place it in data/raw/:
#     Liu et al., 2018, Cell 173(2):400-416  (PMID: 29625055)
#     File: TCGA-CDR-SupplementalTableS1.xlsx  (sheet "TCGA-CDR")
################################################################################


# ---- 0. Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(TCGAbiolinks)          # GDC query/download/prepare
  library(SummarizedExperiment)  # SE accessors: assay(), colData(), rowData()
  library(dplyr)
  library(readxl)                # read the CDR .xlsx
  library(here)                  # project-root-relative paths (reproducible)
})


# ---- 1. Paths & config ------------------------------------------------------
project      <- "TCGA-BRCA"

dir_raw      <- here("data", "raw")
dir_proc     <- here("data", "processed")
dir_gdc      <- here("data", "raw", "GDCdata")
dir_tables   <- here("results", "tables")

for (d in c(dir_raw, dir_proc, dir_gdc, dir_tables)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cdr_file     <- file.path(dir_raw, "TCGA-CDR-SupplementalTableS1.xlsx")


# ---- 2. Query, download & prepare STAR raw counts ---------------------------
message("Querying GDC for ", project, " STAR - Counts ...")
query_exp <- GDCquery(
  project        = project,
  data.category  = "Transcriptome Profiling",
  data.type      = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type  = "STAR - Counts"
)

GDCdownload(
  query          = query_exp,
  method         = "api",
  directory      = dir_gdc,
  files.per.chunk = 20
)

message("Assembling SummarizedExperiment (this can take a few minutes) ...")
se <- GDCprepare(
  query     = query_exp,
  directory = dir_gdc,
  summarizedExperiment = TRUE
)

stopifnot("unstranded" %in% assayNames(se))
message("SE dims: ", nrow(se), " genes x ", ncol(se), " samples")
message("Assays present: ", paste(assayNames(se), collapse = ", "))


# ---- 3. Derive sample-level keys --------------------------------------------
# Barcode: TCGA-A7-A0CE-01A-...  -> chars 1-12 patient, chars 14-15 sample type.
# Sample type codes: 01 primary tumour, 06 metastatic, 11 normal.
cd <- as.data.frame(colData(se))
cd$barcode          <- colnames(se)
cd$patient_barcode  <- substr(cd$barcode, 1, 12)
cd$sample_type_code <- substr(cd$barcode, 14, 15)

message("Sample types (14-15): ",
        paste(names(table(cd$sample_type_code)),
              table(cd$sample_type_code), sep = "=", collapse = "  "))


# ---- 4. PAM50 subtype labels ------------------------------------------------
# Current PAM50 column is 'BRCA_Subtype_PAM50' (older 'PAM50.mRNA' is gone) —
# resolve defensively rather than hard-coding.
message("Fetching PAM50 subtypes ...")
subt <- TCGAquery_subtype(tumor = "brca")

pam50_col <- intersect(
  c("BRCA_Subtype_PAM50", "PAM50.mRNA", "Subtype_mRNA"),
  colnames(subt)
)[1]
if (is.na(pam50_col)) {
  stop("Could not find a PAM50 column in TCGAquery_subtype output. ",
       "Inspect colnames(subt) and update pam50_col.")
}

subt_patient_col <- intersect(c("patient", "sample", "Patient_ID"),
                              colnames(subt))[1]

subt_slim <- subt %>%
  transmute(
    patient_barcode = .data[[subt_patient_col]],
    PAM50           = .data[[pam50_col]]
  ) %>%
  filter(!is.na(PAM50), PAM50 != "") %>%
  distinct(patient_barcode, .keep_all = TRUE)

message("PAM50 labels retrieved for ", nrow(subt_slim), " patients.")
print(table(subt_slim$PAM50, useNA = "ifany"))


# ---- 5. Survival endpoints + clinical covariates (TCGA-CDR) -----------------  [UPDATED]
# Endpoints:
#   PFI / PFI.time -> PRIMARY: time to first new-tumour event OR death; else
#                     censored. Chosen over OS (sparse deaths) and DFI (needs a
#                     curated disease-free start).
#   OS  / OS.time  -> SECONDARY (descriptive only): event = death (any cause).
#                     Carried to summarise event sparsity, not to model.
# Covariates:
#   age    - age at diagnosis (continuous; modelled as the Cox adjustment)
#   stage  - AJCC stage (kept for documentation; not modelled — high missingness)
#   race   - kept for documentation
#   gender - kept to document cohort; ~99% female so not a useful model term
#   vital_status - Alive/Dead descriptor (no time attached)
if (!file.exists(cdr_file)) {
  stop("TCGA-CDR file not found at:\n  ", cdr_file,
       "\nDownload TCGA-CDR-SupplementalTableS1.xlsx (Liu et al. 2018, ",
       "PMID 29625055) and place it in data/raw/. See header NOTE.")
}

message("Reading curated survival endpoints + covariates from TCGA-CDR ...")
cdr <- read_excel(cdr_file, sheet = "TCGA-CDR")

# TCGA encodes missingness with several text tokens; map them all to NA BEFORE
# any numeric coercion (also removes benign 'NAs introduced by coercion' warnings).
na_tokens <- c("#N/A", "[Not Available]", "[Unknown]", "[Not Applicable]",
               "[Not Evaluated]", "[Discrepancy]", "NA", "")
clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% na_tokens] <- NA
  x
}

cdr_brca <- cdr %>%
  filter(type == "BRCA") %>%
  transmute(
    patient_barcode = bcr_patient_barcode,
    # --- primary endpoint ---
    PFI          = as.integer(clean_chr(PFI)),
    PFI.time     = as.numeric(clean_chr(PFI.time)),
    # --- secondary endpoint (descriptive) ---
    OS           = as.integer(clean_chr(OS)),
    OS.time      = as.numeric(clean_chr(OS.time)),
    vital_status = clean_chr(vital_status),
    # --- covariates ---
    age          = as.numeric(clean_chr(age_at_initial_pathologic_diagnosis)),
    gender       = clean_chr(gender),
    race         = clean_chr(race),
    stage        = clean_chr(ajcc_pathologic_tumor_stage)
  ) %>%
  distinct(patient_barcode, .keep_all = TRUE)

message("PFI available for ", sum(!is.na(cdr_brca$PFI)), " BRCA patients; ",
        "PFI events (progression/death): ", sum(cdr_brca$PFI == 1, na.rm = TRUE))
message("OS available for ", sum(!is.na(cdr_brca$OS)), " patients; ",
        "OS events (deaths): ", sum(cdr_brca$OS == 1, na.rm = TRUE),
        "  <- note the sparsity: this is why PFI is the primary endpoint")
message("Age recorded for ", sum(!is.na(cdr_brca$age)), " patients; ",
        "stage recorded for ", sum(!is.na(cdr_brca$stage)), ".")


# ---- 6. Merge subtype + endpoints + covariates into colData -----------------
# GOTCHA: GDCprepare already attaches the GDC's own indexed clinical fields to
# colData, some named identically to our CDR columns (race, gender,
# vital_status, ...). Joining on top would trigger dplyr's .x/.y suffixing and
# leave no clean column name. We PREFER the curated CDR versions, so drop the
# GDC duplicates from cd before joining — the CDR columns then keep clean names.
cdr_value_cols <- c("PFI", "PFI.time", "OS", "OS.time", "vital_status",
                    "age", "gender", "race", "stage")
dup_cols <- intersect(cdr_value_cols, colnames(cd))
if (length(dup_cols) > 0) {
  message("Dropping GDC colData columns superseded by curated CDR: ",
          paste(dup_cols, collapse = ", "))
  cd <- cd %>% select(-all_of(dup_cols))
}

cd_annot <- cd %>%
  left_join(subt_slim, by = "patient_barcode") %>%
  left_join(cdr_brca,  by = "patient_barcode")

stopifnot(identical(cd_annot$barcode, colnames(se)))
rownames(cd_annot) <- cd_annot$barcode
colData(se) <- DataFrame(cd_annot)

message("PAM50 x sample_type coverage:")
print(table(PAM50 = colData(se)$PAM50,
            type  = colData(se)$sample_type_code, useNA = "ifany"))

message("Age (years) summary:")
print(summary(colData(se)$age))
message("Vital status:")
print(table(colData(se)$vital_status, useNA = "ifany"))
message("Stage distribution (raw AJCC):")
print(table(colData(se)$stage, useNA = "ifany"))


# ---- 7. Save outputs --------------------------------------------------------  [UPDATED select]
saveRDS(se, file = file.path(dir_proc, "brca_se_raw_annotated.rds"))

meta_out <- as.data.frame(colData(se)) %>%
  select(barcode, patient_barcode, sample_type_code, PAM50,
         PFI, PFI.time, OS, OS.time, vital_status,
         age, gender, race, stage)
write.csv(meta_out,
          file = file.path(dir_proc, "sample_metadata.csv"),
          row.names = FALSE)

message("Saved: brca_se_raw_annotated.rds  and  sample_metadata.csv")


# ---- 8. Reproducibility record ----------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(dir_tables, "acquisition_session_info.txt"))

message("Done. Next: 02_preprocess.R (subset to LumA/Basal primary tumours, ",
        "filter low-count genes, build the DESeqDataSet).")
