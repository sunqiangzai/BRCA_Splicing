############################################################
## Detect DEAS from PSI matrix (paired tumor vs normal)
##
## Input:
##   1) PSI matrix file:
##      rows = AS events
##      cols = samples
##      first column = event ID
##
##   2) Sample annotation file:
##      must contain:
##         Sample     sample ID matching PSI matrix columns
##         Group      Tumor / Normal
##         PairID     matched pair ID
##
## Output:
##   DEAS_results.tsv
##   DEAS_significant.tsv
############################################################

## =========================
## 0. Package setup
## =========================
cran_pkgs <- c("data.table")

install_if_missing_cran <- function(pkgs) {
  to_install <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(to_install) > 0) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
}

install_if_missing_cran(cran_pkgs)

suppressPackageStartupMessages({
  library(data.table)
})

## =========================
## 1. User settings
## =========================
psi_file <- "PSI_matrix.txt"        # 第一列为event ID，后面为样本
meta_file <- "sample_annotation.txt" # 列: Sample, Group, PairID

out_all <- "DEAS_results.tsv"
out_sig <- "DEAS_significant.tsv"

## filtering thresholds
min_mean_psi <- 0.05
min_detect_rate <- 0.75

## DEAS thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1

## =========================
## 2. Read PSI matrix
## =========================
cat("Reading PSI matrix...\n")
psi_df <- fread(psi_file, data.table = FALSE, check.names = FALSE)

if (ncol(psi_df) < 3) {
  stop("PSI matrix format error: need 1 event column + >=2 sample columns.")
}

event_col <- colnames(psi_df)[1]
event_ids <- psi_df[[event_col]]

psi_mat <- as.matrix(psi_df[, -1, drop = FALSE])
rownames(psi_mat) <- event_ids
mode(psi_mat) <- "numeric"

cat("PSI matrix dimension: ",
    nrow(psi_mat), " events x ", ncol(psi_mat), " samples\n", sep = "")

## =========================
## 3. Read metadata
## =========================
cat("Reading sample annotation...\n")
meta <- fread(meta_file, data.table = FALSE)

required_cols <- c("Sample", "Group", "PairID")
if (!all(required_cols %in% colnames(meta))) {
  stop("sample_annotation.txt must contain columns: Sample, Group, PairID")
}

meta$Sample <- as.character(meta$Sample)
meta$Group <- as.character(meta$Group)
meta$PairID <- as.character(meta$PairID)

## standardize group labels
meta$Group <- ifelse(tolower(meta$Group) %in% c("tumor", "tumour"), "Tumor",
                     ifelse(tolower(meta$Group) %in% c("normal", "nat", "adjacent_normal"), "Normal", meta$Group))

if (!all(meta$Group %in% c("Tumor", "Normal"))) {
  stop("Group column must contain only Tumor / Normal")
}

## keep only samples in PSI matrix
meta <- meta[meta$Sample %in% colnames(psi_mat), , drop = FALSE]

cat("Samples in metadata matched to PSI matrix: ", nrow(meta), "\n", sep = "")

## =========================
## 4. Keep complete pairs
## =========================
pair_tab <- table(meta$PairID, meta$Group)

valid_pairs <- rownames(pair_tab)[pair_tab[, "Tumor"] >= 1 & pair_tab[, "Normal"] >= 1]

meta2 <- meta[meta$PairID %in% valid_pairs, , drop = FALSE]

## if multiple tumor or normal per pair, keep first of each
meta2 <- meta2[order(meta2$PairID, meta2$Group, meta2$Sample), , drop = FALSE]

tumor_meta <- meta2[meta2$Group == "Tumor", , drop = FALSE]
normal_meta <- meta2[meta2$Group == "Normal", , drop = FALSE]

tumor_meta <- tumor_meta[!duplicated(tumor_meta$PairID), , drop = FALSE]
normal_meta <- normal_meta[!duplicated(normal_meta$PairID), , drop = FALSE]

common_pairs <- intersect(tumor_meta$PairID, normal_meta$PairID)

tumor_meta <- tumor_meta[match(common_pairs, tumor_meta$PairID), , drop = FALSE]
normal_meta <- normal_meta[match(common_pairs, normal_meta$PairID), , drop = FALSE]

if (nrow(tumor_meta) != nrow(normal_meta)) {
  stop("Tumor/Normal pairing failed.")
}

cat("Number of valid matched pairs: ", length(common_pairs), "\n", sep = "")

tumor_samples <- tumor_meta$Sample
normal_samples <- normal_meta$Sample

## =========================
## 5. Filter PSI events
## =========================
cat("Filtering events...\n")

pair_samples <- c(tumor_samples, normal_samples)
psi_sub <- psi_mat[, pair_samples, drop = FALSE]

## average PSI across all paired samples
mean_psi <- rowMeans(psi_sub, na.rm = TRUE)

## detection rate = proportion of non-NA values
detect_rate <- rowMeans(!is.na(psi_sub))

keep_idx <- (mean_psi >= min_mean_psi) & (detect_rate >= min_detect_rate)
psi_filt <- psi_sub[keep_idx, , drop = FALSE]

cat("Events retained after filtering: ", nrow(psi_filt), "\n", sep = "")

if (nrow(psi_filt) == 0) {
  stop("No events passed filtering thresholds.")
}

## =========================
## 6. Detect DEAS
## =========================
cat("Running paired Wilcoxon tests...\n")

res_list <- vector("list", nrow(psi_filt))

for (i in seq_len(nrow(psi_filt))) {
  x_tumor <- as.numeric(psi_filt[i, tumor_samples])
  x_normal <- as.numeric(psi_filt[i, normal_samples])
  
  ## keep complete paired observations
  complete_idx <- complete.cases(x_tumor, x_normal)
  x_tumor2 <- x_tumor[complete_idx]
  x_normal2 <- x_normal[complete_idx]
  
  n_pairs_used <- length(x_tumor2)
  
  if (n_pairs_used < 3) {
    pval <- NA
    mean_tumor <- mean(x_tumor2, na.rm = TRUE)
    mean_normal <- mean(x_normal2, na.rm = TRUE)
    log2fc <- NA
  } else {
    wt <- tryCatch(
      wilcox.test(x_tumor2, x_normal2, paired = TRUE, exact = FALSE),
      error = function(e) NULL
    )
    
    pval <- if (is.null(wt)) NA else wt$p.value
    
    mean_tumor <- mean(x_tumor2, na.rm = TRUE)
    mean_normal <- mean(x_normal2, na.rm = TRUE)
    
    ## avoid log2(0)
    log2fc <- log2((mean_tumor + 1e-6) / (mean_normal + 1e-6))
  }
  
  res_list[[i]] <- data.frame(
    EventID = rownames(psi_filt)[i],
    mean_tumor = mean_tumor,
    mean_normal = mean_normal,
    delta_PSI = mean_tumor - mean_normal,
    log2FC = log2fc,
    P.Value = pval,
    n_pairs_used = n_pairs_used,
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, res_list)

## =========================
## 7. Multiple testing correction and significance
## =========================
res$adj.P.Val <- p.adjust(res$P.Value, method = "BH")

res$change <- ifelse(
  !is.na(res$adj.P.Val) & res$adj.P.Val < padj_cutoff & !is.na(res$log2FC) & res$log2FC > log2fc_cutoff, "UP",
  ifelse(
    !is.na(res$adj.P.Val) & res$adj.P.Val < padj_cutoff & !is.na(res$log2FC) & res$log2FC < -log2fc_cutoff, "DOWN",
    "NS"
  )
)

sig_res <- res[res$change != "NS", , drop = FALSE]

cat("Significant DEAS events: ", nrow(sig_res), "\n", sep = "")

## =========================
## 8. Save results
## =========================
write.table(
  res,
  file = out_all,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sig_res,
  file = out_sig,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Saved all results to: ", out_all, "\n", sep = "")
cat("Saved significant DEAS to: ", out_sig, "\n", sep = "")
cat("Done.\n")