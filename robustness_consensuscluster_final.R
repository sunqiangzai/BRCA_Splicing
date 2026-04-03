############################################################
## Robustness test for 3-cluster solution (PSI matrix)
## Final integrated script
##
## Input files:
##   1) BRCA_as_sig.txt   : rows = AS events, cols = samples
##                          first column = AS event ID
##   2) ann_col_2.txt     : reference cluster assignment
##                          two columns, e.g. sub / Group
##
## Output folder:
##   robustness_3clusters_final/
############################################################

## =========================
## 0. Package setup
## =========================
cran_pkgs <- c("data.table", "clue", "mclust", "ggplot2", "pheatmap")
bioc_pkgs <- c("ConsensusClusterPlus", "impute")

install_if_missing_cran <- function(pkgs) {
  to_install <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(to_install) > 0) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
}

install_if_missing_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      BiocManager::install(p, ask = FALSE, update = FALSE)
    }
  }
}

install_if_missing_cran(cran_pkgs)
install_if_missing_bioc(bioc_pkgs)

suppressPackageStartupMessages({
  library(data.table)
  library(ConsensusClusterPlus)
  library(clue)
  library(mclust)
  library(ggplot2)
  library(pheatmap)
  library(impute)
})

## =========================
## 1. User settings
## =========================
expr_file <- "BRCA_as_sig.txt"
ann_file  <- "ann_col_2.txt"

outdir <- "robustness_3clusters_final"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## outer resampling settings
n_resamples <- 200
sample_fraction_outer <- 0.80

## inner ConsensusClusterPlus settings
K <- 3
cc_reps <- 100
pItem_inner <- 0.80
pFeature_inner <- 1
clusterAlg <- "hc"
distance_method <- "pearson"
inner_linkage <- "average"
seed_base <- 12345

## PSI filtering settings
max_na_ratio <- 0.30     # keep events with <=30% NA
use_top_variable_features <- TRUE
top_n_features <- 1500   # set to NULL or FALSE above if not needed

## =========================
## 2. Read expression matrix
## =========================
cat("Reading expression matrix...\n")
expr_df <- fread(expr_file, data.table = FALSE)

if (ncol(expr_df) < 3) {
  stop("Expression matrix format seems wrong: need at least 1 feature column + 2 sample columns.")
}

feature_ids <- expr_df[[1]]
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
rownames(expr_mat) <- feature_ids
mode(expr_mat) <- "numeric"

cat("Raw expression matrix dimension: ",
    nrow(expr_mat), " features x ", ncol(expr_mat), " samples\n", sep = "")

## =========================
## 3. Read annotation / reference labels
## =========================
cat("Reading reference annotation...\n")
ann <- fread(ann_file, data.table = FALSE)

if (ncol(ann) < 2) {
  stop("ann_col_2.txt must contain at least 2 columns.")
}

colnames(ann)[1:2] <- c("Sample", "Group")
ann$Sample <- as.character(ann$Sample)
ann$Group  <- as.character(ann$Group)

## Parse C1/C2/C3 -> 1/2/3
ann$Cluster <- gsub("^C", "", ann$Group, ignore.case = TRUE)
ann$Cluster <- as.integer(ann$Cluster)

if (any(is.na(ann$Cluster))) {
  stop("Cannot parse cluster labels in ann_col_2.txt. Expected values like C1/C2/C3.")
}

## remove duplicated samples if any
ann <- ann[!duplicated(ann$Sample), , drop = FALSE]

ref_labels <- ann$Cluster
names(ref_labels) <- ann$Sample

## =========================
## 4. Match samples between matrix and annotation
## =========================
common_samples <- intersect(colnames(expr_mat), names(ref_labels))

cat("Samples in expression matrix: ", ncol(expr_mat), "\n", sep = "")
cat("Samples in annotation file: ", length(ref_labels), "\n", sep = "")
cat("Common samples: ", length(common_samples), "\n", sep = "")

if (length(common_samples) < 10) {
  stop("Too few overlapping samples between BRCA_as_sig.txt and ann_col_2.txt.")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
ref_labels <- ref_labels[common_samples]

cat("Reference cluster sizes after matching:\n")
print(table(ref_labels))

## =========================
## 5. PSI-specific preprocessing
## =========================
cat("PSI preprocessing...\n")

## 5.1 remove rows with too many NA
na_ratio <- rowMeans(is.na(expr_mat))
expr_mat <- expr_mat[na_ratio <= max_na_ratio, , drop = FALSE]
cat("Features retained after NA filtering: ", nrow(expr_mat), "\n", sep = "")

if (nrow(expr_mat) < 10) {
  stop("Too few features remain after NA filtering. Consider relaxing max_na_ratio.")
}

## 5.2 KNN imputation
cat("Performing KNN imputation...\n")
expr_mat <- impute::impute.knn(expr_mat)$data

## 5.3 remove zero-variance rows
row_sd <- apply(expr_mat, 1, sd)
expr_mat <- expr_mat[row_sd > 0, , drop = FALSE]
cat("Features retained after variance filtering: ", nrow(expr_mat), "\n", sep = "")

if (nrow(expr_mat) < 10) {
  stop("Too few features remain after variance filtering.")
}

## 5.4 row-wise z-score scaling
cat("Performing row-wise z-score scaling...\n")
expr_mat <- t(scale(t(expr_mat)))
expr_mat[is.na(expr_mat)] <- 0

## 5.5 optional MAD selection
if (isTRUE(use_top_variable_features)) {
  cat("Selecting top variable features by MAD...\n")
  mad_vec <- apply(expr_mat, 1, mad, na.rm = TRUE)
  ord <- order(mad_vec, decreasing = TRUE)
  top_n <- min(top_n_features, nrow(expr_mat))
  expr_mat <- expr_mat[ord[1:top_n], , drop = FALSE]
  cat("Features retained after MAD selection: ", nrow(expr_mat), "\n", sep = "")
}

## =========================
## 6. Helper functions
## =========================
run_ccp_k3 <- function(mat_sub, sub_samples, K = 3, reps = 100,
                       pItem = 0.8, pFeature = 1,
                       clusterAlg = "hc", distance = "pearson",
                       inner_linkage = "average",
                       seed = 1, title = tempdir()) {
  
  if (ncol(mat_sub) != length(sub_samples)) {
    stop("ncol(mat_sub) is not equal to length(sub_samples).")
  }
  
  res <- ConsensusClusterPlus(
    d = mat_sub,
    maxK = K,
    reps = reps,
    pItem = pItem,
    pFeature = pFeature,
    clusterAlg = clusterAlg,
    distance = distance,
    innerLinkage = inner_linkage,
    finalLinkage = inner_linkage,
    seed = seed,
    plot = NULL,
    title = title,
    verbose = FALSE
  )
  
  cl <- res[[K]]$consensusClass
  
  if (length(cl) != length(sub_samples)) {
    stop(
      paste0("ConsensusClusterPlus returned ", length(cl),
             " labels, but sub_samples has ", length(sub_samples), " samples.")
    )
  }
  
  cl <- as.integer(cl)
  names(cl) <- sub_samples
  return(cl)
}

align_labels_to_reference <- function(sub_labels, ref_labels, K = 3) {
  if (is.null(names(sub_labels))) stop("sub_labels has no names.")
  if (is.null(names(ref_labels))) stop("ref_labels has no names.")
  
  common_ids <- intersect(names(sub_labels), names(ref_labels))
  if (length(common_ids) == 0) {
    stop("No overlapping samples between sub_labels and ref_labels.")
  }
  
  sub2 <- sub_labels[common_ids]
  ref2 <- ref_labels[common_ids]
  
  sub2 <- sub2[order(names(sub2))]
  ref2 <- ref2[order(names(ref2))]
  
  if (!identical(names(sub2), names(ref2))) {
    stop("Sample order mismatch after matching.")
  }
  
  if (length(sub2) != length(ref2)) {
    stop("sub2 and ref2 still have different lengths.")
  }
  
  tab <- table(
    factor(as.integer(sub2), levels = 1:K),
    factor(as.integer(ref2), levels = 1:K)
  )
  
  cost <- max(tab) - tab
  perm <- clue::solve_LSAP(cost)
  perm <- as.integer(perm)
  
  mapped <- perm[as.integer(sub2)]
  names(mapped) <- names(sub2)
  
  return(mapped)
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- sort(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

## =========================
## 7. Outer resampling
## =========================
samples <- colnames(expr_mat)
n_samples <- length(samples)

assignment_mat <- matrix(
  NA_integer_,
  nrow = n_samples,
  ncol = n_resamples,
  dimnames = list(samples, paste0("Resample_", seq_len(n_resamples)))
)

resampling_summary <- data.frame(
  Resample = seq_len(n_resamples),
  N_subsamples = NA_integer_,
  N_clusters_found = NA_integer_,
  ARI_vs_reference = NA_real_,
  stringsAsFactors = FALSE
)

cat("Starting outer resampling...\n")

for (i in seq_len(n_resamples)) {
  set.seed(seed_base + i)
  
  sub_samples <- sort(sample(samples, size = floor(sample_fraction_outer * n_samples), replace = FALSE))
  sub_mat <- expr_mat[, sub_samples, drop = FALSE]
  
  sub_labels_raw <- run_ccp_k3(
    mat_sub = sub_mat,
    sub_samples = sub_samples,
    K = K,
    reps = cc_reps,
    pItem = pItem_inner,
    pFeature = pFeature_inner,
    clusterAlg = clusterAlg,
    distance = distance_method,
    inner_linkage = inner_linkage,
    seed = seed_base + i,
    title = file.path(outdir, paste0("tmp_run_", i))
  )
  
  sub_labels_aligned <- align_labels_to_reference(
    sub_labels = sub_labels_raw,
    ref_labels = ref_labels,
    K = K
  )
  
  common_ids2 <- intersect(names(sub_labels_aligned), names(ref_labels))
  common_ids2 <- sort(common_ids2)
  
  sub_labels_aligned <- sub_labels_aligned[common_ids2]
  ref_sub2 <- ref_labels[common_ids2]
  
  assignment_mat[names(sub_labels_aligned), i] <- sub_labels_aligned
  
  resampling_summary$N_subsamples[i] <- length(sub_samples)
  resampling_summary$N_clusters_found[i] <- length(unique(sub_labels_aligned))
  resampling_summary$ARI_vs_reference[i] <- mclust::adjustedRandIndex(
    as.integer(sub_labels_aligned),
    as.integer(ref_sub2)
  )
  
  if (i %% 10 == 0) {
    cat("Finished ", i, "/", n_resamples, "\n", sep = "")
  }
}

## =========================
## 8. Patient-level stability
## =========================
cat("Calculating patient-level stability...\n")

patient_stability <- data.frame(
  Sample = samples,
  ReferenceCluster = ref_labels[samples],
  TimesSampled = apply(!is.na(assignment_mat), 1, sum),
  MajorityCluster = apply(assignment_mat, 1, mode_value),
  MajorityProportion = NA_real_,
  AgreementWithReference = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(patient_stability))) {
  x <- assignment_mat[i, ]
  x <- x[!is.na(x)]
  
  if (length(x) == 0) next
  
  tb <- table(x)
  patient_stability$MajorityProportion[i] <- max(tb) / sum(tb)
  patient_stability$AgreementWithReference[i] <- mean(x == patient_stability$ReferenceCluster[i])
}

patient_stability$Stable80 <- patient_stability$MajorityProportion >= 0.80
patient_stability$Stable90 <- patient_stability$MajorityProportion >= 0.90

## =========================
## 9. Pairwise co-clustering matrix
## =========================
cat("Calculating pairwise co-clustering matrix...\n")

co_occurrence <- matrix(0, n_samples, n_samples, dimnames = list(samples, samples))
co_cluster    <- matrix(0, n_samples, n_samples, dimnames = list(samples, samples))

for (i in seq_len(n_resamples)) {
  labs <- assignment_mat[, i]
  idx <- which(!is.na(labs))
  
  if (length(idx) < 2) next
  
  sub_names <- names(labs)[idx]
  sub_labs  <- labs[idx]
  
  co_occurrence[sub_names, sub_names] <- co_occurrence[sub_names, sub_names] + 1
  
  same_cluster <- outer(sub_labs, sub_labs, FUN = "==") * 1
  co_cluster[sub_names, sub_names] <- co_cluster[sub_names, sub_names] + same_cluster
}

co_clustering_freq <- co_cluster / ifelse(co_occurrence == 0, NA, co_occurrence)

## order heatmap by reference cluster
ord_idx <- order(ref_labels[samples], samples)
ord_samples <- samples[ord_idx]
co_clustering_ord <- co_clustering_freq[ord_samples, ord_samples]

ann_heatmap <- data.frame(ReferenceCluster = factor(ref_labels[ord_samples]))
rownames(ann_heatmap) <- ord_samples

## =========================
## 10. Global summary
## =========================
global_summary <- list(
  total_samples = n_samples,
  total_features = nrow(expr_mat),
  n_resamples = n_resamples,
  sample_fraction_outer = sample_fraction_outer,
  mean_ARI = mean(resampling_summary$ARI_vs_reference, na.rm = TRUE),
  median_ARI = median(resampling_summary$ARI_vs_reference, na.rm = TRUE),
  sd_ARI = sd(resampling_summary$ARI_vs_reference, na.rm = TRUE),
  proportion_runs_with_3_clusters = mean(resampling_summary$N_clusters_found == 3),
  mean_patient_stability = mean(patient_stability$MajorityProportion, na.rm = TRUE),
  median_patient_stability = median(patient_stability$MajorityProportion, na.rm = TRUE),
  proportion_patients_stable80 = mean(patient_stability$Stable80, na.rm = TRUE),
  proportion_patients_stable90 = mean(patient_stability$Stable90, na.rm = TRUE)
)

## =========================
## 11. Save tables
## =========================
cat("Saving tables...\n")

write.table(
  resampling_summary,
  file = file.path(outdir, "resampling_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  patient_stability,
  file = file.path(outdir, "patient_stability.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  assignment_mat,
  file = file.path(outdir, "patient_assignments.tsv"),
  sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE
)

write.table(
  co_clustering_freq,
  file = file.path(outdir, "pairwise_co_clustering_frequency.tsv"),
  sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE
)

sink(file.path(outdir, "global_summary.txt"))
cat("Global summary of robustness analysis\n")
cat("=====================================\n")
for (nm in names(global_summary)) {
  cat(nm, ": ", global_summary[[nm]], "\n", sep = "")
}
sink()

## =========================
## 12. Save plots
## =========================
cat("Saving plots...\n")

## 12.1 ARI distribution
p1 <- ggplot(resampling_summary, aes(x = ARI_vs_reference)) +
  geom_histogram(bins = 30, fill = "#2c7fb8", color = "black") +
  theme_classic(base_size = 12) +
  labs(
    title = "Robustness of 3-cluster solution",
    x = "Adjusted Rand Index vs reference clustering",
    y = "Frequency"
  )

ggsave(
  filename = file.path(outdir, "ARI_distribution.pdf"),
  plot = p1, width = 6, height = 4.5
)

## 12.2 Patient stability histogram
p2 <- ggplot(patient_stability, aes(x = MajorityProportion)) +
  geom_histogram(bins = 30, fill = "#41ab5d", color = "black") +
  theme_classic(base_size = 12) +
  labs(
    title = "Patient-level assignment stability",
    x = "Proportion of repeated assignments to majority cluster",
    y = "Number of patients"
  )

ggsave(
  filename = file.path(outdir, "patient_stability_histogram.pdf"),
  plot = p2, width = 6, height = 4.5
)

## 12.3 Patient stability by reference cluster
p3 <- ggplot(patient_stability, aes(x = factor(ReferenceCluster), y = MajorityProportion)) +
  geom_boxplot(fill = "grey85", color = "black") +
  theme_classic(base_size = 12) +
  labs(
    title = "Patient stability by reference cluster",
    x = "Reference cluster",
    y = "Majority assignment proportion"
  )

ggsave(
  filename = file.path(outdir, "patient_stability_by_cluster.pdf"),
  plot = p3, width = 5.5, height = 4.5
)

## 12.4 Co-clustering heatmap
pdf(file.path(outdir, "co_clustering_heatmap.pdf"), width = 7, height = 7)
pheatmap(
  co_clustering_ord,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = ann_heatmap,
  annotation_col = ann_heatmap,
  main = "Pairwise co-clustering frequency across resamples",
  na_col = "white"
)
dev.off()

## =========================
## 13. Console summary
## =========================
cat("\n===== GLOBAL SUMMARY =====\n")
print(global_summary)

cat("\nOutput directory: ", outdir, "\n", sep = "")
cat("Done.\n")