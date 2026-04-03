############################################################
## Centroid-based validation for DEAS-based 3-cluster solution
##
## Input files:
##   1) BRCA_as_sig.txt   : rows = AS events, cols = samples
##                          first column = AS event ID
##   2) ann_col_2.txt     : reference cluster assignment
##                          two columns, e.g. sub / Group
##
## Output folder:
##   centroid_validation_3clusters/
############################################################

## =========================
## 0. Package setup
## =========================
cran_pkgs <- c("data.table", "ggplot2", "pheatmap")
bioc_pkgs <- c("impute")

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
  library(ggplot2)
  library(pheatmap)
  library(impute)
})

## =========================
## 1. User settings
## =========================
expr_file <- "BRCA_as_sig.txt"
ann_file  <- "ann_col_2.txt"

outdir <- "centroid_validation_3clusters"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## split settings
n_splits <- 200
train_fraction <- 0.80
seed_base <- 12345

## PSI preprocessing
max_na_ratio <- 0.30
use_top_variable_features <- TRUE
top_n_features <- 1500

## centroid assignment
similarity_method <- "pearson"   # "pearson" or "spearman"

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

ann$Cluster <- gsub("^C", "", ann$Group, ignore.case = TRUE)
ann$Cluster <- as.integer(ann$Cluster)

if (any(is.na(ann$Cluster))) {
  stop("Cannot parse cluster labels in ann_col_2.txt. Expected values like C1/C2/C3.")
}

ann <- ann[!duplicated(ann$Sample), , drop = FALSE]

ref_labels <- ann$Cluster
names(ref_labels) <- ann$Sample

## =========================
## 4. Match samples
## =========================
common_samples <- intersect(colnames(expr_mat), names(ref_labels))

cat("Samples in expression matrix: ", ncol(expr_mat), "\n", sep = "")
cat("Samples in annotation file: ", length(ref_labels), "\n", sep = "")
cat("Common samples: ", length(common_samples), "\n", sep = "")

if (length(common_samples) < 30) {
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
  stop("Too few features remain after NA filtering.")
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
mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- sort(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

calc_centroids <- function(mat_train, labels_train, K = 3) {
  centroids <- matrix(NA_real_, nrow = nrow(mat_train), ncol = K,
                      dimnames = list(rownames(mat_train), paste0("C", 1:K)))
  for (k in 1:K) {
    idx <- which(labels_train == k)
    if (length(idx) == 0) stop(paste("No samples for cluster", k, "in training set."))
    centroids[, k] <- rowMeans(mat_train[, idx, drop = FALSE])
  }
  return(centroids)
}

assign_by_centroid <- function(mat_test, centroids, method = "pearson") {
  samples_test <- colnames(mat_test)
  K <- ncol(centroids)
  
  cor_mat <- matrix(NA_real_, nrow = length(samples_test), ncol = K,
                    dimnames = list(samples_test, colnames(centroids)))
  
  for (i in seq_along(samples_test)) {
    s <- samples_test[i]
    for (k in 1:K) {
      cor_mat[i, k] <- suppressWarnings(
        cor(mat_test[, s], centroids[, k], method = method)
      )
    }
  }
  
  pred <- apply(cor_mat, 1, function(x) {
    if (all(is.na(x))) return(NA_integer_)
    which.max(x)
  })
  
  pred <- as.integer(pred)
  names(pred) <- samples_test
  
  max_corr <- apply(cor_mat, 1, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })
  
  second_corr <- apply(cor_mat, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sort(x, decreasing = TRUE)[2]
  })
  
  margin <- max_corr - second_corr
  
  return(list(pred = pred, cor_mat = cor_mat, max_corr = max_corr, margin = margin))
}

calc_metrics <- function(true_labels, pred_labels, K = 3) {
  keep <- !(is.na(true_labels) | is.na(pred_labels))
  true_labels <- true_labels[keep]
  pred_labels <- pred_labels[keep]
  
  acc <- mean(true_labels == pred_labels)
  
  conf <- table(
    factor(true_labels, levels = 1:K),
    factor(pred_labels, levels = 1:K)
  )
  
  precision <- rep(NA_real_, K)
  recall <- rep(NA_real_, K)
  f1 <- rep(NA_real_, K)
  
  for (k in 1:K) {
    tp <- conf[k, k]
    fp <- sum(conf[, k]) - tp
    fn <- sum(conf[k, ]) - tp
    
    precision[k] <- ifelse((tp + fp) == 0, NA, tp / (tp + fp))
    recall[k]    <- ifelse((tp + fn) == 0, NA, tp / (tp + fn))
    f1[k] <- ifelse(is.na(precision[k]) | is.na(recall[k]) | (precision[k] + recall[k]) == 0,
                    NA, 2 * precision[k] * recall[k] / (precision[k] + recall[k]))
  }
  
  macro_f1 <- mean(f1, na.rm = TRUE)
  
  return(list(
    accuracy = acc,
    macro_f1 = macro_f1,
    confusion = conf,
    precision = precision,
    recall = recall,
    f1 = f1
  ))
}

## =========================
## 7. Repeated train-test split validation
## =========================
samples <- colnames(expr_mat)
n_samples <- length(samples)

prediction_mat <- matrix(
  NA_integer_,
  nrow = n_samples,
  ncol = n_splits,
  dimnames = list(samples, paste0("Split_", seq_len(n_splits)))
)

tested_mat <- matrix(
  FALSE,
  nrow = n_samples,
  ncol = n_splits,
  dimnames = list(samples, paste0("Split_", seq_len(n_splits)))
)

confidence_mat <- matrix(
  NA_real_,
  nrow = n_samples,
  ncol = n_splits,
  dimnames = list(samples, paste0("Split_", seq_len(n_splits)))
)

margin_mat <- matrix(
  NA_real_,
  nrow = n_samples,
  ncol = n_splits,
  dimnames = list(samples, paste0("Split_", seq_len(n_splits)))
)

split_summary <- data.frame(
  Split = seq_len(n_splits),
  N_train = NA_integer_,
  N_test = NA_integer_,
  Accuracy = NA_real_,
  MacroF1 = NA_real_,
  stringsAsFactors = FALSE
)

cat("Starting repeated train-test centroid validation...\n")

for (i in seq_len(n_splits)) {
  set.seed(seed_base + i)
  
  train_samples <- sort(sample(samples, size = floor(train_fraction * n_samples), replace = FALSE))
  test_samples  <- setdiff(samples, train_samples)
  test_samples  <- sort(test_samples)
  
  train_labels <- ref_labels[train_samples]
  test_labels  <- ref_labels[test_samples]
  
  ## training matrix and test matrix
  mat_train <- expr_mat[, train_samples, drop = FALSE]
  mat_test  <- expr_mat[, test_samples, drop = FALSE]
  
  ## derive centroids from training reference labels
  centroids <- calc_centroids(mat_train, train_labels, K = 3)
  
  ## assign test samples
  pred_res <- assign_by_centroid(mat_test, centroids, method = similarity_method)
  pred_labels <- pred_res$pred
  
  ## save predictions
  prediction_mat[names(pred_labels), i] <- pred_labels
  tested_mat[names(pred_labels), i] <- TRUE
  confidence_mat[names(pred_labels), i] <- pred_res$max_corr
  margin_mat[names(pred_labels), i] <- pred_res$margin
  
  ## evaluate
  met <- calc_metrics(test_labels[names(pred_labels)], pred_labels, K = 3)
  
  split_summary$N_train[i] <- length(train_samples)
  split_summary$N_test[i]  <- length(test_samples)
  split_summary$Accuracy[i] <- met$accuracy
  split_summary$MacroF1[i]  <- met$macro_f1
  
  if (i %% 10 == 0) {
    cat("Finished ", i, "/", n_splits, "\n", sep = "")
  }
}

## =========================
## 8. Patient-level prediction stability
## =========================
cat("Calculating patient-level prediction stability...\n")

patient_stability <- data.frame(
  Sample = samples,
  ReferenceCluster = ref_labels[samples],
  TimesTested = rowSums(tested_mat),
  MajorityPredictedCluster = apply(prediction_mat, 1, mode_value),
  MajorityProportion = NA_real_,
  AgreementWithReference = NA_real_,
  MeanConfidence = apply(confidence_mat, 1, function(x) mean(x, na.rm = TRUE)),
  MeanMargin = apply(margin_mat, 1, function(x) mean(x, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(patient_stability))) {
  x <- prediction_mat[i, ]
  x <- x[!is.na(x)]
  
  if (length(x) == 0) next
  
  tb <- table(x)
  patient_stability$MajorityProportion[i] <- max(tb) / sum(tb)
  patient_stability$AgreementWithReference[i] <- mean(x == patient_stability$ReferenceCluster[i])
}

patient_stability$Stable80 <- patient_stability$MajorityProportion >= 0.80
patient_stability$Stable90 <- patient_stability$MajorityProportion >= 0.90

## =========================
## 9. Cluster-wise accuracy summary
## =========================
cluster_accuracy <- aggregate(
  AgreementWithReference ~ ReferenceCluster,
  data = patient_stability,
  FUN = function(x) mean(x, na.rm = TRUE)
)
colnames(cluster_accuracy)[2] <- "MeanAgreement"

## =========================
## 10. Global summary
## =========================
global_summary <- list(
  total_samples = n_samples,
  total_features = nrow(expr_mat),
  n_splits = n_splits,
  train_fraction = train_fraction,
  mean_accuracy = mean(split_summary$Accuracy, na.rm = TRUE),
  median_accuracy = median(split_summary$Accuracy, na.rm = TRUE),
  sd_accuracy = sd(split_summary$Accuracy, na.rm = TRUE),
  mean_macroF1 = mean(split_summary$MacroF1, na.rm = TRUE),
  median_macroF1 = median(split_summary$MacroF1, na.rm = TRUE),
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
  split_summary,
  file = file.path(outdir, "split_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  patient_stability,
  file = file.path(outdir, "patient_prediction_stability.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  prediction_mat,
  file = file.path(outdir, "patient_predictions.tsv"),
  sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE
)

write.table(
  cluster_accuracy,
  file = file.path(outdir, "cluster_accuracy.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

sink(file.path(outdir, "global_summary.txt"))
cat("Global summary of centroid-based validation\n")
cat("==========================================\n")
for (nm in names(global_summary)) {
  cat(nm, ": ", global_summary[[nm]], "\n", sep = "")
}
sink()

## =========================
## 12. Save plots
## =========================
cat("Saving plots...\n")

## 12.1 Accuracy distribution
p1 <- ggplot(split_summary, aes(x = Accuracy)) +
  geom_histogram(bins = 30, fill = "#2c7fb8", color = "black") +
  theme_classic(base_size = 12) +
  labs(
    title = "Centroid-based train-test validation accuracy",
    x = "Prediction accuracy in test set",
    y = "Frequency"
  )

ggsave(
  filename = file.path(outdir, "accuracy_distribution.pdf"),
  plot = p1, width = 6, height = 4.5
)

## 12.2 Macro-F1 distribution
p2 <- ggplot(split_summary, aes(x = MacroF1)) +
  geom_histogram(bins = 30, fill = "#41ab5d", color = "black") +
  theme_classic(base_size = 12) +
  labs(
    title = "Centroid-based train-test macro-F1",
    x = "Macro-F1 in test set",
    y = "Frequency"
  )

ggsave(
  filename = file.path(outdir, "macroF1_distribution.pdf"),
  plot = p2, width = 6, height = 4.5
)

## 12.3 Patient-level stability histogram
p3 <- ggplot(patient_stability, aes(x = MajorityProportion)) +
  geom_histogram(bins = 30, fill = "#756bb1", color = "black") +
  theme_classic(base_size = 12) +
  labs(
    title = "Patient-level prediction stability",
    x = "Proportion of repeated predictions to majority cluster",
    y = "Number of patients"
  )

ggsave(
  filename = file.path(outdir, "patient_prediction_stability_histogram.pdf"),
  plot = p3, width = 6, height = 4.5
)

## 12.4 Stability by reference cluster
p4 <- ggplot(patient_stability, aes(x = factor(ReferenceCluster), y = MajorityProportion)) +
  geom_boxplot(fill = "grey85", color = "black") +
  theme_classic(base_size = 12) +
  labs(
    title = "Prediction stability by reference cluster",
    x = "Reference cluster",
    y = "Majority predicted cluster proportion"
  )

ggsave(
  filename = file.path(outdir, "prediction_stability_by_cluster.pdf"),
  plot = p4, width = 5.5, height = 4.5
)

## =========================
## 13. Console summary
## =========================
cat("\n===== GLOBAL SUMMARY =====\n")
print(global_summary)

cat("\nOutput directory: ", outdir, "\n", sep = "")
cat("Done.\n")