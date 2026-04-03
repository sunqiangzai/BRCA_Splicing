############################################################
## BRCA splicing cluster robustness + train/validation model
## Complete final script
##
## Input:
##   1) BRCA_as_sig.txt
##      - first column: AS event ID
##      - remaining columns: TCGA samples
##      - values: PSI
##
##   2) ann_col_2.txt
##      - first column: sample ID
##      - second column: Group (C1/C2/C3)
##
## Output folder:
##   paper_figures_and_results/
############################################################

############################
## 0. package setup
############################
cran_pkgs <- c(
  "data.table", "clue", "mclust", "ggplot2", "pheatmap",
  "caret", "glmnet", "xgboost", "pROC", "reshape2", "gridExtra"
)

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
  library(caret)
  library(glmnet)
  library(xgboost)
  library(pROC)
  library(reshape2)
  library(gridExtra)
})

############################
## 1. user settings
############################
expr_file <- "BRCA_as_sig.txt"
ann_file  <- "ann_col_2.txt"

outdir <- "paper_figures_and_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## general
seed_base <- 12345
set.seed(seed_base)

## PSI preprocessing
max_na_ratio <- 0.30
use_top_variable_features <- TRUE
top_n_features <- 1500

## robustness
K <- 3
n_resamples <- 200
sample_fraction_outer <- 0.80

## ConsensusClusterPlus
cc_reps <- 100
pItem_inner <- 0.80
pFeature_inner <- 1
clusterAlg <- "hc"
distance_method <- "pearson"
inner_linkage <- "average"

## train/validation split
train_fraction <- 0.80

############################
## 2. helper functions
############################
theme_paper <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- sort(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

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
  cl
}

align_labels_to_reference <- function(sub_labels, ref_labels, K = 3) {
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
  
  tab <- table(
    factor(as.integer(sub2), levels = 1:K),
    factor(as.integer(ref2), levels = 1:K)
  )
  
  cost <- max(tab) - tab
  perm <- clue::solve_LSAP(cost)
  perm <- as.integer(perm)
  
  mapped <- perm[as.integer(sub2)]
  names(mapped) <- names(sub2)
  mapped
}

save_plot_pdf_png <- function(plot_obj, filename_base, width = 6, height = 4.5, dpi = 300) {
  ggsave(
    filename = file.path(outdir, paste0(filename_base, ".pdf")),
    plot = plot_obj, width = width, height = height
  )
  ggsave(
    filename = file.path(outdir, paste0(filename_base, ".png")),
    plot = plot_obj, width = width, height = height, dpi = dpi
  )
}

############################
## 3. read expression matrix
############################
cat("Reading expression matrix...\n")
expr_df <- fread(expr_file, data.table = FALSE)

if (ncol(expr_df) < 3) {
  stop("Expression matrix format seems wrong.")
}

feature_ids <- expr_df[[1]]
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
rownames(expr_mat) <- feature_ids
mode(expr_mat) <- "numeric"

cat("Raw matrix: ", nrow(expr_mat), " features x ", ncol(expr_mat), " samples\n", sep = "")

############################
## 4. read annotation
############################
cat("Reading reference annotation...\n")
ann <- fread(ann_file, data.table = FALSE)

if (ncol(ann) < 2) {
  stop("ann_col_2.txt must have at least 2 columns.")
}

colnames(ann)[1:2] <- c("Sample", "Group")
ann$Sample <- as.character(ann$Sample)
ann$Group  <- as.character(ann$Group)
ann$Cluster <- as.integer(gsub("^C", "", ann$Group, ignore.case = TRUE))

if (any(is.na(ann$Cluster))) {
  stop("Cannot parse cluster labels. Expected C1/C2/C3.")
}

ann <- ann[!duplicated(ann$Sample), , drop = FALSE]
ref_labels <- ann$Cluster
names(ref_labels) <- ann$Sample

############################
## 5. match samples
############################
common_samples <- intersect(colnames(expr_mat), names(ref_labels))
cat("Samples in matrix: ", ncol(expr_mat), "\n", sep = "")
cat("Samples in annotation: ", length(ref_labels), "\n", sep = "")
cat("Common samples: ", length(common_samples), "\n", sep = "")

if (length(common_samples) < 10) {
  stop("Too few overlapping samples.")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
ref_labels <- ref_labels[common_samples]

############################
## 6. PSI preprocessing
############################
cat("PSI preprocessing...\n")

## keep rows with <=30% NA
na_ratio <- rowMeans(is.na(expr_mat))
expr_mat <- expr_mat[na_ratio <= max_na_ratio, , drop = FALSE]
cat("After NA filtering: ", nrow(expr_mat), " features\n", sep = "")

## KNN imputation
cat("Running KNN imputation...\n")
expr_mat <- impute::impute.knn(expr_mat)$data

## remove zero variance rows
row_sd <- apply(expr_mat, 1, sd)
expr_mat <- expr_mat[row_sd > 0, , drop = FALSE]
cat("After variance filtering: ", nrow(expr_mat), " features\n", sep = "")

## row-wise z-score
expr_mat <- t(scale(t(expr_mat)))
expr_mat[is.na(expr_mat)] <- 0

## select variable features
if (isTRUE(use_top_variable_features)) {
  mad_vec <- apply(expr_mat, 1, mad, na.rm = TRUE)
  ord <- order(mad_vec, decreasing = TRUE)
  top_n <- min(top_n_features, nrow(expr_mat))
  expr_mat <- expr_mat[ord[1:top_n], , drop = FALSE]
  cat("After MAD selection: ", nrow(expr_mat), " features\n", sep = "")
}

############################
## 7. save processed data
############################
write.table(
  expr_mat,
  file = file.path(outdir, "processed_expression_matrix.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

write.table(
  data.frame(Sample = names(ref_labels), ReferenceCluster = ref_labels),
  file = file.path(outdir, "reference_labels.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

############################
## 8. Figure 1: reference cluster sizes
############################
df_ref <- data.frame(Cluster = factor(ref_labels, levels = 1:3))
p_ref <- ggplot(df_ref, aes(x = Cluster)) +
  geom_bar(fill = "#4C78A8", color = "black") +
  labs(title = "Figure 1. Reference cluster sizes",
       x = "Reference cluster", y = "Number of samples") +
  theme_paper()

save_plot_pdf_png(p_ref, "Figure1_reference_cluster_sizes", 5, 4.2)

############################
## 9. robustness analysis
############################
cat("Starting robustness analysis...\n")

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
    title = file.path(outdir, paste0("tmp_ccp_run_", i))
  )
  
  sub_labels_aligned <- align_labels_to_reference(
    sub_labels = sub_labels_raw,
    ref_labels = ref_labels,
    K = K
  )
  
  common_ids2 <- sort(intersect(names(sub_labels_aligned), names(ref_labels)))
  sub_labels_aligned <- sub_labels_aligned[common_ids2]
  ref_sub2 <- ref_labels[common_ids2]
  
  assignment_mat[names(sub_labels_aligned), i] <- sub_labels_aligned
  
  resampling_summary$N_subsamples[i] <- length(sub_samples)
  resampling_summary$N_clusters_found[i] <- length(unique(sub_labels_aligned))
  resampling_summary$ARI_vs_reference[i] <- adjustedRandIndex(
    as.integer(sub_labels_aligned),
    as.integer(ref_sub2)
  )
  
  if (i %% 10 == 0) {
    cat("Finished ", i, "/", n_resamples, "\n", sep = "")
  }
}

write.table(
  resampling_summary,
  file = file.path(outdir, "resampling_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  assignment_mat,
  file = file.path(outdir, "patient_assignments.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

############################
## 10. patient stability
############################
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

write.table(
  patient_stability,
  file = file.path(outdir, "patient_stability.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

############################
## 11. pairwise co-clustering
############################
cat("Calculating co-clustering matrix...\n")

co_occurrence <- matrix(0, n_samples, n_samples, dimnames = list(samples, samples))
co_cluster    <- matrix(0, n_samples, n_samples, dimnames = list(samples, samples))

for (i in seq_len(n_resamples)) {
  labs <- assignment_mat[, i]
  idx <- which(!is.na(labs))
  if (length(idx) < 2) next
  
  sub_names <- names(labs)[idx]
  sub_labs <- labs[idx]
  
  co_occurrence[sub_names, sub_names] <- co_occurrence[sub_names, sub_names] + 1
  same_cluster <- outer(sub_labs, sub_labs, FUN = "==") * 1
  co_cluster[sub_names, sub_names] <- co_cluster[sub_names, sub_names] + same_cluster
}

co_clustering_freq <- co_cluster / ifelse(co_occurrence == 0, NA, co_occurrence)

write.table(
  co_clustering_freq,
  file = file.path(outdir, "pairwise_co_clustering_frequency.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

############################
## 12. Figure 2: ARI distribution
############################
p_ari <- ggplot(resampling_summary, aes(x = ARI_vs_reference)) +
  geom_histogram(bins = 30, fill = "#2C7FB8", color = "black") +
  labs(title = "Figure 2. ARI distribution across resampling runs",
       x = "Adjusted Rand Index vs reference",
       y = "Frequency") +
  theme_paper()

save_plot_pdf_png(p_ari, "Figure2_ARI_distribution", 6, 4.5)

############################
## 13. Figure 3: patient stability histogram
############################
p_stab_hist <- ggplot(patient_stability, aes(x = MajorityProportion)) +
  geom_histogram(bins = 30, fill = "#41AB5D", color = "black") +
  labs(title = "Figure 3. Patient-level assignment stability",
       x = "Proportion assigned to majority cluster",
       y = "Number of patients") +
  theme_paper()

save_plot_pdf_png(p_stab_hist, "Figure3_patient_stability_histogram", 6, 4.5)

############################
## 14. Figure 4: stability by cluster
############################
p_stab_box <- ggplot(patient_stability,
                     aes(x = factor(ReferenceCluster), y = MajorityProportion)) +
  geom_boxplot(fill = "grey85", color = "black") +
  labs(title = "Figure 4. Patient stability by reference cluster",
       x = "Reference cluster", y = "Majority assignment proportion") +
  theme_paper()

save_plot_pdf_png(p_stab_box, "Figure4_patient_stability_by_cluster", 5.5, 4.5)

############################
## 15. Figure 5: co-clustering heatmap
############################
ord_idx <- order(ref_labels[samples], samples)
ord_samples <- samples[ord_idx]
co_clustering_ord <- co_clustering_freq[ord_samples, ord_samples]

ann_heatmap <- data.frame(ReferenceCluster = factor(ref_labels[ord_samples]))
rownames(ann_heatmap) <- ord_samples

pdf(file.path(outdir, "Figure5_co_clustering_heatmap.pdf"), width = 7.5, height = 7.5)
pheatmap(
  co_clustering_ord,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = ann_heatmap,
  annotation_col = ann_heatmap,
  main = "Figure 5. Pairwise co-clustering frequency",
  na_col = "white"
)
dev.off()

png(file.path(outdir, "Figure5_co_clustering_heatmap.png"), width = 2200, height = 2200, res = 300)
pheatmap(
  co_clustering_ord,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = ann_heatmap,
  annotation_col = ann_heatmap,
  main = "Figure 5. Pairwise co-clustering frequency",
  na_col = "white"
)
dev.off()

############################
## 16. train/validation model
############################
cat("Building train/validation model...\n")

X <- t(expr_mat)
y <- factor(ref_labels[colnames(expr_mat)], levels = c(1, 2, 3))

set.seed(seed_base)
train_idx <- createDataPartition(y, p = train_fraction, list = FALSE)

X_train <- X[train_idx, , drop = FALSE]
X_valid <- X[-train_idx, , drop = FALSE]
y_train <- y[train_idx]
y_valid <- y[-train_idx]

write.table(
  data.frame(Sample = rownames(X_train), Cluster = y_train),
  file = file.path(outdir, "training_samples.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  data.frame(Sample = rownames(X_valid), Cluster = y_valid),
  file = file.path(outdir, "validation_samples.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

############################
## 17. LASSO feature selection
############################
cat("Running LASSO feature selection...\n")

cvfit <- cv.glmnet(
  x = as.matrix(X_train),
  y = y_train,
  family = "multinomial",
  alpha = 1,
  nfolds = 5,
  type.measure = "class"
)

coef_list <- coef(cvfit, s = "lambda.min")

selected_features <- unique(unlist(lapply(coef_list, function(x) {
  rownames(x)[which(as.matrix(x) != 0)]
})))
selected_features <- setdiff(selected_features, "(Intercept)")

cat("Selected features by LASSO: ", length(selected_features), "\n", sep = "")

if (length(selected_features) < 2) {
  stop("Too few features selected by LASSO. Check data or use lambda.1se/lambda.min carefully.")
}

write.table(
  data.frame(Feature = selected_features),
  file = file.path(outdir, "lasso_selected_features.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

X_train_sel <- X_train[, selected_features, drop = FALSE]
X_valid_sel <- X_valid[, selected_features, drop = FALSE]

############################
## 18. XGBoost model
############################
cat("Training XGBoost model...\n")

y_train_num <- as.numeric(y_train) - 1
y_valid_num <- as.numeric(y_valid) - 1

dtrain <- xgb.DMatrix(data = as.matrix(X_train_sel), label = y_train_num)
dvalid <- xgb.DMatrix(data = as.matrix(X_valid_sel), label = y_valid_num)

params <- list(
  objective = "multi:softprob",
  num_class = 3,
  eval_metric = "mlogloss",
  eta = 0.1,
  max_depth = 4,
  subsample = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 1
)

xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 200,
  watchlist = list(train = dtrain, valid = dvalid),
  early_stopping_rounds = 20,
  verbose = 0
)

xgb.save(xgb_model, fname = file.path(outdir, "xgboost_multiclass_model.model"))

############################
## 19. validation predictions
############################
cat("Evaluating validation set...\n")

pred_prob <- predict(xgb_model, dvalid)
pred_matrix <- matrix(pred_prob, ncol = 3, byrow = TRUE)
colnames(pred_matrix) <- c("C1", "C2", "C3")
rownames(pred_matrix) <- rownames(X_valid_sel)

pred_class <- max.col(pred_matrix)
pred_factor <- factor(pred_class, levels = c(1, 2, 3))

validation_results <- data.frame(
  Sample = rownames(X_valid_sel),
  TrueCluster = y_valid,
  PredictedCluster = pred_factor,
  Prob_C1 = pred_matrix[, 1],
  Prob_C2 = pred_matrix[, 2],
  Prob_C3 = pred_matrix[, 3],
  stringsAsFactors = FALSE
)

write.table(
  validation_results,
  file = file.path(outdir, "validation_predictions.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cm <- confusionMatrix(pred_factor, y_valid)
accuracy <- mean(pred_factor == y_valid)

############################
## 20. Figure 6: confusion matrix
############################
cm_df <- as.data.frame(cm$table)
colnames(cm_df) <- c("Prediction", "Reference", "Freq")

p_cm <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "white", high = "#2166AC") +
  labs(title = "Figure 6. Validation confusion matrix",
       x = "Reference cluster", y = "Predicted cluster") +
  theme_paper()

save_plot_pdf_png(p_cm, "Figure6_validation_confusion_matrix", 5.5, 4.8)

############################
## 21. Figure 7: multiclass ROC
############################
roc_list <- list()
auc_values <- numeric(3)

pdf(file.path(outdir, "Figure7_multiclass_ROC.pdf"), width = 6, height = 6)
plot(
  0, 0,
  type = "n", xlim = c(1, 0), ylim = c(0, 1),
  xlab = "False Positive Rate", ylab = "True Positive Rate",
  main = "Figure 7. One-vs-rest ROC curves"
)
abline(0, 1, lty = 2, col = "grey50")

line_types <- c(1, 1, 1)
line_cols  <- c("red", "blue", "darkgreen")

for (i in 1:3) {
  truth_i <- as.numeric(y_valid == levels(y_valid)[i])
  roc_obj <- roc(truth_i, pred_matrix[, i], quiet = TRUE)
  roc_list[[i]] <- roc_obj
  auc_values[i] <- as.numeric(auc(roc_obj))
  lines(1 - roc_obj$specificities, roc_obj$sensitivities,
        lwd = 2, lty = line_types[i], col = line_cols[i])
}

legend(
  "bottomright",
  legend = c(
    paste0("C1 AUC = ", sprintf("%.3f", auc_values[1])),
    paste0("C2 AUC = ", sprintf("%.3f", auc_values[2])),
    paste0("C3 AUC = ", sprintf("%.3f", auc_values[3]))
  ),
  col = line_cols, lwd = 2, bty = "n"
)
dev.off()

png(file.path(outdir, "Figure7_multiclass_ROC.png"), width = 1800, height = 1800, res = 300)
plot(
  0, 0,
  type = "n", xlim = c(1, 0), ylim = c(0, 1),
  xlab = "False Positive Rate", ylab = "True Positive Rate",
  main = "Figure 7. One-vs-rest ROC curves"
)
abline(0, 1, lty = 2, col = "grey50")
for (i in 1:3) {
  lines(1 - roc_list[[i]]$specificities, roc_list[[i]]$sensitivities,
        lwd = 2, lty = line_types[i], col = line_cols[i])
}
legend(
  "bottomright",
  legend = c(
    paste0("C1 AUC = ", sprintf("%.3f", auc_values[1])),
    paste0("C2 AUC = ", sprintf("%.3f", auc_values[2])),
    paste0("C3 AUC = ", sprintf("%.3f", auc_values[3]))
  ),
  col = line_cols, lwd = 2, bty = "n"
)
dev.off()

############################
## 22. Figure 8: feature importance
############################
imp <- xgb.importance(feature_names = selected_features, model = xgb_model)
imp_top <- head(imp, 20)

write.table(
  imp,
  file = file.path(outdir, "xgboost_feature_importance.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

p_imp <- ggplot(imp_top, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "#B2182B", color = "black") +
  coord_flip() +
  labs(title = "Figure 8. Top XGBoost features",
       x = "Feature", y = "Gain") +
  theme_paper()

save_plot_pdf_png(p_imp, "Figure8_top_xgboost_features", 6.5, 5.5)

############################
## 23. summary stats
############################
global_summary <- list(
  total_samples = n_samples,
  total_features_after_preprocessing = nrow(expr_mat),
  n_resamples = n_resamples,
  mean_ARI = mean(resampling_summary$ARI_vs_reference, na.rm = TRUE),
  median_ARI = median(resampling_summary$ARI_vs_reference, na.rm = TRUE),
  sd_ARI = sd(resampling_summary$ARI_vs_reference, na.rm = TRUE),
  proportion_runs_with_3_clusters = mean(resampling_summary$N_clusters_found == 3, na.rm = TRUE),
  mean_patient_stability = mean(patient_stability$MajorityProportion, na.rm = TRUE),
  median_patient_stability = median(patient_stability$MajorityProportion, na.rm = TRUE),
  proportion_patients_stable80 = mean(patient_stability$Stable80, na.rm = TRUE),
  proportion_patients_stable90 = mean(patient_stability$Stable90, na.rm = TRUE),
  training_samples = nrow(X_train_sel),
  validation_samples = nrow(X_valid_sel),
  lasso_selected_features = length(selected_features),
  validation_accuracy = accuracy,
  auc_C1 = auc_values[1],
  auc_C2 = auc_values[2],
  auc_C3 = auc_values[3],
  mean_auc = mean(auc_values)
)

sink(file.path(outdir, "global_summary.txt"))
cat("Global summary\n")
cat("==============\n")
for (nm in names(global_summary)) {
  cat(nm, ": ", global_summary[[nm]], "\n", sep = "")
}
sink()

############################
## 24. figure legends draft
############################
figure_legends <- c(
  "Figure 1. Distribution of reference clusters in the TCGA cohort.",
  "Figure 2. Distribution of adjusted Rand index (ARI) values across repeated 80% resampling runs. Higher ARI indicates better recovery of the original three-cluster solution.",
  "Figure 3. Histogram of patient-level cluster assignment stability across repeated resampling.",
  "Figure 4. Patient-level stability stratified by reference cluster.",
  "Figure 5. Pairwise co-clustering frequency heatmap across repeated resampling. Strong block structure supports the robustness of the three-cluster solution.",
  "Figure 6. Confusion matrix of the XGBoost classifier in the 20% validation set.",
  "Figure 7. One-vs-rest ROC curves for the three cluster classes in the validation set.",
  "Figure 8. Top-ranked features contributing to the XGBoost classifier after LASSO feature selection."
)

writeLines(figure_legends, con = file.path(outdir, "figure_legends.txt"))

############################
## 25. console summary
############################
cat("\n===== ANALYSIS FINISHED =====\n")
print(global_summary)
cat("\nResults saved in: ", outdir, "\n", sep = "")