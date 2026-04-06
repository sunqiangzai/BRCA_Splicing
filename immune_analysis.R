############################################################
## Immune infiltration analysis using CIBERSORT + ESTIMATE
##
## Input:
##   1) expression_matrix.txt
##      first column = GeneSymbol
##      remaining columns = samples
##
##   2) sample_subtype.txt
##      columns: Sample, Cluster
##
##   3) CIBERSORT_results.txt
##      downloaded from CIBERSORT web portal
##
## Output:
##   immune_analysis/
############################################################

## =========================
## 0. Package setup
## =========================
cran_pkgs <- c("data.table", "ggplot2", "reshape2", "dplyr")
bioc_pkgs <- c("estimate")

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
  library(reshape2)
  library(dplyr)
  library(estimate)
})

## =========================
## 1. User settings
## =========================
expr_file <- "expression_matrix.txt"
subtype_file <- "sample_subtype.txt"
cibersort_file <- "CIBERSORT_results.txt"

outdir <- "immune_analysis"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cluster_levels <- c("C1", "C2", "C3")

## =========================
## 2. Helper
## =========================
normalize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]", "", x)
  x
}

## =========================
## 3. Read expression matrix
## =========================
cat("Reading expression matrix...\n")
expr_df <- fread(expr_file, data.table = FALSE, check.names = FALSE)

if (ncol(expr_df) < 3) {
  stop("Expression matrix must contain 1 gene column + >=2 sample columns.")
}

gene_col <- colnames(expr_df)[1]
expr_df[[gene_col]] <- as.character(expr_df[[gene_col]])

## remove duplicated gene symbols by keeping max mean expression
expr_mat0 <- expr_df[, -1, drop = FALSE]
mode(expr_mat0) <- "numeric"
expr_df$mean_expr_tmp <- rowMeans(expr_mat0, na.rm = TRUE)
expr_df <- expr_df[order(expr_df[[gene_col]], -expr_df$mean_expr_tmp), ]
expr_df <- expr_df[!duplicated(expr_df[[gene_col]]), ]
expr_df$mean_expr_tmp <- NULL

## =========================
## 4. Prepare CIBERSORT input
## =========================
cat("Writing CIBERSORT input file...\n")
write.table(
  expr_df,
  file = file.path(outdir, "CIBERSORT_input.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## =========================
## 5. Run ESTIMATE
## =========================
cat("Running ESTIMATE...\n")

estimate_input <- file.path(outdir, "estimate_input.gct")
estimate_common <- file.path(outdir, "estimate_common.gct")
estimate_score_file <- file.path(outdir, "estimate_scores.gct")

## ESTIMATE requires GCT-like input
## format first two columns: NAME, Description
estimate_df <- expr_df
colnames(estimate_df)[1] <- "NAME"
estimate_df$Description <- estimate_df$NAME
estimate_df <- estimate_df[, c("NAME", "Description", setdiff(colnames(estimate_df), c("NAME", "Description")))]

## write GCT manually
write_gct <- function(df, file) {
  con <- file(file, open = "wt")
  on.exit(close(con))
  writeLines("#1.2", con)
  writeLines(paste(nrow(df), ncol(df) - 2, sep = "\t"), con)
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_gct(estimate_df, estimate_input)

filterCommonGenes(
  input.f = estimate_input,
  output.f = estimate_common,
  id = "GeneSymbol"
)

estimateScore(
  input.ds = estimate_common,
  output.ds = estimate_score_file,
  platform = "illumina"
)

## read ESTIMATE output
read_estimate_gct <- function(file) {
  x <- read.table(file, skip = 2, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  return(x)
}

estimate_scores <- read_estimate_gct(estimate_score_file)

## rows are scores
score_df <- estimate_scores[, -2, drop = FALSE]   # remove Description
rownames(score_df) <- score_df[, 1]
score_df <- score_df[, -1, drop = FALSE]

score_df <- as.data.frame(t(score_df), stringsAsFactors = FALSE)
score_df$Sample <- rownames(score_df)
rownames(score_df) <- NULL

for (cc in c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")) {
  if (cc %in% colnames(score_df)) {
    score_df[[cc]] <- suppressWarnings(as.numeric(score_df[[cc]]))
  }
}

## =========================
## 6. Read subtype info
## =========================
cat("Reading subtype file...\n")
subtype_df <- fread(subtype_file, data.table = FALSE)

if (!all(c("Sample", "Cluster") %in% colnames(subtype_df))) {
  stop("sample_subtype.txt must contain columns: Sample, Cluster")
}

subtype_df$Sample_raw <- as.character(subtype_df$Sample)
subtype_df$Sample_norm <- normalize_id(subtype_df$Sample_raw)
subtype_df$Cluster <- factor(as.character(subtype_df$Cluster), levels = cluster_levels)

## =========================
## 7. Read CIBERSORT output
## =========================
cat("Reading CIBERSORT results...\n")
ciber_df <- fread(cibersort_file, data.table = FALSE, check.names = FALSE)

## first column usually sample ID
sample_id_col <- colnames(ciber_df)[1]
ciber_df$Sample_raw <- as.character(ciber_df[[sample_id_col]])
ciber_df$Sample_norm <- normalize_id(ciber_df$Sample_raw)

## common immune cell columns in LM22 output
immune_cols <- setdiff(
  colnames(ciber_df),
  c(sample_id_col, "P-value", "Correlation", "RMSE", "Sample_raw", "Sample_norm")
)

## =========================
## 8. Merge data
## =========================
score_df$Sample_norm <- normalize_id(score_df$Sample)

score_merge <- merge(
  subtype_df[, c("Sample_raw", "Sample_norm", "Cluster"), drop = FALSE],
  score_df,
  by = "Sample_norm",
  all.x = FALSE,
  all.y = FALSE
)

ciber_merge <- merge(
  subtype_df[, c("Sample_raw", "Sample_norm", "Cluster"), drop = FALSE],
  ciber_df,
  by = "Sample_norm",
  all.x = FALSE,
  all.y = FALSE
)

cat("Samples matched for ESTIMATE:", nrow(score_merge), "\n")
cat("Samples matched for CIBERSORT:", nrow(ciber_merge), "\n")

## save merged files
write.table(score_merge,
            file = file.path(outdir, "ESTIMATE_scores_merged.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(ciber_merge,
            file = file.path(outdir, "CIBERSORT_merged.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## =========================
## 9. Statistics for ESTIMATE scores
## =========================
cat("Running Kruskal-Wallis tests for ESTIMATE scores...\n")

estimate_stats <- data.frame(
  Variable = character(),
  P.Value = numeric(),
  stringsAsFactors = FALSE
)

estimate_vars <- intersect(c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"), colnames(score_merge))

for (v in estimate_vars) {
  tmp <- score_merge[!is.na(score_merge[[v]]) & !is.na(score_merge$Cluster), ]
  if (nrow(tmp) > 0) {
    kw <- kruskal.test(tmp[[v]] ~ tmp$Cluster)
    estimate_stats <- rbind(estimate_stats,
                            data.frame(Variable = v, P.Value = kw$p.value, stringsAsFactors = FALSE))
  }
}

write.table(estimate_stats,
            file = file.path(outdir, "ESTIMATE_KruskalWallis.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## =========================
## 10. Boxplots for ESTIMATE scores
## =========================
cluster_colors <- c("C1" = "#ff2d20", "C2" = "#4c67af", "C3" = "#16a34a")

for (v in estimate_vars) {
  tmp <- score_merge[!is.na(score_merge[[v]]) & !is.na(score_merge$Cluster), ]
  if (nrow(tmp) == 0) next
  
  pval_show <- estimate_stats$P.Value[match(v, estimate_stats$Variable)]
  
  p <- ggplot(tmp, aes(x = Cluster, y = .data[[v]], fill = Cluster)) +
    geom_boxplot(width = 0.65, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
    scale_fill_manual(values = cluster_colors) +
    theme_classic(base_size = 12) +
    labs(title = v, x = "Subtype", y = v) +
    annotate("text",
             x = 2,
             y = max(tmp[[v]], na.rm = TRUE) * 1.05,
             label = paste0("Kruskal-Wallis P = ", signif(pval_show, 3)),
             size = 4) +
    theme(legend.position = "none")
  
  ggsave(file.path(outdir, paste0(v, "_boxplot.pdf")), p, width = 4.5, height = 4.2)
}

## =========================
## 11. Statistics for CIBERSORT cell fractions
## =========================
cat("Running Kruskal-Wallis tests for CIBERSORT fractions...\n")

ciber_stats <- data.frame(
  CellType = character(),
  P.Value = numeric(),
  stringsAsFactors = FALSE
)

for (cell in immune_cols) {
  tmp <- ciber_merge[!is.na(ciber_merge[[cell]]) & !is.na(ciber_merge$Cluster), ]
  if (nrow(tmp) > 0) {
    kw <- kruskal.test(tmp[[cell]] ~ tmp$Cluster)
    ciber_stats <- rbind(ciber_stats,
                         data.frame(CellType = cell, P.Value = kw$p.value, stringsAsFactors = FALSE))
  }
}

ciber_stats$adj.P.Val <- p.adjust(ciber_stats$P.Value, method = "BH")

write.table(ciber_stats,
            file = file.path(outdir, "CIBERSORT_KruskalWallis.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## =========================
## 12. Boxplot for all immune cell fractions
## =========================
plot_df <- ciber_merge[, c("Cluster", immune_cols), drop = FALSE]
plot_df_long <- melt(plot_df, id.vars = "Cluster", variable.name = "CellType", value.name = "Fraction")

## add significance labels
plot_df_long$CellType <- factor(plot_df_long$CellType, levels = immune_cols)

sig_map <- ciber_stats[, c("CellType", "adj.P.Val"), drop = FALSE]
sig_map$label <- ifelse(sig_map$adj.P.Val < 0.001, "***",
                        ifelse(sig_map$adj.P.Val < 0.01, "**",
                               ifelse(sig_map$adj.P.Val < 0.05, "*", "ns")))

ymax_df <- plot_df_long %>%
  group_by(CellType) %>%
  summarise(ymax = max(Fraction, na.rm = TRUE), .groups = "drop")

sig_map <- merge(sig_map, ymax_df, by = "CellType", all.x = TRUE)

p_cells <- ggplot(plot_df_long, aes(x = CellType, y = Fraction, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = cluster_colors) +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "CIBERSORT Estimated Cell Fraction") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

## add significance above each cell type
for (i in seq_len(nrow(sig_map))) {
  p_cells <- p_cells + annotate(
    "text",
    x = which(levels(plot_df_long$CellType) == as.character(sig_map$CellType[i])),
    y = sig_map$ymax[i] + 0.03,
    label = sig_map$label[i],
    size = 3
  )
}

ggsave(file.path(outdir, "CIBERSORT_all_celltypes_boxplot.pdf"),
       p_cells, width = 12, height = 4.8)

## =========================
## 13. Optional stacked barplot of mean fractions
## =========================
mean_frac <- plot_df_long %>%
  group_by(Cluster, CellType) %>%
  summarise(mean_fraction = mean(Fraction, na.rm = TRUE), .groups = "drop")

p_bar <- ggplot(mean_frac, aes(x = Cluster, y = mean_fraction, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 12) +
  labs(y = "Mean estimated fraction", x = "Subtype")

ggsave(file.path(outdir, "CIBERSORT_mean_fraction_stackedbar.pdf"),
       p_bar, width = 6.5, height = 5)

## =========================
## 14. Console summary
## =========================
cat("Done.\n")
cat("Results saved in folder:", outdir, "\n")
cat("CIBERSORT input file written to:", file.path(outdir, "CIBERSORT_input.txt"), "\n")