#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript filter_psi.R <input.psi> <output.psi>")
}

input_file <- args[1]
output_file <- args[2]

cat("Reading PSI file...\n")

psi <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Total events:", nrow(psi), "\n")
cat("Total samples:", ncol(psi)-1, "\n")

# 第一列是 event_id
event_id <- psi[,1]
psi_mat <- psi[,-1]

# 处理 nan
psi_mat[psi_mat == "nan"] <- NA

# 转 numeric
psi_mat <- as.data.frame(lapply(psi_mat, function(x) as.numeric(x)))

# ----------------------------
# 1. 检测率（非NA比例）
# ----------------------------
detect_ratio <- rowSums(!is.na(psi_mat)) / ncol(psi_mat)

# ----------------------------
# 2. 平均 PSI
# ----------------------------
mean_psi <- rowMeans(psi_mat, na.rm = TRUE)

# ----------------------------
# 3. 过滤条件
# ----------------------------
keep <- detect_ratio >= 0.75 & mean_psi >= 0.05

psi_filtered <- psi[keep, ]

cat("After filtering:\n")
cat("Remaining events:", nrow(psi_filtered), "\n")
cat("Removed events:", nrow(psi) - nrow(psi_filtered), "\n")

# 写出
psi_filtered_out <- data.frame(
  AS_ID = rownames(psi_filtered),
  psi_filtered,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 写出
write.table(
  psi_filtered_out,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


cat("Output saved to:", output_file, "\n")