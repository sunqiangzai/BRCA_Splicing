############################################################
## Correlation analysis between splicing factor expression
## and DEAS PSI values
##
## Input:
##   1) DEAS_PSI_matrix.txt
##      rows = DEAS events
##      cols = samples
##      first column = event ID
##
##   2) SF_expression_matrix.txt
##      rows = splicing factors
##      cols = samples
##      first column = gene symbol
##
## Output:
##   SF_DEAS_correlation_all.tsv
##   SF_DEAS_correlation_significant.tsv
##   SF_summary_by_event_count.tsv
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
psi_file <- "DEAS_PSI_matrix.txt"
sf_expr_file <- "SF_expression_matrix.txt"

out_all <- "SF_DEAS_correlation_all.tsv"
out_sig <- "SF_DEAS_correlation_significant.tsv"
out_sf_summary <- "SF_summary_by_event_count.tsv"

r_cutoff <- 0.4
p_cutoff <- 0.05

## =========================
## 2. Read DEAS PSI matrix
## =========================
cat("Reading DEAS PSI matrix...\n")
psi_df <- fread(psi_file, data.table = FALSE, check.names = FALSE)

if (ncol(psi_df) < 3) {
  stop("DEAS_PSI_matrix.txt must contain 1 event column + >=2 sample columns.")
}

psi_event_col <- colnames(psi_df)[1]
psi_ids <- as.character(psi_df[[psi_event_col]])

psi_mat <- as.matrix(psi_df[, -1, drop = FALSE])
rownames(psi_mat) <- psi_ids
mode(psi_mat) <- "numeric"

cat("PSI matrix dimension: ",
    nrow(psi_mat), " events x ", ncol(psi_mat), " samples\n", sep = "")

## =========================
## 3. Read SF expression matrix
## =========================
cat("Reading SF expression matrix...\n")
sf_df <- fread(sf_expr_file, data.table = FALSE, check.names = FALSE)

if (ncol(sf_df) < 3) {
  stop("SF_expression_matrix.txt must contain 1 gene column + >=2 sample columns.")
}

sf_gene_col <- colnames(sf_df)[1]
sf_names <- as.character(sf_df[[sf_gene_col]])

sf_mat <- as.matrix(sf_df[, -1, drop = FALSE])
rownames(sf_mat) <- sf_names
mode(sf_mat) <- "numeric"

cat("SF expression matrix dimension: ",
    nrow(sf_mat), " SFs x ", ncol(sf_mat), " samples\n", sep = "")

## =========================
## 4. Match common samples
## =========================
common_samples <- intersect(colnames(psi_mat), colnames(sf_mat))

cat("Common samples: ", length(common_samples), "\n", sep = "")

if (length(common_samples) < 10) {
  stop("Too few common samples between PSI and SF expression matrices.")
}

psi_mat <- psi_mat[, common_samples, drop = FALSE]
sf_mat <- sf_mat[, common_samples, drop = FALSE]

## =========================
## 5. Correlation analysis
## =========================
cat("Running Pearson correlation analysis...\n")

res_list <- vector("list", nrow(sf_mat) * nrow(psi_mat))
idx <- 1

for (i in seq_len(nrow(sf_mat))) {
  sf_name <- rownames(sf_mat)[i]
  sf_vec <- as.numeric(sf_mat[i, ])
  
  for (j in seq_len(nrow(psi_mat))) {
    event_id <- rownames(psi_mat)[j]
    psi_vec <- as.numeric(psi_mat[j, ])
    
    ## keep complete cases
    keep <- complete.cases(sf_vec, psi_vec)
    x <- sf_vec[keep]
    y <- psi_vec[keep]
    
    n_used <- length(x)
    
    if (n_used < 5 || sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
      r_val <- NA
      p_val <- NA
    } else {
      ct <- tryCatch(
        cor.test(x, y, method = "pearson"),
        error = function(e) NULL
      )
      
      if (is.null(ct)) {
        r_val <- NA
        p_val <- NA
      } else {
        r_val <- unname(ct$estimate)
        p_val <- ct$p.value
      }
    }
    
    res_list[[idx]] <- data.frame(
      SF = sf_name,
      EventID = event_id,
      R = r_val,
      P.Value = p_val,
      N = n_used,
      stringsAsFactors = FALSE
    )
    
    idx <- idx + 1
  }
  
  if (i %% 10 == 0) {
    cat("Processed ", i, "/", nrow(sf_mat), " splicing factors\n", sep = "")
  }
}

res <- do.call(rbind, res_list)

## =========================
## 6. Filter significant pairs
## =========================
sig_res <- res[
  !is.na(res$R) &
    !is.na(res$P.Value) &
    abs(res$R) > r_cutoff &
    res$P.Value < p_cutoff,
  ,
  drop = FALSE
]

sig_res$Direction <- ifelse(sig_res$R > 0, "Positive", "Negative")

cat("Significant SF-DEAS pairs: ", nrow(sig_res), "\n", sep = "")

## =========================
## 7. Summarize by SF
## =========================
if (nrow(sig_res) > 0) {
  sf_summary <- aggregate(
    EventID ~ SF,
    data = sig_res,
    FUN = function(x) length(unique(x))
  )
  colnames(sf_summary)[2] <- "N_significant_DEAS"
  
  sf_summary$N_positive <- sapply(sf_summary$SF, function(sf) {
    sum(sig_res$SF == sf & sig_res$R > 0, na.rm = TRUE)
  })
  
  sf_summary$N_negative <- sapply(sf_summary$SF, function(sf) {
    sum(sig_res$SF == sf & sig_res$R < 0, na.rm = TRUE)
  })
  
  sf_summary <- sf_summary[order(sf_summary$N_significant_DEAS, decreasing = TRUE), ]
} else {
  sf_summary <- data.frame(
    SF = character(),
    N_significant_DEAS = integer(),
    N_positive = integer(),
    N_negative = integer(),
    stringsAsFactors = FALSE
  )
}

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

write.table(
  sf_summary,
  file = out_sf_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Saved all results to: ", out_all, "\n", sep = "")
cat("Saved significant results to: ", out_sig, "\n", sep = "")
cat("Saved SF summary to: ", out_sf_summary, "\n", sep = "")
cat("Done.\n")