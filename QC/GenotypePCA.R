# 1. Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript generate_pca_plink2.R <plink_prefix> <phenotype_file.csv> <color_column>\n")
  cat("Example: Rscript generate_pca_plink2.R my_data sample_info.csv Group\n")
  quit(status = 1)
}

prefix     <- args[1]
pheno_file <- args[2]
color_col  <- args[3]

# 2. Check for PLINK2 and install required R packages
if (Sys.which("plink2") == "") {
  stop("Error: 'plink2' command not found. Please ensure plink2 is installed and in your system PATH.")
}

cran_packages <- c("ggplot2", "dplyr")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# 3. Perform LD Pruning using plink2
cat("\nPerforming LD pruning with plink2 (r^2 threshold = 0.2)...\n")
prune_out <- paste0(prefix, "_pruned")
prune_cmd_args <- c(
  "--bfile", prefix,
  "--indep-pairwise", "50", "5", "0.2",
  "--allow-extra-chr",
  "--out", prune_out
)
system2("plink2", args = prune_cmd_args)

prune_in_file <- paste0(prune_out, ".prune.in")
if (!file.exists(prune_in_file)) {
  stop("Error: LD pruning failed. .prune.in file was not generated.")
}

# 4. Run PCA on the pruned SNPs using plink2
cat("\nRunning PCA with plink2...\n")
pca_out <- paste0(prefix, "_pca")
pca_cmd_args <- c(
  "--bfile", prefix,
  "--extract", prune_in_file,
  "--pca", "10", 
  "--allow-extra-chr",
  "--out", pca_out
)
system2("plink2", args = pca_cmd_args)

eigenvec_file <- paste0(pca_out, ".eigenvec")
eigenval_file <- paste0(pca_out, ".eigenval")

if (!file.exists(eigenvec_file) || !file.exists(eigenval_file)) {
  stop("Error: PCA calculation failed. .eigenvec or .eigenval files are missing.")
}

# 5. Read PLINK2 PCA outputs
cat("\nProcessing PCA results...\n")
# Read eigenvectors
pca_data <- read.table(eigenvec_file, header = TRUE, comment.char = "", stringsAsFactors = FALSE)

# PLINK2 headers often start with #IID or #FID. We need to standardize the ID column.
colnames(pca_data)[colnames(pca_data) %in% c("#IID", "IID")] <- "IID"
if (!"IID" %in% colnames(pca_data)) {
  # If PLINK2 was run with a custom ID format, fallback to the first column
  colnames(pca_data)[1] <- "IID" 
}
pca_data$IID <- as.character(pca_data$IID)

# Read eigenvalues to calculate percentage of variance
eigenval <- read.table(eigenval_file, header = FALSE)$V1
# In PLINK's GRM, the sum of all eigenvalues equals the number of samples
N_samples <- nrow(pca_data)
pc.percent <- (eigenval / N_samples) * 100

# 6. Load and merge phenotype data
cat("\nLoading phenotype data and merging...\n")
pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)

# Force the first column to be named 'IID' and ensure it's a character string
first_col <- colnames(pheno)[1]
pheno <- pheno %>% 
  rename(IID = all_of(first_col)) %>%
  mutate(IID = as.character(IID))

# Merge PCA results with Phenotype data
merged_data <- inner_join(pca_data, pheno, by = "IID")

# Check if the requested color column exists in the merged data
if (!(color_col %in% colnames(merged_data))) {
  stop(paste("\nError: The column", color_col, "was not found in the phenotype file.\n"))
}

# 7. Generate Plots
cat("\nGenerating PCA plots...\n")

# Create dynamic axis labels with variance percentage
xlab1 <- paste0("PC1 (", round(pc.percent[1], 2), "%)")
ylab2 <- paste0("PC2 (", round(pc.percent[2], 2), "%)")
ylab3 <- paste0("PC3 (", round(pc.percent[3], 2), "%)")

# Determine if color_col is character. If so, convert to factor for discrete colors.
if (is.character(merged_data[[color_col]])) {
  merged_data[[color_col]] <- as.factor(merged_data[[color_col]])
}

# PC1 vs PC2
p1 <- ggplot(merged_data, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = xlab1, y = ylab2, title = "PCA: PC1 vs PC2", color = color_col) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# PC1 vs PC3
p2 <- ggplot(merged_data, aes(x = PC1, y = PC3, color = .data[[color_col]])) +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = xlab1, y = ylab3, title = "PCA: PC1 vs PC3", color = color_col) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# PC2 vs PC3
p3 <- ggplot(merged_data, aes(x = PC2, y = PC3, color = .data[[color_col]])) +
  geom_point(alpha = 0.8, size = 2) +
  labs(x = ylab2, y = ylab3, title = "PCA: PC2 vs PC3", color = color_col) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 8. Save to PDF
out_pdf <- paste0(prefix, "_PCA_plots.pdf")
pdf(out_pdf, width = 8, height = 6)
print(p1)
print(p2)
print(p3)
invisible(dev.off())

cat(paste("\nSuccess! Saved 3 plots to", out_pdf, "\n"))
