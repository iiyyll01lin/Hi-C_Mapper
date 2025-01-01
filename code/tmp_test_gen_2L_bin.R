# Install and load pheatmap package
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)

# Read data
bin_data <- read.table("GSE99106_nm_none_10000.bins.txt", header = TRUE, sep = "\t")

# Check data structure
str(bin_data)
head(bin_data)

# Filter data for specific chromosome (e.g. "chr2L")
specific_chr <- "2L"
filtered_data <- bin_data[bin_data$chr == specific_chr, ]

# Check filtered data
if (nrow(filtered_data) == 0) {
  stop("No data found for the specified chromosome.")
}

# Extract unique bin IDs
bin_ids <- unique(filtered_data$cbin)

# Check number of bin IDs
if (length(bin_ids) <= 1) {
  stop("Not enough bins to generate a meaningful heatmap.")
}

# Initialize a zero matrix
hic_matrix <- matrix(0, nrow = length(bin_ids), ncol = length(bin_ids))

# Add row and column names to matrix
rownames(hic_matrix) <- bin_ids
colnames(hic_matrix) <- bin_ids

# Fill the matrix
for (i in 1:nrow(filtered_data)) {
  bin_row <- as.character(filtered_data$cbin[i])
  bin_col <- as.character(filtered_data$cbin[i])  # Assuming data is symmetric
  hic_matrix[bin_row, bin_col] <- filtered_data$count[i]
}

# Check if matrix is empty
if (all(hic_matrix == 0)) {
  stop("The Hi-C matrix is empty or contains only zeros.")
}

# Log transform the matrix
log_hic_matrix <- log2(hic_matrix + 1)

# Add example annotation tracks (simulated functional domains)
annotation_data <- data.frame(
  Region = c("Domain1", "Domain2", "Domain3"),
  Type = c("Red", "Blue", "Black")
)
rownames(annotation_data) <- bin_ids[1:3]  # Assume first 3 bins correspond to domains

# Set annotation colors
annotation_colors <- list(
  Type = c(Red = "red", Blue = "blue", Black = "black")
)

# Plot heatmap
pheatmap(log_hic_matrix, 
         color = colorRampPalette(c("white", "blue", "red"))(100), 
         main = paste("Hi-C Heatmap for", specific_chr),
         annotation_row = annotation_data,
         annotation_colors = annotation_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

# Save heatmap to file
output_file <- "hic_heatmap_chr2L.png"
png(output_file, width = 800, height = 800)
pheatmap(log_hic_matrix, 
         color = colorRampPalette(c("white", "blue", "red"))(100), 
         main = paste("Hi-C Heatmap for", specific_chr),
         annotation_row = annotation_data,
         annotation_colors = annotation_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

cat("Hi-C heatmap saved to:", output_file, "\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("HiContacts")

library(HiContacts)
plotMatrix(hic)

# Load required libraries
library(ggplot2)
library(reshape2)

# Load the data
# Replace "your_file.csv" with the path to your data file
data <- read.table("GSE99104_nm_none_160000.n_contact.txt", header = TRUE)

# Create a matrix from observed counts
# Assuming cbin1 and cbin2 define rows and columns, and observed_count is the value
n_bins <- max(data$cbin1, data$cbin2)
matrix_observed <- matrix(0, nrow = n_bins, ncol = n_bins)

for (i in 1:nrow(data)) {
  row <- data$cbin1[i]
  col <- data$cbin2[i]
  matrix_observed[row, col] <- data$observed_count[i]
  matrix_observed[col, row] <- data$observed_count[i]  # Symmetry for Hi-C data
}

# Convert matrix to a data frame for ggplot2
heatmap_data <- melt(matrix_observed)
colnames(heatmap_data) <- c("x", "y", "value")

# Plot heatmap
ggplot(heatmap_data, aes(x = x, y = y, fill = log2(value + 1))) +
  geom_tile() +
  scale_fill_gradient(low = "skyblue", high = "yellow", na.value = "white") +
  theme_minimal() +
  labs(title = "Hi-C Contact Heatmap", x = "Genomic Bins", y = "Genomic Bins", fill = "log2(Counts)") +
  coord_fixed()
