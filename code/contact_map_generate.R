# Load required libraries
library(ggplot2)
library(reshape2)

# Load the data
# Replace "your_file.csv" with the path to your data file
data <- read.table("GSE99104_nm_none_40000.n_contact.txt", header = TRUE)

# Create a matrix for log-transformed observed/expected ratios
n_bins <- max(data$cbin1, data$cbin2)
matrix_log_ratio <- matrix(NA, nrow = n_bins, ncol = n_bins)  # Initialize with NA for symmetry handling

# Fill the matrix with log-transformed observed/expected ratios
for (i in 1:nrow(data)) {
  row <- data$cbin1[i]
  col <- data$cbin2[i]
  obs <- data$observed_count[i]
  exp <- data$expected_count[i]
  
  if (exp > 0) {  # Avoid division by zero
    ratio <- log2(obs / exp + 1)  # Log2 transformation with a +1 offset to handle very small ratios
  } else {
    ratio <- NA
  }
  
  matrix_log_ratio[row, col] <- ratio
  matrix_log_ratio[col, row] <- ratio  # Symmetry for Hi-C data
}

# Convert the log-ratio matrix to a data frame for ggplot2
heatmap_data <- melt(matrix_log_ratio)
colnames(heatmap_data) <- c("x", "y", "value")

# Plot heatmap with three to five colors (如果跑不到，用下面的)
ggplot(heatmap_data, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#ffffff", "#44a2fe", "#653056", "#d26500", "#fec14f", "#faeeda", "#3434a3"), 
    values = scales::rescale(c(0, 2, 7, 11, 14, 16)), # Adjust positions of colors
    na.value = "white",
    name = "Contact Enrichment (log2 scale)"
  ) +
  scale_x_continuous(
    breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000),  # Specify where the ticks should be
    labels = c("10,000kb", "10,500kb", "11,000kb", "11,500kb", "12,000kb", "12,500kb", "13,000kb")  # Custom labels for those ticks
  ) +
  scale_y_continuous(
    breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000),  # Specify where the ticks should be
    labels = c("10,000kb", "10,500kb", "11,000kb", "11,500kb", "12,000kb", "12,500kb", "13,000kb")  # Custom labels for those ticks
  ) +
  labs(
    title = "Hi-C Contact Heatmap"
  ) +
  coord_fixed()

# Plot heatmap with three to five colors
ggplot(heatmap_data, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#ffffff", "#44a2fe", "#653056", "#d26500", "#fec14f", "#faeeda", "#3434a3"), 
    values = scales::rescale(c(0, 2, 7, 11, 14, 16)), # Adjust positions of colors
    na.value = "white",
    name = "Contact Enrichment (log2 scale)"
  ) +
  theme_minimal() +
  labs(title = "Hi-C Log-transformed Observed/Expected Ratio Heatmap") +
  coord_fixed()
