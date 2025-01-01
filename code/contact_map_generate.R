# Load required libraries
library(ggplot2)
library(reshape2)

# Load the data
# Replace ".txt" with the path to your data file
data <- read.table("GSE99104_nm_none_10000.n_contact.txt", header = TRUE)

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

library(tidyverse)
# Filter data for the range of 2500~3000 (select range@interval of 500)
filtered_data <- heatmap_data %>%
  filter(x >= 2500 & x <= 3000, y >= 2500 & y <= 3000)

ggplot(filtered_data, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#ffffff", "#2793fe", "#8c0703", "#f99d01", "#f7f4f2", "#2b2b9f"),
    values = scales::rescale(c(0, 2, 7, 11, 14, 16)),
    breaks = c(0, 1, 6, 12, 16),
    labels = c("0", "2", "7", "14", "16"),
    name = "Contact Enrichment (log"[2]~" scale)"
  ) +
  theme_minimal() +
  scale_x_continuous(
    breaks = c(2500, 2600, 2700, 2800, 2900, 3000),  # Specify where the ticks should be
    labels = c("10,000kb", "10,500kb", "11,000kb", "11,500kb", "12,000kb", "12,500kb")  # Custom labels for those ticks
  ) +
  scale_y_continuous(
    breaks = c(2500, 2600, 2700, 2800, 2900, 3000),  # Specify where the ticks should be
    labels = c("10,000kb", "10,500kb", "11,000kb", "11,500kb", "12,000kb", "12,500kb")  # Custom labels for those ticks
  ) +
  labs(
    title = "Hi-C Contact Heatmap",
    x = "2L",
    y = "2L"
  ) +
  coord_fixed()
