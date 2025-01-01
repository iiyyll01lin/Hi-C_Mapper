# Load Bin file
bins <- read.table("GSE99104_nm_none_160000.bins.txt", header = TRUE, sep = "\t")

# Load Pair file
pairs <- read.table("pair_file.txt", header = TRUE, sep = "\t", fill=TRUE)

# Filter only chr=2L
bins <- subset(bins, chr == "2L")
pairs <- subset(pairs, chrom1 == "2L" & chrom2 == "2L")

# Initialize dataframe to store results
result <- data.frame(cbin1 = integer(), cbin2 = integer(), expected_count = numeric(), observed_count = numeric())

# Calculate observed count for each cbin pair
observed <- matrix(0, nrow = nrow(bins), ncol = nrow(bins))
for (i in 1:nrow(pairs)) {
  # Find corresponding cbin1 and cbin2 in pairs
  cbin1 <- with(bins, cbin[chr == pairs$chrom1[i] & from.coord <= pairs$pos1[i] & to.coord > pairs$pos1[i]])
  cbin2 <- with(bins, cbin[chr == pairs$chrom2[i] & from.coord <= pairs$pos2[i] & to.coord > pairs$pos2[i]])
  
  # Ensure cbin exists
  if (length(cbin1) == 1 & length(cbin2) == 1) {
    cbin1_index <- which(bins$cbin == cbin1)
    cbin2_index <- which(bins$cbin == cbin2)
    
    # Increment observation count
    observed[cbin1_index, cbin2_index] <- observed[cbin1_index, cbin2_index] + 1
    if (cbin1_index != cbin2_index) {
      observed[cbin2_index, cbin1_index] <- observed[cbin1_index, cbin2_index]  # Symmetric matrix
    }
  }
}

# Calculate bin center positions
bins$center <- (bins$from.coord + bins$to.coord) / 2

# Calculate expected count (distance-dependent model)
expected <- matrix(0, nrow = nrow(bins), ncol = nrow(bins))
for (i in 1:nrow(bins)) {
  for (j in i:nrow(bins)) {
    # Calculate distance between bin centers
    distance <- abs(bins$center[i] - bins$center[j])
    distance_bin <- floor(distance / 160000) * 160000  # Distance binning (e.g., every 160kb)
    
    # Find observed counts in this distance range
    observed_at_distance <- observed[abs(outer(bins$center, bins$center, "-")) == distance_bin]
    expected_count <- mean(observed_at_distance, na.rm = TRUE)
    
    # Fill the matrix
    expected[i, j] <- expected_count
    expected[j, i] <- expected[i, j]  # Symmetric matrix
    
    # Save to result dataframe
    result <- rbind(result, data.frame(cbin1 = bins$cbin[i], cbin2 = bins$cbin[j], 
                                       expected_count = expected[i, j], observed_count = observed[i, j]))
  }
}

# Save results to file
write.table(result, "cbin_observed_expected_2L_final_with_pairs.txt", sep = "\t", row.names = FALSE, quote = FALSE)


