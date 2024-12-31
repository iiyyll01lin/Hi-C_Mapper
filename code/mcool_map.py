import cooler
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

mcool_file = "output.mcool"

resolution = 10000
clr = cooler.Cooler(f"{mcool_file}::resolutions/{resolution}")

# Extract matrix (careful! can be large!)
matrix = clr.matrix(balance=True)[:]

# or, extract a subregion (eg, chromosome 2L, 0-1 Mb)
chrom = "2L"
start, end = 0, 1_000_000
bin_start = clr.offset(chrom) + start // resolution
bin_end = clr.offset(chrom) + end // resolution
sub_matrix = clr.matrix(balance=True)[bin_start:bin_end, bin_start:bin_end]

plt.figure(figsize=(10, 8))
sns.heatmap(
    sub_matrix,
    cmap="Reds",
    square=True,
    xticklabels=False,
    yticklabels=False,
    cbar_kws={"label": "Normalized Contact Frequency"},
)
plt.title(f"Contact Map for {chrom}: {start}-{end}")

plt.savefig("mcool_contact_map.png")

plt.show()
