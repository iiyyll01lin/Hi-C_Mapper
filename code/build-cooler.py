import cooler
import pandas as pd

contacts = pd.read_csv("GSE99104_nm_none_10000.n_contact.txt", sep="\t")
bins = pd.read_csv("GSE99104_nm_none_10000.bins.txt", sep="\t")

# Convert bin data to cooler format
bins_cooler = bins[["chr", "from.coord", "to.coord"]].copy()
bins_cooler["start"] = bins_cooler["from.coord"]
bins_cooler["end"] = bins_cooler["to.coord"]
bins_cooler = bins_cooler.rename(columns={"chr": "chrom"})

# Map cbin IDs to 0-based indices
bin_id_map = bins.set_index("cbin").index - 1
bin_id_map = bin_id_map.to_series().to_dict()
contacts["bin1_id"] = contacts["cbin1"].map(bin_id_map)
contacts["bin2_id"] = contacts["cbin2"].map(bin_id_map)

# Ensure bin1_id is always less than or equal to bin2_id
contacts["bin1_id"], contacts["bin2_id"] = (
    contacts[["bin1_id", "bin2_id"]].min(axis=1),
    contacts[["bin1_id", "bin2_id"]].max(axis=1),
)

# Combine duplicate pixels by summing their counts
cooler_contacts = (
    contacts.groupby(["bin1_id", "bin2_id"])
    .agg({"observed_count": "sum"})
    .reset_index()
)
cooler_contacts = cooler_contacts.rename(columns={"observed_count": "count"})

# Save
cooler.create_cooler(
    "output.cool",
    bins_cooler,
    cooler_contacts,
)
