import cooler
import matplotlib.pyplot as plt
import numpy as np

c = cooler.Cooler("output.cool")

chunk_size = 1000

# total number of bins
n_bins = c.info["nbins"]

fig, ax = plt.subplots()
cax = ax.matshow(np.zeros((chunk_size, chunk_size)), cmap="coolwarm", aspect="auto")

for i in range(0, n_bins, chunk_size):
    for j in range(0, n_bins, chunk_size):
        # actual size of the chunk
        i_end = min(i + chunk_size, n_bins)
        j_end = min(j + chunk_size, n_bins)

        matrix_chunk = c.matrix(balance=True)[i:i_end, j:j_end]

        cax.set_data(matrix_chunk)
        ax.set_xlim(j, j_end)
        ax.set_ylim(i_end, i)
        plt.draw()
        plt.pause(0.01)

plt.colorbar(cax, label="Contact Frequency")
plt.title("Contact Map")
plt.xlabel("Bin 1")
plt.ylabel("Bin 2")

plt.savefig("contact_map.png")

plt.show()
