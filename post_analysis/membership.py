import re

import matplotlib.pyplot as plt
import numpy as np

# >>> Read stratified data and classification results
patients = []
with open("../dataset/sample_matrices/sample_matrix_Pan.txt") as f:
    for line in f:
        regex_match = re.search(r"(P-\d+)", line)
        if regex_match:
            patients.append(regex_match.group(1))

data = np.genfromtxt(
    "../dataset/sample_matrices/sample_matrix_Pan.txt",
    skip_header=1,
    dtype=int,
)
data = data[:, 2:]

predictions_MHN = np.loadtxt(
    f"../classification_results/sample_Ps_fastMHN_13groups.dat"
).T
predictions_CBN = np.loadtxt(
    f"../classification_results/sample_Ps_CBN_13groups.dat"
).T
predictions_baserate = np.loadtxt(
    f"../classification_results/sample_Ps_baserate_13groups.dat"
).T
group_names = [f"Group {i+1}" for i in range(predictions_MHN.shape[1])]
# <<< Read sample matrix and classification


# >>> sort samples according to mutational burden
def sort(data, predictions):
    mutational_burdens = np.sum(data, axis=1)

    sorted_indices = np.lexsort(
        (-np.max(predictions, axis=1), mutational_burdens)
    )
    return sorted_indices, mutational_burdens[sorted_indices]


sorted_indices, mutational_burdens = sort(data, predictions_MHN)
predictions_MHN = predictions_MHN[sorted_indices]
predictions_CBN = predictions_CBN[sorted_indices]
predictions_baserate = predictions_baserate[sorted_indices]
# <<< sort samples according to mutational burden

# >>> get y ticks and labels
unique_mbs, indices = np.unique(mutational_burdens, return_index=True)
indices = np.append(indices, len(patients))
ytick_poss = [(i1 + i2) / 2 for i1, i2 in zip(indices[:-1], indices[1:])]
# <<< get y ticks and labels

# >>> plot results
fig, ax = plt.subplots(1, 3)
ax[0].imshow(
    predictions_MHN, cmap="Grays", interpolation="none", vmin=0, vmax=1
)
ax[1].imshow(
    predictions_CBN, cmap="Grays", interpolation="none", vmin=0, vmax=1
)
im = ax[2].imshow(
    predictions_baserate, cmap="Grays", interpolation="none", vmin=0, vmax=1
)
for a in ax:
    for ytick_pos, index in zip(ytick_poss, indices):
        a.text(
            -0.002,
            ytick_pos,
            mutational_burdens[index],
            transform=a.get_yaxis_transform(),
            ha="right",
            va="center",
        )
    a.set_yticks(indices, ["" for _ in indices])
    a.set_ylabel("Patients sorted by mutational burden")
    a.grid(axis="y", color="black", linewidth=0.5)
    a.axis("auto")
    a.set_xticks(
        list(range(predictions_MHN.shape[1])), group_names, rotation=90
    )
ax[0].set_title("fastMHN")
ax[1].set_title("CBN")
ax[2].set_title("base-rate")
fig.colorbar(im)
fig.suptitle(f"Membership probabilities for different classification methods")
# <<< plot results

plt.show()
