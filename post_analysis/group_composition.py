import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# >>> helping variables
tissue_names = [
    "Breast",
    "Colorectal",
    "Non-Small_Cell_Lung",
    "Pancreatic",
    "Prostate",
]
with open("../dataset/gene_panel.txt") as f:
    gene_panel = f.readline().strip().split(",")
# <<< helping variables

# >>> Read tissue and group classifications
df = pd.read_csv("../dataset/sample_matrices/sample_matrix_Pan.txt", sep=" ")
patient_from_sample_index = (
    df["studyID:sampleId"].str.split("-").str[1].tolist()
)
data = np.genfromtxt(
    "../dataset/sample_matrices/sample_matrix_Pan.txt",
    skip_header=1,
    dtype=int,
)
data = data[:, 2:]

# >>> get tissue classification
dfs = [
    pd.read_csv(
        f"../dataset/sample_matrices/sample_matrix_{tissue}.txt", sep=" "
    )
    for tissue in tissue_names
]
patient_from_sample_index_types = [
    dfs_t["studyID:sampleId"].str.split("-").str[1].tolist() for dfs_t in dfs
]

tissues = np.zeros((5, len(patient_from_sample_index)))
for i in range(len(patient_from_sample_index)):
    for t in range(5):
        if patient_from_sample_index[i] in patient_from_sample_index_types[t]:
            tissues[t, i] = 1
tissues = np.argmax(tissues, axis=0).astype(int)
# <<< get tissue classification

classifications_MHN = (
    np.loadtxt(f"../classification_results/classification_fastMHN_13groups.dat")
    - 1
).astype(int)
classifications_CBN = (
    np.loadtxt("../classification_results/classification_CBN_13groups.dat") - 1
).astype(int)
group_names = [
    f"Group {i+1}" for i in range(int(np.max(classifications_MHN) + 1))
]
# <<< Read sample matrix and classification

# >>> calculate classification matrix
classification_matrix_MHN = np.zeros((len(tissue_names), len(group_names)))
classification_matrix_CBN = np.zeros((len(tissue_names), len(group_names)))
for i in range(len(tissues)):
    classification_matrix_MHN[tissues[i], classifications_MHN[i]] += 1
    classification_matrix_CBN[tissues[i], classifications_CBN[i]] += 1
# <<< calculate classification matrix

# >>> plot group composition
fig, ax = plt.subplots(1, 2)

# >>> fastMHN classification
bottom = np.zeros(len(group_names))
for tissue_index in range(len(tissue_names)):
    ax[0].bar(
        list(range(len(group_names))),
        classification_matrix_MHN[tissue_index],
        0.8,
        label=tissue_names[tissue_index],
        bottom=bottom,
    )
    bottom += classification_matrix_MHN[tissue_index]
ax[0].set_title("fastMHN classification")
ax[0].legend(loc="upper right")
ax[0].set_ylim(0, 3300)
ax[0].set_xticks(list(range(len(group_names))), labels=group_names, rotation=90)
# <<< fastMHN classification

# >>> CBN classification
bottom = np.zeros(len(group_names))
for tissue_index in range(len(tissue_names)):
    ax[1].bar(
        list(range(len(group_names))),
        classification_matrix_CBN[tissue_index],
        width=0.8,
        label=tissue_names[tissue_index],
        bottom=bottom,
    )
    bottom += classification_matrix_CBN[tissue_index]
ax[1].set_title("CBN classification")
ax[1].legend(loc="upper right")
ax[1].set_ylim(0, 3300)
ax[1].set_xticks(list(range(len(group_names))), labels=group_names, rotation=90)
# <<< CBN classification
# <<< plot group composition

plt.show()
