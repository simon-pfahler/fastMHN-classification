import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, statistics

# >>> helping variables
data = np.genfromtxt(
    "../dataset/sample_matrices/sample_matrix_Pan.txt",
    skip_header=1,
)
data = data[:, 2:]
with open("../dataset/gene_panel.txt") as f:
    gene_panel = f.readline().strip().split(",")
# <<< helping variables

# Read the KM data file
df = pd.read_csv("../dataset/KM_Data.txt", sep=" ")
statuss = df["OS_STATUS"].str.split(":").str[0].astype(int).to_numpy()
times = df["OS_MONTHS"].to_numpy() / 12

max_time = max(times)
ordered_indices = np.argsort(times)

# Read the patient list
classification_groups = np.loadtxt(
    "../classification_results/classification_fastMHN_13groups.dat"
)
group_names = [f"Group {i+1}" for i in range(len(classification_groups))]

# >>> get tissue classification
df = pd.read_csv("../dataset/sample_matrices/sample_matrix_Pan.txt", sep=" ")
patient_from_sample_index = (
    df["studyID:sampleId"].str.split("-").str[1].tolist()
)
tissue_names = [
    "Breast",
    "Colorectal",
    "Non-Small_Cell_Lung",
    "Pancreatic",
    "Prostate",
]
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
classification_tissues = np.argmax(tissues, axis=0).astype(int)
# <<< get tissue classification

fig, ax = plt.subplots(2, 2)

for x in range(2):
    for y in range(2):
        times_clustered = [list() for _ in range(2)]
        statuss_clustered = [list() for _ in range(2)]

        names = [
            "STK11=0",
            "STK11=1",
        ]

        for i, group in enumerate(classification_groups):

            def in_group(x, y, i):
                match x:
                    case 0:
                        if classification_groups[i] != 11:
                            return False
                    case 1:
                        if classification_groups[i] == 11:
                            return False
                match y:
                    case 0:
                        if classification_tissues[i] != 2:
                            return False
                    case 1:
                        if classification_tissues[i] == 2:
                            return False
                return True

            if in_group(x, y, i):
                times_clustered[int(data[i, gene_panel.index("STK11")])].append(
                    times[i]
                )
                statuss_clustered[
                    int(data[i, gene_panel.index("STK11")])
                ].append(statuss[i])

        print(len(times_clustered[0]), len(times_clustered[1]))

        kmfs = [KaplanMeierFitter() for _ in range(len(times_clustered))]

        for group in range(len(times_clustered)):
            kmfs[group].fit(
                durations=times_clustered[group],
                event_observed=statuss_clustered[group],
                label=names[group],
            )

        # >>> quantitative analysis
        results = statistics.logrank_test(
            times_clustered[0],
            times_clustered[1],
            event_observed_A=statuss_clustered[0],
            event_observed_B=statuss_clustered[1],
        )

        chi2 = results.test_statistic
        p_value = results.p_value
        # <<< quantitative analysis

        # >>> plotting
        for group in range(len(times_clustered)):
            kmfs[group].plot(ax=ax[x][y])
        ax[x][y].set_xlabel("Survival time (Years)")
        ax[x][y].set_ylabel("Survival probability")

        if x == y == 0:
            ax[x][y].set_title(
                f"Survival of lung cancers with and without STK11 in group 11\n"
                f"LR={chi2:.4f}, p={p_value:.4e}"
            )
        if x == 0 and y == 1:
            ax[x][y].set_title(
                f"Survival of non-lung cancers with and without STK11 in group 11\n"
                f"LR={chi2:.4f}, p={p_value:.4e}"
            )
        if y == 0 and x == 1:
            ax[x][y].set_title(
                f"Survival of lung cancers with and without STK11 in other groups\n"
                f"LR={chi2:.4f}, p={p_value:.4e}"
            )
        if x == y == 1:
            ax[x][y].set_title(
                f"Survival of non-lung cancers with and without STK11 in other groups\n"
                f"LR={chi2:.4f}, p={p_value:.4e}"
            )

        ax[x][y].set_ylim(0, 1)

        # <<< plotting

plt.show()
