#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def pca_plot(eigenvec, eigenval, PED, super_populations, colors, color_mapping):
    # Calculate the proportion of variance explained by each principal component
    total_variance = eigenval.sum()
    proportion_variance_explained = (eigenval / total_variance * 100).round(2)

    # Ensure the PED and eigenvec dataframes are in the same order
    PED = PED.loc[eigenvec.index]
    assert all(PED.index == eigenvec.index)

    # Map each population to a superpopulation
    pop_to_superpop = {
        pop: superpop for superpop, pops in super_populations.items() for pop in pops
    }
    PED["Superpopulation"] = PED["Population"].map(pop_to_superpop)

    # Assign colors based on the 'Superpopulation' column
    PED["Color"] = PED["Superpopulation"].map(color_mapping)

    # Separate data for 'Severe COVID-19' population and the others
    covid_data = PED[PED["Superpopulation"] == "Severe COVID-19"]
    other_data = PED[PED["Superpopulation"] != "Severe COVID-19"]

    # Create plots
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))

    # First plot
    axs[0].scatter(
        other_data.index.map(eigenvec.iloc[:, 0]),
        other_data.index.map(eigenvec.iloc[:, 1]),
        c=other_data["Color"],
        s=53.6,
    )
    axs[0].scatter(
        covid_data.index.map(eigenvec.iloc[:, 0]),
        covid_data.index.map(eigenvec.iloc[:, 1]),
        facecolors="white",
        edgecolors=covid_data["Color"],
        s=53.6,
    )
    axs[0].set_xlabel(
        f"PC1: variance explained = {proportion_variance_explained.iloc[0, 0]}%"
    )
    axs[0].set_ylabel(
        f"PC2: variance explained = {proportion_variance_explained.iloc[1, 0]}%"
    )

    # Second plot
    axs[1].scatter(
        other_data.index.map(eigenvec.iloc[:, 0]),
        other_data.index.map(eigenvec.iloc[:, 2]),
        c=other_data["Color"],
        s=53.6,
    )
    axs[1].scatter(
        covid_data.index.map(eigenvec.iloc[:, 0]),
        covid_data.index.map(eigenvec.iloc[:, 2]),
        facecolors="white",
        edgecolors=covid_data["Color"],
        s=53.6,
    )
    axs[1].set_xlabel(
        f"PC1: variance explained = {proportion_variance_explained.iloc[0, 0]}%"
    )
    axs[1].set_ylabel(
        f"PC3: variance explained = {proportion_variance_explained.iloc[2, 0]}%"
    )

    # Create legend elements
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markersize=10)
        for color in colors[:-1]
    ] + [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="white",
            markeredgecolor=colors[-1],
            markersize=10,
        )
    ]

    # Add legend to the first plot only
    axs[0].legend(
        legend_elements,
        super_populations.keys(),
        title="Populations",
        loc="lower right",
    )

    # Adjust and save figure
    fig.set_size_inches(12, 6)
    fig.tight_layout()
    return fig, axs


# Load eigenvectors and eigenvalues
eigenvec = pd.read_csv(
    "pca/1000G_pca/plink.eigenvec", header=None, delim_whitespace=True
)
eigenval = pd.read_csv("pca/1000G_pca/plink.eigenval", header=None)

# Set index and update column names for the eigenvectors
eigenvec.index = eigenvec[1]
eigenvec = eigenvec.iloc[:, 2:]
eigenvec.columns = ["Principal Component " + str(i) for i in range(1, 21)]

# Load and process PED data
PED = pd.read_csv("pca/raw/20130606_g1k.ped", sep="\t")
PED = PED.set_index("Individual ID")

# Identify and assign missing individuals to 'Severe COVID-19' population
missing_individuals = set(eigenvec.index) - set(PED.index)
for individual in missing_individuals:
    PED.loc[individual, "Population"] = "Severe COVID-19"

# Fill any other missing values in the 'Population' column
PED["Population"].fillna("Severe COVID-19", inplace=True)

# Load new set of eigenvectors and eigenvalues for the IBS PCA
eigenvec_IBS = pd.read_csv(
    "pca/IBS_pca/plink.eigenvec", header=None, delim_whitespace=True
)
eigenval_IBS = pd.read_csv("pca/IBS_pca/plink.eigenval", header=None)

# Set index and update column names for the new eigenvectors
eigenvec_IBS.index = eigenvec_IBS[1]
eigenvec_IBS = eigenvec_IBS.iloc[:, 2:]
eigenvec_IBS.columns = ["Principal Component " + str(i) for i in range(1, 21)]

# Define superpopulations for the 1000G plot
super_populations = {
    "Africans (AFR)": ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"],
    "Admixed Americans (AMR)": ["CLM", "MXL", "PEL", "PUR"],
    "East Asians (EAS)": ["CDX", "CHB", "CHS", "JPT", "KHV"],
    "Europeans (EUR)": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "South Asians (SAS)": ["BEB", "GIH", "ITU", "PJL", "STU"],
    "Severe COVID-19": ["Severe COVID-19"],
}

# Define colors for each superpopulation for the 1000G plot
colors = ["#b2df8a", "#ff7f00", "#33a02c", "#1f78b4", "#e31a1c", "#4d4d4d"]

# Define color mapping for the 1000G plot
color_mapping = dict(zip(super_populations.keys(), colors))

# Create PCA plot for 1000G data
fig, axs = pca_plot(
    eigenvec,
    eigenval,
    PED,
    super_populations,
    colors,
    color_mapping,
)

# Save the 1000G PCA plot
fig.savefig("pca/1000G_pca/1000G_PCA_plot.png", dpi=300)

# Define populations for the IBS plot
populations_IBS = {
    "Iberian (IBS)": ["IBS"],
    "Severe COVID-19": ["Severe COVID-19"],
}

# Define colors for each population for the IBS plot
colors_IBS = ["#1f78b4", "#4d4d4d"]

# Define color mapping for the IBS plot
color_mapping_IBS = dict(zip(populations_IBS.keys(), colors_IBS))

# Create PCA plot for IBS data
fig_IBS, axs_IBS = pca_plot(
    eigenvec_IBS,
    eigenval_IBS,
    PED,
    populations_IBS,
    colors_IBS,
    color_mapping_IBS,
)

# Save the IBS PCA plot
fig_IBS.savefig("pca/IBS_pca/IBS_PCA_plot.png", dpi=300)

# Combine the plots
fig_combined, axs_combined = plt.subplots(nrows=2, figsize=(13, 13))

# Plot the first PCA plot on the top row
axs_combined[0].imshow(plt.imread("pca/1000G_pca/1000G_PCA_plot.png"))
axs_combined[0].axis("off")
axs_combined[0].annotate(
    "A",
    xy=(0.00, 0.98),
    xycoords="axes fraction",
    fontsize=16,
    fontweight="bold",
    va="top",
)

# Plot the second PCA plot on the bottom row
axs_combined[1].imshow(plt.imread("pca/IBS_pca/IBS_PCA_plot.png"))
axs_combined[1].axis("off")
axs_combined[1].annotate(
    "B",
    xy=(0.00, 0.98),
    xycoords="axes fraction",
    fontsize=16,
    fontweight="bold",
    va="top",
)

# Adjust and save the combined figure
fig_combined.tight_layout()
fig_combined.savefig("pca/combined_plot.png", dpi=300)
