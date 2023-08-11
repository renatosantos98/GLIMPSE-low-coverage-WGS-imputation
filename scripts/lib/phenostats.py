#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
df = pd.read_csv("phenostats/phenostats.csv")

# Set color scheme
colors = [
    "#1f78b4",
    "#b2df8a",
    "#a6cee3",
    "#33a02c",
    "#fb9a99",
    "#e31a1c",
    "#fdbf6f",
    "#ff7f00",
    "#cab2d6",
    "#6a3d9a",
    "#ffff99",
    "#b15928",
    "#4d4d4d",
]

# Convert 'ICU_Days' to a boolean column indicating ICU admission
df["ICU"] = df["ICU_Days"] > 0
df["ICU admittance"] = df["ICU"].map({True: "Admitted", False: "Not admitted"})

# Ignore patients that did not go to the ICU for the ICU days analysis
icu_df = df[df["ICU_Days"] > 0]

# Plot 1: Patient Characterisation
sns.set_style("whitegrid")
fig, axs = plt.subplots(2, 2, figsize=(9, 9))

sns.histplot(
    ax=axs[0, 0],
    data=df,
    x="Age",
    binwidth=5,
    binrange=(20, 69),
    kde=True,
    color=colors[0],
    fill=True,
    alpha=0.5,
)
axs[0, 0].set_xlabel("Age (years)")
axs[0, 0].set_ylabel("Number of severe COVID-19 patients")
axs[0, 0].text(
    -0.15, 1.05, "A", transform=axs[0, 0].transAxes, fontsize=16, fontweight="bold"
)

sns.barplot(
    ax=axs[0, 1],
    x=df["Sex"].value_counts(),
    y=df["Sex"].value_counts().index,
    palette=colors,
)
for index, value in enumerate(df["Sex"].value_counts()):
    axs[0, 1].text(value + 1, index, str(value), ha="left", va="center")
axs[0, 1].set_xlabel("Number of severe COVID-19 patients")
axs[0, 1].set_ylabel("Sex")
axs[0, 1].text(
    -0.15, 1.05, "B", transform=axs[0, 1].transAxes, fontsize=16, fontweight="bold"
)

sns.boxplot(ax=axs[1, 0], x="Sex", y="Age", data=df, palette=colors)
axs[1, 0].set_ylabel("Age (years)")
axs[1, 0].text(
    -0.15, 1.05, "C", transform=axs[1, 0].transAxes, fontsize=16, fontweight="bold"
)

sns.barplot(
    ax=axs[1, 1],
    x=df["Country_of_origin"].value_counts(),
    y=df["Country_of_origin"].value_counts().index,
    palette=colors,
)
for index, value in enumerate(df["Country_of_origin"].value_counts()):
    axs[1, 1].text(value, index, str(value), ha="left", va="center")
axs[1, 1].set_xlabel("Number of severe COVID-19 patients")
axs[1, 1].set_ylabel("Country of origin")
axs[1, 1].text(
    -0.15, 1.05, "D", transform=axs[1, 1].transAxes, fontsize=16, fontweight="bold"
)

plt.tight_layout()
plt.savefig("phenostats/plot_1.png", dpi=300)

# Plot 2: Hospital Stays
fig, axs = plt.subplots(3, 2, figsize=(10, 12))

sns.histplot(
    ax=axs[0, 0],
    data=df,
    x="Days_in_hospital",
    bins=20,
    kde=True,
    color=colors[0],
    fill=True,
    alpha=0.5,
)
axs[0, 0].text(
    -0.15, 1.05, "A", transform=axs[0, 0].transAxes, fontsize=16, fontweight="bold"
)

sns.boxplot(ax=axs[0, 1], x="Sex", y="Days_in_hospital", data=df, palette=colors)
axs[0, 1].text(
    -0.15, 1.05, "B", transform=axs[0, 1].transAxes, fontsize=16, fontweight="bold"
)

sns.countplot(ax=axs[1, 0], x="ICU admittance", data=df, palette=colors)
axs[1, 0].text(
    -0.15, 1.05, "C", transform=axs[1, 0].transAxes, fontsize=16, fontweight="bold"
)

sns.histplot(
    ax=axs[1, 1],
    data=df,
    x="Age",
    hue="ICU admittance",
    multiple="stack",
    kde=True,
    palette=colors,
    fill=True,
    alpha=0.5,
)
axs[1, 1].text(
    -0.15, 1.05, "D", transform=axs[1, 1].transAxes, fontsize=16, fontweight="bold"
)

sns.countplot(ax=axs[2, 0], x="Sex", hue="ICU admittance", data=df, palette=colors)
axs[2, 0].text(
    -0.15, 1.05, "E", transform=axs[2, 0].transAxes, fontsize=16, fontweight="bold"
)

sns.histplot(
    ax=axs[2, 1],
    data=icu_df,
    x="ICU_Days",
    bins=20,
    kde=True,
    color=colors[0],
    fill=True,
    alpha=0.5,
)
axs[2, 1].text(
    -0.15, 1.05, "F", transform=axs[2, 1].transAxes, fontsize=16, fontweight="bold"
)

plt.tight_layout()
plt.savefig("phenostats/plot_2.png", dpi=300)

# Plot 3: Phenotype Analysis
sns.set_style("ticks")
phenotype_df = df.loc[:, "V-1":"V-28"].applymap(lambda x: 1 if x == "YES" else 0)
correlation_matrix = phenotype_df.corr(method="spearman")
plt.figure(figsize=(16, 12))
sns.heatmap(correlation_matrix, cmap="coolwarm", center=0, annot=True, fmt=".2f")
plt.title("Heatmap of Phenotype Correlations")
plt.tight_layout()
plt.savefig("phenostats/plot_3.png", dpi=300)

# Plot 4: Distribution of SNPs
file_path = "filtered_vcf/snps.txt"
with open(file_path, "r") as file:
    sample_lines = [file.readline().strip() for _ in range(10)]

sample_ids = []
number_of_snps = []

with open(file_path, "r") as file:
    for i in range(79):
        line = file.readline().strip()
        sample_id = line.split("/")[1].split(".")[0]
        snps = int(line.split(": ")[-1])
        sample_ids.append(sample_id)
        number_of_snps.append(snps)

sns.set_style("whitegrid")
plt.figure(figsize=(12, 5))
sns.barplot(x=sample_ids, y=[snps / 1e6 for snps in number_of_snps], color="#1f78b4")
plt.xticks(rotation=90)
plt.xlabel("Sample ID")
plt.ylabel("Number of SNVs (in millions)")
plt.tight_layout()
plt.savefig("filtered_vcf/plot_4.png", dpi=300)
