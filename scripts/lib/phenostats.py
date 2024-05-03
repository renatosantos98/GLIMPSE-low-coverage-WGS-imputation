#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

df = pd.read_csv("phenostats/phenostats.csv")

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
    stat="percent",
    binwidth=5,
    binrange=(20, 69),
    kde=True,
    color=colors[0],
    fill=True,
    alpha=0.7,
)
axs[0, 0].set_xlabel("Age (years)")
axs[0, 0].set_ylabel("Percentage of severe COVID-19 patients")
axs[0, 0].text(
    -0.15, 1.05, "A", transform=axs[0, 0].transAxes, fontsize=14, fontweight="bold"
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
    -0.15, 1.05, "B", transform=axs[0, 1].transAxes, fontsize=14, fontweight="bold"
)

sns.boxplot(ax=axs[1, 0], x="Sex", y="Age", data=df, palette=colors)
axs[1, 0].set_ylabel("Age (years)")
axs[1, 0].text(
    -0.15, 1.05, "C", transform=axs[1, 0].transAxes, fontsize=14, fontweight="bold"
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
    -0.15, 1.05, "D", transform=axs[1, 1].transAxes, fontsize=14, fontweight="bold"
)

plt.tight_layout()
plt.savefig("phenostats/plot_1.png", dpi=300)

# Plot 2: Hospital Stays
fig, axs = plt.subplots(3, 2, figsize=(10, 12))

sns.histplot(
    ax=axs[0, 0],
    data=df,
    x="Days_in_hospital",
    stat="percent",
    binwidth=5,
    binrange=(0, 204),
    kde=True,
    color=colors[0],
    fill=True,
    alpha=0.7,
)
axs[0, 0].text(
    -0.15, 1.05, "A", transform=axs[0, 0].transAxes, fontsize=14, fontweight="bold"
)
axs[0, 0].set_xlabel("Days in hospital")
axs[0, 0].set_ylabel("Percentage of severe COVID-19 patients")

sns.boxplot(ax=axs[0, 1], x="Days_in_hospital", y="Sex", data=df, palette=colors)
axs[0, 1].text(
    -0.15, 1.05, "B", transform=axs[0, 1].transAxes, fontsize=14, fontweight="bold"
)
axs[0, 1].set_xlabel("Days in hospital")
axs[0, 1].set_ylabel("Sex")

sns.countplot(ax=axs[1, 0], x="ICU admittance", data=df, palette=colors)
for p in axs[1, 0].patches:
    axs[1, 0].annotate(
        f"{p.get_height():.0f}",
        (p.get_x() + p.get_width() / 2.0, p.get_height()),
        ha="center",
        va="center",
        fontsize=12,
        color="black",
        xytext=(0, 10),
        textcoords="offset points",
    )
axs[1, 0].text(
    -0.15, 1.05, "C", transform=axs[1, 0].transAxes, fontsize=14, fontweight="bold"
)
axs[1, 0].set_xlabel("ICU admittance")
axs[1, 0].set_ylabel("Number of severe COVID-19 patients")

sns.countplot(ax=axs[1, 1], x="Sex", hue="ICU admittance", data=df, palette=colors)
for p in axs[1, 1].patches:
    axs[1, 1].annotate(
        f"{p.get_height():.0f}",
        (p.get_x() + p.get_width() / 2.0, p.get_height()),
        ha="center",
        va="center",
        fontsize=12,
        color="black",
        xytext=(0, 10),
        textcoords="offset points",
    )
axs[1, 1].text(
    -0.15, 1.05, "D", transform=axs[1, 1].transAxes, fontsize=14, fontweight="bold"
)
axs[1, 1].set_xlabel("Sex")
axs[1, 1].set_ylabel("Number of severe COVID-19 patients")

sns.histplot(
    ax=axs[2, 0],
    data=df,
    x="Age",
    stat="percent",
    binwidth=5,
    binrange=(20, 69),
    hue="ICU admittance",
    multiple="stack",
    kde=True,
    palette=colors,
    fill=True,
    alpha=0.7,
)
axs[2, 0].text(
    -0.15, 1.05, "E", transform=axs[2, 0].transAxes, fontsize=14, fontweight="bold"
)
axs[2, 0].set_xlabel("Age (years)")
axs[2, 0].set_ylabel("Percentage of severe COVID-19 patients")

sns.histplot(
    ax=axs[2, 1],
    data=icu_df,
    x="ICU_Days",
    stat="percent",
    binwidth=5,
    binrange=(1, 175),
    kde=True,
    color=colors[0],
    fill=True,
    alpha=0.7,
)
axs[2, 1].text(
    -0.15, 1.05, "F", transform=axs[2, 1].transAxes, fontsize=14, fontweight="bold"
)
axs[2, 1].set_xlabel("Days in ICU")
axs[2, 1].set_ylabel("Percentage of severe COVID-19 patients")

plt.tight_layout()
plt.savefig("phenostats/plot_2.png", dpi=300)

# Plot 3: Phenotype Analysis
phenotype_labels = [
    "Pneumonia",
    "ARDS",
    "ARDS & ICU",
    "Exanthema",
    "Myocarditis",
    "Arrhythmia",
    "Hepatitis",
    "Glomerulonephritis",
    "Tubulopathy",
    "Encephalitis/encephalopathy",
    "Psychiatric",
    "Polyneuropathy",
    "Myelitis",
    "Seizure",
    "Diarrhoea",
    "Nausea/vomiting",
    "Endocrine dysfunction",
    "Myopathy",
    "Arthritis",
    "Bone marrow",
    "Pulmonary embolism",
    "Deep venous thrombosis",
    "Peripheral arterial thrombosis",
    "Stroke",
    "Ischemic heart event",
    "Disseminated intravascular coagulation",
    "Persistent fever",
    "Fatigue, malaise, headache",
]

sns.set_style("ticks")
phenotype_df = df.loc[:, "V-1":"V-28"].applymap(lambda x: 1 if x == "YES" else 0)
correlation_matrix, p_values = spearmanr(phenotype_df)

mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))

# Significance level (alpha)
alpha = 0.05

annot_matrix = np.where(p_values < alpha, correlation_matrix.round(2).astype(str), "")

plt.figure(figsize=(16, 12))
sns.heatmap(
    correlation_matrix,
    cmap="coolwarm",
    center=0,
    annot=annot_matrix,
    fmt="",
    xticklabels=phenotype_labels,
    yticklabels=phenotype_labels,
    mask=mask,
    cbar_kws={"label": "Spearman Correlation Coefficient"},
)

plt.title("Heatmap of Phenotype Correlations")
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig("phenostats/plot_3.png", dpi=300)
