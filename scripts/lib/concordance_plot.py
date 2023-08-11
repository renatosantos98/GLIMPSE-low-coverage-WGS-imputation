#!/usr/bin/env python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import csv
import gzip

matplotlib.use("Agg")

plt.rc("axes", labelsize=15)


# Function to read data from a CSV file
def read_csv_data(filename):
    x_maf = []
    y_aggr = []
    with gzip.open(filename, "rt") as csvfile:
        plots = csv.reader(csvfile, delimiter=" ")
        for row in plots:
            y_aggr.append(float(row[4]))
            x_maf.append(float(row[2]))
    return x_maf, y_aggr


x_maf1, y_aggr1 = read_csv_data("IBS001/bgi1_bgi30/output.rsquare.grp.txt.gz")
x_maf2, y_aggr2 = read_csv_data("IBS001/illumina1_illumina30/output.rsquare.grp.txt.gz")
x_maf3, y_aggr3 = read_csv_data("IBS001/bgi1_illumina30/output.rsquare.grp.txt.gz")
x_maf4, y_aggr4 = read_csv_data("IBS001/illumina1_bgi30/output.rsquare.grp.txt.gz")
x_maf5, y_aggr5 = read_csv_data("IBS001_filtered/bgi1_bgi30/output.rsquare.grp.txt.gz")
x_maf6, y_aggr6 = read_csv_data(
    "IBS001_filtered/illumina1_illumina30/output.rsquare.grp.txt.gz"
)
x_maf7, y_aggr7 = read_csv_data(
    "IBS001_filtered/bgi1_illumina30/output.rsquare.grp.txt.gz"
)
x_maf8, y_aggr8 = read_csv_data(
    "IBS001_filtered/illumina1_bgi30/output.rsquare.grp.txt.gz"
)

fig, axs = plt.subplots(2, 1, sharex=True)  # Create two subplots in the same figure

axs[0].semilogx(
    x_maf1,
    y_aggr1,
    "-",
    marker="o",
    lw=1,
    label="Aggregate $r^2$ MGI 1x vs MGI 40x",
    markersize=10,
    color="#ff7f00",
    alpha=0.9,
)
axs[0].semilogx(
    x_maf2,
    y_aggr2,
    "--",
    marker="o",
    lw=1,
    label="Aggregate $r^2$ Illumina 1x vs Illumina 40x",
    markersize=10,
    color="#33a02c",
    alpha=0.9,
)
axs[0].semilogx(
    x_maf3,
    y_aggr3,
    "-.",
    marker="o",
    lw=1,
    label="Aggregate $r^2$ MGI 1x vs Illumina 40x",
    markersize=10,
    color="#e31a1c",
    alpha=0.9,
)
axs[0].semilogx(
    x_maf4,
    y_aggr4,
    ":",
    marker="o",
    lw=1,
    label="Aggregate $r^2$ Illumina 1x vs MGI 40x",
    markersize=10,
    color="#1f78b4",
    alpha=0.9,
)
axs[0].set_ylabel("$r^2$ imputed vs true genotypes")
axs[0].grid(linestyle="--")
axs[0].legend(loc="lower right")
axs[0].set_ylim([0.1, 1.05])
axs[0].set_xlim([0.0005, 0.6])
axs[0].minorticks_off()


axs[1].semilogx(
    x_maf5,
    y_aggr5,
    "-",
    marker="o",
    lw=1,
    label="Aggregate $r^2$ MGI 1x vs MGI 40x",
    markersize=10,
    color="#ff7f00",
    alpha=0.9,
)
axs[1].semilogx(
    x_maf6,
    y_aggr6,
    "--",
    marker="o",
    lw=1,
    label="Aggregate $r^2$ Illumina 1x vs Illumina 40x",
    markersize=10,
    color="#33a02c",
    alpha=0.9,
)
axs[1].semilogx(
    x_maf7,
    y_aggr7,
    "-.",
    marker="o",
    lw=1,
    label="Aggregate $r^2$ MGI 1x vs Illumina 40x",
    markersize=10,
    color="#e31a1c",
    alpha=0.9,
)
axs[1].semilogx(
    x_maf8,
    y_aggr8,
    ":",
    marker="o",
    lw=1,
    label="Aggregate $r^2$ Illumina 1x vs MGI 40x",
    markersize=10,
    color="#1f78b4",
    alpha=0.9,
)
axs[1].set_xlabel("1000 Genomes Project Minor allele frequency (%)")
axs[1].set_ylabel("$r^2$ imputed vs true genotypes")
axs[1].grid(linestyle="--")
axs[1].legend(loc="lower right")
axs[1].set_ylim([0.1, 1.05])
axs[1].set_xlim([0.0005, 0.6])
mybins = [
    0.0002,
    0.0005,
    0.001,
    0.002,
    0.005,
    0.01,
    0.02,
    0.05,
    0.10,
    0.15,
    0.20,
    0.3,
    0.4,
    0.50,
    0.60,
]
labels = [
    "0.02",
    "0.05",
    "0.1",
    "0.2",
    "0.5",
    "1",
    "2",
    "5",
    "10",
    "15",
    "20",
    "30",
    "40",
    "50",
    "60",
]
axs[0].set_xticks(mybins)
axs[1].set_xticks(mybins)
axs[0].set_xticklabels(labels)
axs[1].set_xticklabels(labels)
axs[0].minorticks_off()
axs[1].minorticks_off()

# Add labels (letters) to each plot
axs[0].text(
    -0.10,
    0.98,
    "A",
    transform=axs[0].transAxes,
    fontsize=18,
    va="top",
    ha="right",
    fontweight="bold",
)
axs[1].text(
    -0.10,
    0.98,
    "B",
    transform=axs[1].transAxes,
    fontsize=18,
    va="top",
    ha="right",
    fontweight="bold",
)

fig.tight_layout()
fig.set_size_inches(10, 13)
fig.savefig("accplot_whole_genome.png", dpi=300)
