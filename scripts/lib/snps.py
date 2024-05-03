#!/usr/bin/env python
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import pysam
from multiprocessing import Pool


def plot_number_of_snps(ax, file_path="filtered_vcf/snps.txt"):
    sample_ids = []
    number_of_snps = []

    with open(file_path, "r") as file:
        for i in range(79):
            line = file.readline().strip()
            sample_id = line.split("/")[1].split(".")[0]
            snps = int(line.split(": ")[-1])
            sample_ids.append(sample_id)
            number_of_snps.append(snps)

    sns.barplot(
        x=sample_ids, y=[snps / 1e6 for snps in number_of_snps], color="#1f78b4", ax=ax
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_xlabel("Sample ID")
    ax.set_ylabel("Number of SNVs (in millions)")
    ax.text(-0.05, 1.05, "A", transform=ax.transAxes, fontsize=14, fontweight="bold")

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


def extract_snps(args):
    vcf_file, output_dir = args
    output_dir_path = Path(output_dir)
    output_path = output_dir_path / f"{Path(vcf_file).stem}_snps.csv"

    if output_path.exists():
        print(f"File {output_path} already exists. Skipping processing for {vcf_file}.")
        return

    print(f"Processing file {vcf_file}...")
    snp_data = []
    with pysam.VariantFile(vcf_file) as vcf:
        for rec in vcf.fetch():
            snp_data.append((rec.chrom, rec.pos, rec.ref, ",".join(rec.alts)))

    if not output_dir_path.exists():
        output_dir_path.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(snp_data, columns=["Chromosome", "Position", "Ref", "Alt"]).to_csv(
        output_path, index=False
    )
    print(f"File {vcf_file} processed and saved to {output_path}")


def calculate_positions(
    vcf_files,
    snp_data_dir="vcf_stats/snp_data",
    positions_file="vcf_stats/snp_positions_count.csv",
):
    positions_file_path = Path(positions_file)
    if positions_file_path.exists():
        print(
            f"Positions file {positions_file_path} already exists. Skipping position calculation."
        )
        return

    print("Calculating SNP positions and counts...")
    all_snp_data = []
    for vcf_file in vcf_files:
        snp_csv_path = Path(snp_data_dir) / f"{Path(vcf_file).stem}_snps.csv"
        if snp_csv_path.exists():
            snp_df = pd.read_csv(
                snp_csv_path,
                dtype={"Chromosome": str, "Position": int, "Ref": str, "Alt": str},
            )
            all_snp_data.append(snp_df)

    if all_snp_data:
        all_snp_data_df = pd.concat(all_snp_data, ignore_index=True)
        positions_count_df = (
            all_snp_data_df.groupby(["Chromosome", "Position"])
            .size()
            .reset_index(name="Count")
        )
        positions_count_df.to_csv(positions_file, index=False)
        print(f"SNP positions and counts saved to {positions_file}")
    else:
        print("No SNP data found.")


def plot_density(ax, positions_file="vcf_stats/snp_positions_count.csv", bin_size_mb=1):
    chrom_lengths_bp = {
        "1": 249250621,
        "2": 243199373,
        "3": 198022430,
        "4": 191154276,
        "5": 180915260,
        "6": 171115067,
        "7": 159138663,
        "8": 146364022,
        "9": 141213431,
        "10": 135534747,
        "11": 135006516,
        "12": 133851895,
        "13": 115169878,
        "14": 107349540,
        "15": 102531392,
        "16": 90354753,
        "17": 81195210,
        "18": 78077248,
        "19": 59128983,
        "20": 63025520,
        "21": 48129895,
        "22": 51304566,
        "X": 155270560,
    }
    chrom_order = [str(i) for i in range(1, 23)] + ["X"]
    chrom_lengths_mb = {
        chrom: length / 1e6 for chrom, length in chrom_lengths_bp.items()
    }

    positions_df = pd.read_csv(
        positions_file, dtype={"Chromosome": str, "Position": int}
    )
    positions_df["Bin"] = positions_df.apply(
        lambda row: int(row["Position"] / (1e6 * bin_size_mb)), axis=1
    )

    weighted_counts = (
        positions_df.groupby(["Chromosome", "Bin"]).size().reset_index(name="Count")
    )

    heatmap_df = pd.DataFrame()
    for chrom in chrom_order:
        max_bin = int(chrom_lengths_mb[chrom] // bin_size_mb)
        bins_df = pd.DataFrame({"Bin": range(max_bin + 1), "Chromosome": chrom})
        heatmap_df = pd.concat([heatmap_df, bins_df], ignore_index=True)

    heatmap_df = heatmap_df.merge(
        weighted_counts, on=["Chromosome", "Bin"], how="left"
    ).fillna(0)

    pivot_table = heatmap_df.pivot(index="Chromosome", columns="Bin", values="Count")
    pivot_table = pivot_table.loc[chrom_order]

    sns.heatmap(
        pivot_table,
        cmap="viridis",
        norm=PowerNorm(gamma=0.5),
        ax=ax,
        cbar_kws={"shrink": 0.5},
    )
    ax.set_title("SNP Density across Chromosomes")
    ax.set_xlabel("Position (Mb)")
    ax.set_ylabel("Chromosome")
    ax.text(-0.05, 1.05, "B", transform=ax.transAxes, fontsize=14, fontweight="bold")


def calculate_overlap(vcf_files, snp_data_dir, overlap_matrix_file):
    overlap_matrix_file_path = Path(overlap_matrix_file)
    if overlap_matrix_file_path.exists():
        print(
            f"Overlap matrix file {overlap_matrix_file_path} already exists. Skipping overlap calculation."
        )
        return

    snp_sets = {}
    total_snps_counts = {}

    for vcf_file in vcf_files:
        snp_csv_path = Path(snp_data_dir) / f"{Path(vcf_file).stem}_snps.csv"
        if snp_csv_path.exists():
            snp_df = pd.read_csv(
                snp_csv_path,
                dtype={"Chromosome": str, "Position": int, "Ref": str, "Alt": str},
                usecols=["Chromosome", "Position", "Ref", "Alt"],
            )
            snp_sets[Path(vcf_file).stem] = set(map(tuple, snp_df.to_numpy()))
            total_snps_counts[Path(vcf_file).stem] = len(snp_df)

    overlap_matrix_df = (
        pd.DataFrame(columns=snp_sets.keys(), index=snp_sets.keys())
        .fillna(0.0)
        .astype(float)
    )

    for file1, snps1 in snp_sets.items():
        for file2, snps2 in snp_sets.items():
            if file1 != file2:
                overlap_count = len(snps1.intersection(snps2))
                overlap_percentage = (
                    (overlap_count / total_snps_counts[file1]) * 100
                    if total_snps_counts[file1] > 0
                    else 0
                )
                overlap_matrix_df.loc[file1, file2] = overlap_percentage

    overlap_matrix_df.to_csv(overlap_matrix_file_path)
    print(f"Overlap matrix saved to {overlap_matrix_file_path}")


def plot_overlap(ax, overlap_matrix_file="vcf_stats/overlap_matrix.csv"):
    overlap_matrix_df = pd.read_csv(overlap_matrix_file, index_col=0).transpose()

    colors = ["#deebf7", "#9ecae1", "#3182bd"]
    n_bins = 100  # Define number of bins for granularity in color mapping
    cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=n_bins)

    sns.heatmap(
        overlap_matrix_df,
        cmap=cmap,
        ax=ax,
        cbar_kws={"shrink": 0.5},
        vmax=100,
    )
    ax.set_title("Percentage of Overlap of SNPs Between Samples")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.text(-0.05, 1.05, "C", transform=ax.transAxes, fontsize=14, fontweight="bold")


if __name__ == "__main__":
    vcf_dir = "filtered_vcf/"
    vcf_files = [str(path) for path in Path(vcf_dir).glob("*.vcf.gz")]

    with Pool() as pool:
        pool.map(
            extract_snps, [(vcf_file, "vcf_stats/snp_data") for vcf_file in vcf_files]
        )

    calculate_positions(
        vcf_files, "vcf_stats/snp_data", "vcf_stats/snp_positions_count.csv"
    )
    calculate_overlap(vcf_files, "vcf_stats/snp_data", "vcf_stats/overlap_matrix.csv")

    sns.set_style("ticks")
    fig, axs = plt.subplots(3, 1, figsize=(19, 26))

    plot_number_of_snps(axs[0])
    axs[0].set_title("Number of SNVs")

    plot_density(axs[1])
    axs[1].set_title("SNV Density")

    plot_overlap(axs[2])
    axs[2].set_title("SNV Overlap Percentage")

    plt.tight_layout()
    plt.savefig("vcf_stats/plot_4.png", dpi=300)
