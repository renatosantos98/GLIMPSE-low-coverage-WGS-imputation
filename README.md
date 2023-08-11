# Low-coverage whole genome sequencing for a highly selective cohort of severe COVID-19 patients

#### by Renato Santos, Víctor Moreno-Torres, Ilduara Pintos, Carmen de Mendoza, Octavio Corral, Vicente Soriano, and Manuel Corpas

This repository contains the scripts used in Santos et al. **Low-coverage whole genome sequencing for a highly selective cohort of severe COVID-19 patients**. Supporting data, including validation files and patient clinical histories, are archived in our [Figshare collection](https://doi.org/10.25452/figshare.plus.c.6347534). Genetic data for the patient cohort can be found in our EGA study. This repository will remain open for support with reproducibility and issues.

# Abstract

### Background

With the symptomatic diversity of severe COVID-19 cases, understanding the genetic risk factors influencing the progression of disease severity has become vital. Imputation of low-coverage whole genome sequencing has emerged as an economical method to study such disease-related genetic markers. This study aimed to explore hospitalisation trends, patient phenotypes, and the potential use of imputation in the context of severe COVID-19.

### Findings

From a well-characterised cohort of patients who exhibited severe COVID-19 symptoms during the first wave of the SARS-CoV-2 pandemic in Madrid, Spain, we generated an imputed dataset of 79 VCF files using the GLIMPSE1 tool, each containing, on average, 9.5 million single nucleotide variants. The imputation concordance assessment yielded a squared Pearson correlation of approximately 0.97 across sequencing platforms, showing that GLIMPSE1 can be used to confidently impute variants with minor allele frequency up to approximately 2% in Spanish ancestry individuals. We conducted a comprehensive analysis on the patient cohort, examining hospitalisation and intensive care utilisation, sex and age-based differences, and clinical phenotypes using a standardised set of medical terms we developed.

### Conclusion

This dataset not only highlights the utility and accuracy of low-coverage whole genome sequencing imputation in the study of COVID-19 severity, but also sets a precedent for other applications, especially in resource-constrained environments. The methodology validation, coupled with detailed patient characterisation, paves the way for comprehensive analyses of genetic components linked to various complex diseases. The methods and findings may be leveraged in future genomic projects, providing vital insights for health challenges like COVID-19.

![Principal component analysis of genetic variation in the severe COVID-19 patient cohort against the 1000 Genomes Project global superpopulations ](pca/1000G_pca/1000G_PCA_plot.png)
_Projection of imputed low-coverage whole-genome sequencing (lcWGS) data from severe COVID-19 patients against the backdrop of global superpopulations from the 1000 Genomes Project. Each point represents an individual, colour-coded according to their superpopulation. Severe COVID-19 patients are distinguished by points with a white fill and coloured border. The x-axis and y-axis on the two subplots represent the first and second, and first and third principal components, respectively, with the percentage of variance explained by each component indicated in the axis label._

# Software implementation

All the source code used to generate the results and figures in the paper are in the [`scripts`](scripts) folder. See the `README.md` files in each directory for a full description.

# Setup

## Getting the code

You can download a copy of all the files in this repository by cloning this git repository.

```
git clone https://github.com/renatosantos98/GLIMPSE-low-coverage-WGS-imputation.git
```

A copy of the repository is also archived at [doi.org/10.25452/figshare.plus.21679799](https://doi.org/10.25452/figshare.plus.21679799).

## Dependencies

You'll need a working Python environment to run the code. We recommend you set up your environment through [Anaconda](https://www.anaconda.com/download/), which provides the `conda` package manager.

Run the following command in the main repository folder (where ([`environment.yml`](environment.yml)) is located) to create a conda environment and install all required dependencies in it.

```
conda env create -f environment.yml
conda activate glimpse
```

## Running the code

All scripts were designed to be run from the main repository folder. To reproduce the data generated in the paper, run the scripts in the following order and syntax:

```
scripts/1_setup.sh
scripts/2_gl_calling.sh
scripts/3_glimpse_impute_parallel.sh
scripts/4_vcf_filtering.sh
scripts/5_glimpse_concordance.sh
scripts/6_pca.sh
```

See the `README.md` files in the [`scripts`](scripts) directory for a full description of each script and required files.

# License

All source code is made available under an MIT license. You can freely use and modify the code, without warranty. See [`LICENSE.md`](LICENSE.md) for the full license text. The authors reserve the rights to the
article content, which is currently submitted for publication.
