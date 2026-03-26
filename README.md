# ASV distribution and depth plots

This repository contains an R script used to analyze the distribution of selected 18S V4 amplicons/ASVs related to _Lemnamoeba sinaloensis_ 

## Input files

The script uses the following input files from EukBank:

* `eukbank_18S_V4_samples.tsv.gz`: sample metadata
* `eukbank_18S_V4_counts.tsv.gz`: ASV read counts per sample

It also requires a text file containing the ASVs of interest:

* `M6MM_clade_in_tree.txt`

## Data availability

The EukBank dataset files used in this analysis can be downloaded from Zenodo:

* https://doi.org/10.5281/zenodo.7804946

## Required R packages

The script requires:

* `readr`
* `dplyr`
* `ggplot2`
* `sf`
* `rnaturalearth`
* `rnaturalearthdata`

## What the script does

The script:

1. loads sample metadata and ASV count data
2. filters the count table to retain only ASVs of interest
3. merges ASV counts with sample metadata
4. calculates relative abundance
5. generates:

   * a global distribution map of ASVs
   * a depth distribution plot of ASVs

## Output files

The script saves:

* `asv_distribution.pdf`
* `depth_plot.pdf`

## Notes

* File paths in the script should be updated by the user before running.
* ASV colors are manually assigned in the script.
* ASV display names were manually adjusted later for the final figure version; see the Supplementary Table for details.
* In the map, point size represents the maximum relative abundance per site.
* In the depth plot, point size represents relative abundance.

## Usage

Open the script in R or RStudio, update the input and output paths, and run the script.
