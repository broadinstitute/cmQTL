# Cell Morphology Quantitative Trait Loci (cmQTL)

Defining the genetic architecture of cell morphology

## Data

We are collecting cell painting data from induced pluripotent stem cells (iPSC) derived from ~1,000 patients.
We also have matched genotype data from these patients.
Our goal is to determine if specific genotypes can be mapped to cell morphology.

## Determining Cell Culture Conditions

Our first step is to determine optimal conditions at which to plate samples in order to generate morphological profiles.
We describe our approach in [0.pilot-determine-conditions](0.pilot-determine-conditions/).

Based on these analyses, we determined that the optimal conditions included plating cells for 6 hours at 10,000 plate density.

## Reproducibility

To reproduce the computational environment, first fork and clone the repository, then create and activate the conda environment.

```bash
# Using conda version > 4.6
conda env create --force --file environment.yml

# Activate environment
conda activate cmqtl
```

## Documents

**[GDrive folder](https://drive.google.com/drive/folders/1bBXGhk8Ic5KPD7W-BDYHjf-Gf3Q-RNohA)** 
