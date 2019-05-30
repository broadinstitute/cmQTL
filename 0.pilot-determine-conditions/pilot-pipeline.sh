#!/bin/bash
#
# Cellular Morphology Resistance Mechanisms
#
# Pilot Conditions Analysis Pipeline

# Step 0 - Convert all notebooks to scripts
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb

# Step 1 - Determine optimal conditions and output all figures
Rscript --vanilla scripts/nbconverted/0.determine-optimal-conditions.r

# Step 2 - Extract single cell profiles into single files
# Note that this notebook is executed in AWS and assumes absolute paths pointing to single cell databases
jupyter nbconvert --to=html \
    --FilesWriter.build_directory=scripts/html \
    --ExecutePreprocessor.kernel_name=python3 \
    --ExecutePreprocessor.timeout=10000000 \
    --execute 1.extract-single-cell-profiles.ipynb

# Step 3 - Visualize the number of cells per well
Rscript --vanilla scripts/nbconverted/2.visualize-cell-count.r
