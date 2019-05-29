# Determine Optimal Conditions

Our goal with this module is to determine the optimal cell painting conditions before scaling up data collection

## Conditions

We tested primary cell lines from 6 patients under the following conditions:

* Two plates: `BR00103267` and `BR00103268`
* Eight plating densities: `1000`, `2000`, `3000`, `4000`, `5000`, `10000`, `15000`, and `20000`
* Six cell lines: `A`, `B`, `C`, `D`, `E`, `F`
* Two time points: `6` and `24`
* Three profile types: `All Cells`, `Isolated Cells`, and `Cells in Colonies`

There are a total of 144 different combinations.

We measured each condition 16 times.
Our goal is to determine which of these conditions is the most reproducible across replicates.

## Profile Types

We extracted three different profiles.

1. All cells that passed quality control
2. Isolated Cells (Cells that were not touching any neighbor)
3. Colony Cells (Cells that were in contact with more than 3 neighbors)

## Pilot Conclusions

We determined that the most consistently reproducible profiles were acquired at **6 hours** and at **10,000 plating density**

![conditions](https://raw.githubusercontent.com/broadinstitute/cmQTL/master/0.pilot-determine-conditions/figures/condition_pilot_kstest.png?token=ABYADBJWGUZE6CZKPOSF2LK47BKFK)

## Reproduce Analysis

The first step to reproduce the analysis is to generate the three profile types listed above.

```bash
# Generate all profiles
bash generate-pilot-profiles.sh

# Generate single cell profiles
bash generate-sc-profiles.sh
```

**Note:** These profiles were created on AWS using security credentials.
The processed data are included in this repository to ensure reproducibility.

Next, perform the remaining analyses to determine optimal conditions.

```bash
# Ensure active conda environment
conda activate cmqtl

# Perform pilot analyses
bash pilot-pipeline.sh
```
