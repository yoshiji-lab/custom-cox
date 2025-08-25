# Survival Plots (Q1–Q5) Toy Example

This repo contains a compact R script that simulates Q1–Q5 groups, fits a Cox model adjusted for **age**, **sex**, and **centre**, and outputs an **adjusted cumulative incidence** plot, an **HR table**, and an **RDS bundle** for easy re-plotting.

## Script
- `toy_km_targetpattern_100y_Q5_adj.R`

The script auto-installs required packages: `survival`, `survminer`, `dplyr`, `ggplot2`.

## Output naming
All results are written to `output/` using the pattern:
- `<cluster>_<outcome>.pdf` (figure)
- `<cluster>_<outcome>.csv` (HR table)
- `<cluster>_<outcome>.rds` (bundle)

Example already set in the script:
- `betacell1_type2diabetes.pdf`
- `betacell1_type2diabetes.csv`
- `betacell1_type2diabetes.rds`

## Quick start
```bash
Rscript toy_km_targetpattern_100y_Q5_adj.R
