# Post-Hoc Multiple Comparison Procedures: A Monte Carlo Simulation Study

This repository contains the R code and results for the simulation study:

**"An Empirical Comparison of Seven Post-Hoc Procedures under Variance 
Heterogeneity: A Monte Carlo Investigation in Small-Sample Designs"**

## Repository Contents

| File | Description |
|------|-------------|
| `Post Hoc kodlar.R` | Main Monte Carlo simulation code |
| `figures.R` | Code for reproducing Figure 1 (Power) and Figure 2 (Type I Error) |
| `results_final.csv` | Simulation output used in the manuscript |
| `results_final.rds` | Simulation output in R format |

## Requirements

```r
install.packages(c("PMCMRplus", "ggplot2", "dplyr"))
```

## How to Reproduce

1. Run `simulation.R` to generate the results 
2. Run `figures.R` to reproduce the figures

Alternatively, load `results_final.csv` directly and run only `figures.R`.

## Simulation Design

- **Procedures compared:** Tukey's HSD, Scheffé, Fisher's LSD, 
  Games–Howell, Dunnett T3, Tamhane T2, Šidák
- **Groups (k):** 3, 4, 5
- **Sample size per group (n):** 6, 8, 10
- **Variance conditions:** Homogeneous, Moderate heterogeneity, 
  High heterogeneity
- **Replications:** 10,000 per condition
- **Effect size:** Cohen's f = 0.40
- **Distribution:** Normal

## Reproducibility

All results are fully reproducible via `set.seed(12345)`.

## License

This code is released under the MIT License.

## Contact

For questions, please open an issue or contact via GitHub.
