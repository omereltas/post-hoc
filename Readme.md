# Monte Carlo Simulation of Post-Hoc Tests under Variance Heterogeneity

This repository contains the R code and results for the study:

"An Empirical Comparison of Seven Post-Hoc Procedures under Variance Heterogeneity: A Monte Carlo Investigation in Small-Sample Designs"

## Overview

This study evaluates seven post hoc multiple comparison procedures:

- Tukey HSD  
- Scheffé  
- Fisher LSD (protected)  
- Games–Howell  
- Dunnett T3  
- Tamhane T2  
- Šidák correction  

The evaluation is based on:

- Type I error rate  
- Statistical power  

across a fully crossed factorial design:

- Number of groups: k = 3, 4, 5  
- Sample size per group: n = 6, 8, 10  
- Variance structure: homogeneous, moderate, and high heterogeneity  

## Files

- `R/Post Hoc kodlar.R` → main simulation script  
- `results/results_final.csv` → summarized simulation results  
- `results/results_final.rds` → full simulation output  

## Reproducibility

All simulations were conducted in R (version 4.5).

To reproduce the results:

1. Open the project in RStudio  
2. Run `Post Hoc kodlar.R`  
3. Results will be generated and saved in the `/results` folder  

## Author

Ömer Eltas  
Atatürk University  

## License

MIT License
