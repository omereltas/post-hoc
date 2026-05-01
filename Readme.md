# Monte Carlo Simulation of Post-Hoc Tests in Small Samples

This repository contains the R code and simulation results for the study:

Eltas, Ö. (2026). *A Monte Carlo Comparison of Seven Post Hoc Procedures under Homogeneous and Heterogeneous Variance Conditions in Small Samples*.
## Overview

This study evaluates the performance of seven post hoc multiple comparison procedures:

- Tukey HSD  
- Scheffé  
- Fisher LSD (protected)  
- Games–Howell  
- Dunnett T3  
- Tamhane T2  
- Šidák correction  

The evaluation is based on:

- Familywise Type I error rate  
- Statistical power  

across a fully crossed factorial design:

- Number of groups: k = 3, 4, 5  
- Sample size per group: n = 6, 8, 10  
- Variance structure: homogeneous, moderate, and high heterogeneity  
- Replications: 10,000 per condition  

## Files

- `R/Post Hoc kodlar.R` → main simulation script  
- `results/results_final.csv` → full simulation results  
- `results/results_final.rds` → full simulation results (R object)  

## Reproducibility

All simulations were conducted in R (version 4.5).

To reproduce the results:

1. Open the project in RStudio  
2. Run `Post Hoc kodlar.R`  
3. Results will be generated and saved in the `/results` folder  

## Notes

This repository is provided for transparency and reproducibility.  

## Author

Ömer Eltas  
Department of Biometry  
Atatürk University, Turkey  

## License

MIT License

MIT License
