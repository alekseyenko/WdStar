

# WdStar: robust distance-based multivariate analysis of variance for microbiome data

## Introduction 

`WdStar` is an R package designed to improve the statistical robustness of microbiome data analyses. 
- Have you ever questioned the reliability of community-wide microbiome analysis results? 
- Are you concerned about how heteroscedasticity and unbalanced sample sizes may affect your analyses? 
>`WdStar` addresses these concerns by providing a more robust method for multivariate analysis of variance. 

## Background
Microbiome data is incredibly rich and complex, making its analysis a challenging task. Commonly, community-wide analysis is used to evaluate the effects of interventions on microbiome composition. However,  methods such as PERMANOVA have limitations, prone to misleading results stemming from uneven sample sizes or difference in variance between groups (heteroscedasticity). 

## Our Solution
`WdStar` introduces a more robust approach based on Welch's MANOVA, specially designed to handle the intricacies of microbiome data. Our method: - Is robust to heteroscedasticity in the data. - Can handle multi-level factors and stratification. - Allows for multiple post hoc testing scenarios. 

## Quick Start
Install the package directly from GitHub: 
```R
remotes::install_github("alekseyenko/WdStar/pkg")
```
