

# $W_d^*$: distance-based multivariate analysis of variance for multivariate data

## Introduction 

`WdStar` is an R package for $W_d^*$, a method for distance-based multivariate
ANOVA designed to address challenges of existing methods including PERMANOVA.
- Are you using PERMANOVA and concerned about how heteroscedasticity and unbalanced sample sizes may affect your analyses?
- Have you ever questioned the reliability of community-wide microbiome analysis results?
- Are you looking for a global test for your multivariate dataset?
- Do you need to remove the effects of confounders and technical variables in your data?

>**$W_d^*$ addresses these concerns and is**
>- robust to heteroscedasticity; 
>- handles multi-level factors and stratification;
>- allows for multiple post hoc testing scenarios;
>- allows for adjustment of covariates; and
>- can be used with any distance or dissimilarity matrix.


## Peer-Reviewed Publications on the $W_d^*$-test Family  
**Preprint: Covariate-adjusted multivariate analysis for omics data**  
- Hamidi B, Fanning L, Wallace K, & Alekseyenko AV. *DOI forthcoming*
- [Code repository](https://github.com/alekseyenko/WdStar/tree/master/publications/Hamidi%20et%20al.%20FORTHCOMING%20(2026))


**$W_d^*$-test: robust distance-based multivariate analysis of variance**  
- [Hamidi B, Wallace K, Vasu C, & Alekseyenko AV. *Microbiome.* 2019.](https://doi.org/10.1186/s40168-019-0659-9)
- [Code repository](https://github.com/alekseyenko/WdStar/tree/master/publications/Hamidi%20et%20al.%20Microbiome%20(2019))


**$T_w^2$: multivariate Welch t-test on distances**  
- [Alekseyenko AV. Multivariate Welch t-test on distances. *Bioinformatics*. 2016.](https://doi.org/10.1093/bioinformatics/btw524) 
- [Code repository](https://github.com/alekseyenko/Tw2)


## Installation  
Source installation of the `WdStar` R package is available directly from GitHub using `remotes` for R 3.6 or later:
```R
install.packages("remotes")
remotes::install_github("alekseyenko/WdStar", force = TRUE)
library(WdStar)
packageVersion("WdStar")
```

Alternatively, the (latest released packaged source archives can be downloaded)[https://github.com/bashirhamidi/WdStar/releases]
and installed using:
```R
install.packages("[directory location]/WdStar_2.5.3.tar.gz", repos = NULL, type = "source")
```

## Quick Start  

For detailed and complex examples please refer to [our publication repositories](https://github.com/alekseyenko/WdStar/tree/master/publications), which contain Markdown files with application datasets and code.


The following examples use the `phyloseq::enterotype` microbiome dataset to test whether community composition differs by enterotype using Bray-Curtis distances. 

```R
## Load packages and example microbiome data
library(WdStar)
library(phyloseq) # we use the enterotype dataset for demonstration purposes
data("enterotype", package = "phyloseq")

## Use all samples with observed Enterotype
ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
dm <- phyloseq::distance(ent, method = "bray")

###########################################################
###########################################################
## (Eg 1) Basic multivariate test:
(res_unadjusted <- 
      WdS.test(dm = dm,
               group = Enterotype,
               data = ent))
### Data parameter can be given as a data frame, list, or phyloseq object.

###########################################################
###########################################################
## (Eg 2) Stratified permutation example:
### This restricts permutations within sequencing technology.
### Stratification does not change statistics, estimate, or goodness-of-fit.
(res_stratified <- 
      WdS.test(dm = dm,
               group = Enterotype,
               strata = SeqTech,
               data = ent))

###########################################################
###########################################################
## (Eg 3) Covariate adjustment by a categorical covariate (sequencing technology):
(res_adjusted <-
      WdS.test(dm = dm,
               group = Enterotype,
               adjustment = ~SeqTech,
               data = ent))

###########################################################
###########################################################
## Goodness-of-fit pseudo-R^2 components
### By default, goodness.of.fit reports the R^2 value of the tested group's for 
###    an unadjusted test and semi-partial R^2 of the tested group for an adjusted test.

### Request all available goodness-of-fit components:
###    group-only, adjustment-only, full-model, semi-partial, and partial
WdS.test(..., goodness = "all")

###########################################################
###########################################################
## Get distance diagnostics
### Eigenvalue and tolerance diagnostics for residual distance matrices are stored in the distance.diagnostics slot of the WdS.test output object.
res_adjusted$distance.diagnostics
```

Further examples are provided in the package documentation and may be accessed by running the following commands:
```R
?WdS.test
?a.dist
?dist.goodness.of.fit
```

## Feature Requests and Bugs
We welcome feature requests and bug reports and kindly ask you to submit them via [our GitHub issue tracker.](https://github.com/alekseyenko/WdStar/issues)
