

# $W_d^*$: distance-based multivariate analysis of variance for multivariate data

## Introduction 

`WdStar` is an R package for $W_d^*$, a method for multivariate analysis of variance based on Welch's MANOVA designed to address challenges of existing methods including PERMANOVA.
- Are you using PERMANOVA and concerned about how heteroscedasticity and unbalanced sample sizes may affect your analyses? 
- Have you ever questioned the reliability of community-wide microbiome analysis results?
- Are you looking for a global test for your multivariate dataset?

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

## Methods Knowledge Base

See the [goodness-of-fit and covariate-adjustment methods inventory](https://github.com/alekseyenko/WdStar/blob/main/docs/goodness-of-fit-methods-knowledge-base.md) for implementation notes and comparisons with vegan, fast.adonis, DISTLM, MultANOVA, and microbiome-analysis packages.

## Installation  
Source installation of the `WdStar` R package is available directly from GitHub using `remotes` for R 3.6 or later:
```R
install.packages("remotes")
remotes::install_github("alekseyenko/WdStar", force = TRUE)
library(WdStar)
packageVersion("WdStar")
```

Until a CRAN release is available, install `WdStar` from GitHub as shown above.


## Quick Start  

For detailed and complex examples please refer to [our publication repositories](https://github.com/alekseyenko?tab=readme-ov-file#peer-reviewed-publications-on-the-w_d-test-family), which contain Markdown files with application datasets and code.


The following example uses the `phyloseq::enterotype` microbiome dataset to test whether community composition differs by enterotype using Bray-Curtis distances. The main examples use the full enterotype analysis set after removing samples with missing `Enterotype`.

```R
# Load packages and example microbiome data.
library(WdStar)
library(phyloseq)
data("enterotype", package = "phyloseq")

# Use all samples with observed Enterotype.
ent <- subset_samples(enterotype, !is.na(Enterotype))
dm <- distance(ent, method = "bray")
meta <- data.frame(sample_data(ent))
group <- factor(meta$Enterotype)

# Basic multivariate test.
WdS.test(
  dm = dm,
  group = Enterotype,
  data = sample_data(ent)
)

## By default, the unadjusted test's goodness.of.fit reports the distance-based
## pseudo-R-squared for the tested group.
unadjusted_res <- WdS.test(
  dm = dm,
  group = "Enterotype",
  data = sample_data(ent)
)
unadjusted_res$goodness.of.fit

## Use goodness = "none" to skip goodness-of-fit calculations.
WdS.test(
  dm = dm,
  group = Enterotype,
  data = sample_data(ent),
  goodness = "none"
)

# Stratified permutation example.
## This restricts permutations within sequencing technology. It is shown as a
## separate design choice from adjustment by SeqTech below.
WdS.test(
  dm = dm,
  group = Enterotype,
  strata = SeqTech,
  data = sample_data(ent)
)

# Covariate adjustment/elimination examples #
#############################################
## Adjustment example 1: pass unadjusted `dm` and an adjustment formula to
## WdS.test(). Here SeqTech is a categorical covariate.
WdS.test(
  dm = dm,
  group = Enterotype,
  adjustment = ~ SeqTech,
  data = sample_data(ent)
)

## By default, goodness.of.fit reports the distance-based semi-partial
## pseudo-R-squared for the tested group after adjustment.
## Additional components computed along the way are stored but not printed.
res <- WdS.test(
  dm = dm,
  group = Enterotype,
  adjustment = ~ SeqTech,
  data = sample_data(ent)
)
res$goodness.components

## Eigenvalue and tolerance diagnostics for residual distance matrices are
## stored separately from goodness-of-fit values.
res$distance.diagnostics

## Request all available goodness-of-fit components, including group-only,
## adjustment-only, full-model, semi-partial, and partial pseudo-R-squared.
WdS.test(
  dm = dm,
  group = Enterotype,
  adjustment = ~ SeqTech,
  data = sample_data(ent),
  goodness = "all"
)

## Interpretation, in default component order:
## - group.only: variation explained by the group alone, ignoring covariates.
## - adjustment.only: variation explained by covariates alone, ignoring the group.
## - full: variation explained by both together.
## - semi.partial = full - adjustment.only: the group's additional fraction of
##   ORIGINAL TOTAL distance variation.
## - partial = semi.partial / (1 - adjustment.only): the fraction of variation
##   REMAINING AFTER ADJUSTMENT explained by adding the group.
## Example: adjustment.only = 0.40 and full = 0.55 give semi.partial = 0.15
## (15 percentage points of total variation) and partial = 0.25
## (25% of the remaining 60%).
## These distance-based diagnostics generalize ordinary R2 and match it in
## the Euclidean setting described below. They do not establish causal effects. Components are not an additive partition:
## group.only + adjustment.only generally differs from full.
## Negative eigenvalue removal for non-Euclidean distances can disrupt usual
## nested-model inequalities; negative semi.partial means full < adjustment.only.
## Partial is NA when adjustment leaves essentially no residual variation.
## strata restricts permutations, not the goodness-of-fit calculations.
## goodness.of.fit contains selected/printed values in the requested order;
## goodness.components includes all calculated components and a displayed flag.

## Continuous covariates can also be used. Age has many missing values in
## enterotype, so this example builds a matching Age-complete distance matrix.
ent_age <- subset_samples(enterotype, !is.na(Enterotype) & !is.na(Age))
dm_age <- distance(ent_age, method = "bray")

WdS.test(
  dm = dm_age,
  group = Enterotype,
  adjustment = ~ Age,
  data = sample_data(ent_age)
)

## Adjustment example 2: create the adjusted distance matrix `a.dm` outside the
## function, then pass it to WdS.test().
a.dm <- a.dist(
  dm = dm,
  adjustment = ~ SeqTech,
  data = sample_data(ent)
)
WdS.test(
  dm = a.dm,
  group = Enterotype,
  data = sample_data(ent)
)
attr(a.dm, "distance.diagnostics")

## Store raw eigenvalues only when deeper diagnostics are needed.
a.dm.with.eigenvalues <- a.dist(
  dm = dm,
  adjustment = ~ SeqTech,
  data = sample_data(ent),
  keep.eigenvalues = TRUE
)
length(attr(a.dm.with.eigenvalues, "distance.diagnostics")$eigenvalues[[1]])

## Distance-based pseudo-R-squared can also be computed directly
dist.goodness.of.fit(dm = dm, dm_residual = a.dm)

```

Further examples are provided in the package documentation and may be accessed by running the following commands:
```R
?WdS.test
?a.dist
?dist.goodness.of.fit
```

## Component and diagnostic tables

`res$goodness.components` contains `component`, `pseudo.R2`, `model`, and
`displayed`. Definitions and interpretations are documented in `?WdS.test`.
Actual grouping names replace generic placeholders, for example:

| component | model |
|---|---|
| group.only | `~ Enterotype` |
| adjustment.only | `~ SeqTech` |
| full | `~ SeqTech + Enterotype` |
| semi.partial | `(Enterotype | SeqTech) / total` |
| partial | `(Enterotype | SeqTech) / residual` |

The conditioning bar denotes additional variation explained by adding the group
following adjustment. `total` means original total distance variation; `residual`
means variation remaining after adjustment. Multiple adjustment terms are grouped,
for example `(Enterotype | (Age + SeqTech)) / total`. These are descriptive labels,
not necessarily executable R formulas. Unadjusted labels use `none` after the bar;
both reference variations then equal original total variation, so the existing
equal component values are retained.

`res$distance.diagnostics` also uses `component` and `model`, retained even with
`goodness = "none"`. It contains only matrices actually reconstructed, ordered
as group.only, adjustment.only, full when present; inherited `input` diagnostics
precede them. Direct `a.dist()` diagnostics use adjustment.only. No extra matrices
are constructed for semi.partial or partial. Raw eigenvalues precede tolerance
handling and negative removal; final values follow both steps. Small negative
values zeroed by tolerance are counted separately from remaining negative values
removed afterward. `prop.negative.raw` measures absolute eigenvalue magnitude,
not eigenvalue counts. See `?a.dist` for each diagnostic field.

## Reading printed test results

`WdS.test()` prints the data description, parameters, WdS statistic, p-value,
omega-squared effect size estimate, and selected goodness-of-fit components.
The data description includes the tested variable name and number of levels.
Names longer than 30 characters are shortened only for display; full labels
remain in `res$data.info` and `res$data.name`. Variable names preserve their case.

- `res$statistic` contains WdS; `res$p.value` contains the permutation p-value.
- `res$estimate` contains omega-squared, estimated from WdS as
  `(k - 1) * (WdS - 1) / ((k - 1) * (WdS - 1) + N)`.
  This uses the ordinary one-way ANOVA omega-squared transformation with WdS
  in place of the classical F statistic.
- `res$goodness.of.fit` contains the selected distance-based pseudo-R2 values;
  `res$goodness.components` includes the calculated components and model descriptions.
  Their values depend on squared original and residual distances. For univariate
  Euclidean data and a corresponding least-squares model with an intercept,
  the basic measure equals ordinary R2, apart from numerical reconstruction.
- Stratification restricts permutations for the p-value; it does not adjust
  WdS, omega-squared, or R2. Use the `adjustment` formula for covariate adjustment.

## Feature Requests and Bugs
We welcome feature requests and bug reports and kindly ask you to submit them via [our GitHub issue tracker.](https://github.com/alekseyenko/WdStar/issues)
