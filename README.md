

# $W_d^*$: distance-based multivariate analysis of variance for multivariate data

## 2025 Updates
In version 2.3.0 we update exisiting funtions and introduce new ones including:
  - allows installation of package using `devtools::install_github("alekseyenko/WdStar", ref="sandbox", subdir = "pkg", force=T)`
  - Package versioning  
  - Addition of `aWdS.test()` to peform covariate-adjusted tests 
  - Addition of omega squared estimate in output
  - Addition of between degrees of freedom
  - Changes of parameters within the `a.dist()` function to not request duplicate objects
  - error and datatype handling improvements
  - improvements to test outputs including `WdS.test()`
  - Updated examples for `a.dist()`, `WdS.test()`, and `aWdS.test()`


## Introduction 

`WdStar` is an R package for $W_d^*$, a more robust approach based on Welch's MANOVA designed to improve the statistical robustness of multivariate data analyses. 
- Are you using PERMANOVA and concerned about how heteroscedasticity and unbalanced sample sizes may affect your analyses? 
- Have you ever questioned the reliability of community-wide microbiome analysis results?
- Are you looking for a global test for your multivariate dataset?

>`WdStar` addresses these concerns by providing a more robust method for multivariate analysis of variance.  

**Our methods:** 
- are robust to heteroscedasticity; 
- can handle multi-level factors and stratification;
- allow for multiple post hoc testing scenarios; and
- allow for adjustment of covariates.

## Peer-Reviewed Articles on $W_d^*$
$aW_d^*$: preprint title
- Publication authors: [preprint doi link](https://doi.org/10.1093/bioinformatics/btw524)
- [Code repository](https://github.com/alekseyenko/WdStar/tree/master/publications/Bioinformatics%20(2025))

$W_d^*$-test: robust distance-based multivariate analysis of variance 
- [Hamidi B, Wallace K, Vasu C, & Alekseyenko AV. *Microbiome.* 2019.](https://doi.org/10.1186/s40168-019-0659-9)
- [Code repository](https://github.com/alekseyenko/WdStar/tree/master/publications/Hamidi%20et%20al.%20Microbiome%20(2019))


$T_w^2$: multivariate Welch t-test on distances
- [Alekseyenko AV. Multivariate Welch t-test on distances. *Bioinformatics*. 2016.](https://doi.org/10.1093/bioinformatics/btw524) 
- [Code repository](https://github.com/alekseyenko/Tw2)

## Installation  
Source installation of R package is available directly from GitHub using devtools:
```R
install.packages("devtools")
library("devtools")
devtools::install_github("alekseyenko/WdStar/pkg", force=T)
library(WdStar)
packageVersion("WdStar")
```


## Quick Start

We will use the `mtcars` dataset and look at the the effect of `gears` on `mpg`, `cyl`, and `disp` (first three variables of the dataset):   

```R
# Load dataset
data(mtcars)

# The outcome could be a single variable or multiple variables (such as multidimensional omics data).  

### This is an example with a single variable (`mpg`):
dm <- dist(mtcars$mpg, method="euclidean")

### This is an example with multiple variables (`mpg`, `cyl`, and `disp`):
dm <- dist(mtcars[1:3], method="euclidean") 

# Independent variable
f <- factor(mtcars$gear)

# Perform multivariate test  
WdS.test(dm=dm, f=f)

# Perform stratified test by variable 'vs'
strata <- factor(mtcars$vs)
WdS.test(dm=dm, f=f, strata=strata)
```

Further example codes are provided in the package documentation and may be accessed by running the following commands:
```R
?WdS.test #OR
?aWdS.test
```

Additionally [our publication repositories](./README.md#peer-reviewed-articles-on) contain Markdown files with application datasets and code.


## Feature Requests and Bugs
We welcome feature requests and bug reports and kindly ask you to submit them via [our GitHub issue tracker.](https://github.com/alekseyenko/WdStar/issues)
