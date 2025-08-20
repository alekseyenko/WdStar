# WdStar
Distance-based Multivariate Welch ANOVA

## Updates
In this version we update exisiting funtions and introduce new ones including:
  - allows installation of package using `devtools::install_github("alekseyenko/WdStar", ref="sandbox", subdir = "pkg", force=T)`  
  - aWdS.test() to peform covariate-adjusted tests 
  - error and datatype handling improvements
  - improvements to test outputs including WdS.test()
  - Addition of omega squared estimate in output
  - Addition of between degrees of freedom
  - Changes of parameters within the a.dist() function to not request duplicate objects
