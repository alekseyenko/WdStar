#' Test if an object is of class 'dist'
#'
#' This function checks whether the given object belongs to class 'dist'.
#'
#' @param x Object to test.
#'
#' @return Logical indicating whether the object is a distance matrix.
#'
#' @examples
#' is.dist(as.dist(matrix(1:4, nrow = 2)))
#' is.dist(matrix(1:4, nrow = 2))
is.dist = function(x) any(class(x) == 'dist')

#' Calculate Sigma Squared for Distance Matrix
#'
#' This function computes the sigma squared statistic for a given distance matrix.
#'
#' @param dm Distance matrix.
#'
#' @return Sigma squared value.
#'
#' @examples
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' dist.sigma2(dm)
dist.sigma2 = function(dm) {
  dd = as.matrix(dm)
  dd[upper.tri(dd)] = 0
  sum(dd^2) / nrow(dd) / (nrow(dd) - 1)
}

#' Calculate Sum of Squares for Distance Matrix with Factor
#'
#' This function calculates the sum of squares for a given distance matrix,
#' weighted by a factor variable.
#'
#' @param dm2 Square distance matrix.
#' @param f Factor variable for weighting.
#'
#' @return Sum of squares matrix.
#'
#' @examples
#' dm2 <- matrix(runif(100), nrow = 10)
#' f <- factor(c(rep("A", 5), rep("B", 5)))
#' dist.ss2(dm2, f)
dist.ss2 = function(dm2, f) {
  K = sapply(levels(f), function(lev) f == lev)
  t(K) %*% dm2 %*% K / 2
}

#' Calculate Group-wise Sigma Squared
#'
#' This function calculates the sigma squared statistic for each group defined
#' by a factor variable in a given distance matrix.
#'
#' @param dm Distance matrix.
#' @param f Factor variable for group definition.
#'
#' @return A diagonal matrix of group-wise sigma squared values.
#'
#' @examples
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' f <- factor(c(rep("A", 5), rep("B", 5)))
#' dist.group.sigma2(dm, f)
dist.group.sigma2 = function(dm, f) {
  diag(dist.ss2(as.matrix(dm)^2, f)) / table(f) / (table(f) - 1)
}

#' Calculate Cohen's d for Distance Matrix
#'
#' This function calculates the Cohen's d statistic for a distance matrix,
#' assuming two levels in the factor variable.
#'
#' @param dm Distance matrix.
#' @param f Factor variable with two levels.
#'
#' @return Cohen's d value if factor has exactly two levels; NULL otherwise.
#'
#' @examples
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' f <- factor(c(rep("A", 5), rep("B", 5)))
#' dist.cohen.d(dm, f)
dist.cohen.d = function(dm, f){
  if(nlevels(f) != 2) return(NULL)
  SS2 = dist.ss2(as.matrix(dm^2),f)
  ns = summary(f)
  N = sum(ns)

  SST = sum(SS2)/N
  SSW = sum(diag(SS2)/ns)

  mean.diff = (sqrt((ns[1]+ns[2])/(ns[1]*ns[2])))*sqrt(SST-SSW)

  sigmas = diag(SS2)/ns/(ns-1)
  s1 = sigmas[1]
  s2 = sigmas[2]

  mean.diff/sqrt(((ns[1]-1)*s1 + (ns[2]-1)*s2)/(sum(ns)-2))
}
