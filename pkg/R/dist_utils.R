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
is.dist <- function(x) any(class(x) == "dist")

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
dist.sigma2 <- function(dm) {
  dd <- as.matrix(dm)
  dd[upper.tri(dd)] <- 0
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
dist.ss2 <- function(dm2, f) {
  K <- sapply(levels(f), function(lev) f == lev)
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
dist.group.sigma2 <- function(dm, f) {
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
dist.cohen.d <- function(dm, f) {
  if (nlevels(f) != 2) {
    return(NULL)
  }
  SS2 <- dist.ss2(as.matrix(dm^2), f)
  ns <- summary(f)
  N <- sum(ns)

  SST <- sum(SS2) / N
  SSW <- sum(diag(SS2) / ns)

  mean.diff <- (sqrt((ns[1] + ns[2]) / (ns[1] * ns[2]))) * sqrt(SST - SSW)

  sigmas <- diag(SS2) / ns / (ns - 1)
  s1 <- sigmas[1]
  s2 <- sigmas[2]

  mean.diff / sqrt(((ns[1] - 1) * s1 + (ns[2] - 1) * s2) / (sum(ns) - 2))
}

#' Covariate Adjusted Principal Coordinates Analysis
#'
#' This function takes in a formula and distance matrix, adjusts for covariates,
#' performs PCoA, and returns the resuting corrected matrix
#'
#' @param formula A typical formula such as Y~ A, but here Y is a dissimilarity distance. 
#'                The formula has the same requirements as in adonis function of the vegan package.
#'                
#' @param data A dataset with the rownames the same as the rownames in distance.
#'             This dataset should include both the confounding covariate and the primary covariate.
#'
#' @return Returns a distance matrix of class `"dist"` representing the Euclidean distances 
#'         between the aPCoA coordinates.
#'
#' @details The 'aPCoA.dist' function first evaluates the left-hand side (LHS) of the formula 
#'          within the data frame. Then, it constructs a model matrix from the right-hand side 
#'          (RHS) of the formula. After performing necessary matrix operations and eigen 
#'          decomposition, it calculates the Euclidean distances between the aPCoA coordinates.
#'
#' @examples
#' data(iris)
#' formula <- Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width
#' dist_matrix <- apcoa.dist(formula, iris)
#' print(dist_matrix)
#'
aPCoA.dist = function (formula, data){
  Terms <- terms(formula, data = data)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  rhs <- model.matrix(formula, rhs.frame)
  grps <- attr(rhs, "assign")
  qrhs <- qr(rhs)
  rhs <- rhs[, qrhs$pivot, drop = FALSE]
  rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
  grps <- grps[qrhs$pivot][1:qrhs$rank]
  u.grps <- unique(grps)
  nterms <- length(u.grps) - 1
  if (nterms < 1) 
    stop("right-hand-side of formula has no usable terms")
  dmat <- as.matrix(lhs^2)
  X <- rhs
  y <- lhs
  X <- X[rownames(dmat), ]
  X <- as.matrix(X[, -1], nrow = nrow(X))
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  A <- -1/2 * as.matrix(y)^2
  J <- diag(nrow(X)) - matrix(rep(1/(nrow(X)), length(A)), 
                              nrow = nrow(A))
  E <- (diag(nrow(H)) - H) %*% J %*% A %*% J %*% (diag(nrow(H)) - 
                                                    H)
  rownames(E) <- rownames(data)
  colnames(E) <- rownames(data)
  
  eigen_stat = eigen(E)
  eigenE <- eigen_stat$vectors
  eigenvalue <- eigen_stat$values
  
  rownames(eigenE) <- rownames(data)
  plotMatrix <- eigenE * matrix(rep(eigenvalue^(1/2), each = nrow(eigenE)), 
                                nrow = nrow(eigenE))
  plotMatrix <- plotMatrix[, !is.na(apply(plotMatrix, 2, sum))]
  
  dist(plotMatrix, method="euclidean")
}