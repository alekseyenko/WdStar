#' Test if an object is of class 'dist'
#'
#' This function checks whether the given object belongs to class 'dist'.
#'
#' @param x Object to test.
#'
#' @return Logical indicating whether the object is a distance matrix.
#'
#' @examples
#' 
#' is.dist(as.dist(matrix(1:4, nrow = 2)))
#' is.dist(matrix(1:4, nrow = 2))
#' 
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
#' 
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' dist.sigma2(dm)
#' 
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
#' 
#' dm <- matrix(runif(100), nrow = 10)
#' f <- factor(c(rep("A", 5), rep("B", 5)))
#' dist.ss2(dm, f)
#' 
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
#' \dontrun{
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' f <- factor(c(rep("A", 5), rep("B", 5)))
#' dist.group.sigma2(dm, f)
#' }
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
#' 
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' f <- factor(c(rep("A", 5), rep("B", 5)))
#' dist.cohen.d(dm, f)
#' 
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
#' performs PCoA, and returns the resulting corrected matrix
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#'
#' @param formula Only the right hand side of a typical formula such as Y~ A is necessary.
#'                The formula has the same requirements as in vegan::adonis() function.
#'
#' @param formula_data A dataset which contains the variables specified in formula. 
#'            It must be in a data.frame format with the row names the same as the 
#'            row names in distance matrix dm. This dataset should include both the 
#'            confounding covariate and the primary covariate.
#'             If not provided, the parent data.frame will be used.
#' @param tol Tolerance for eigenvalues. This is the cutoff for the eigenvalues 
#'            to be considered zero. Default is 10^-8.
#'
#' @return Returns a distance matrix of class `"dist"` representing the Euclidean distances
#'
#' @details The 'a.dist' function only requires a right-hand side of the formula. 
#'          Instead of the left-hand side, it uses the dissimilarity distance matrix dm.
#'          The function constructs a model matrix from the right-hand side (RHS) 
#'          of the formula. After performing necessary matrix operations and 
#'          eigen-decomposition, it calculates the Euclidean distances. It preserves the 
#'          labels of the input dm. 
#' @export
#' @examples
#' data(mtcars)
#' 
#' # The outcome could be a single variable or multiple variables (such as multidimensional omics data). 

#' ## This is an example with a single variable:
#' dm <- dist(mtcars$mpg, method="euclidean")
#' 
#' ## This is an example with multiple variables:
#' dm <- dist(mtcars[1:3], method="euclidean") 
#' 
#' # Right-hand side adjustment formula. Note that you may use any data type 
#' # including factor, character, integer, and numeric.
#' formula <- ~ as.factor(gear) + as.integer(hp) + wt
#' 
#' #' # Create the adjusted distance matrix 'a.dm'
#' a.dm <- a.dist(dm=dm, formula=formula, formula_data=mtcars)
#' a.dm
#' 
a.dist = function(dm, formula, formula_data=parent.frame(), tol=10^-8) 
{
  data <- formula_data
  Terms <- stats::terms(formula, data = data)
  # lhs <- formula[[2]]
  # lhs <- eval(lhs, data, parent.frame())
  # formula[[2]] <- NULL
  lhs <- dm
  rhs.frame <- stats::model.frame(formula, data, drop.unused.levels = TRUE) ## extract data.frame from a formula
  rhs <- stats::model.matrix(formula, rhs.frame) ## expand and create dummy variables for the terms
  grps <- attr(rhs, "assign") ## tells which column from the data set each dummy variable came from
  qrhs <- qr(rhs) ##computes the QR decomposition of rhs matrix
  rhs <- rhs[, qrhs$pivot, drop = FALSE]
  rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
  grps <- grps[qrhs$pivot][1:qrhs$rank]
  u.grps <- unique(grps)
  nterms <- length(u.grps) - 1
  if (nterms < 1) 
    stop("right-hand-side of formula has no usable terms")
  dmat <- as.matrix(lhs^2) #still has the labels of lhs
  X <- rhs
  X <- as.matrix(X)
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  A <- -1/2 * dmat
  J <- diag(nrow(X)) - matrix(rep(1/(nrow(X)), length(A)), nrow = nrow(A))
  E <- (diag(nrow(H)) - H) %*% J %*% A %*% J %*% (diag(nrow(H)) - H)

  # rownames(E) <- rownames(data)
  # colnames(E) <- rownames(data)
  
  # Correct rounding errors that cause E to be not symmetric
  E <- (E + t(E))/2
  
  eig <- eigen(E)
  eig$values[abs(eig$values) < tol] = 0
  lambda <- eig$values
  
  if (any(lambda < 0)) {
    warning(paste(sum(lambda < 0), "out of", length(lambda), "eigenvalues are negative!"))
  }
  
  w <- t(t(eig$vectors) * sqrt(lambda))
  w <- stats::dist(w)
  attr(w, "Labels") <- attr(lhs, "Labels")
  return(w)
}