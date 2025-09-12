#' Calculate the Tw2 Statistic for Heteroscedastic Test
#'
#' This function calculates the Tw2 statistic to compare means of k populations
#' with potentially unequal observations. It is suitable for microbiome data
#' and utilizes permutation testing for significance estimation.
#'
#' @param dm A distance matrix, representing dissimilarity between observations.
#' @param f A factor variable indicating the group for each observation.
#'
#' @return The calculated Tw2 statistic.
#'
#' @details
#' The function first checks if the factor variable has exactly two levels.
#' The Tw2 statistic is a modification of Hotelling's T-square statistic
#' adapted for heteroscedasticity and specifically suitable for microbiome data.
#' It calculates the sum of squares within each group and then computes the
#' Tw2 statistic based on these sum of squares.
#'
#' @examples
#' \dontrun{
#' # Generate a synthetic distance matrix and a factor variable
#' dm <- matrix(runif(100), nrow = 10)
#' f <- factor(rep(1:2, each = 5))
#' # Calculate the Tw2 statistic
#' Tw2_stat <- Tw2(dm, f)
#' }
#'
#' @references
#' Hamidi, Bashir, et al. "$ W_ {d}^{*} $-test: robust distance-based multivariate analysis of variance." Microbiome 7.1 (2019): 1-9.
#'
#' @export
Tw2 <- function(dm, f) {
  if (nlevels(f) != 2) {
    return(NULL)
  }
  SS2 <- dist.ss2(as.matrix(dm^2), f)
  ns <- summary(f)
  N <- sum(ns)

  SST <- sum(SS2) / N
  SSW1 <- SS2[1, 1] / ns[1]
  SSW2 <- SS2[2, 2] / ns[2]
  SSW <- SSW1 + SSW2

  s1 <- SSW1 / (ns[1] - 1)
  s2 <- SSW2 / (ns[2] - 1)

  t.stat <- (ns[1] + ns[2]) / (ns[1] * ns[2]) * (SST - SSW) / (s1 / ns[1] + s2 / ns[2])
  t.stat
}

#' Calculate the Wd* Statistic for Heteroscedastic Test
#'
#' This function calculates the Wd* statistic to compare means of k populations
#' with potentially unequal variances and observations. It is suitable for microbiome data
#' and utilizes permutation testing for significance estimation.
#'
#' @param dm A distance matrix, representing dissimilarity between observations.
#' @param f A factor variable indicating the group for each observation.
#'
#' @return The calculated Wd* statistic.
#'
#' @details
#' This function is an extension of Welch's ANOVA statistic suitable for multivariate
#' data and specifically for microbiome data. The Wd* statistic is computed based on
#' pairwise square differences between group means and variances. It explicitly accounts
#' for potentially unbalanced number of observations and differences in multivariate
#' spread in the two samples.
#'
#' @examples
#' \dontrun{
#' # Generate a synthetic distance matrix and a factor variable
#' dm <- matrix(runif(100), nrow = 10)
#' f <- factor(rep(1:2, each = 5))
#' # Calculate the Wd* statistic
#' WdS_stat <- WdS(dm, f)
#' }
#'
#' @references
#' Hamidi, Bashir, et al. "$ W_ {d}^{*} $-test: robust distance-based multivariate analysis of variance." Microbiome 7.1 (2019): 1-9.
#'
#' @seealso \url{https://github.com/alekseyenko/WdStar}
#'
#' @export
WdS <- function(dm, f) {
  ns <- table(f)
  SS2 <- dist.ss2(as.matrix(dm)^2, f)
  s2 <- diag(SS2) / ns / (ns - 1)
  W <- sum(ns / s2)

  idxs <- apply(utils::combn(levels(f), 2), 2, function(idx) levels(f) %in% idx)

  Ws <- sum(apply(
    idxs, 2,
    function(idx) {
      sum(ns[idx]) / prod(s2[idx]) *
        (sum(SS2[idx, idx]) / sum(ns[idx]) - sum(diag(SS2[idx, idx]) / ns[idx]))
    }
  ))
  k <- nlevels(f)
  h <- sum((1 - ns / s2 / W)^2 / (ns - 1))
  Ws / W / (k - 1) / (1 + (2 * (k - 2) / (k^2 - 1)) * h)
}


