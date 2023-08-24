#' Perform posthoc tests for Tw2 statistic
#'
#' This function performs posthoc tests for the Tw2 statistic using a permutation approach.
#'
#' @param dm A distance matrix.
#' @param f A factor variable.
#' @param nrep Number of permutations (default: 999).
#' @param strata A factor variable for stratified permutation (default: NULL).
#'
#' @return A matrix containing posthoc test results.
#'
#' @examples
#' data <- ...  # TODO add data here
#' Tw2.posthoc.tests(data$dm, data$f)
#'
#' @export
Tw2.posthoc.tests <- function(dm, f, nrep = 999, strata = NULL) {
  dd <- as.matrix(dm)

  Tw2.subset.test <- function(include.levels) {
    subs <- f %in% include.levels
    c(include.levels,
      table(f[subs, drop = TRUE]),
      Tw2.test(dd[subs, subs], f[subs, drop = TRUE], nrep = nrep, strata = strata[subs]))
  }

  res <- t(utils::combn(levels(f), 2, Tw2.subset.test))
  colnames(res) <- c("Level1", "Level2", "N1", "N2", "p.value", "tw2.stat", "nrep")
  res
}

#' Perform 1-vs-All posthoc tests for Tw2 statistic
#'
#' This function performs 1-vs-All posthoc tests for the Tw2 statistic using a permutation approach.
#'
#' @param dm A distance matrix.
#' @param f A factor variable.
#' @param nrep Number of permutations (default: 999).
#' @param strata A factor variable for stratified permutation (default: NULL).
#'
#' @return A matrix containing 1-vs-All posthoc test results.
#'
#' @examples
#' data <- ...  # TODO add data here
#' Tw2.posthoc.1vsAll.tests(data$dm, data$f)
#'
#' @export
Tw2.posthoc.1vsAll.tests <- function(dm, f, nrep = 999, strata = NULL) {
  Tw2.subset.test <- function(level) {
    fs <- factor(f == level)
    c(table(fs), Tw2.test(dm, fs, nrep = nrep, strata = strata))
  }

  res <- t(sapply(levels(f), Tw2.subset.test))
  colnames(res) <- c("N1", "N2", "p.value", "tw2.stat", "nrep")
  res
}
