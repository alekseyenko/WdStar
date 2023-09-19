#' Perform posthoc tests for Tw2 statistic
#'
#' This function performs posthoc tests for the Tw2 statistic using a permutation approach.
#' The Tw2 statistic is an extension of the Welch ANOVA test for multivariate distances.
#' It is especially useful for analyzing microbiome data.
#'
#' @param dm A distance matrix representing the pairwise distances between observations.
#' @param f A factor variable indicating the group membership of each observation.
#' @param nrep Number of permutations to perform (default: 999).
#' @param strata A factor variable for stratified permutation. This allows for controlling confounders
#'               in repeated measures or other hierarchical designs (default: NULL).
#'
#' @return A matrix containing posthoc test results, with columns for level combinations,
#'         sample sizes, p-values, Tw2 statistics, and number of permutations.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- ... # TODO add data here
#' result <- Tw2.posthoc.tests(data$dm, data$f)
#' print(result)
#' }
#'
#' @seealso \url{https://github.com/alekseyenko/WdStar}
#' @export

Tw2.posthoc.tests <- function(dm, f, nrep = 999, strata = NULL) {
  dd <- as.matrix(dm)

  Tw2.subset.test <- function(include.levels) {
    subs <- f %in% include.levels
    c(
      include.levels,
      table(f[subs, drop = TRUE]),
      Tw2.test(dd[subs, subs], f[subs, drop = TRUE], nrep = nrep, strata = strata[subs])
    )
  }

  res <- t(utils::combn(levels(f), 2, Tw2.subset.test))
  colnames(res) <- c("Level1", "Level2", "N1", "N2", "p.value", "tw2.stat", "nrep")
  res
}

#' Perform 1-vs-All posthoc tests for Tw2 statistic
#'
#' This function performs 1-vs-All posthoc tests for the Tw2 statistic using a permutation approach.
#' This is useful for comparing each group against all other groups collectively.
#' This method is particularly useful in settings where multiple testing can lead to power loss,
#' such as microbiome data analyses.
#'
#' @param dm A distance matrix representing the pairwise distances between observations.
#' @param f A factor variable indicating the group membership of each observation.
#' @param nrep Number of permutations to perform (default: 999).
#' @param strata A factor variable for stratified permutation. This allows for controlling confounders
#'               in repeated measures or other hierarchical designs (default: NULL).
#'
#' @return A matrix containing 1-vs-All posthoc test results, with columns for sample sizes,
#'         p-values, Tw2 statistics, and number of permutations.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- ... # TODO add data here
#' result <- Tw2.posthoc.1vsAll.tests(data$dm, data$f)
#' print(result)
#' }
#'
#' @seealso \url{https://github.com/alekseyenko/WdStar}
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
