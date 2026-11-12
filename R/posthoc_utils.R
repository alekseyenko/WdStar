#' Perform Pairwise Post Hoc Tests for the Tw2 Statistic
#'
#' This function performs pairwise post hoc tests for the \eqn{T_w^2} statistic
#' using permutation tests. The \eqn{T_w^2} statistic is a two-group
#' distance-based extension of Welch's test for multivariate data.
#'
#' @param dm A distance matrix representing the pairwise distances between
#'   observations.
#' @param f A factor variable indicating the group membership of each
#'   observation.
#' @param nrep Number of permutations to perform (default: 999).
#'
#' @return A matrix containing post hoc test results, with columns for level
#'   combinations, sample sizes, p-values, Tw2 statistics, and number of
#'   permutations.
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   Tw2.posthoc.tests(dm, f, nrep = 9)
#' }
#'
#' @seealso \url{https://github.com/alekseyenko/WdStar}
#' @export

Tw2.posthoc.tests <- function(dm, f, nrep = 999) {
  dd <- as.matrix(dm)

  Tw2.subset.test <- function(include.levels) {
    subs <- f %in% include.levels
    c(
      include.levels,
      table(f[subs, drop = TRUE]),
      Tw2.test(dd[subs, subs], f[subs, drop = TRUE], nrep = nrep)
    )
  }

  res <- t(utils::combn(levels(f), 2, Tw2.subset.test))
  colnames(res) <- c("Level1", "Level2", "N1", "N2", "p.value", "tw2.stat", "nrep")
  res
}

#' Perform 1-vs-All Post Hoc Tests for the Tw2 Statistic
#'
#' This function performs 1-vs-all post hoc tests for the \eqn{T_w^2} statistic
#' using permutation tests. This is useful for comparing each group against all
#' other groups collectively, especially when exhaustive pairwise testing may
#' reduce power.
#'
#' @param dm A distance matrix representing the pairwise distances between
#'   observations.
#' @param f A factor variable indicating the group membership of each
#'   observation.
#' @param nrep Number of permutations to perform (default: 999).
#'
#' @return A matrix containing 1-vs-all post hoc test results, with columns for
#'   sample sizes, p-values, Tw2 statistics, and number of permutations.
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   Tw2.posthoc.1vsAll.tests(dm, f, nrep = 9)
#' }
#'
#' @seealso \url{https://github.com/alekseyenko/WdStar}
#' @export

Tw2.posthoc.1vsAll.tests <- function(dm, f, nrep = 999) {
  Tw2.subset.test <- function(level) {
    fs <- factor(f == level)
    c(table(fs), Tw2.test(dm, fs, nrep = nrep))
  }

  res <- t(sapply(levels(f), Tw2.subset.test))
  colnames(res) <- c("N1", "N2", "p.value", "tw2.stat", "nrep")
  res
}
