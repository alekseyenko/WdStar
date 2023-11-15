#' Conducts a generic distance-based permutation test for k-group differences
#'
#' This function performs a distance-based permutation test using an arbitrary test statistic.
#' It is built on the principle of Welch's ANOVA and extends it to handle multivariate distance data.
#' It is particularly useful for analyzing microbiome data.
#'
#' @param test.statistic A function to calculate the test statistic.
#' @param dm A distance metric (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param nrep The number of permutations to conduct. Default is 999.
#' @param strata An optional stratifying variable for restricted permutation.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{p.value}: The p-value of the test
#'   \item \code{statistic}: The observed test statistic
#'   \item \code{nrep}: The number of permutations performed
#' }
#'
#' @references Hamidi, Bashir, et al. "$ W_ {d}^{*} $-test: robust distance-based multivariate analysis of variance." Microbiome 7.1 (2019): 1-9.
#' @seealso \code{\link{Tw2.test}}, \code{\link{WdS.test}}
#' @export
#' @examples
#' # TODO: Add examples
generic.distance.permutation.test <-
  function(test.statistic, dm, f, nrep = 999, strata = NULL) {
    N <- length(f)
    generate.permutation <- function() {
      f[sample(N)]
    }

    if (!is.null(strata)) {
      # map elements of each strata back to their positions in the factor variable
      strata.map <- order(unlist(tapply(seq_along(f), strata, identity)))
      generate.permutation <- function() {
        p <- unlist(tapply(f, strata, sample)) # permute within strata
        p[strata.map]
      }
    }

    stats <- c(
      test.statistic(dm, f),
      replicate(
        nrep,
        test.statistic(dm, generate.permutation())
      )
    )

    p.value <- sum(stats >= stats[1]) / (nrep + 1)
    statistic <- stats[1]
    list(p.value = p.value, statistic = statistic, nrep = nrep)
  }

#' Conducts a Tw2 distance-based permutation test for k-group differences
#'
#' This function is a specialized version of the \code{\link{generic.distance.permutation.test}},
#' specifically designed to use the Tw2 test statistic for k-group comparison.
#'
#' @param dm A distance metric (any arbitrary distance or dissimilarity metric)
#' @param f A factor variable representing the groups.
#' @param nrep The number of permutations to conduct. Default is 999.
#' @param strata An optional stratifying variable for restricted permutation.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{p.value}: The p-value of the test
#'   \item \code{statistic}: The observed test statistic
#'   \item \code{nrep}: The number of permutations performed
#' }
#'
#' @seealso \code{\link{generic.distance.permutation.test}}, \code{\link{WdS.test}}
#' @export
#' @examples
#' # TODO: Add examples
Tw2.test <- function(dm, f, nrep = 999, strata = NULL) {
  generic.distance.permutation.test(Tw2, dm = dm, f = f, nrep = nrep, strata = strata)
}

#' Conducts a WdS Distance-Based Permutation Test for k-Group Differences with Optional aPCoA Preprocessing
#'
#' This function performs the WdS test statistic for k-group comparison using a given distance matrix.
#' Optionally, it can preprocess the provided data through the \code{aPCoA.dist} function using a specified
#' formula before performing the test, if 'data', 'dm', and 'formula' are provided.
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param nrep The number of permutations to conduct. Default is 999.
#' @param strata An optional stratifying variable for restricted permutation.
#' @param data (optional) Data frame to be used in conjunction with 'formula' and 'dm' for aPCoA preprocessing.
#' @param formula (optional) A formula to be used with 'data' for aPCoA preprocessing.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{p.value}: The p-value of the test.
#'   \item \code{statistic}: The observed test statistic.
#'   \item \code{nrep}: The number of permutations performed.
#' }
#'
#' @seealso \code{\link{generic.distance.permutation.test}}, \code{\link{Tw2.test}},
#'          \code{\link{aPCoA.dist}}
#' @export
#' @examples
#' # Example with provided distance matrix
#' \dontrun{
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' f <- factor(c(rep("A", 5), rep("B", 5)))
#' WdS.test(dm, f)
#'
#' # Example with optional aPCoA preprocessing
#' data(iris)
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' formula <- Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width
#' WdS.test(dm = dm, f = iris$Species, data = iris, formula = formula)
#' }
#'
WdS.test <- function(dm, f, nrep = 999, strata = NULL, data = NULL, formula = NULL) {
  if (!is.null(data) != !is.null(formula)) {
    stop("Both 'data' and 'formula' must be provided together for aPCoA processing.")
  }
  
  if (!is.null(data) && !is.null(formula)) {
    dm <- aPCoA.dist(formula, data)
  }
  
  generic.distance.permutation.test(WdS, dm = dm, f = f, nrep = nrep, strata = strata)
}
