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
#' @param strata A factor variable representing the strata. Default is NULL.
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
generic.distance.permutation.test =
  function(test.statistic, dm, f, nrep=999, strata = NULL){
    N = length(f)
    generate.permutation=function(){
      f[sample(N)]
    }

    if(!is.null(strata)){
      # map elements of each strata back to their positions in the factor variable
      strata.map = order(unlist(tapply(seq_along(f), strata, identity)))
      generate.permutation=function(){
        p = unlist(tapply(f,strata,sample)) # permute within strata
        p[strata.map]
      }
    }

    stats = c(test.statistic(dm, f),
              replicate(nrep,
                        test.statistic(dm, generate.permutation())))

    p.value = sum(stats>=stats[1])/(nrep+1)
    statistic = stats[1]
    list(p.value = p.value, statistic = statistic, nrep=nrep)
  }

#' Conducts a Tw2 distance-based permutation test for k-group differences
#'
#' This function is a specialized version of the \code{\link{generic.distance.permutation.test}},
#' specifically designed to use the Tw2 test statistic for k-group comparison.
#'
#' @param dm A distance metric (any arbitrary distance or dissimilarity metric)
#' @param f A factor variable representing the groups.
#' @param nrep The number of permutations to conduct. Default is 999.
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
Tw2.test <- function(dm, f, nrep = 999) {
  generic.distance.permutation.test(Tw2, dm = dm, f = f, nrep = nrep)
}

#' Conducts an adjusted WdS Distance-Based Permutation Test for k-Group differences with optional confounder elimination
#'
#' This function performs the WdS test statistic for k-group comparison using a given distance matrix.
#' Optionally, it can preprocess the provided data through the \code{a.dist} function using a specified
#' formula before performing the test, if 'data', 'dm', and 'formula' are provided.
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param nrep The number of permutations to conduct. Default is 999.
#' @param data (optional) Data frame to be used in conjunction with 'formula' and 'dm' for confounder elimination.
#' @param formula (optional) A formula to be used with 'data' for confounder elimination
#' @param strata (optional) A factor variable to be used for stratification in the permutation test.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{p.value}: The p-value of the test.
#'   \item \code{statistic}: The observed test statistic.
#'   \item \code{nrep}: The number of permutations performed.
#' }
#'
#' @seealso \code{\link{generic.distance.permutation.test}}, \code{\link{Tw2.test}},
#'          \code{\link{a.dist}}
#' @export
#' @examples
#' # Example with provided distance matrix
#' \dontrun{
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' f <- factor(c(rep("A", 5), rep("B", 5)))
#' aWdS.test(dm, f)
#'
#' # Example with optional confounder elimination with a.dist()
#' data(iris)
#' formula <- Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width
#'
#' distance_matrix = dist(iris[2:4], method="euclidean") # Numerical columns only
#'
#' aWdS.test(dm = distance_matrix, f = iris$Species)}
aWdS.test <- function(dm, f, nrep = 999, data = NULL, formula = NULL, strata=NULL) {

  if (!is.null(data) != !is.null(formula)) {
    stop("Both 'data' and 'formula' must be provided together for a.dist() processing.")
  }

  if (!is.null(data) && !is.null(formula)) {
    dm <- a.dist(formula, data)
  }

  test.results <- generic.distance.permutation.test(WdS, dm = dm, f = f, nrep = nrep,strata=strata)

  # Name of the hypothesis test.
  method <- "Multivariate Welch ANOVA"

  # Description of the data
  data.name <- paste0("Distance matrix with grouping factor ", deparse(substitute(f)))

  # Value of the parameter under the null hypothesis
  null.value <- 0
  attr(null.value, "names") <- "expected WdS under H0"

  # Direction of the alternative hypothesis relative to the null value
  alternative <- "two.sided"

  # Statistic value
  statistic <- test.results$statistic
  attr(statistic, "names") <- "WdS"

  # P value
  p.value <- test.results$p.value

  # Creating object of class 'htest'
  TEST <- list(method = method, data.name = data.name, null.value = null.value, alternative = alternative,
               statistic = statistic, p.value = p.value)
  class(TEST) <- "htest"

  return(TEST)
}

#' Conducts the original WdS Distance-Based Permutation Test for k-Group differences
#'
#' This function performs the WdS test statistic for k-group comparison using a given distance matrix.
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param nrep The number of permutations to conduct. Default is 999.
#' @param strata (optional) A factor variable to be used for stratified permutation.
#' @return A list containing:
#' \itemize{
#'   \item \code{p.value}: The p-value of the test.
#'   \item \code{statistic}: The observed test statistic.
#'   \item \code{nrep}: The number of permutations performed.
#' }
#'
#' @seealso \code{\link{generic.distance.permutation.test}}, \code{\link{Tw2.test}}
#' @export
#' @examples
#' # Example with provided distance matrix
#' \dontrun{
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' f <- factor(c(rep("A", 5), rep("B", 5)))
#' WdS.test(dm, f)
#'}
WdS.test = function(dm, f, nrep=999, strata=NULL){

  test.results <- generic.distance.permutation.test(WdS, dm = dm, f = f, nrep = nrep, strata=strata)

  # Name of the hypothesis test.
  method <- "Multivariate Welch ANOVA"

  # Description of the data
  data.name <- paste0("Distance matrix with grouping factor ", deparse(substitute(f)))

  # Value of the parameter under the null hypothesis
  null.value <- 0
  attr(null.value, "names") <- "expected WdS under H0"

  # Direction of the alternative hypothesis relative to the null value
  alternative <- "two.sided"

  estimate <- NA
  attr(estimate, "names") <- NA

  # Statistic value
  statistic <- test.results$statistic
  attr(statistic, "names") <- "WdS"

  # P value
  p.value <- test.results$p.value

  # Creating object of class 'htest'
  TEST <- list(method = method, data.name = data.name, null.value = null.value, alternative = alternative,
               statistic = statistic, p.value = p.value, estimate = NA)
  class(TEST) <- "htest"

  return(TEST)
}

