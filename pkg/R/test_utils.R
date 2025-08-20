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

#' Conducts an adjusted WdS distance-based multivariate Welch ANOVA with optional confounder elimination
#'
#' This function performs the WdS test statistic for k-group comparison using a given distance matrix.
#' Optionally, it can preprocess the provided data through the \code{a.dist} function using a specified
#' formula before performing the test, if 'data', 'dm', and 'formula' are provided.  
#' 
#' If optional adjustment parameters are provided, it projects matrices to remove the effect of covariate. 
#'
#'
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param nrep The number of permutations to conduct. Default is 999.
#' @param strata (optional) A factor variable representing strata. If specified,
#'   the test will perform p-value computations using stratified permutation.
#'   Default is NULL.
#' @param formula (optional) A right-hand side ONLY formula to be used with 'data' for confounder adjustment of WdS statistic, Omega squared (ω²) effect size estimate, and p-value. 
#' @param formula_data (optional) A data frame to be used in conjunction with 'formula' for confounder adjustment.
#' @return A list containing:
#' \itemize{
#'   \item \code{method}: The name of the method used.
#'   \item \code{data.name}: A string containing the name of input data.
#'   \item \code{statistic}: The observed WdS test statistic.
#'   \item \code{estimate}: The Omega squared (ω²) effect size estimate.
#'   \item \code{nrep}: The number of permutations performed.
#'   \item \code{p.value}: The p-value of the test.
#'   \item \code{parameter}: a list of 2-4 containing strata (optional), 
#'   adjustment formula (optional), between degrees of freedom (dfb), and number of permutations performed
#'   (nrep).
#' }
#'
#' @seealso \code{\link{generic.distance.permutation.test}}, \code{\link{Tw2.test}},
#'          \code{\link{a.dist}}
#' @export
#' @examples
#' # Example with provided distance matrix
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' f <- factor(c(rep("A", 5), rep("B", 5)))
#' aWdS.test(dm, f)
#'
#' # Example with optional confounder elimination with a.dist()
#' data(iris)
#' formula <- ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width
#'
#' dm = dist(iris[2:4], method="euclidean") 
#'
#' aWdS.test(dm = dm, f = iris$Species, formula=formula)
aWdS.test <- function(dm, f, nrep = 999, strata=NULL, formula = NULL, formula_data = parent.frame()) {
  
#   if (!is.null(formula) != !is.null(formula_data)) { #in certain cases didn't trigger error when formula was provided but data wasn't. 
#    if (xor(is.null(formula), is.null(formula_data))) { # same issue as above. This has to do with default value of formula_data=parent.frame()
#  if (!is.null(formula) && (is.null(formula_data) || !(is.data.frame(formula_data) || is.list(formula_data)))) { #function proceeds even without formula
  has_formula <- !is.null(formula)
  has_data <- !is.null(formula_data) && (is.data.frame(formula_data) || is.list(formula_data))
  if (has_formula != has_data) {
      stop("Both 'formula' and 'formula_data' must be provided together for a.dist() processing.")
   }

  if (#!is.null(data) && 
      !is.null(formula)) {
    dm <- a.dist(dm, formula, formula_data)
  }

  test.results <- generic.distance.permutation.test(WdS, dm = dm, f = f, nrep = nrep, strata=strata)

  # Name of the hypothesis test
  if(!is.null(formula)) {
    method <- "Adjusted Distance-based Multivariate Welch ANOVA"
    adjustment_formula <- deparse(substitute(formula))
    } else {
    method <- "Distance-based Multivariate Welch ANOVA"
    adjustment_formula <- NULL
    }
  
  # Insert strata information
  if(!is.null(strata)){
    strata_selected <- "\n\tP-value is based on stratified permutations. WdS statistic and ω² are not stratified by strata."} 
    else {strata_selected <- NULL}
  
  # Description of the data
  data.name <- paste0(attr(dm, "method"), " distance matrix with ", attr(dm, "Size"), " observations and grouping factor with ", nlevels(f), " levels", strata_selected)
  #paste0("Distance matrix ", deparse(substitute(dm)), " with grouping factor ", deparse(substitute(f)), adjustment_formula, strata_selected)
  
  # Value of the parameter under the null hypothesis
  # null.value <- 0
  # attr(null.value, "names") <- "expected WdS under H0"

  # Direction of the alternative hypothesis relative to the null value
  # alternative <- "two.sided"
  
  # Statistic value
  statistic <- test.results$statistic
  attr(statistic, "names") <- "WdS"

  # Estimate using omega squared
  estimate <- (((length(levels(f))-1)*(unname(statistic)-1))/((length(levels(f))-1)*(unname(statistic) - 1) + attr(dm, "Size")))
  attr(estimate, "names") <- "Omega squared (ω²) effect size"
  
  # P value
  p.value <- test.results$p.value

  # Creating object of class 'htest'
  TEST <- list(method = method, data.name = data.name, # null.value = null.value, # alternative = alternative,
               statistic = statistic, p.value = p.value, 
               estimate = estimate,
               parameter = c(if(!is.null(strata)){c("strata"= deparse(substitute(strata)))},
                             if(!is.null(formula)){c("formula"= paste0(deparse(substitute(formula)),
                                                                       " with ",
                                                                       length(attr(terms(formula), "term.labels")),
                                                                       " terms"))},
                             c("dfb" = (length(levels(f))-1)),
                             c("number of permutations" = nrep)
                             )
               )
  class(TEST) <- c("wdstest", "htest")

  return(TEST)
}

#' Conducts the original WdS distance-based multivariate Welch ANOVA
#'
#' This function performs the WdS test statistic for k-group differences comparison using a given distance matrix.
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param nrep The number of permutations to conduct. Default is 999.
#' @param strata (optional) A factor variable representing strata. If specified,
#'   the test will perform p-value computations using stratified permutation.
#'   Default is NULL.
#' @return A list containing:
#' \itemize{
#'   \item \code{method}: The name of the method used.
#'   \item \code{data.name}: A string containing the name of input data.
#'   \item \code{statistic}: The observed WdS test statistic.
#'   \item \code{estimate}: The Omega squared (ω²) effect size estimate.
#'   \item \code{nrep}: The number of permutations performed.
#'   \item \code{p.value}: The p-value of the test.
#'   \item \code{parameter}: a list of 2 or 3 containing Strata (optional),
#'   between degrees of freedom (dfb), and number of permutations performed
#'   (nrep).
#' }
#'
#' @seealso \code{\link{generic.distance.permutation.test}}, \code{\link{Tw2.test}}
#' @export
#' @examples{
#' # Example with provided distance matrix
#' 
#' dm <- as.dist(matrix(runif(100), nrow = 10))
#' f <- factor(c(rep("A", 5), rep("B", 5)))
#' WdS.test(dm, f)
#' }
WdS.test <- function(dm, f, nrep=999, strata=NULL){

  test.results <- generic.distance.permutation.test(WdS, dm = dm, f = f, nrep = nrep, strata=strata)

  # Name of the hypothesis test
  method <- "Distance-based Multivariate Welch ANOVA" #"Multivariate Welch ANOVA"
  
  # Insert strata information
  if(!is.null(strata)){
    strata_selected <- "\n\tp-value is based on stratified permutations. WdS statistic and ω² are not adjusted for strata"} else {
      strata_selected <- NULL}
  
  # Description of the data
  data.name <- paste0(attr(dm, "method"), " distance matrix with ", attr(dm, "Size"), " observations and grouping factor with ", nlevels(f), " levels", strata_selected)
  # paste0("Distance matrix ", deparse(substitute(dm)), " with grouping factor ", deparse(substitute(f)), strata_selected)
      
  # Value of the parameter under the null hypothesis
  # null.value <- 0
  # attr(null.value, "names") <- "expected WdS under H0"

  # Direction of the alternative hypothesis relative to the null value
  # alternative <- "two.sided" # this is not necessary

  # Statistic value
  statistic <- test.results$statistic
  attr(statistic, "names") <- "WdS"

  # Estimate using omega squared
  estimate <- (((length(levels(f))-1)*(unname(statistic)-1))/((length(levels(f))-1)*(unname(statistic) - 1) + attr(dm, "Size")))
  attr(estimate, "names") <- "Omega squared (ω²) effect size"

  # P value
  p.value <- test.results$p.value

  # Creating object of class 'htest'
  TEST <- list(method = method, data.name = data.name, # null.value = null.value, # alternative = alternative,
               statistic = statistic, p.value = p.value, 
               estimate = estimate,
               parameter = c(if(!is.null(strata)){c("strata"= deparse(substitute(strata)))},
                             c("dfb" = (length(levels(f))-1)),
                             c("number of permutations" = nrep)
                             )
               )
  class(TEST) <- "htest"

  return(TEST)
}