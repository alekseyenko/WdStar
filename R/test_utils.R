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
#' @seealso \code{\link{WdS.test}}, \code{\link{generic.distance.permutation.test}}
#' @export
#' @examples
#' # TODO: Add examples
Tw2.test <- function(dm, f, nrep = 999) {
  generic.distance.permutation.test(Tw2, dm = dm, f = f, nrep = nrep)
}

#' Conducts distance-based multivariate Welch ANOVA
#' 
#' This function performs the \eqn{\mathnormal{W}_d^*} test statistic for k-group comparison
#' using a given distance matrix.
#'
#' This method uses permutation testing to establish the significance by
#' computing \eqn{\mathnormal{W}_d^*} statistic on \eqn{\mathnormal{m}} permutations of the original data and
#' estimate the significance as the fraction of times the permuted statistic is
#' greater than or equal to \eqn{\mathnormal{W}_d}.
#'
#' Our approach can accommodate multi-level factors, stratification (via
#' restricted permutations), multiple post-hoc testing scenarios, and covariate
#' adjustment/ elimination (via projection of residuals).
#'
#' If optional covariate adjustment parameters (\code{`formula`}, and
#' \code{`formula_data`}) are provided as input, they are passed to the \code{WdStar::a.dist()} 
#' function to project residual matrices and remove the effect of specified 
#' covariate(s) from distance matrix \code{dm} before performing the test. 
#' Users may choose to make adjustments within the function or use \code{WdStar::a.dist()} 
#' outside the function and then pass on the adjusted distance matrix to \code{WdStar::WdS.test()}.
#' 
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param nrep (optional) The number of permutations to conduct. Default is 999.
#' @param strata (optional) A factor variable representing strata. If specified,
#'   the test will perform p-value computations using stratified permutation.
#'   Default is NULL.
#' @param formula (optional) A right-hand side ONLY formula to be used with
#'  'formula_data' for confounder adjustment/elimination from `dm`. Results in 
#'  adjusted WdS statistic, Omega squared (ω²) effect size estimate, and p-value. 
#'  Note that you may use any data type including factor, character, integer, and
#'  numeric. Default is NULL
#' @param formula_data (optional) A data frame to be used in conjunction with
#' 'formula' for confounder adjustment. Default is parent.frame()
#' @return A list containing:
#' \itemize{
#'   \item \code{method}: name of the method used.
#'   \item \code{data.name}: a string describing the input data.
#'   \item \code{statistic}: observed WdS test statistic.
#'   \item \code{estimate}: Omega squared (ω²) effect size estimate.
#'   \item \code{p.value}: The p-value of the test.
#'   \item \code{parameter}: A list of 2-4 containing:  
#'          \itemize{
#'          \item \code{strata}: (optional) strata variable used to perform restricted permutations.
#'          \item \code{formula}: (optional) formula for residual projection to perform matrix adjustment.
#'          \item \code{dfb}: between degrees of freedom.
#'          \item \code{nrep}: number of permutations performed
#'          } 
#' }
#'
#' @seealso \code{\link{generic.distance.permutation.test}}, \code{\link{Tw2.test}},
#'          \code{\link{a.dist}}
#' @export
#' @examples
#' # The following is a simple example using the mtcars dataset to assess the  
#' # effect of gears on mpg, cyl, and disp (first three variables of the dataset):
#' data(mtcars)
#'
#' # The outcome could be a single variable or multiple variables (such as  
#' #  multidimensional omics data).
#'
#' ## This is an example of outcome with a single variable:
#' dm <- dist(mtcars$mpg, method="euclidean")
#'
#' ## This is an example of outcome with multiple variables:
#' dm <- dist(mtcars[1:3], method="euclidean")
#'
#' # Grouping/independent variable. You could use multiple variables here too.
#' f <- factor(mtcars$gear)
#'
#' # Basic multivariate test example ###########
#' #############################################
#' WdS.test(dm=dm, f=f)
#'
#' # Stratified example ########################
#' #############################################
#' strata <- factor(mtcars$vs)
#' WdS.test(dm=dm, f=f, strata=strata)
#'
#' # Covariate adjustment/elimination examples #
#' #############################################
#' ## Right-hand side adjustment formula to specify adjustment covariates. 
#' formula <- ~ wt + as.factor(am) 
#'
#' ## Adjustment example 1: pass unadjusted `dm` and formula to WdS.test()
#' WdS.test(dm=dm, f=f, formula=formula, formula_data=mtcars) ## Perform adjusted test
#'
#' ## Adjustment example 2: Create the adjusted distance matrix `a.dm` outside the function
#' a.dm <- a.dist(dm=dm, formula=formula, formula_data=mtcars) 
#' WdS.test(dm=a.dm, f=f) ## Perform adjusted test with `a.dm`
#'
WdS.test <- function(dm, f, nrep=999, strata=NULL, formula=NULL, formula_data=parent.frame()){
#   if (!is.null(formula) != !is.null(formula_data)) { #in certain cases didn't trigger error when formula was provided but data wasn't.
#    if (xor(is.null(formula), is.null(formula_data))) { # same issue as above. This has to do with default value of formula_data=parent.frame()
#  if (!is.null(formula) && (is.null(formula_data) || !(is.data.frame(formula_data) || is.list(formula_data)))) { #function proceeds even without formula
  has_formula <- !is.null(formula)
  has_data <- !is.null(formula_data) && (is.data.frame(formula_data) || is.list(formula_data))
  if (has_formula != has_data) {
      stop("Both 'formula' and 'formula_data' must be provided together for a.dist() adjustment processing.")
   }
  if (!is.null(formula)) {
    dm <- a.dist(dm, formula, formula_data)
  }
  test.results <- generic.distance.permutation.test(WdS, dm=dm, f=f, nrep=nrep, strata=strata)

    # Name of the hypothesis test
  if(!is.null(formula)) {
    method <- "Adjusted Distance-based Multivariate Welch ANOVA"
    adjustment_formula <- deparse(substitute(formula))
    } else {
    method <- "Distance-based Multivariate Welch ANOVA"
    adjustment_formula <- NULL
    }

  # Strata information
  if(!is.null(strata)){
    strata_selected <- "\n\tP-value is based on stratified permutations. WdS statistic and ω² are not computed by strata."}
    else {strata_selected <- NULL}

  # Description of the data
  data.name <- paste0(attr(dm, "method"), " distance matrix with ", attr(dm, "Size"), " observations and grouping factor with ", nlevels(f), " levels", strata_selected)
  
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

  # P-value
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

