#' Conducts a generic distance-based permutation test for k-group differences
#'
#' This function performs a distance-based permutation test using an arbitrary test statistic.
#' It is built on the principle of Welch's ANOVA and extends it to handle multivariate distance data.
#' It is particularly useful for analyzing microbiome data.
#'
#' @param test.statistic A function to calculate the test statistic.
#' @param dm A distance or dissimilarity matrix.
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
#' @param dm A distance or dissimilarity matrix.
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
Tw2.test = function(dm, f, nrep=999, strata=NULL){
  generic.distance.permutation.test(Tw2, dm = dm, f = f, nrep = nrep, strata=strata)
}

#' Conducts a WdS distance-based permutation test for k-group differences
#'
#' This function is a specialized version of the \code{\link{generic.distance.permutation.test}},
#' specifically designed to use the WdS test statistic for k-group comparison.
#'
#' @param dm A distance or dissimilarity matrix.
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
#' @seealso \code{\link{generic.distance.permutation.test}}, \code{\link{Tw2.test}}
#' @export
#' @examples
#' # TODO: Add examples
WdS.test = function(dm, f, nrep=999, strata=NULL){
  generic.distance.permutation.test(WdS, dm = dm, f = f, nrep = nrep, strata=strata)
}
