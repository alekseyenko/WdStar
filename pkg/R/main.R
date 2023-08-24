#' Calculate the Tw2 Statistic for Heteroscedastic Test
#'
#' This function calculates the Tw2 statistic for a heteroscedastic test
#'
#' @param dm A distance matrix.
#' @param f A factor variable.
#'
#' @return The calculated Tw2 statistic.
#'
#' @details
#' The Tw2 statistic is used for comparing means of k populations with potentially unequal observations.
#' The test is suitable for analysis of microbiome data using permutation testing.
#'
#' @examples
#' # Example usage
#' data <- ...  # TODO Provide an example distance matrix and factor variable
#' Tw2(data$dm, data$f)
#'
#' @references
#' Hamidi, Bashir, et al. "$ W_ {d}^{*} $-test: robust distance-based multivariate analysis of variance." Microbiome 7.1 (2019): 1-9.
#'
#' @export
Tw2 = function(dm, f){
  if(nlevels(f) != 2) return(NULL)
  SS2 = dist.ss2(as.matrix(dm^2),f)
  ns = summary(f)
  N = sum(ns)

  SST = sum(SS2)/N
  SSW1 = SS2[1,1]/ns[1]
  SSW2 = SS2[2,2]/ns[2]
  SSW = SSW1 + SSW2

  s1 = SSW1/(ns[1]-1)
  s2 = SSW2/(ns[2]-1)

  t.stat = (ns[1]+ns[2])/(ns[1]*ns[2])*(SST-SSW)/(s1/ns[1] + s2/ns[2])
  t.stat
}

#' Calculate the Wd* Statistic for Heteroscedastic Test
#'
#' This function calculates the Wd* statistic for a heteroscedastic test, extending Welch's solution to multivariate data.
#'
#' @param dm A distance matrix.
#' @param f A factor variable.
#'
#' @return The calculated Wd* statistic.
#'
#' @details
#' The Wd* statistic is used for comparing means of k populations with unequal variances, suitable for microbiome data.
#' The test uses permutation testing for significance estimation.
#'
#' @examples
#' # Example usage
#' data <- ...  # TODO provide an example distance matrix and factor variable
#' WdS(data$dm, data$f)
#'
#' @references
#' Hamidi, Bashir, et al. "$ W_ {d}^{*} $-test: robust distance-based multivariate analysis of variance." Microbiome 7.1 (2019): 1-9.
#'
#' @export
WdS = function(dm, f){
  ns = table(f)
  SS2 = dist.ss2(as.matrix(dm)^2, f)
  s2 = diag(SS2)/ns/(ns-1)
  W = sum(ns/s2)

  idxs = apply(utils::combn(levels(f), 2),2, function(idx) levels(f) %in% idx)

  Ws = sum(apply(idxs, 2,
                 function(idx) sum(ns[idx])/prod(s2[idx]) *
                   (sum(SS2[idx, idx])/sum(ns[idx]) - sum(diag(SS2[idx, idx])/ns[idx]))))
  k=nlevels(f)
  h = sum( (1-ns/s2/W)^2/(ns-1))
  Ws/W/(k-1)/(1+(2*(k-2)/(k^2-1))*h)
}

