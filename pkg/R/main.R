# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' Tw2
#'
#' @param dm Distance matrix
#' @param f Factor
#'
#' @return
#' @export
#'
#' @examples
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

#' This method computes Wd* statistic for distance matrix dm and factor f
#'
#' @param dm Distance matrix
#' @param f Factor
#'
#' @return
#' @export
#'
#' @examples
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

