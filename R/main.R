#' Calculate the Tw2 Statistic for Heteroscedastic Test
#'
#' This function calculates the \eqn{T_w^2} statistic to compare two groups
#' with potentially unequal sample sizes. It is suitable for microbiome and
#' other multivariate distance data.
#'
#' @param dm A distance matrix, representing dissimilarity between observations.
#' @param f A factor variable indicating the group for each observation.
#'
#' @return The calculated Tw2 statistic.
#'
#' @details
#' The function first checks if the factor variable has exactly two levels.
#' The \eqn{T_w^2} statistic is a modification of Hotelling's T-square statistic
#' adapted for heteroscedasticity in multivariate distance data. This function
#' returns the observed statistic only; use \code{\link{Tw2.test}} for
#' permutation-based significance testing.
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   keep <- f %in% c("1", "2")
#'   dm_12 <- as.dist(as.matrix(dm)[keep, keep])
#'   f_12 <- droplevels(f[keep])
#'
#'   Tw2(dm_12, f_12)
#' }
#'
#' @references
#' Alekseyenko AV. Multivariate Welch t-test on distances. \emph{Bioinformatics}.
#' 2016;32(23):3552-3558. doi:10.1093/bioinformatics/btw524
#'
#' @seealso \code{\link{Tw2.test}}
#'
#' @export
Tw2 <- function(dm, f) {
  if (nlevels(f) != 2) {
    return(NULL)
  }
  SS2 <- dist.ss2(as.matrix(dm^2), f)
  ns <- summary(f)
  N <- sum(ns)

  SST <- sum(SS2) / N
  SSW1 <- SS2[1, 1] / ns[1]
  SSW2 <- SS2[2, 2] / ns[2]
  SSW <- SSW1 + SSW2

  s1 <- SSW1 / (ns[1] - 1)
  s2 <- SSW2 / (ns[2] - 1)

  t.stat <- (ns[1] + ns[2]) / (ns[1] * ns[2]) * (SST - SSW) / (s1 / ns[1] + s2 / ns[2])
  t.stat
}

#' Calculate the Wd* Statistic for Heteroscedastic Test
#'
#' This function calculates the \eqn{\mathnormal{W}_d^*} statistic to compare
#' k groups with potentially unequal multivariate dispersions and sample sizes.
#' It is suitable for microbiome and other multivariate distance data.
#'
#' @param dm A distance matrix, representing dissimilarity between observations.
#' @param group A grouping factor with at least two levels, in the same
#'   observation order as \code{dm}. Supply the factor directly; metadata-column
#'   lookup is available in \code{\link{WdS.test}}, not this function.
#' @param ... Deprecated named argument \code{f}, an alias for \code{group}.
#'   It warns when used and cannot be combined with \code{group}. No other
#'   extra arguments are accepted.
#'
#' @return The calculated Wd* statistic.
#'
#' @details
#' This function is an extension of Welch's ANOVA statistic suitable for
#' multivariate distance data. The \eqn{\mathnormal{W}_d^*} statistic is
#' computed from pairwise squared distances among groups while accounting for
#' unequal group sizes and differences in multivariate spread across groups.
#' This function returns the observed statistic only; use
#' \code{\link{WdS.test}} for permutation-based significance testing.
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   group <- factor(phyloseq::sample_data(ent)$Enterotype)
#'   print(WdS(dm, group = group))
#' }
#'
#' @references Hamidi Bashir, Wallace K, Vasu C, Alekseyenko AV. The
#' \eqn{\mathnormal{W}_d^*}-test: robust distance-based multivariate analysis
#' of variance. \emph{Microbiome}. 2019;7:51.
#'
#' @seealso \code{\link{WdS.test}}, \code{\link{Tw2}},
#'   \url{https://github.com/alekseyenko/WdStar}
#'
#' @export
WdS <- function(dm, group, ...) {
  dots <- list(...)
  # Only the old grouping name is supported here. Adjustment and permutation
  # options belong to WdS.test(), not this statistic-only function.
  dot_names <- names(dots)
  if (length(dots) && (is.null(dot_names) || any(!nzchar(dot_names)) ||
                       any(dot_names != "f"))) {
    stop("Only the deprecated named 'f' argument is accepted in '...'.", call. = FALSE)
  }
  if (length(dots)) {
    if (length(dots) > 1L) {
      stop("Argument 'f' was supplied more than once.", call. = FALSE)
    }
    if (!missing(group)) {
      stop("Use only one of 'group' or deprecated 'f'.", call. = FALSE)
    }
    .warn_deprecated_argument("f", "group", "WdS")
    group <- dots$f
  }
  if (missing(group)) stop("'group' must be provided.", call. = FALSE)
  # WdS.test() prepares the factor once before permutations. Keep only cheap
  # checks here rather than repeating metadata lookup or factor conversion.
  if (!is.factor(group) || nlevels(group) < 2L) {
    stop("'group' must be a factor with at least two levels.", call. = FALSE)
  }
  .validate_wdstest_vector_length(group, attr(dm, "Size"), "group")
  ns <- table(group)
  SS2 <- dist.ss2(as.matrix(dm)^2, group)
  s2 <- diag(SS2) / ns / (ns - 1)
  W <- sum(ns / s2)

  idxs <- apply(utils::combn(levels(group), 2), 2, function(idx) levels(group) %in% idx)

  Ws <- sum(apply(
    idxs, 2,
    function(idx) {
      sum(ns[idx]) / prod(s2[idx]) *
        (sum(SS2[idx, idx]) / sum(ns[idx]) - sum(diag(SS2[idx, idx]) / ns[idx]))
    }
  ))
  k <- nlevels(group)
  h <- sum((1 - ns / s2 / W)^2 / (ns - 1))
  Ws / W / (k - 1) / (1 + (2 * (k - 2) / (k^2 - 1)) * h)
}
