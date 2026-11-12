#' Test if an object is of class 'dist'
#'
#' This function checks whether the given object belongs to class 'dist'.
#'
#' @param x Object to test.
#'
#' @return Logical indicating whether the object is a distance matrix.
#' @export
#' @examples
#'
#' is.dist(as.dist(matrix(1:4, nrow = 2)))
#' is.dist(matrix(1:4, nrow = 2))
#'
is.dist <- function(x) any(class(x) == "dist")

.as_formula_data <- function(data, arg = "data", dm = NULL) {
  # Whole phyloseq inputs need extraction and ID matching before coercion.
  # Other input types retain their existing ordering and lookup behavior.
  if (inherits(data, "phyloseq")) { # if input is phyloseq object
    data <- .phyloseq_metadata(data, dm)
  }
  if (is.environment(data)) {
    return(data)
  }

  data <- tryCatch( # Catches conversion errors and stores the error object in data
    data.frame(data, check.names = FALSE), # Converts the object while preserving column names exactly—including spaces or special characters.
    error = function(e) e
  )

  if (inherits(data, "error")) { # If there was a conversion error
    stop(
      "'", arg, "' must be an environment, data frame, list, or coercible to a data frame.",
      call. = FALSE
    )
  }

  data
}

# Extract only explicitly supplied phyloseq metadata, matching distance labels.
## Return the metadata rows selected by sample ID and in distance-matrix order.
## Extra metadata samples are excluded. drop = FALSE ensures the result remains
## a data frame, even if it has only one column.
### For example:
### Input               | Sample order
  # Distance matrix	    | S3, S1, S2
  # Original metadata	  | S1, S2, S3, S4
  # Returned metadata	  | S3, S1, S2

## This matching prevents metadata from being assigned to the wrong samples during analysis.
.phyloseq_metadata <- function(data, dm) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required for phyloseq data inputs.", call. = FALSE)
  }
  # Keep phyloseq optional and report absent metadata without an accessor error.
  metadata <- phyloseq::sample_data(data, errorIfNULL = FALSE)
  if (is.null(metadata)) {
    stop("The phyloseq object has no sample metadata.", call. = FALSE)
  }
  # Row position alone cannot establish correspondence with a separately
  # computed distance matrix, so require unambiguous IDs on both inputs.
  labels <- attr(dm, "Labels")
  valid_ids <- function(ids) { # helper to check for unique, nonempty IDs
    !is.null(ids) &&  # exists
      !anyNA(ids) &&  # no missing
      all(nzchar(ids)) &&  # no empty
      !anyDuplicated(ids)  # no duplicates (are uniques)
  }
  # Will stop if the label count differs from the number of samples in dm
  if (!valid_ids(labels) || length(labels) != attr(dm, "Size")) {
    stop("Phyloseq data requires unique, nonempty sample labels on 'dm'.", call. = FALSE)
  }
  # Convert the phyloseq::sample_data() to standard data.frame()
  metadata <- data.frame(metadata, check.names = FALSE)
  if (!valid_ids(rownames(metadata))) {
    stop("Phyloseq metadata must have unique, nonempty sample IDs.", call. = FALSE)
  }
  # Validate the metadata’s sample IDs, which are stored as row names.
  missing_samples <- setdiff(labels, rownames(metadata))
  # Find IDs present in the distance matrix but absent from the metadata.
  # If any samples are missing, stop and list their IDs.
  # Every distance-matrix sample needs a matching metadata row.
  if (length(missing_samples)) {
    stop("Phyloseq metadata is missing samples found in 'dm': ",
         paste(missing_samples, collapse = ", "), call. = FALSE)
  }
  # The distance matrix defines the analysis sample set and order. Extra
  # metadata rows are safe to discard, but missing samples must fail above.

  # Return the metadata rows selected by sample ID and in distance-matrix order.
  # Extra metadata samples are excluded.
  # drop = FALSE ensures the result remains a data.frame, even if it has only one column.
  metadata[labels, , drop = FALSE]
}

.warn_deprecated_argument <- function(old, new, fn) {
  warning(
    "'", old, "' is deprecated in ", fn, "(); use '", new, "' instead.",
    call. = FALSE
  )
}

#' Calculate Sigma Squared for Distance Matrix
#'
#' This function computes the sigma squared statistic for a given distance matrix.
#'
#' @param dm Distance matrix.
#'
#' @return Sigma squared value.
#' @export
#' @examples
#'
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'
#'   dist.sigma2(dm)
#' }
#'
dist.sigma2 <- function(dm) {
  dd <- as.matrix(dm)
  dd[upper.tri(dd)] <- 0
  sum(dd^2) / nrow(dd) / (nrow(dd) - 1)
}

#' Calculate Distance-Based Pseudo-R-Squared for Residual Distances
#'
#' This function computes a distance-based pseudo-coefficient of determination
#' (pseudo-\eqn{R^2}) by comparing total variation in a raw distance matrix with
#' residual variation in another distance matrix.
#'
#' @param dm The original/raw distance matrix.
#' @param dm_residual A residual distance matrix after removing a model
#'   component, such as a grouping factor, adjustment variables, or their full
#'   model.
#'
#' @return Distance-based pseudo-\eqn{R^2}. Returns \code{NA_real_} when the
#'   original distance matrix has zero total variation.
#' @details Let \eqn{V(D) = [N(N-1)]^{-1}\sum_{i<j} d_{ij}^2}. The returned
#'   diagnostic is
#'   \deqn{R^2_{pseudo} = 1 - V(D_{residual}) / V(D_{original}).}
#'   The shared sample-size factor cancels in the ratio. This is a
#'   distance-based variation-reduction diagnostic. With univariate Euclidean
#'   distances and residual distances from the corresponding least-squares model
#'   with an intercept, it equals ordinary \eqn{R^2} (apart from numerical
#'   reconstruction). For arbitrary dissimilarities this equivalence need not hold.
#'   Values near zero indicate little
#'   reduction in distance variation, larger positive values indicate greater
#'   reduction, and negative values may occur when the residual distance matrix
#'   has more variation than the original distance matrix.
#' @export
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'
#'   dm_residual <- a.dist(
#'     dm = dm,
#'     adjustment = ~ SeqTech,
#'     data = ent
#'   )
#'   dist.goodness.of.fit(dm, dm_residual)
#' }
#'
dist.goodness.of.fit <- function(dm, dm_residual = NULL) {
  if (is.null(dm_residual)) {
    stop("'dm_residual' must be provided.")
  }
  if (!is.dist(dm) || !is.dist(dm_residual)) {
    stop("'dm' and 'dm_residual' must both be distance matrices of class 'dist'.")
  }
  if (attr(dm, "Size") != attr(dm_residual, "Size")) {
    stop("'dm' and 'dm_residual' must contain the same number of observations.")
  }

  # dist.sigma2() is proportional to distance-based total SS; the shared
  # sample-size scaling cancels in SS_residual / SS_total.
  ss_total <- dist.sigma2(dm)
  ss_residual <- dist.sigma2(dm_residual)

  # If all original distances are zero, the total variation denominator is zero,
  # so the pseudo-R2 is undefined rather than 0.
  goodness.of.fit <- if (ss_total == 0) NA_real_ else 1 - (ss_residual / ss_total)
  attr(goodness.of.fit, "names") <- "pseudo-R\u00B2"
  goodness.of.fit
}

#' Calculate Sum of Squares for Distance Matrix with Factor
#'
#' This function calculates the sum of squares for a given distance matrix,
#' weighted by a factor variable.
#'
#' @param dm2 Square distance matrix.
#' @param f Factor variable for weighting.
#'
#' @return Sum of squares matrix.
#' @export
#' @examples
#'
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   dist.ss2(as.matrix(dm)^2, f)
#' }
#'
dist.ss2 <- function(dm2, f) {
  K <- sapply(levels(f), function(lev) f == lev)
  t(K) %*% dm2 %*% K / 2
}

#' Calculate Group-wise Sigma Squared
#'
#' This function calculates the sigma squared statistic for each group defined
#' by a factor variable in a given distance matrix.
#'
#' @param dm Distance matrix.
#' @param f Factor variable for group definition.
#'
#' @return A named numeric vector of group-wise sigma squared values.
#' @export
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   dist.group.sigma2(dm, f)
#' }
dist.group.sigma2 <- function(dm, f) {
  diag(dist.ss2(as.matrix(dm)^2, f)) / table(f) / (table(f) - 1)
}

#' Calculate Cohen's d for Distance Matrix
#'
#' This function calculates the Cohen's d statistic for a distance matrix,
#' assuming two levels in the factor variable.
#'
#' @param dm Distance matrix.
#' @param f Factor variable with two levels.
#'
#' @return Cohen's d value if factor has exactly two levels; NULL otherwise.
#' @export
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
#'   dist.cohen.d(dm_12, f_12)
#' }
#'
dist.cohen.d <- function(dm, f) {
  if (nlevels(f) != 2) {
    return(NULL)
  }
  SS2 <- dist.ss2(as.matrix(dm^2), f)
  ns <- summary(f)
  N <- sum(ns)

  SST <- sum(SS2) / N
  SSW <- sum(diag(SS2) / ns)

  mean.diff <- (sqrt((ns[1] + ns[2]) / (ns[1] * ns[2]))) * sqrt(SST - SSW)

  sigmas <- diag(SS2) / ns / (ns - 1)
  s1 <- sigmas[1]
  s2 <- sigmas[2]

  mean.diff / sqrt(((ns[1] - 1) * s1 + (ns[2] - 1) * s2) / (sum(ns) - 2))
}

#' Compute a Covariate-Adjusted Distance Matrix
#'
#' This function takes a formula and distance matrix, projects out covariate
#' effects through principal-coordinate residualization, and returns the
#' resulting adjusted distance matrix.
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#'
#' @param adjustment A required right-hand-side formula specifying the
#'   covariates to remove, such as \code{~ SeqTech} or \code{~ Age + SeqTech}.
#'
#' @param data A dataset which contains the variables specified in \code{adjustment}.
#'   It may be an environment, data frame, list, or object coercible to a data
#'   frame, such as \code{phyloseq::sample_data()}. Row names should match the
#'   row names in distance matrix \code{dm}, with samples in the same order.
#'   If omitted or \code{NULL}, the calling environment is used. Formula
#'   variables are looked up in \code{data}, then in the formula environment.
#'   A whole \code{phyloseq} object is also accepted: sample metadata is
#'   extracted and matched to unique, nonempty \code{dm} sample labels. Extra
#'   metadata samples are excluded; missing samples cause an error.
#'   Adjustment variables must have one value per sample and no missing values.
#' @param tol Tolerance for eigenvalues. This is the cutoff for the eigenvalues
#'   to be considered zero. Default is 10^-8.
#' @param distance.diagnostics Logical indicating whether to attach a compact
#'   diagnostic summary of the eigendecomposition used to reconstruct the
#'   residual distance matrix. Default is \code{TRUE}.
#' @param keep.eigenvalues Logical indicating whether to store the raw
#'   eigenvalues in the diagnostic summary. Default is \code{FALSE} to keep
#'   repeated simulations and large analyses memory-light.
#' @param ... Deprecated named aliases: \code{formula} for \code{adjustment}
#'   and \code{formula_data} for \code{data}. These warn when used. Supplying
#'   both an alias and its replacement is an error. Other arguments, including
#'   unnamed arguments in \code{...}, are rejected. The first six positional
#'   arguments retain their meaning; the fourth argument is still \code{tol}.
#'
#' @return A covariate-adjusted distance matrix of class \code{dist}. When
#'   \code{distance.diagnostics = TRUE}, the returned object has a
#'   \code{"distance.diagnostics"} attribute containing one row of eigenvalue
#'   diagnostics.
#'
#' @details The \code{a.dist()} function only requires a right-hand side of the
#'   formula. Instead of the left-hand side, it uses the dissimilarity distance
#'   matrix \code{dm}. The function constructs a model matrix from the
#'   right-hand side (RHS) of the formula. After performing necessary matrix
#'   operations and eigen-decomposition, it calculates the Euclidean distances.
#'   It preserves the labels of the input dm.
#'
#'   During reconstruction, eigenvalues with \code{abs(lambda) < tol} are
#'   treated as zero, and remaining negative eigenvalues are set to zero before
#'   distances are recomputed. The \code{"distance.diagnostics"} attribute records
#'   how many eigenvalues were positive, negative, treated as zero by
#'   \code{tol}, and removed as negative after tolerance handling. These
#'   diagnostics describe construction of the residual distance matrix and are
#'   separate from goodness-of-fit statistics such as pseudo-\eqn{R^2}. If
#'   \code{keep.eigenvalues = TRUE}, the diagnostics also include the raw
#'   eigenvalues before tolerance and negative-value zeroing.
#'
#'   The diagnostic fields include:
#'   \itemize{
#'     \item \code{component}: \code{adjustment.only} for a direct
#'     \code{a.dist()} call; \code{WdS.test()} assigns \code{group.only},
#'     \code{adjustment.only}, or \code{full} according to the model, and
#'     \code{input} for inherited diagnostics on the supplied distance matrix.
#'     \item \code{model} and \code{tol}: full model description and tolerance.
#'     The model column is retained even when goodness-of-fit is disabled.
#'     \item \code{n.eigenvalues}: number of eigenvalues from the decomposition.
#'     \item \code{n.positive.raw}, \code{n.negative.raw}, and
#'     \code{n.zero.raw}: signs before tolerance handling.
#'     \item \code{n.zeroed.by.tol}: eigenvalues with
#'     \code{abs(lambda) < tol}, including exact zeros.
#'     \item \code{n.positive.zeroed.by.tol} and
#'     \code{n.negative.zeroed.by.tol}: positive and negative raw eigenvalues
#'     counted in \code{n.zeroed.by.tol}.
#'     \item \code{n.negative.removed}: negative eigenvalues remaining after
#'     tolerance handling, which are set to zero before reconstruction.
#'     \item \code{n.positive.retained} and \code{n.zero.final}: counts after
#'     tolerance and negative-value zeroing.
#'     \item \code{min.eigenvalue}, \code{max.eigenvalue},
#'     \code{sum.positive.raw}, \code{sum.negative.raw},
#'     \code{sum.negative.magnitude.raw}, and \code{prop.negative.raw}: compact
#'     summaries of the raw eigenspectrum. Here raw means before either
#'     zeroing step; final means after tolerance handling and removal of remaining
#'     negative eigenvalues. \code{prop.negative.raw} is the negative absolute
#'     magnitude divided by total absolute eigenvalue magnitude, not the
#'     fraction of eigenvalues that are negative. The counts satisfy
#'     \code{n.negative.raw = n.negative.zeroed.by.tol + n.negative.removed}.
#'     Individual raw eigenvalues are stored only when requested; their zeroing
#'     reason can be reconstructed from their sign and the recorded tolerance.
#'     \item \code{eigenvalues}: optional raw eigenvalue vector, present only
#'     when \code{keep.eigenvalues = TRUE}.
#'   }
#'
#'   This function refactors and generalizes functionality from
#'   \code{aPCoA::aPCoA()} function in the aPCoA package.
#' @references  Shi Y, Zhang L, Do KA, Peterson CB, Jenq RR. aPCoA: covariate
#'   adjusted principal coordinates analysis. Bioinformatics.
#'   2020;36(13):4099-4101. doi:10.1093/bioinformatics/btaa276
#'
#'   Shi Y (2021). aPCoA: Covariate Adjusted PCoA Plot. R package version 1.3,
#'   https://CRAN.R-project.org/package=aPCoA.
#'
#'   Please cite both the package and the paper when using this function.
#' @seealso \code{\link{WdS.test}}, \code{\link{dist.goodness.of.fit}}
#' @export
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'
#'   ## Use the full enterotype data after removing samples without Enterotype.
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'
#'   ## Adjust the distance matrix by sequencing technology, a categorical
#'   ## covariate stored in phyloseq::sample_data(ent).
#'   a.dm <- a.dist(
#'     dm = dm,
#'     adjustment = ~ SeqTech,
#'     data = phyloseq::sample_data(ent)
#'   )
#'   a.dm
#'   attr(a.dm, "distance.diagnostics")
#'
#'   ## Continuous covariates can also be used. Age has many missing values in
#'   ## enterotype, so this example builds a matching Age-complete distance
#'   ## matrix before using Age in the adjustment formula.
#'   ent_age <- phyloseq::subset_samples(
#'     enterotype,
#'     !is.na(Enterotype) & !is.na(Age)
#'   )
#'   dm_age <- phyloseq::distance(ent_age, method = "bray")
#'   a.dm.age <- a.dist(
#'     dm = dm_age,
#'     adjustment = ~ Age,
#'     data = phyloseq::sample_data(ent_age)
#'   )
#'   attr(a.dm.age, "distance.diagnostics")
#' }
#'
a.dist = function(dm, adjustment, data = parent.frame(), tol=10^-8,
                  distance.diagnostics = TRUE, keep.eigenvalues = FALSE,
                  ...)
{
  # Capture aliases before evaluating them so conflicts and misspellings fail
  # clearly. Evaluate their values in the user's frame, not a helper's frame.
  caller <- parent.frame()
  dots <- as.list(match.call(expand.dots = FALSE)$...)
  dot_names <- names(dots)
  if (length(dots) && (is.null(dot_names) || any(!nzchar(dot_names)))) {
    stop("Arguments in '...' must be named.", call. = FALSE)
  }
  unknown <- setdiff(dot_names, c("formula", "formula_data"))
  if (length(unknown)) {
    stop("Unused argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(dot_names)) {
    stop("Deprecated arguments must not be supplied more than once.", call. = FALSE)
  }
  if ("formula" %in% dot_names) {
    if (!missing(adjustment)) {
      stop("Use only one of 'adjustment' or deprecated 'formula'.", call. = FALSE)
    }
    .warn_deprecated_argument("formula", "adjustment", "a.dist")
    adjustment <- eval(dots$formula, envir = caller)
  } else if (missing(adjustment)) {
    stop("'adjustment' must be provided.", call. = FALSE)
  }
  if ("formula_data" %in% dot_names) {
    if (!missing(data)) {
      stop("Use only one of 'data' or deprecated 'formula_data'.", call. = FALSE)
    }
    .warn_deprecated_argument("formula_data", "data", "a.dist")
    data <- eval(dots$formula_data, envir = caller)
  }
  if (is.null(adjustment)) {
    stop("'adjustment' must be provided.", call. = FALSE)
  }
  # Keep the existing calculation's local name and stored diagnostic field.
  # This validates formula syntax without changing intercept or rank handling.
  formula <- .as_wdstest_adjustment(adjustment)
  if (is.null(data)) data <- caller
  .validate_distance_diagnostic_options(distance.diagnostics, keep.eigenvalues)
  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol < 0) {
    stop("'tol' must be a single non-negative numeric value.")
  }

  data <- .as_formula_data(data, arg = "data", dm = dm)
  Terms <- stats::terms(formula, data = data)
  # lhs <- formula[[2]]
  # lhs <- eval(lhs, data, parent.frame())
  # formula[[2]] <- NULL
  lhs <- dm
  rhs.frame <- stats::model.frame(formula, data, drop.unused.levels = TRUE,
                                  na.action = stats::na.fail) ## extract data.frame from a formula
  # Stop when adjustment values are missing
  ## Check that the number of rows in the right-hand-side data frame matches the
  ## number of samples in the distance matrix. This ensures that each sample has
  ## a corresponding row in the adjustment metadata.
  if (nrow(rhs.frame) != attr(dm, "Size")) {
    stop("Adjustment metadata must have one row per distance-matrix sample.", call. = FALSE)
  }
  rhs <- stats::model.matrix(formula, rhs.frame) ## expand and create dummy variables for the terms
  grps <- attr(rhs, "assign") ## tells which column from the data set each dummy variable came from
  qrhs <- qr(rhs) ##computes the QR decomposition of rhs matrix
  rhs <- rhs[, qrhs$pivot, drop = FALSE]
  rhs <- rhs[, seq_len(qrhs$rank), drop = FALSE] # 1:qrhs$rank fails for rank zero
  grps <- grps[qrhs$pivot][seq_len(qrhs$rank)] #same as above, seq_length() is more robust
  u.grps <- unique(grps)
  nterms <- length(setdiff(u.grps, 0L)) # previously length(u.grps) -1 assumed
  # that an intercept was present. This removes the intercept identifier (0) from
  # the count of usable terms, which is more robust to formulas with no intercept.
  # This allows formulas like ~ SeqTech - 1 to be valid, resulting in a 1 where
  # previoulsy it would have been 0.
  if (nterms < 1)
    stop("right-hand-side of formula has no usable terms")
  dmat <- as.matrix(lhs^2) #still has the labels of lhs
  X <- rhs
  X <- as.matrix(X)
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  A <- -1/2 * dmat
  J <- diag(nrow(X)) - matrix(rep(1/(nrow(X)), length(A)), nrow = nrow(A))
  E <- (diag(nrow(H)) - H) %*% J %*% A %*% J %*% (diag(nrow(H)) - H)

  # rownames(E) <- rownames(data)
  # colnames(E) <- rownames(data)

  # Correct rounding errors that cause E to be not symmetric
  E <- (E + t(E))/2

  eig <- eigen(E)
  lambda_raw <- eig$values
  lambda <- lambda_raw

  below_tol <- abs(lambda_raw) < tol
  if (any(below_tol)) {
    message(paste(
      sum(below_tol), "out of", length(lambda),
      "eigenvalues were smaller than 'tol' and were set to zero."
    ))
    if (missing(tol)) {
      message("The default 'tol' is 1e-8; use a smaller 'tol' to retain more small eigenvalues.")
    }
    lambda[below_tol] <- 0
  }

  negative <- lambda < 0
  if (any(negative)) {
    message(paste(
      sum(negative), "out of", length(lambda),
      "eigenvalues were negative and were set to zero."
    ))
    lambda[negative] <- 0
  }

  w <- t(t(eig$vectors) * sqrt(lambda))
  w <- stats::dist(w)
  attr(w, "Labels") <- attr(lhs, "Labels")
  if (isTRUE(distance.diagnostics)) {
    # Diagnostics are stored as metadata on the residual distance object. This
    # keeps a.dist() users informed without changing the dist values themselves.
    attr(w, "distance.diagnostics") <- .make_distance_diagnostics(
      lambda.raw = lambda_raw,
      lambda.final = lambda,
      formula = formula,
      tol = tol,
      keep.eigenvalues = keep.eigenvalues
    )
  }
  return(w)
}

.validate_distance_diagnostic_options <- function(distance.diagnostics, keep.eigenvalues) {
  if (!is.logical(distance.diagnostics) ||
      length(distance.diagnostics) != 1 ||
      is.na(distance.diagnostics)) {
    stop("'distance.diagnostics' must be TRUE or FALSE.")
  }
  if (!is.logical(keep.eigenvalues) ||
      length(keep.eigenvalues) != 1 ||
      is.na(keep.eigenvalues)) {
    stop("'keep.eigenvalues' must be TRUE or FALSE.")
  }
  if (isTRUE(keep.eigenvalues) && !isTRUE(distance.diagnostics)) {
    stop("'keep.eigenvalues = TRUE' requires 'distance.diagnostics = TRUE'.")
  }
}

.make_distance_diagnostics <- function(lambda.raw, lambda.final, formula, tol,
                                       keep.eigenvalues = FALSE) {
  raw.positive <- lambda.raw > 0
  raw.negative <- lambda.raw < 0
  raw.zero <- lambda.raw == 0
  below.tol <- abs(lambda.raw) < tol
  total.magnitude <- sum(abs(lambda.raw))
  negative.magnitude <- sum(abs(lambda.raw[raw.negative]))

  diagnostics <- data.frame(
    component = "adjustment.only",
    model = .display_adjustment_model(formula),
    tol = tol,
    n.eigenvalues = length(lambda.raw),
    n.positive.raw = sum(raw.positive),
    n.negative.raw = sum(raw.negative),
    n.zero.raw = sum(raw.zero),
    n.zeroed.by.tol = sum(below.tol),
    n.positive.zeroed.by.tol = sum(raw.positive & below.tol),
    n.negative.zeroed.by.tol = sum(raw.negative & below.tol),
    n.negative.removed = sum(raw.negative & !below.tol),
    n.positive.retained = sum(lambda.final > 0),
    n.zero.final = sum(lambda.final == 0),
    min.eigenvalue = min(lambda.raw),
    max.eigenvalue = max(lambda.raw),
    sum.positive.raw = sum(lambda.raw[raw.positive]),
    sum.negative.raw = sum(lambda.raw[raw.negative]),
    sum.negative.magnitude.raw = negative.magnitude,
    prop.negative.raw = if (total.magnitude == 0) NA_real_ else negative.magnitude / total.magnitude,
    keep.eigenvalues = isTRUE(keep.eigenvalues),
    stringsAsFactors = FALSE
  )

  if (isTRUE(keep.eigenvalues)) {
    # Store raw eigenvalues only when requested. They can be useful for deep
    # diagnosis but add memory cost when many tests are retained.
    diagnostics$eigenvalues <- I(list(lambda.raw))
  }

  diagnostics
}
