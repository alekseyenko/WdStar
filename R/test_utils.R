#' Conduct a Generic Distance-Based Permutation Test
#'
#' This function performs a distance-based permutation test using an arbitrary
#' test statistic. It is useful when a statistic such as \code{\link{WdS}} or
#' \code{\link{Tw2}} has been defined separately and permutation-based
#' significance testing is needed for multivariate distance data.
#'
#' @param test.statistic A function that takes \code{dm} and \code{f} and
#'   returns the observed test statistic.
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param nrep A single non-negative integer giving the number of permutations
#'   to conduct. Default is 999.
#' @param strata A factor variable representing permutation strata. If supplied,
#'   permutations are restricted within strata. Supply one non-missing value
#'   per observation. Singleton strata remain fixed. Default is \code{NULL}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{p.value}: The p-value of the test
#'   \item \code{statistic}: The observed test statistic
#'   \item \code{nrep}: The number of permutations performed
#' }
#'
#' @references Hamidi B, Wallace K, Vasu C, Alekseyenko AV. The
#' \eqn{\mathnormal{W}_d^*}-test: robust distance-based multivariate analysis
#' of variance. \emph{Microbiome}. 2019;7:51.
#' @seealso \code{\link{Tw2.test}}, \code{\link{WdS.test}}
#' @export
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   generic.distance.permutation.test(WdS, dm = dm, f = f, nrep = 9)
#' }
generic.distance.permutation.test =
  function(test.statistic, dm, f, nrep = 999, strata = NULL){
    N <- length(f)
    # Check that the number of permutations is a single non-negative integer.
    if (!is.numeric(nrep) || length(nrep) != 1L || !is.finite(nrep) ||
        nrep < 0 || nrep != floor(nrep)) {
      stop("'nrep' must be a single non-negative integer.", call. = FALSE)
    }
    # Reject missing group labels because they cannot be permuted. Users need to
    # subset the distance matrix to exclude missing labels before calling this
    # function.
    if (anyNA(f)) stop("'f' must not contain missing values.", call. = FALSE)
    generate.permutation <- function() f[sample(N)]
    # Restricted permutations require one non-missing stratum per observation.
    if (!is.null(strata)) {
      if (length(strata) != N || anyNA(strata)) {
        stop("'strata' must have one non-missing value per observation.", call. = FALSE)
      }
      # map elements of each strata back to their positions in the factor variable
      # tapply() groups its output by stratum. Map the concatenated results
      # back to the original observation order after shuffling within strata.
      # using seq_len(N) instead of seq_along(f) since we already pulled the
      # length of f into N.
      strata.map <- order(unlist(tapply(seq_len(N), strata, identity)))
      generate.permutation <- function() {
        # Shuffle observation indices rather than labels. This avoids both
        # sample()'s numeric-singleton special case and factor coercion when
        # tapply()/unlist() combine their results.
        p <- unlist(tapply(seq_len(N), strata, # splits those positions into strata
                           # sample))
        # This implementation passed each stratum's labels to sample(). This worked
        # except for singleton numeric labels caused issue, for example,
        # sample(5) shuffled 1:5 instead of returning only 5.
                           function(idx) { #Receive position belonging to one strata
          idx[sample.int(length(idx))] # Generates a random ordering within that stratum
        # and idx[] returns that stratum's shuffled observation position
        }), use.names = FALSE) # Prevents unlist() from adding stratum-generated names
        f[p[strata.map]] #  Restores the original observation order after shuffling within strata
        # f[] returns the group labels at the shuffled positions
      }
    }

    observed <- test.statistic(dm, f)
    # replicate(0, ...) returns a list; keep the zero-permutation case numeric.
    permuted_stats <- if (nrep > 0) {
      replicate(nrep, test.statistic(dm, generate.permutation()))
    } else {
      numeric(0)
    }
    stats <- c(observed, permuted_stats)

    p.value = sum(stats>=stats[1])/(nrep+1)
    statistic = stats[1]
    list(p.value = p.value, statistic = statistic, nrep=nrep)
  }

#' Conduct a Tw2 Distance-Based Permutation Test
#'
#' This function is a specialized version of
#' \code{\link{generic.distance.permutation.test}} that uses the
#' \code{\link{Tw2}} statistic for a two-group comparison.
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
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
#' @seealso \code{\link{Tw2}}, \code{\link{WdS.test}},
#'   \code{\link{generic.distance.permutation.test}}
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
#'   Tw2.test(dm_12, f_12, nrep = 9)
#' }
Tw2.test <- function(dm, f, nrep = 999) {
  generic.distance.permutation.test(Tw2, dm = dm, f = f, nrep = nrep)
}

#' Conduct Distance-Based Multivariate Welch ANOVA
#'
#' This function performs the multivariate Welch ANOVA \eqn{\mathnormal{W}_d^*} test using
#' optional stratified restricted permutations and confounder elimination for
#' k-group comparison using a given distance matrix.
#'
#' This method uses permutation testing to establish the significance by
#' computing \eqn{\mathnormal{W}_d^*} statistic on \eqn{\mathnormal{m}}
#' permutations of the original data and estimate the significance as the
#' fraction of times the permuted statistic is greater than or equal to
#' \eqn{\mathnormal{W}_d^*} (see \code{\link{WdS}} for details).
#'
#' Our approach can accommodate multi-level factors, stratification (via
#' restricted permutations), multiple post-hoc testing scenarios, and covariate
#' adjustment/ elimination by projection of residuals from the distance matrix
#' (see \code{\link{a.dist}} for details).
#'
#' The preferred interface is \code{group} for the tested grouping variable,
#' \code{adjustment} for optional covariates to remove, and \code{data} for
#' sample metadata. Deprecated function parameter names \code{f}, \code{formula},
#' and \code{formula_data} may still be supplied through \code{...}.
#'
#' Printed output shows parameters, WdS, p-value, effect size estimate, and
#' selected goodness-of-fit components, in that order. The tested group label is
#' limited to 30 characters for display; full labels are retained in the result.
#' Variable names and conventional mathematical notation preserve their case.
#'
#' The effect size \eqn{R^2}omega-squared is estimated from the observed WdS statistic:
#' \deqn{\widehat{\omega}^2 =
#' \frac{(k-1)(W_d^*-1)}{(k-1)(W_d^*-1)+N},}
#' where \eqn{k} is the number of tested group levels and \eqn{N} the sample size.
#' This is the ordinary one-way ANOVA omega-squared transformation with WdS
#' substituted for the classical F statistic. In ordinary ANOVA the equivalent
#' sums-of-squares expression is
#' \eqn{[SS_{between}-(k-1)MS_{within}]/[SS_{total}+MS_{within}]}.
#' In adjusted tests WdS is calculated on the adjustment-residual distances.
#'
#'
#'
#' Goodness-of-fit is reported as a distance-based pseudo-\eqn{R^2}.
#' For original distances \eqn{D} and model-residual distances \eqn{D_{Res}},
#' the basic goodness-of-fit measure is
#' \deqn{R^2 = 1 - \frac{\sum_{i<j} d_{Res,ij}^2}{\sum_{i<j} d_{ij}^2}.}
#' For univariate Euclidean distances and the corresponding least-squares model
#' with an intercept, this equals ordinary regression \eqn{R^2} when numerical
#' eigenvalue truncation does not alter the residual reconstruction. For
#' multivariate Euclidean data the analogous total sums squared variation across
#' response dimensions.
#'
#' The \code{goodness.components} table stores all components computed directly
#' and intermediate quantities derived from the residual distance matrices; For
#' adjusted \code{"auto"} tests, this includes \code{adjustment.only}, \code{full},
#' \code{semi.partial}, and \code{partial}, but not \code{group.only}, which
#' requires an additional residual distance matrix.
#'
#' The \code{goodness.components$displayed} column identifies selected values
#' printed in \code{goodness.of.fit}. With \code{goodness = "auto"}, unadjusted
#' tests report \code{group.only}, the variation explained by the tested group
#' alone. Adjusted tests report \code{semi.partial}, the additional fraction of
#' total distance variation explained by the tested group after the adjustment
#' variables.
#'
#' Distance diagnostics are stored separately from goodness-of-fit components,
#' in table \code{distance.diagnostics}. When \code{distance.diagnostics = TRUE},
#' \code{WdS.test()} collects the diagnostic summaries attached by \code{\link{a.dist}()}
#' for each residual distance matrix it constructs. These diagnostics record the
#' tolerance used, how many eigenvalues were treated as zero by \code{tol}, and
#' how many negative eigenvalues were removed during reconstruction. Use
#' \code{keep.eigenvalues = TRUE} to store the raw eigenvalue vectors for deeper
#' diagnosis; the default keeps only the compact summary.
#'
#' Tables \code{goodness.components} and \code{distance.diagnostics} use
#' \code{component} and \code{model} to identify their rows. For adjustment
#' variables \eqn{Z} and tested group \eqn{X}, the components and their
#' interpretations are as follows:
#' \itemize{
#'   \item \code{group.only}: \eqn{R^2_X}, variation explained by \eqn{X}
#'   alone as a fraction of original total distance variation, ignoring
#'   adjustment variables.
#'   \item \code{adjustment.only}: \eqn{R^2_Z}, variation explained by \eqn{Z}
#'   alone as a fraction of original total distance variation, ignoring \eqn{X}.
#'   \item \code{full}: \eqn{R^2_{X+Z}}, variation explained by \eqn{X} and
#'   \eqn{Z} together as a fraction of original total distance variation. The
#'   full-model residual matrix is computed from the original \code{dm} with a
#'   combined formula. It is not computed by applying \code{a.dist()} a second
#'   time to an already residualized distance matrix.
#'   \item \code{semi.partial}: \eqn{R^2_{X+Z} - R^2_Z}, additional variation
#'   explained by adding the tested group after adjustment, relative to original
#'   total distance variation. A negative semi.partial means full is less than
#'   adjustment.only.
#'   \item \code{partial}: \eqn{(R^2_{X+Z} - R^2_Z)/(1 - R^2_Z)}, additional
#'   variation explained by adding the tested group after adjustment, relative
#'   to distance variation remaining after adjustment. Partial is returned as
#'   \code{NA} when adjustment leaves essentially no residual variation.
#' }
#'
#' The components are not an additive partition: group.only plus adjustment.only
#' generally differs from full because their explained variation can overlap.
#' For non-Euclidean distances, removal of negative eigenvalues during
#' reconstruction can disrupt usual nested-model inequalities. For group.only,
#' adjustment.only, and full, values near zero indicate little distance-variation
#' reduction; negative values mean the reconstructed residual distances have
#' more variation than the original distances. Strata restrict permutations for
#' the p-value and do not enter these goodness-of-fit calculations.
#'
#' The rows semi-partial and partial are only included in goodness table as
#' those are derivations of other computations, and not residual distance
#' matrices to have associated eigendecomposition. These derived components have
#' different denominators: semi-partial uses original total variation as its
#' denominator; partial uses variation remaining after adjustment. For example,
#' if adjustment.only is 0.40 and full is 0.55, semi.partial is 0.15 (15
#' percentage points of original total variation), whereas partial is
#' 0.15 / 0.60 = 0.25 (25 percent of the variation remaining after adjustment).
#'
#' The diagnostic rows correspond to residual distance matrices actually
#' constructed by the call. For example, an unadjusted test with the default
#' \code{goodness = "auto"} stores a \code{"group.only"} row. When no adjustment
#' formula is supplied, \code{adjustment.only} is unavailable; if requested, \code{full},
#' \code{semi.partial}, and \code{partial} reduce to the same value as \code{group.only}.
#' If \code{goodness = "none"} and no adjustment formula is supplied, no residual
#' distance matrix is constructed and \code{distance.diagnostics} is \code{NULL},
#' unless the input \code{dm} already carries diagnostics from a previous
#' \code{a.dist()} call. In that case, the input diagnostics are preserved with
#' component \code{"input"}.
#' In case of an adjusted test with \code{goodness =
#' "auto"} stores rows for \code{"adjustment.only"} and \code{"full"};
#' \code{"group.only"} is added only when requested, such as with \code{goodness
#' = "all"}.
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#' @param group The tested grouping variable. This may be a vector/factor, an
#'   unquoted column name in \code{data}, or a length-one character column name
#'   in \code{data}.
#' @param nrep (optional) The number of permutations to conduct. Default is 999.
#' @param strata (optional) A factor variable representing strata. If specified,
#'   the test will perform p-value computations using stratified permutation.
#'   This may be a vector/factor, an unquoted column name in \code{data}, or a
#'   length-one character column name in \code{data}. Default is \code{NULL}.
#' @param adjustment (optional) A right-hand-side formula specifying covariates
#'   to remove from \code{dm} before testing \code{group}. Formula is passed to
#'   \code{\link{a.dist}} to project residual matrices before the test is
#'   performed. Users may also make adjustments outside the function with
#'   \code{\link{a.dist}} and then pass the adjusted distance matrix to
#'   \code{WdS.test()}. Note that you may use any data type supported by
#'   \code{model.matrix()}, including factor, character, integer, and numeric
#'   variables. Default is \code{NULL}.
#' @param data (optional) An environment, data frame, list, or object coercible
#'   to a data frame, such as \code{phyloseq::sample_data()}. Variables in
#'   \code{group}, \code{strata}, and \code{adjustment} are first looked up
#'   here. If omitted or \code{NULL}, the calling environment is used. Metadata
#'   rows and external vectors must follow the sample order in \code{dm}; this
#'   function does not reorder these inputs by row names. Adjustment variables
#'   absent from \code{data} are looked up in the formula environment. A whole
#'   \code{phyloseq} object is also accepted: sample metadata is extracted and
#'   matched to unique, nonempty \code{dm} sample labels. Extra metadata samples
#'   are excluded; missing samples cause an error.
#' @param goodness Which distance-based pseudo-\eqn{R^2} goodness-of-fit
#'   component to report. Options are \code{"auto"}, \code{"group.only"},
#'   \code{"adjustment.only"}, \code{"full"}, \code{"semi.partial"},
#'   \code{"partial"}, \code{"all"}, and \code{"none"}. The default
#'   \code{"auto"} reports \code{"group.only"} for unadjusted tests and
#'   \code{"semi.partial"} for adjusted tests. Multiple explicit components may
#'   be requested with a character vector, except for \code{"auto"},
#'   \code{"all"}, and \code{"none"}. Use \code{"none"} to skip goodness-of-fit
#'   calculations.
#' @param distance.diagnostics Logical indicating whether to store compact
#'   eigendecomposition diagnostics for residual distance matrices constructed
#'   inside \code{WdS.test()}. Default is \code{TRUE}. \code{goodness} and
#'   \code{distance.diagnostics} are placed after \code{...} so older calls that
#'   pass \code{a.dist()} arguments through \code{...} keep their positional
#'   behavior; users should pass these newer controls by name.
#' @param keep.eigenvalues Logical indicating whether to include raw eigenvalue
#'   vectors in \code{distance.diagnostics}. Default is \code{FALSE} to avoid
#'   extra memory use in large analyses and simulations.
#' @param ... Additional arguments passed to \code{\link{a.dist}()}, such as
#'   \code{tol} (tolerance for eigenvalues; default \code{1e-8}). These
#'   arguments are used when residual distance matrices are computed for
#'   adjustment or goodness-of-fit components. Deprecated named arguments
#'   \code{f}, \code{formula}, and \code{formula_data} are also accepted here
#'   for backward compatibility. Unexpected arguments are rejected even when
#'   \code{goodness = "none"}; a single unnamed extra argument is treated as
#'   \code{tol} for positional compatibility.
#' @return A list containing:
#' \itemize{
#'   \item \code{method}: name of the method used.
#'   \item \code{data.name}: a string describing the input data with the full group label.
#'   \item \code{data.info}: distance method, observation count, full group label,
#'   and number of group levels, used to format the data description.
#'   \item \code{statistic}: observed WdS test statistic.
#'   \item \code{estimate}: omega-squared effect size estimate calculated from WdS (see Details).
#'   \item \code{goodness.of.fit}: requested/default distance-based
#'   pseudo-\eqn{R^2} component.
#'   \item \code{goodness.components}: data frame of the goodness-of-fit
#'   components computed while evaluating \code{goodness.of.fit}. This may
#'   include intermediate values that are not printed by default.
#'   \item \code{distance.diagnostics}: data frame of eigendecomposition
#'   diagnostics for residual distance matrices constructed by \code{a.dist()},
#'   or \code{NULL} when no diagnostics were collected.
#'   \item \code{p.value}: The p-value of the test.
#'   \item \code{parameter}: A list of 2-4 containing:
#'          \itemize{
#'          \item \code{strata}: (optional) strata variable used for restricted permutations.
#'          \item \code{adjustment}: (optional) formula for residual projection to perform distance adjustment.
#'          \item \code{dfb}: between degrees of freedom.
#'          \item \code{nrep}: number of permutations performed
#'          }
#' }
#'
#' @seealso \code{\link{WdS}}, \code{\link{generic.distance.permutation.test}},
#'   \code{\link{Tw2.test}}, \code{\link{a.dist}},
#'   \code{\link{dist.goodness.of.fit}}
#' @references Hamidi Bashir, Wallace K, Vasu C, Alekseyenko AV. The
#' \eqn{\mathnormal{W}_d^*}-test: robust distance-based multivariate analysis
#' of variance. \emph{Microbiome}. 2019;7:51.
#' @importFrom stats terms
#' @export
#' @examples
#' # We use the optional package phyloseq in these examples for their enterotype dataset;
#' # Examples run when phyloseq is already installed.
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq") # load phyloseq's enterotype data
#'   ## Use the full enterotype data after removing samples without Enterotype.
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'
#'   ###########################################################
#'   ###########################################################
#'   ## (Eg 1) Basic multivariate test:
#'   res_unadjusted <-
#'     WdS.test(dm = dm,
#'              group = Enterotype,
#'              data = ent,
#'              nrep = 1) # Change nrep to 999 or more for a real analysis
#'   ### Data parameter can be given as a data frame, list, or phyloseq object.
#'
#'   ###########################################################
#'   ###########################################################
#'   ## (Eg 2) Stratified permutation example:
#'   ### This restricts permutations within sequencing technology.
#'   ### Stratification does not change statistics, estimate, or goodness-of-fit.
#'   res_stratified <-
#'     WdS.test(dm = dm,
#'              group = Enterotype,
#'              strata = SeqTech,
#'              data = phyloseq::sample_data(ent),
#'              nrep = 1)
#'
#'   ###########################################################
#'   ###########################################################
#'   ## (Eg 3) Covariate adjustment by a categorical covariate (sequencing technology):
#'   res_adj_categorical <-
#'     WdS.test(dm = dm,
#'              group = Enterotype,
#'              adjustment = ~SeqTech,
#'              data = ent,
#'              goodness = "all", ## goodness = "none"
#'              nrep = 1)
#'
#'   ### Goodness-of-Fit:
#'   #### By default, the unadjusted test's goodness.of.fit reports the
#'   ####    distance-based pseudo-R-squared for the tested group.
#'   #### Use goodness = "all" to request all available components, or
#'   ####   goodness = "none" to skip the computations.
#'   #### goodness.components contains all components calculated for this test.
#'   res_adj_categorical$goodness.components
#'
#'   ### Get distance diagnostics
#'   res_adj_categorical$distance.diagnostics
#'
#'   ###########################################################
#'   ###########################################################
#'   ## (Eg 4) Covariate adjustment by a continuous covariate (age):
#'   ### In the enterotype dataset, age has many missing values so this example
#'   ### builds a matching age-complete distance matrix before using Age in the
#'   ### adjustment formula.
#'   ent_age <-
#'     phyloseq::subset_samples(enterotype,
#'                              !is.na(Enterotype) & !is.na(Age))
#'   dm_age <- phyloseq::distance(ent_age, method = "bray")
#'
#'   res_adj_continuous <-
#'     WdS.test(dm = dm_age,
#'              group = Enterotype,
#'              adjustment = ~Age,
#'              data = ent_age,
#'              nrep = 1)
#'
#'   ###########################################################
#'   ###########################################################
#'   ## (Eg 5) Covariate adjustment by building distance matrix outside the
#'   ## function, then pass dm to WdS.test():
#'   ### First build adjusted dm
#'   a.dm <-
#'     a.dist(dm = dm,
#'            adjustment = ~SeqTech,
#'            data = phyloseq::sample_data(ent))
#'
#'   ### Then use adjusted dm in testing
#'   res_eg5 <-
#'     WdS.test(dm = a.dm,
#'              group = Enterotype,
#'              data = phyloseq::sample_data(ent),
#'              nrep = 1)
#'
#'   ### You can get the diagnostics of adjusted dm
#'   # attr(a.dm, "distance.diagnostics")
#'
#'   print(res_unadjusted)
#'   print(res_stratified)
#'   print(res_adj_categorical)
#'   print(res_adj_continuous)
#' }
#'
#' # goodness and diagnostic controls are placed after ... so older calls that pass
#' # a.dist() arguments through ... keep their positional behavior; users should
#' # pass these newer controls by name.
WdS.test <- function(dm, group, nrep=999, strata=NULL, adjustment=NULL, data=NULL, ...,
                     goodness="auto", distance.diagnostics = TRUE,
                     keep.eigenvalues = FALSE){
  # Capture the call before evaluating group/strata. This is what lets users
  # write group = Enterotype when Enterotype lives inside data.
  call <- match.call(expand.dots = FALSE)
  parent <- parent.frame()
  dot_exprs <- as.list(call$...)

  # Helper function to validate inputs: normalize the new names and the old
  # aliases once at the top. The below code should only have to think about
  # resolved vectors, one adjustment formula, one metadata source, and a clean
  # list for a.dist().
  args <- .resolve_wdstest_arguments(
    group = if (!missing(group)) substitute(group),
    strata = if (!missing(strata)) substitute(strata),
    adjustment = if (!missing(adjustment)) substitute(adjustment),
    data = if (!missing(data)) substitute(data),
    dots = dot_exprs,
    envir = parent,
    dm_size = attr(dm, "Size"), dm = dm
  )
  group <- args$group
  strata <- args$strata
  adjustment <- args$adjustment
  data <- args$data
  a.dist.args <- args$a.dist.args
  strata_label <- args$strata_label

  .validate_distance_diagnostic_options(distance.diagnostics, keep.eigenvalues)
  has_adjustment <- !is.null(adjustment)

  # Goodness-of-fit always compares a residual distance matrix back to the
  # original dm, even when the test statistic itself uses an adjusted dm.
  original_dm <- dm

  if (has_adjustment) {
    dm <- .call_a_dist(
      dm = dm,
      formula = adjustment,
      data = data,
      a.dist.args = a.dist.args,
      distance.diagnostics = distance.diagnostics,
      keep.eigenvalues = keep.eigenvalues,
      quiet = FALSE
    )
  }
  test.results <- generic.distance.permutation.test(WdS, dm=dm, f=group, nrep=nrep, strata=strata)

  # Store full labels; shortening and wrapping belong only to the printer.
  method <- paste0(
    if (has_adjustment && !is.null(strata)) "Adjusted and stratified " else
      if (has_adjustment) "Adjusted " else if (!is.null(strata)) "Stratified " else "",
    if (has_adjustment || !is.null(strata)) "distance-based" else "Distance-based",
    " multivariate Welch ANOVA"
  )
  adjustment_formula <- if (has_adjustment) paste(deparse(adjustment), collapse = "") else NULL
  data.info <- list(
    distance = if (is.null(attr(original_dm, "method"))) "Distance matrix" else attr(original_dm, "method"),
    observations = attr(original_dm, "Size"),
    group = args$group_label,
    levels = nlevels(group)
  )
  data.name <- .wdstest_data_description(data.info)

  # Statistic value
  statistic <- test.results$statistic
  attr(statistic, "names") <- "WdS"

  # Effect size estimator of variance, omega-squared
  estimate <- (((length(levels(group))-1)*(unname(statistic)-1))/((length(levels(group))-1)*(unname(statistic) - 1) + attr(dm, "Size")))
  attr(estimate, "names") <- "effect size, omega-squared (\u03C9\u00B2)"

  # Build the requested pseudo-R2 output after the WdS statistic is computed.
  # For adjusted tests, dm is now the adjustment-residual distance matrix.
  goodness <- .compute_wdstest_goodness(
    dm = original_dm,
    group_label = args$group_label,
    f = group,
    formula = adjustment,
    formula_data = data,
    dm_adjustment_residual = if (has_adjustment) dm else NULL,
    goodness = goodness,
    distance.diagnostics = distance.diagnostics,
    keep.eigenvalues = keep.eigenvalues,
    a.dist.args = a.dist.args
  )

  # P-value
  p.value <- test.results$p.value

  # Creating object of class 'htest'
  TEST <- list(method = method, data.name = data.name, data.info = data.info, # null.value = null.value, # alternative = alternative,
               statistic = statistic, p.value = p.value,
               estimate = estimate,
               goodness.of.fit = goodness$goodness.of.fit,
               goodness.components = goodness$goodness.components,
               distance.diagnostics = goodness$distance.diagnostics,
               parameter = c(if(!is.null(strata)){c("strata"= strata_label)},
                             if(!is.null(adjustment)){c("adjustment"= paste0(adjustment_formula,
                                                                             " with ",
                                                                             length(attr(terms(adjustment), "term.labels")),
                                                                             " terms"))},
                             c("dfb" = (length(levels(group))-1)),
                             c("number of permutations" = nrep)
                             )
               )

  class(TEST) <- c("wdstest", "htest")
  return(TEST)
}

#' @export
print.wdstest <- function(x, digits = 4, ...) {
  # Own the layout instead of relying on print.htest()'s combined summary line.
  # Formatting operates on labels and copies, never on the stored result.
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 1 || digits > 22 || digits != floor(digits)) {
    stop("'digits' must be an integer between 1 and 22.", call. = FALSE)
  }
  width <- max(20L, getOption("width", 80L))
  write_line <- function(text) cat(strwrap(text, width = width, exdent = 2L), sep = "\n")
  number <- function(value) format(unname(value), digits = digits, trim = TRUE)

  cat("\n")
  # Match the traditional htest heading indentation.
  cat(strwrap(x$method, width = width, prefix = "\t"), sep = "\n")
  cat("\n")
  description <- if (is.null(x$data.info)) x$data.name else
    .wdstest_data_description(x$data.info, truncate = TRUE)
  write_line(paste0("Data: ", description))

  if (length(x$parameter)) {
    # Show the formula without the historical term-count suffix; retain the
    # underlying parameter vector for programmatic access and htest tools.
    labels <- names(x$parameter)
    values <- as.character(x$parameter)
    values[labels == "adjustment"] <- sub(" with [0-9]+ terms$", "", values[labels == "adjustment"])
    labels[labels == "number of permutations"] <- "permutations"
    order <- order(match(labels, c("adjustment", "strata", "dfb", "permutations")))
    write_line(paste0("Parameters: ", paste(paste0(labels[order], " = ", values[order]), collapse = ", ")))
  }
  if ("strata" %in% names(x$parameter)) {
    note <- "Stratification affects only the p-value (not WdS, \u03C9\u00B2, or R\u00B2)."
    write_line(.format_r2_name(.format_estimate_name(note)))
  }
  if (length(x$statistic)) write_line(paste0("WdS = ", number(x$statistic)))
  if (length(x$p.value)) {
    value <- format.pval(x$p.value, digits = digits)
    write_line(paste0("P-value ", if (startsWith(value, "<")) value else paste0("= ", value)))
  }
  if (length(x$estimate)) {
    label <- .format_estimate_name("Effect size estimate (\u03C9\u00B2)")
    write_line(paste0(label, " = ", number(x$estimate)))
  }
  if (length(x$goodness.of.fit)) {
    labels <- .wdstest_goodness_labels(x)
    cat(.format_r2_name("Goodness-of-fit (pseudo-R\u00B2):"), "\n", sep = "")
    for (i in seq_along(labels)) {
      label <- paste0(toupper(substr(labels[i], 1L, 1L)), substring(labels[i], 2L))
      cat(paste0("  ", strwrap(paste0(label, " = ", number(x$goodness.of.fit[i])),
                               width = width - 2L, exdent = 2L)), sep = "\n")
    }
  }
  invisible(x)
}

.wdstest_data_description <- function(info, truncate = FALSE) {
  label <- info$group
  if (truncate && nchar(label, type = "chars") > 30L) {
    label <- paste0(substr(label, 1L, 27L), "...")
  }
  distance <- info$distance
  # phyloseq/vegan label Bray-Curtis distances as 'bray'.
  if (tolower(distance) == "bray") distance <- "Bray-Curtis"
  distance <- if (identical(distance, "Distance matrix")) distance else paste(distance, "distances")
  paste0(distance, ", ", info$observations, " observations, ", label,
         " (", info$levels, " levels)")
}

.wdstest_goodness_labels <- function(x) {
  components <- c("group.only", "adjustment.only", "full", "semi.partial", "partial")
  adjusted <- "adjustment" %in% names(x$parameter)
  # Match stored names to component IDs, preserving explicit selection order.
  ids <- components[match(names(x$goodness.of.fit),
                          .goodness_component_label(components, adjusted))]
  labels <- c("group-only", "adjustment-only", "full-model", "semi-partial", "partial")[match(ids, components)]
  if (adjusted) {
    labels[which(ids == "semi.partial")] <-
      "semi-partial (group | adjustment) / total)"
    labels[which(ids == "partial")] <-
      "partial (group | adjustment) / residual)"
  }
  # Unknown labels from externally constructed objects remain printable.
  labels[is.na(labels)] <- names(x$goodness.of.fit)[is.na(labels)]
  .format_r2_name(labels)
}

.resolve_wdstest_arguments <- function(group, strata, adjustment, data, dots,
                                       envir, dm_size, dm = NULL) {
  if (is.null(dots)) {
    dots <- list()
  }

  # Old named arguments arrive through ... because the new parameters now use
  # group/adjustment/data. We keep them working, but we do not let users supply
  # both old and new names for the same job.
  has_legacy_group <- .has_named_dot_arg(dots, "f")
  has_legacy_adjustment <- .has_named_dot_arg(dots, "formula")
  has_legacy_data <- .has_named_dot_arg(dots, "formula_data")

  if (!is.null(group) && has_legacy_group) {
    stop("Use only one of 'group' or deprecated 'f'.", call. = FALSE)
  }
  if (!is.null(adjustment) && has_legacy_adjustment) {
    stop("Use only one of 'adjustment' or deprecated 'formula'.", call. = FALSE)
  }
  if (!is.null(data) && has_legacy_data) {
    stop("Use only one of 'data' or deprecated 'formula_data'.", call. = FALSE)
  }

  if (is.null(group)) {
    if (!has_legacy_group) {
      stop("'group' must be provided.", call. = FALSE)
    }
    # Named f = ... is the old WdS.test() API. Positional second arguments still
    # work naturally because they now match group.
    .warn_deprecated_argument("f", "group", "WdS.test")
    group <- .get_named_dot_expr(dots, "f")
  }

  if (is.null(adjustment) && has_legacy_adjustment) {
    legacy_adjustment <- eval(.get_named_dot_expr(dots, "formula"), envir = envir)
    if (!is.null(legacy_adjustment)) {
      .warn_deprecated_argument("formula", "adjustment", "WdS.test")
      adjustment <- legacy_adjustment
    }
  } else if (!is.null(adjustment)) {
    adjustment <- eval(adjustment, envir = envir)
  }
  adjustment <- .as_wdstest_adjustment(adjustment)

  metadata <- .resolve_wdstest_data(
    data = data,
    formula_data = if (has_legacy_data) .get_named_dot_expr(dots, "formula_data"),
    envir = envir, dm = dm
  )

  # Resolve group and strata after metadata is known. That gives data columns
  # first chance, then falls back to ordinary objects in the caller's frame.
  group_resolved <- .resolve_wdstest_variable(
    expr = group,
    data = metadata$data,
    envir = envir,
    arg = "group",
    allow_string_column = metadata$explicit
  )
  group_value <- .as_wdstest_group(group_resolved$value, arg = "group")

  if (is.null(strata)) {
    strata_value <- NULL
    strata_label <- NULL
  } else {
    strata_resolved <- .resolve_wdstest_variable(
      expr = strata,
      data = metadata$data,
      envir = envir,
      arg = "strata",
      allow_string_column = metadata$explicit
    )
    if (is.null(strata_resolved$value)) {
      # Internal callers often forward strata = strata where the local strata
      # object is NULL. Treat that the same as not supplying strata at all.
      strata_value <- NULL
      strata_label <- NULL
    } else {
      strata_value <- .as_wdstest_vector(strata_resolved$value, arg = "strata")
      .validate_wdstest_vector_length(strata_value, dm_size = dm_size, arg = "strata")
      strata_label <- strata_resolved$label
    }
  }

  .validate_wdstest_vector_length(group_value, dm_size = dm_size, arg = "group")

  list(
    group = group_value,
    group_label = group_resolved$label,
    strata = strata_value,
    strata_label = strata_label,
    adjustment = adjustment,
    data = metadata$data,
    # These are true a.dist() options such as tol. Legacy aliases have already
    # been consumed above, so they cannot leak into a.dist() by mistake.
    a.dist.args = .validate_wdstest_distance_args(.eval_dot_exprs(
      .drop_named_dot_args(dots, c("f", "formula", "formula_data")),
      envir = envir
    ))
  )
}

.resolve_wdstest_data <- function(data, formula_data, envir, dm = NULL) {
  has_data <- !is.null(data)
  has_formula_data <- !is.null(formula_data)

  # data/formula_data can themselves be expressions, for example
  # data = phyloseq::sample_data(ent), so evaluate only after conflicts are
  # checked at the expression level in .resolve_wdstest_arguments().
  data_value <- if (has_data) eval(data, envir = envir) else NULL
  formula_data_value <- if (has_formula_data) eval(formula_data, envir = envir) else NULL
  has_data <- !is.null(data_value)
  has_formula_data <- !is.null(formula_data_value)

  if (has_data && has_formula_data) {
    stop("Use only one of 'data' or deprecated 'formula_data'.", call. = FALSE)
  }
  if (has_formula_data) {
    .warn_deprecated_argument("formula_data", "data", "WdS.test")
    data_value <- formula_data_value
  }

  # Without a user-supplied data object, formulas and unquoted group/strata
  # expressions retain the existing parent.frame() lookup behavior.
  if (!has_data && !has_formula_data) {
    return(list(data = envir, explicit = FALSE))
  }

  # Extract and align whole phyloseq inputs before resolving group or strata
  # columns; subsequent model fits then receive the same ordered metadata.
  list(data = .as_formula_data(data_value, arg = "data", dm = dm), explicit = TRUE)
}

.resolve_wdstest_variable <- function(expr, data, envir, arg,
                                      allow_string_column = FALSE) {
  expr_label <- paste(deparse(expr), collapse = "")

  # Data columns intentionally take precedence, matching formula/data APIs such
  # as lm() and adonis2() while still allowing external vectors as a fallback.
  if (allow_string_column && is.symbol(expr)) {
    column <- as.character(expr)
    if (.data_has_name(data, column)) {
      return(list(value = .data_get_name(data, column), label = column))
    }
  }

  value <- tryCatch(
    eval(expr, envir = envir),
    error = function(e) e
  )
  if (inherits(value, "error")) {
    stop(
      "Could not find '", expr_label, "' for '", arg,
      "'. Supply a vector directly, or supply a column name with 'data'.",
      call. = FALSE
    )
  }
  # This makes programmatic use convenient:
  #   group_var <- "Enterotype"
  #   WdS.test(dm, group = group_var, data = meta)
  if (allow_string_column && is.character(value) && length(value) == 1 &&
      .data_has_name(data, value)) {
    return(list(value = .data_get_name(data, value), label = value))
  }

  list(value = value, label = expr_label)
}

.as_wdstest_group <- function(x, arg) {
  x <- .as_wdstest_vector(x, arg = arg)
  group <- factor(x)
  if (nlevels(group) < 2) {
    stop("'group' must contain at least two levels.", call. = FALSE)
  }
  group
}

.as_wdstest_vector <- function(x, arg) {
  if (is.data.frame(x) || is.matrix(x) || (is.list(x) && !is.factor(x))) {
    stop("'", arg, "' must resolve to a vector or factor, not a data frame, matrix, or list.",
         call. = FALSE)
  }
  if (length(x) < 1) {
    stop("'", arg, "' must not be empty.", call. = FALSE)
  }
  x
}

.validate_wdstest_vector_length <- function(x, dm_size, arg) {
  if (!is.null(dm_size) && length(x) != dm_size) {
    stop(
      "'", arg, "' must have length ", dm_size,
      " to match the distance matrix, but has length ", length(x), ".",
      call. = FALSE
    )
  }
  invisible(x)
}

.as_wdstest_adjustment <- function(adjustment) {
  if (is.null(adjustment)) {
    return(NULL)
  }
  if (!inherits(adjustment, "formula") || length(adjustment) != 2L) {
    stop("'adjustment' must be a right-hand-side formula, such as ~ Age + SeqTech.",
         call. = FALSE)
  }
  adjustment
}

.has_named_dot_arg <- function(dots, arg) {
  dot_names <- names(dots)
  if (is.null(dot_names)) {
    return(FALSE)
  }
  any(dot_names == arg)
}

.get_named_dot_expr <- function(dots, arg) {
  index <- which(names(dots) == arg)
  if (length(index) > 1) {
    stop("Argument '", arg, "' was supplied more than once.", call. = FALSE)
  }
  dots[[index]]
}

.drop_named_dot_args <- function(dots, args) {
  dot_names <- names(dots)
  if (is.null(dot_names)) {
    return(dots)
  }
  dots[!(dot_names %in% args)]
}

.eval_dot_exprs <- function(dots, envir) {
  values <- lapply(dots, eval, envir = envir)
  names(values) <- names(dots)
  values
}

.validate_wdstest_distance_args <- function(args) {
  # Validate routing even if no residual matrix will be computed. Otherwise a
  # misspelled option can appear to work only when goodness is disabled.
  arg_names <- names(args)
  unknown <- setdiff(arg_names, c("", "tol"))
  if (length(unknown)) {
    stop("Unused argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if (length(args) > 1L) {
    stop("Supply 'tol' only once; no other extra distance options are supported.",
         call. = FALSE)
  }
  # The historical seventh positional WdS.test() argument is tol.
  if (length(args)) names(args) <- "tol"
  args
}

.data_has_name <- function(data, name) {
  if (is.environment(data)) {
    exists(name, envir = data, inherits = FALSE)
  } else {
    name %in% names(data)
  }
}

.data_get_name <- function(data, name) {
  if (is.environment(data)) {
    get(name, envir = data, inherits = FALSE)
  } else {
    data[[name]]
  }
}

.call_a_dist <- function(dm, formula, data, a.dist.args,
                         distance.diagnostics = TRUE,
                         keep.eigenvalues = FALSE,
                         quiet = FALSE) {
  # Use the public adjustment/data names so internal calls do not trigger
  # deprecation warnings. Keep extra options such as tol in the same order.
  args <- c(
    list(dm = dm, adjustment = formula, data = data),
    a.dist.args,
    list(
      distance.diagnostics = distance.diagnostics,
      keep.eigenvalues = keep.eigenvalues
    )
  )

  if (isTRUE(quiet)) {
    return(suppressMessages(do.call(a.dist, args)))
  }
  do.call(a.dist, args)
}

.compute_wdstest_goodness <- function(dm, f, formula = NULL, formula_data = NULL,
                                      dm_adjustment_residual = NULL,
                                      goodness = "auto",
                                      distance.diagnostics = TRUE,
                                      keep.eigenvalues = FALSE,
                                      a.dist.args = list(), group_label = "group") {
  has_formula <- !is.null(formula)
  group_text <- .display_group_label(group_label, formula_data)
  group_model <- paste("~", group_text)
  adjustment_model <- if (has_formula) .display_adjustment_model(formula) else NULL
  conditional_model <- .display_conditional_model(group_text, formula, formula_data)
  goodness <- .resolve_goodness(goodness, has_formula = has_formula)
  diagnostics <- list()

  if (has_formula && isTRUE(distance.diagnostics)) {
    diagnostics[["adjustment.only"]] <- .collect_distance_diagnostics(
      dm_adjustment_residual,
      component = "adjustment.only",
      model = adjustment_model
    )
  }
  if (!has_formula && isTRUE(distance.diagnostics) &&
      !is.null(attr(dm, "distance.diagnostics", exact = TRUE))) {
    diagnostics[["input"]] <- .collect_distance_diagnostics(
      dm,
      component = "input"
    )
  }

  if (identical(goodness$requested, "none")) {
    return(list(
      goodness.of.fit = NULL,
      goodness.components = NULL,
      distance.diagnostics = .combine_distance_diagnostics(diagnostics)
    ))
  }

  # Add the tested group as an internal formula variable so users do not need to
  # add it to data themselves. Name generation avoids collisions
  # with existing columns or environment bindings.
  factor_name <- .make_internal_factor_name(formula_data)
  factor_data <- .add_internal_factor(formula_data, f, factor_name)
  factor_formula <- stats::as.formula(paste("~", factor_name))
  components <- list()

  if (has_formula) {
    # The adjustment residual matrix is already produced for the adjusted test.
    r2_adjustment <- unname(dist.goodness.of.fit(dm, dm_adjustment_residual))
    components[["adjustment.only"]] <- .goodness_component(
      component = "adjustment.only",
      model = adjustment_model,
      pseudo.R2 = r2_adjustment
    )

    # group.only requires its own F-only residual distance matrix. It is
    # skipped for adjusted goodness="auto" to avoid unnecessary computation.
    needs_factor_only <- any(goodness$computed %in% "group.only")
    if (needs_factor_only) {
      dm_factor_residual <- .a_dist_for_goodness(
        dm,
        factor_formula,
        factor_data,
        a.dist.args = a.dist.args,
        distance.diagnostics = distance.diagnostics,
        keep.eigenvalues = keep.eigenvalues
      )
      diagnostics[["group.only"]] <- .collect_distance_diagnostics(
        dm_factor_residual,
        component = "group.only",
        model = group_model
      )
      r2_factor <- unname(dist.goodness.of.fit(dm, dm_factor_residual))
      components[["group.only"]] <- .goodness_component(
        component = "group.only",
        model = group_model,
        pseudo.R2 = r2_factor
      )
    }

    needs_full <- any(goodness$computed %in% c("full", "semi.partial", "partial"))
    if (needs_full) {
      # Fit Z + F from the original distances, avoiding sequential residualization.
      full_formula <- .add_factor_to_formula(formula, factor_name, formula_data)
      dm_full_residual <- .a_dist_for_goodness(
        dm,
        full_formula,
        factor_data,
        a.dist.args = a.dist.args,
        distance.diagnostics = distance.diagnostics,
        keep.eigenvalues = keep.eigenvalues
      )
      diagnostics[["full"]] <- .collect_distance_diagnostics(
        dm_full_residual,
        component = "full",
        model = .display_formula_with_factor(formula, formula_data, factor_label = group_text)
      )
      r2_full <- unname(dist.goodness.of.fit(dm, dm_full_residual))
      r2_semi_partial <- r2_full - r2_adjustment
      adjustment_residual_fraction <- 1 - r2_adjustment

      # partial R2 is undefined if the adjustment model leaves no residual
      # variation to be explained by the tested group.
      r2_partial <- if (is.na(r2_adjustment) ||
                        abs(adjustment_residual_fraction) < .Machine$double.eps) {
        NA_real_
      } else {
        r2_semi_partial / adjustment_residual_fraction
      }

      components[["full"]] <- .goodness_component(
        component = "full",
        model = .display_formula_with_factor(formula, formula_data, factor_label = group_text),
        pseudo.R2 = r2_full
      )
      components[["semi.partial"]] <- .goodness_component(
        component = "semi.partial",
        model = paste(conditional_model, "/ total"),
        pseudo.R2 = r2_semi_partial
      )
      components[["partial"]] <- .goodness_component(
        component = "partial",
        model = paste(conditional_model, "/ residual"),
        pseudo.R2 = r2_partial
      )
    }
  } else {
    if (any(goodness$computed %in% "adjustment.only")) {
      stop("'adjustment.only' goodness-of-fit requires an 'adjustment' formula.")
    }

    dm_factor_residual <- .a_dist_for_goodness(
      dm,
      factor_formula,
      factor_data,
      a.dist.args = a.dist.args,
      distance.diagnostics = distance.diagnostics,
      keep.eigenvalues = keep.eigenvalues
    )
    diagnostics[["group.only"]] <- .collect_distance_diagnostics(
      dm_factor_residual,
      component = "group.only",
      model = group_model
    )
    r2_factor <- unname(dist.goodness.of.fit(dm, dm_factor_residual))

    # With no adjustment set, full, semi-partial, and partial reduce to the
    # tested group's model pseudo-R-squared.
    for (component in goodness$computed) {
      components[[component]] <- .goodness_component(
        component = component,
        model = switch(component,
                       semi.partial = paste(conditional_model, "/ total"),
                       partial = paste(conditional_model, "/ residual"),
                       group_model),
        pseudo.R2 = r2_factor
      )
    }
  }

  # Use one presentation order for stored components, independent of the
  # order in which residual matrices were needed for the calculations.
  component_order <- c("group.only", "adjustment.only", "full", "semi.partial", "partial")
  components <- components[component_order[component_order %in% names(components)]]
  goodness.components <- do.call(rbind, components)
  goodness.components$displayed <- goodness.components$component %in% goodness$displayed
  rownames(goodness.components) <- NULL

  displayed_index <- match(goodness$displayed, goodness.components$component)
  goodness.of.fit <- goodness.components$pseudo.R2[displayed_index]
  names(goodness.of.fit) <- .goodness_component_label(goodness$displayed, has_formula)

  list(
    goodness.of.fit = goodness.of.fit,
    goodness.components = goodness.components,
    distance.diagnostics = .combine_distance_diagnostics(diagnostics)
  )
}

.resolve_goodness <- function(goodness, has_formula) {
  allowed <- c("auto", "group.only", "adjustment.only", "full",
               "semi.partial", "partial", "all", "none")
  goodness <- unname(as.character(goodness))
  if (length(goodness) < 1) {
    stop("'goodness' must contain at least one value.")
  }
  if (anyNA(goodness)) {
    stop("'goodness' cannot contain missing values.")
  }
  unknown <- setdiff(goodness, allowed)
  if (length(unknown) > 0) {
    stop("Unknown 'goodness' value: ", paste(unknown, collapse = ", "), ".")
  }
  if (length(intersect(goodness, c("auto", "all", "none"))) > 0 && length(goodness) > 1) {
    stop("'auto', 'all', and 'none' cannot be combined with other 'goodness' values.")
  }

  if (identical(goodness, "none")) {
    return(list(requested = "none", displayed = "none", computed = "none"))
  }
  if (identical(goodness, "auto")) {
    if (has_formula) {
      # Adjusted auto shows the group's semi-partial pseudo-R2. The stored
      # components include partial because it is derived from the same scalars.
      return(list(
        requested = "auto",
        displayed = "semi.partial",
        computed = c("adjustment.only", "full", "semi.partial", "partial")
      ))
    }
    # With no adjustment variables, the natural default is the group-only
    # pseudo-R2 rather than a semi-partial label.
    return(list(
      requested = "auto",
      displayed = "group.only",
      computed = "group.only"
    ))
  }
  if (identical(goodness, "all")) {
    if (has_formula) {
      all_components <- c("group.only", "adjustment.only", "full", "semi.partial", "partial")
    } else {
      all_components <- c("group.only", "full", "semi.partial", "partial")
    }
    return(list(requested = "all", displayed = all_components, computed = all_components))
  }

  computed <- unique(goodness)
  if (has_formula) {
    # semi.partial and partial both require the adjustment-only and full-model
    # pseudo-R2 values, so store those intermediate components too.
    if (any(computed %in% c("semi.partial", "partial"))) {
      computed <- unique(c("adjustment.only", "full", "semi.partial", "partial", computed))
    }
    if (any(computed %in% "full")) {
      computed <- unique(c("adjustment.only", computed))
    }
    if (any(computed %in% "group.only")) {
      computed <- unique(c("adjustment.only", computed))
    }
  }

  list(requested = goodness, displayed = goodness, computed = computed)
}

.a_dist_for_goodness <- function(dm, formula, formula_data, a.dist.args = list(),
                                 distance.diagnostics = TRUE,
                                 keep.eigenvalues = FALSE) {
  # Keep WdS.test output focused; direct a.dist() calls still report eigenvalue
  # handling messages to users.
  .call_a_dist(
    dm = dm,
    formula = formula,
    data = formula_data,
    a.dist.args = a.dist.args,
    distance.diagnostics = distance.diagnostics,
    keep.eigenvalues = keep.eigenvalues,
    quiet = TRUE
  )
}

.collect_distance_diagnostics <- function(dm, component, model = NULL) {
  diagnostics <- attr(dm, "distance.diagnostics", exact = TRUE)
  if (is.null(diagnostics)) return(NULL)
  diagnostics$component <- component
  if (!is.null(model)) diagnostics$model <- model
  diagnostics
}

# Combine a list of distance diagnostics data frames into one, removing NULLs and
#  aligning the optional eigenvalues column. The resulting data frame is sorted
#  by component in the order input, group.only, adjustment.only, full.
# If using a.dist(..., keep.eigenvalues = TRUE), the eigenvalues column will
#  contain lists of numeric vectors. The eigenvalues column will be NULL if
#  keep.eigenvalues = FALSE. This results in a list column of NULLs, which is
#  still a valid data frame column.

.combine_distance_diagnostics <- function(diagnostics) {
  diagnostics <- Filter(Negate(is.null), diagnostics)
  if (length(diagnostics) == 0) return(NULL)
  # Inherited and newly computed diagnostics can differ in whether raw
  # eigenvalues were saved. Align the optional list column before binding.
  if (any(                   # Check whether at least one table has that column
    vapply(diagnostics,      # Check each table for an eigenvalues column
           function(x) "eigenvalues" %in% names(x), #
           logical(1)))) {   # Specifies that each check must return a single logical value
    # If the eigenvalues column is missing from any table, add it as a list column of NULLs.
    diagnostics <- lapply(diagnostics, function(x) {
      if (!"eigenvalues" %in% names(x)) {
        x$eigenvalues <-
          I(                 # Keeps those entries together as a list col inside data.frame
            rep(             # Creates one entry per diagnostic row
              list(NULL),    # Creates one empty list entry
                nrow(x)))    # x returns the updated table to lapply()
      }
      x
    })
  }
  combined <- do.call(rbind, diagnostics)
  # Inherited input diagnostics precede newly constructed model components.
  component_order <- c("input", "group.only", "adjustment.only", "full")
  combined <- combined[order(match(combined$component, component_order)), , drop = FALSE]
  rownames(combined) <- NULL
  combined
}

.make_internal_factor_name <- function(data, base_name = ".wdstar_f") {
  if (is.null(data)) {
    return(base_name)
  }

  name_exists <- function(candidate) {
    if (is.environment(data)) {
      exists(candidate, envir = data, inherits = FALSE)
    } else {
      candidate %in% names(data)
    }
  }

  candidate <- base_name
  counter <- 1L
  while (name_exists(candidate)) {
    candidate <- paste0(base_name, counter)
    counter <- counter + 1L
  }
  candidate
}

.add_internal_factor <- function(data, f, factor_name) {
  if (is.null(data)) {
    data <- data.frame(f, check.names = FALSE)
    names(data) <- factor_name
    return(data)
  }
  if (is.environment(data)) {
    factor_env <- new.env(parent = data)
    assign(factor_name, f, envir = factor_env)
    return(factor_env)
  }

  data[[factor_name]] <- f
  data
}

.add_factor_to_formula <- function(formula, factor_name, data) {
  formula_terms <- stats::terms(formula, data = data)
  term_labels <- attr(formula_terms, "term.labels")
  rhs <- paste(c(term_labels, factor_name), collapse = " + ")

  # Preserve explicit no-intercept formulas. Without this, ~ x - 1 would become
  # ~ x + group when building the full model, silently reintroducing an intercept
  # and changing full, semi.partial, and partial pseudo-R2 values.
  if (identical(attr(formula_terms, "intercept"), 0L)) {
    rhs <- paste(rhs, "- 1")
  }
  stats::as.formula(paste("~", rhs), env = environment(formula))
}

.display_formula_with_factor <- function(formula, data, factor_label = "group") {
  formula_terms <- stats::terms(formula, data = data)
  term_labels <- attr(formula_terms, "term.labels")
  rhs <- paste(c(term_labels, factor_label), collapse = " + ")

  # Match .add_factor_to_formula() so goodness.components displays the same
  # intercept convention used in the internal full-model calculation.
  if (identical(attr(formula_terms, "intercept"), 0L)) {
    rhs <- paste(rhs, "- 1")
  }
  paste("~", rhs)
}

.display_group_label <- function(label, data = NULL) {
  # Metadata columns may have nonsyntactic names; expressions are retained.
  expr <- if (.data_has_name(data, label)) as.name(label) else
    tryCatch({
      parsed <- parse(text = label)
      if (length(parsed) != 1L) stop("Expected one expression")
      parsed[[1L]]
    }, error = function(e) as.name(label))
  text <- paste(deparse(expr, backtick = TRUE, width.cutoff = 500L), collapse = " ")
  if (is.call(expr) && as.character(expr[[1L]])[1L] %in% c("+", "-", "*", "/", "|", "~")) {
    text <- paste0("(", text, ")")
  }
  text
}

.display_adjustment_model <- function(formula) {
  paste("~", paste(deparse(formula[[2L]], width.cutoff = 500L), collapse = " "))
}

.display_conditional_model <- function(group_text, formula, data) {
  adjustment_text <- if (is.null(formula)) "none" else
    paste(deparse(formula[[2L]], width.cutoff = 500L), collapse = " ")
  # Parenthesize compound RHS expressions, preserving interactions, transforms,
  # and no-intercept conventions rather than reducing them to variable names.
  if (!is.null(formula) && is.call(formula[[2L]]) &&
      as.character(formula[[2L]][[1L]])[1L] %in% c("+", "-", "*", "/", "|")) {
    adjustment_text <- paste0("(", adjustment_text, ")")
  }
  paste0("(", group_text, " | ", adjustment_text, ")")
}

.goodness_component <- function(component, model, pseudo.R2) {
  data.frame(component = component, model = model, pseudo.R2 = pseudo.R2,
             stringsAsFactors = FALSE)
}

.goodness_component_label <- function(component, has_formula) {
  vapply(component, function(x) {
    switch(
      x,
      "adjustment.only" = "pseudo-R\u00B2 for adjustment variables only",
      "group.only" = "pseudo-R\u00B2 for group only",
      "full" = "pseudo-R\u00B2 for full model",
      "semi.partial" = if (has_formula) {
        "semi-partial pseudo-R\u00B2 for group given adjustment"
      } else {
        "semi-partial pseudo-R\u00B2 for group"
      },
      "partial" = if (has_formula) {
        "partial pseudo-R\u00B2 for group given adjustment"
      } else {
        "partial pseudo-R\u00B2 for group"
      },
      x
    )
  }, character(1))
}

.format_estimate_name <- function(name,
                                  utf8 = isTRUE(l10n_info()[["UTF-8"]]),
                                  html = .is_html_output(),
                                  latex = .is_latex_output()) {
  if (.use_ascii_math_labels(utf8 = utf8, html = html, latex = latex)) {
    return(sub("\u03C9\u00B2", "omega^2", name, fixed = TRUE))
  }
  name
}

.format_r2_name <- function(name,
                            utf8 = isTRUE(l10n_info()[["UTF-8"]]),
                            html = .is_html_output(),
                            latex = .is_latex_output()) {
  if (.use_ascii_math_labels(utf8 = utf8, html = html, latex = latex)) {
    return(sub("R\u00B2", "R^2", name, fixed = TRUE))
  }
  name
}

.use_ascii_math_labels <- function(
    utf8 = isTRUE(l10n_info()[["UTF-8"]]),
    html = .is_html_output(),
    latex = .is_latex_output()) {
  isTRUE(latex) || (!isTRUE(utf8) && isTRUE(html))
}

.is_html_output <- function() {
  if (requireNamespace("knitr", quietly = TRUE)) {
    isTRUE(knitr::is_html_output())
  } else {
    FALSE
  }
}

.is_latex_output <- function() {
  if (requireNamespace("knitr", quietly = TRUE)) {
    isTRUE(knitr::is_latex_output())
  } else {
    FALSE
  }
}
