set.seed(44)
n_rows <- 100

example_df <- data.frame(
  X1 = rnorm(n_rows, mean = 0, sd = 1),
  X2 = rnorm(n_rows, mean = 0, sd = 1),
  X3 = rnorm(n_rows, mean = 0, sd = 1),
  X4 = rnorm(n_rows, mean = 0, sd = 1)
)

# Generate factors
factor_var <- factor(rep(c("A", "B"), each = n_rows / 2))

# Make distance matrix (type: dist)
distance_matrix <- dist(example_df)

# Make square matrix form (type: matrix, array)
distance_matrix_sq <- as.matrix(distance_matrix)

ss2 <- dist.ss2(distance_matrix_sq, factor_var)
cat("Sum of Squares Matrix:\n")
print(ss2)

test_that("'distance' objects are successfully identified", {
  expect_true(is.dist(distance_matrix))
  expect_false(is.dist(distance_matrix_sq))
})

test_that("correct sigma squared value is calculated", {
  expect_equal(dist.sigma2(distance_matrix), 4.2206406)
})

test_that("dist.ss2 correctly calcualtes sum of squares matrix", {
  expect_equal(dist.ss2(distance_matrix_sq, factor_var)[1], 3323.9046)
})

test_that("dist.group.sigma2 correctly calcualtes grou-wise sigma squared stats", {
  expect_equal(dist.group.sigma2(distance_matrix, factor_var)[1], setNames(4.1662956, "A"))
  expect_equal(dist.group.sigma2(distance_matrix, factor_var)[2], setNames(4.2468826, "B"))
})

test_that("dist.cohen.d correctly calculates Cohen's d for the distance matrix", {
  expect_equal(dist.cohen.d(distance_matrix, factor_var), setNames(0.23071167, "A"),
               tolerance = 0.001)
})

test_that("WdS.test validates extra options even without goodness calculations", {
  expect_error(WdS.test(distance_matrix, factor_var, nrep = 1,
                       goodness = "none", typo = 1), "Unused argument")
  expect_error(WdS.test(distance_matrix, factor_var, nrep = 1,
                       goodness = "none", tol = 1e-8, tol = 1e-9),
               "Supply 'tol' only once", fixed = TRUE)
  expect_error(WdS.test(distance_matrix, factor_var, 1, NULL, NULL, NULL,
                       1e-8, tol = 1e-9, goodness = "none"),
               "Supply 'tol' only once", fixed = TRUE)
})

test_that("a.dist uses formula variables from the parent frame by default", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)

  adjusted_dm <- suppressMessages(a.dist(distance_matrix, adjustment = ~ adjustment))

  expect_s3_class(adjusted_dm, "dist")
  expect_equal(attr(adjusted_dm, "Size"), n_rows)
})

test_that("a.dist accepts data and warns for deprecated formula_data", {
  metadata <- data.frame(adjustment = rep(c("low", "high"), length.out = n_rows))

  adjusted_dm <- suppressMessages(a.dist(
    dm = distance_matrix,
    adjustment = ~ adjustment,
    data = metadata
  ))
  expect_warning(
    legacy_adjusted_dm <- suppressMessages(a.dist(
      dm = distance_matrix,
      adjustment = ~ adjustment,
      formula_data = metadata
    )),
    "'formula_data' is deprecated"
  )

  expect_s3_class(adjusted_dm, "dist")
  expect_s3_class(legacy_adjusted_dm, "dist")
  expect_equal(as.vector(adjusted_dm), as.vector(legacy_adjusted_dm))
})

test_that("a.dist resolves legacy formula and caller metadata consistently", {
  covariate <- seq_len(n_rows)
  metadata <- data.frame(covariate = covariate)
  expected <- suppressMessages(a.dist(distance_matrix, ~ covariate, metadata))
  expect_warning(
    legacy <- suppressMessages(a.dist(distance_matrix, formula = ~ covariate,
                                      data = metadata)),
    "'formula' is deprecated", fixed = TRUE
  )
  expect_equal(legacy, expected)
  # NULL and omitted data must both find covariate in this calling frame.
  expect_equal(suppressMessages(a.dist(distance_matrix, ~ covariate)), expected)
  expect_equal(suppressMessages(a.dist(distance_matrix, ~ covariate, data = NULL)), expected)
  expect_equal(suppressMessages(a.dist(distance_matrix, ~ covariate,
                                      data = as.list(metadata))), expected)
  expect_equal(suppressMessages(a.dist(distance_matrix, ~ covariate,
                                      data = list2env(metadata))), expected)
  # Variables absent from supplied metadata still use the formula environment.
  expect_equal(suppressMessages(a.dist(distance_matrix, ~ covariate,
                                      data = data.frame(other = covariate))), expected)
})

test_that("a.dist rejects ambiguous arguments and invalid adjustment syntax", {
  covariate <- seq_len(n_rows)
  expect_error(a.dist(distance_matrix), "'adjustment' must be provided", fixed = TRUE)
  expect_error(a.dist(distance_matrix, NULL), "'adjustment' must be provided", fixed = TRUE)
  expect_error(a.dist(distance_matrix, "covariate"), "right-hand-side formula")
  expect_error(a.dist(distance_matrix, covariate ~ covariate), "right-hand-side formula")
  expect_error(a.dist(distance_matrix, ~ covariate, formula = ~ covariate),
               "Use only one of 'adjustment'", fixed = TRUE)
  expect_error(a.dist(distance_matrix, ~ covariate, data = NULL, formula_data = NULL),
               "Use only one of 'data'", fixed = TRUE)
  expect_error(a.dist(distance_matrix, ~ covariate, typo = 1), "Unused argument")
  expect_error(a.dist(distance_matrix, formula = ~ covariate, formula = ~ covariate),
               "more than once")
  expect_error(a.dist(distance_matrix, ~ covariate, NULL, 1e-8, TRUE, FALSE, list()),
               "must be named")
})

test_that("a.dist accepts phyloseq sample metadata with the preferred names", {
  skip_if_not_installed("phyloseq")
  metadata <- data.frame(covariate = seq_len(n_rows))
  expected <- suppressMessages(a.dist(distance_matrix, ~ covariate, metadata))
  actual <- suppressMessages(a.dist(distance_matrix, adjustment = ~ covariate,
                                    data = phyloseq::sample_data(metadata)))
  expect_equal(actual, expected)
})

test_that("a.dist and WdS.test keep old positional tol behavior", {
  metadata <- data.frame(adjustment = rep(c("low", "high"), length.out = n_rows))

  adjusted_dm <- suppressMessages(a.dist(
    distance_matrix,
    ~ adjustment,
    metadata,
    0.1
  ))
  test_result <- suppressMessages(WdS.test(
    distance_matrix,
    factor_var,
    1,
    NULL,
    ~ adjustment,
    metadata,
    0.1,
    goodness = "adjustment.only"
  ))

  expect_equal(attr(adjusted_dm, "distance.diagnostics")$tol, 0.1)
  expect_equal(test_result$distance.diagnostics$tol, 0.1)
})

test_that("WdS.test uses formula variables from the parent frame by default", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)

  result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9,
    adjustment = ~ adjustment
  ))

  expect_s3_class(result, "htest")
  expect_match(result$method, "Adjusted")
})

test_that("WdS.test accepts phyloseq sample_data as data", {
  skip_if_not_installed("phyloseq")

  sample_df <- data.frame(
    adjustment = rep(c("low", "high"), length.out = n_rows),
    row.names = attr(distance_matrix, "Labels")
  )
  sample_df <- phyloseq::sample_data(sample_df)

  result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9,
    adjustment = ~ adjustment,
    data = sample_df
  ))

  expect_s3_class(result, "htest")
  expect_match(result$method, "Adjusted")
})

test_that("WdS.test resolves group, strata, and adjustment variables from data", {
  metadata <- data.frame(
    group_col = factor_var,
    strata_col = rep(c("site1", "site2"), length.out = n_rows),
    adjustment = rep(c("low", "high"), length.out = n_rows),
    row.names = attr(distance_matrix, "Labels")
  )

  vector_result <- WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 1,
    goodness = "none"
  )
  unquoted_result <- WdS.test(
    dm = distance_matrix,
    group = group_col,
    data = metadata,
    nrep = 1,
    goodness = "none"
  )
  string_result <- WdS.test(
    dm = distance_matrix,
    group = "group_col",
    data = metadata,
    nrep = 1,
    goodness = "none"
  )
  stratified_result <- WdS.test(
    dm = distance_matrix,
    group = group_col,
    strata = "strata_col",
    data = metadata,
    nrep = 1,
    goodness = "none"
  )
  adjusted_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    group = group_col,
    adjustment = ~ adjustment,
    data = metadata,
    nrep = 1,
    goodness = "adjustment.only"
  ))

  expect_equal(unname(unquoted_result$statistic), unname(vector_result$statistic))
  expect_equal(unname(string_result$statistic), unname(vector_result$statistic))
  expect_equal(unname(stratified_result$parameter["strata"]), "strata_col")
  expect_match(adjusted_result$method, "Adjusted")
  expect_named(adjusted_result$goodness.of.fit,
               "pseudo-R\u00B2 for adjustment variables only")
})

test_that("WdS.test keeps legacy argument names working with deprecation warnings", {
  metadata <- data.frame(group_col = factor_var)
  adjustment <- rep(c("low", "high"), length.out = n_rows)

  expect_warning(
    legacy_group <- WdS.test(
      dm = distance_matrix,
      f = factor_var,
      nrep = 1,
      goodness = "none"
    ),
    "'f' is deprecated"
  )
  expect_warning(
    legacy_formula <- suppressMessages(WdS.test(
      dm = distance_matrix,
      group = factor_var,
      formula = ~ adjustment,
      nrep = 1,
      goodness = "adjustment.only"
    )),
    "'formula' is deprecated"
  )
  expect_warning(
    legacy_data <- WdS.test(
      dm = distance_matrix,
      group = "group_col",
      formula_data = metadata,
      nrep = 1,
      goodness = "none"
    ),
    "'formula_data' is deprecated"
  )

  expect_s3_class(legacy_group, "htest")
  expect_s3_class(legacy_formula, "htest")
  expect_s3_class(legacy_data, "htest")
})

test_that("WdS.test rejects conflicting old and new argument names", {
  metadata <- data.frame(group_col = factor_var)
  adjustment <- rep(c("low", "high"), length.out = n_rows)

  expect_error(
    WdS.test(
      dm = distance_matrix,
      group = factor_var,
      f = factor_var,
      nrep = 1
    ),
    "Use only one of 'group' or deprecated 'f'"
  )
  expect_error(
    WdS.test(
      dm = distance_matrix,
      group = factor_var,
      adjustment = ~ adjustment,
      formula = ~ adjustment,
      nrep = 1
    ),
    "Use only one of 'adjustment' or deprecated 'formula'"
  )
  expect_error(
    WdS.test(
      dm = distance_matrix,
      group = "group_col",
      data = metadata,
      formula_data = metadata,
      nrep = 1
    ),
    "Use only one of 'data' or deprecated 'formula_data'"
  )
})

test_that("WdS.test reports auto goodness of fit for adjusted and unadjusted tests", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(distance_matrix, adjustment = ~ adjustment))
  factor_data <- data.frame(.wdstar_f = factor_var)
  factor_dm <- suppressMessages(a.dist(distance_matrix, adjustment = ~ .wdstar_f, data = factor_data))
  full_data <- data.frame(adjustment = adjustment, .wdstar_f = factor_var)
  full_dm <- suppressMessages(a.dist(distance_matrix, adjustment = ~ adjustment + .wdstar_f, data = full_data))
  expected_factor_r_squared <- dist.goodness.of.fit(distance_matrix, factor_dm)
  expected_adjustment_r_squared <- dist.goodness.of.fit(distance_matrix, adjusted_dm)
  expected_full_r_squared <- dist.goodness.of.fit(distance_matrix, full_dm)
  expected_semi_partial <- unname(expected_full_r_squared) - unname(expected_adjustment_r_squared)
  expected_partial <- expected_semi_partial / (1 - unname(expected_adjustment_r_squared))

  unadjusted_result <- WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9
  )
  adjusted_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9,
    adjustment = ~ adjustment
  ))

  expect_named(unadjusted_result$estimate, "effect size, omega-squared (\u03C9\u00B2)")
  expect_named(adjusted_result$estimate, "effect size, omega-squared (\u03C9\u00B2)")
  expect_named(unadjusted_result$goodness.of.fit, "pseudo-R\u00B2 for group only")
  expect_named(adjusted_result$goodness.of.fit, "semi-partial pseudo-R\u00B2 for group given adjustment")
  expect_equal(unname(unadjusted_result$goodness.of.fit), unname(expected_factor_r_squared))
  expect_equal(unname(adjusted_result$goodness.of.fit), expected_semi_partial)
  expect_setequal(
    adjusted_result$goodness.components$component,
    c("adjustment.only", "full", "semi.partial", "partial")
  )
  expect_equal(
    adjusted_result$goodness.components$pseudo.R2[
      adjusted_result$goodness.components$component == "adjustment.only"
    ],
    unname(expected_adjustment_r_squared)
  )
  expect_equal(
    adjusted_result$goodness.components$pseudo.R2[
      adjusted_result$goodness.components$component == "full"
    ],
    unname(expected_full_r_squared)
  )
  expect_equal(
    adjusted_result$goodness.components$pseudo.R2[
      adjusted_result$goodness.components$component == "partial"
    ],
    expected_partial
  )
  expect_equal(
    adjusted_result$goodness.components$component[adjusted_result$goodness.components$displayed],
    "semi.partial"
  )
  expect_equal(
    adjusted_result$goodness.components$model[
      adjusted_result$goodness.components$component == "full"
    ],
    "~ adjustment + factor_var"
  )
  expect_equal(adjusted_result$distance.diagnostics$component, c("adjustment.only", "full"))
  expect_false("eigenvalues" %in% names(adjusted_result$distance.diagnostics))
})

test_that("WdS.test can report selected goodness-of-fit components", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)

  adjustment_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9,
    adjustment = ~ adjustment,
    goodness = "adjustment.only"
  ))
  all_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9,
    adjustment = ~ adjustment,
    goodness = "all"
  ))
  none_result <- WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9,
    goodness = "none"
  )

  expect_named(adjustment_result$goodness.of.fit, "pseudo-R\u00B2 for adjustment variables only")
  expect_setequal(
    all_result$goodness.components$component,
    c("group.only", "adjustment.only", "full", "semi.partial", "partial")
  )
  expect_equal(
    all_result$goodness.components$component[all_result$goodness.components$displayed],
    c("group.only", "adjustment.only", "full", "semi.partial", "partial")
  )
  expect_null(none_result$goodness.of.fit)
  expect_null(none_result$goodness.components)
  expect_null(none_result$distance.diagnostics)
  expect_named(
    WdS.test(
      dm = distance_matrix,
      group = factor_var,
      nrep = 9,
      goodness = stats::setNames("auto", "user_supplied_name")
    )$goodness.of.fit,
    "pseudo-R\u00B2 for group only"
  )
  expect_error(
    WdS.test(dm = distance_matrix, group = factor_var, nrep = 9, goodness = "adjustment.only"),
    "requires an 'adjustment' formula"
  )
  expect_error(
    WdS.test(dm = distance_matrix, group = factor_var, nrep = 9, goodness = NA_character_),
    "cannot contain missing values"
  )
})

test_that("WdS.test collects distance diagnostics separately from goodness components", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  pre_adjusted_dm <- suppressMessages(a.dist(distance_matrix, adjustment = ~ adjustment))

  adjustment_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9,
    adjustment = ~ adjustment,
    goodness = "adjustment.only"
  ))
  all_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9,
    adjustment = ~ adjustment,
    goodness = "all",
    keep.eigenvalues = TRUE
  ))
  no_diagnostics_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9,
    adjustment = ~ adjustment,
    distance.diagnostics = FALSE
  ))
  pre_adjusted_result <- WdS.test(
    dm = pre_adjusted_dm,
    group = factor_var,
    nrep = 9,
    goodness = "none"
  )

  expect_equal(adjustment_result$distance.diagnostics$component, "adjustment.only")
  expect_equal(
    all_result$distance.diagnostics$component,
    c("group.only", "adjustment.only", "full")
  )
  expect_true("eigenvalues" %in% names(all_result$distance.diagnostics))
  expect_length(all_result$distance.diagnostics$eigenvalues[[1]], n_rows)
  expect_null(no_diagnostics_result$distance.diagnostics)
  expect_equal(pre_adjusted_result$distance.diagnostics$component, "input")
  expect_equal(pre_adjusted_result$distance.diagnostics$model, "~ adjustment")
})

test_that("WdS.test goodness formula construction preserves no-intercept formulas", {
  formula_data <- data.frame(adjustment = seq_len(n_rows))
  full_formula <- WdStar:::.add_factor_to_formula(~ adjustment - 1, ".wdstar_f", formula_data)

  expect_equal(deparse(full_formula), "~adjustment + .wdstar_f - 1")
  expect_equal(
    WdStar:::.display_formula_with_factor(~ adjustment - 1, formula_data),
    "~ adjustment + group - 1"
  )
})

test_that("WdS.test suppresses internal goodness residualization messages", {
  expect_silent(WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9
  ))
})

test_that("estimate labels fall back safely for non-UTF-8 HTML output", {
  x <- WdS.test(
    dm = distance_matrix,
    group = factor_var,
    nrep = 9
  )

  expect_match(names(x$estimate), "\u03C9\u00B2", fixed = TRUE)
  expect_equal(
    WdStar:::.format_estimate_name(
      names(x$estimate), utf8 = FALSE, html = TRUE, latex = FALSE
    ),
    "effect size, omega-squared (omega^2)"
  )
  expect_equal(
    WdStar:::.format_estimate_name(
      names(x$estimate), utf8 = TRUE, html = TRUE, latex = FALSE
    ),
    names(x$estimate)
  )
  expect_equal(
    WdStar:::.format_estimate_name(
      names(x$estimate), utf8 = TRUE, html = FALSE, latex = TRUE
    ),
    "effect size, omega-squared (omega^2)"
  )
})

test_that("the unreleased factor.only name is not accepted", {
  # Only the documented group.only name is supported.
  for (adjusted in c(FALSE, TRUE)) {
    expect_error(
      WdStar:::.resolve_goodness("factor.only", adjusted),
      "Unknown 'goodness' value", fixed = TRUE
    )
  }
})

test_that("dist.goodness.of.fit computes distance-based pseudo-R squared", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  dm_residual <- suppressMessages(a.dist(distance_matrix, adjustment = ~ adjustment))
  expected_r_squared <- 1 - (dist.sigma2(dm_residual) / dist.sigma2(distance_matrix))

  result <- dist.goodness.of.fit(distance_matrix, dm_residual = dm_residual)

  expect_named(result, "pseudo-R\u00B2")
  expect_equal(unname(result), expected_r_squared)
})

test_that("distance-based pseudo-R squared has expected boundary cases", {
  expect_equal(unname(dist.goodness.of.fit(distance_matrix, distance_matrix)), 0)
  expect_equal(
    unname(dist.goodness.of.fit(distance_matrix, distance_matrix / 2)),
    0.75
  )

  zero_dm <- stats::as.dist(matrix(0, nrow = 4, ncol = 4))
  expect_true(is.na(unname(dist.goodness.of.fit(zero_dm, zero_dm))))
})

test_that("dist.goodness.of.fit validates distance matrix inputs", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  dm_residual <- suppressMessages(a.dist(distance_matrix, adjustment = ~ adjustment))

  expect_error(
    dist.goodness.of.fit(as.matrix(distance_matrix), dm_residual),
    "must both be distance matrices"
  )
  expect_error(
    dist.goodness.of.fit(distance_matrix, as.dist(matrix(0, nrow = 3, ncol = 3))),
    "same number of observations"
  )
  expect_error(
    dist.goodness.of.fit(distance_matrix),
    "'dm_residual' must be provided"
  )
})

test_that("a.dist stores compact distance diagnostics by default", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(distance_matrix, adjustment = ~ adjustment, tol = 0.1))
  diagnostics <- attr(adjusted_dm, "distance.diagnostics")

  expect_s3_class(diagnostics, "data.frame")
  expect_equal(nrow(diagnostics), 1)
  expect_equal(diagnostics$tol, 0.1)
  expect_equal(diagnostics$n.eigenvalues, n_rows)
  expect_gt(diagnostics$n.zeroed.by.tol, 0)
  expect_equal(
    diagnostics$n.negative.removed,
    diagnostics$n.negative.raw - diagnostics$n.negative.zeroed.by.tol
  )
  expect_true(all(c(
    "sum.positive.raw",
    "sum.negative.raw",
    "sum.negative.magnitude.raw",
    "prop.negative.raw"
  ) %in% names(diagnostics)))
  expect_false("eigenvalues" %in% names(diagnostics))
})

test_that("a.dist can keep raw eigenvalues when requested", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(
    distance_matrix,
    adjustment = ~ adjustment,
    keep.eigenvalues = TRUE
  ))
  diagnostics <- attr(adjusted_dm, "distance.diagnostics")

  expect_true(diagnostics$keep.eigenvalues)
  expect_true("eigenvalues" %in% names(diagnostics))
  expect_type(diagnostics$eigenvalues[[1]], "double")
  expect_length(diagnostics$eigenvalues[[1]], n_rows)
})

test_that("a.dist can skip distance diagnostics", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(
    distance_matrix,
    adjustment = ~ adjustment,
    distance.diagnostics = FALSE
  ))

  expect_null(attr(adjusted_dm, "distance.diagnostics", exact = TRUE))
  expect_error(
    a.dist(
      distance_matrix,
      adjustment = ~ adjustment,
      distance.diagnostics = FALSE,
      keep.eigenvalues = TRUE
    ),
    "requires 'distance.diagnostics = TRUE'"
  )
})

test_that("a.dist reports and zeroes eigenvalues below tolerance", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)

  expect_message(
    adjusted_dm <- a.dist(distance_matrix, adjustment = ~ adjustment, tol = 0.1),
    "eigenvalues were smaller than 'tol'"
  )

  expect_s3_class(adjusted_dm, "dist")
  expect_false(any(is.nan(adjusted_dm)))
})

test_that("a.dist reports and zeroes negative eigenvalues", {
  non_euclidean_dm <- as.dist(matrix(
    c(
      0.000, 2.114, 1.444, 1.293,
      2.114, 0.000, 1.133, 1.278,
      1.444, 1.133, 0.000, 2.059,
      1.293, 1.278, 2.059, 0.000
    ),
    nrow = 4
  ))
  attr(non_euclidean_dm, "method") <- "non-euclidean"
  adjustment <- c(1, 2, 1, 2)

  expect_message(
    adjusted_dm <- a.dist(non_euclidean_dm, adjustment = ~ adjustment, tol = 0),
    "eigenvalues were negative and were set to zero"
  )

  expect_s3_class(adjusted_dm, "dist")
  expect_false(any(is.nan(adjusted_dm)))
})

test_that("whole phyloseq data aligns samples for tests and adjusted distances", {
  skip_if_not_installed("phyloseq")
  ids <- paste0("sample", 1:8)
  meta <- data.frame(group = factor(rep(c("A", "B"), each = 4)),
                     covariate = c(1, 3, 2, 5, 4, 8, 6, 7), row.names = ids)
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(matrix(seq_len(24), nrow = 3,
                               dimnames = list(paste0("taxon", 1:3), rev(ids))),
                        taxa_are_rows = TRUE),
    phyloseq::sample_data(meta[rev(ids), , drop = FALSE]))
  dm <- dist(matrix(seq_len(14)^2, nrow = 7, dimnames = list(ids[1:7], NULL)))
  aligned <- meta[ids[1:7], , drop = FALSE]
  set.seed(201)
  expected <- suppressMessages(WdS.test(dm, group = group, adjustment = ~ covariate,
                                      data = aligned, nrep = 19))
  set.seed(201)
  actual <- suppressMessages(WdS.test(dm, group = group, adjustment = ~ covariate,
                                    data = ps, nrep = 19))
  expect_equal(actual$statistic, expected$statistic)
  expect_equal(actual$p.value, expected$p.value)
  expect_equal(actual$goodness.components, expected$goodness.components)
  expect_equal(suppressMessages(a.dist(dm, adjustment = ~ covariate, data = ps)),
               suppressMessages(a.dist(dm, adjustment = ~ covariate, data = aligned)))
  set.seed(201)
  quoted <- suppressMessages(WdS.test(dm, group = "group", data = ps, nrep = 19))
  expect_s3_class(quoted, "wdstest")
  bad <- dm
  attr(bad, "Labels") <- NULL
  expect_error(WdS.test(bad, group = group, data = ps), "sample labels")
  attr(bad, "Labels") <- rep(ids[1], 7)
  expect_error(a.dist(bad, adjustment = ~ covariate, data = ps), "sample labels")
  attr(bad, "Labels") <- c(ids[1:6], "missing")
  expect_error(WdS.test(bad, group = group, data = ps), "missing samples")
  no_meta <- phyloseq::phyloseq(phyloseq::otu_table(ps),
    phyloseq::tax_table(matrix("Bacteria", nrow = 3,
      dimnames = list(paste0("taxon", 1:3), "Kingdom"))))
  expect_error(WdS.test(dm, group = group, data = no_meta), "no sample metadata")
})
