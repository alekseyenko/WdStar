# WdStar 2.5.3

* Standardized goodness-of-fit and distance-diagnostic tables on `component`
  and `model`, using actual grouping labels. Removed the interpretation column
  in favor of consistent documented definitions; model notation identifies
  total and residual denominators. Numerical results and printed summaries are
  unchanged. Diagnostic components now share goodness-of-fit names and order.

* Moved pseudo-R-squared to the goodness-of-fit print heading and shortened
  component labels. Adjusted semi-partial and partial labels now specify their
  total and adjustment-residual reference variation.

* Clarified the pseudo-R-squared and WdS-derived omega-squared definitions and
  equations. Printed results now separate parameters, WdS, p-value, effect size
  estimate, and goodness-of-fit; show the grouping label (up to 30 characters);
  and use sentence-case labels. Help-index titles use plain text.

* Renamed the goodness component `adjustment` to `adjustment.only` without an
  alias. The adjustment formula argument is unchanged. Components now use
  group.only, adjustment.only, full, semi.partial, partial order; explicit
  selections retain their requested display order. Clarified denominators,
  interpretation limits, and selected versus calculated output fields.

# WdStar 2.5.2

* Print omega-squared on one line for consistency.
* Distinguish stratified and adjusted tests in the method heading.
* Fixed a bug that caused issues with singleton strata. Now we validate 
  permutation counts and strata; and we preserve singleton numeric strata
  when generating restricted permutations.
* Reject missing or mismatched adjustment metadata explicitly and support
  single-term adjustment formulas without an intercept.
* Refresh examples for clarity, 
* Updated the automated checks to include checks for new changes, including 
  outputs of omega-squared, sample estimates, and pseudo-R^2.
* Updated package-build and git ignore file exclusions in preparation for 
  CRAN submission. 

# WdStar 2.5.1 2026-09-09

* `WdS.test()` and `a.dist()` accept whole phyloseq objects through `data`,
  extracting sample metadata and aligning it to distance-matrix sample IDs.

* Moved feature and taxa importance functions to the
  `feature-ranking-selection` development branch.

* Updated `WdS()` to use `group` for a directly supplied grouping factor.
  Named `f` remains a deprecated alias; positional calls still work. Metadata
  lookup stays in `WdS.test()` to keep repeated statistic calculations light.
  `WdS()` continues to return only the observed statistic.

* Reject unexpected `WdS.test()` options even when goodness-of-fit is disabled.
  Clarified metadata lookup/order.
  Help examples skip optional phyloseq usage when it is unavailable rather
  than installing dependencies, and explicitly print the main test results.

* Goodness-of-fit labels now refer to the tested group. Use `group.only` for
  the group-only component in `goodness`, returned components, and distance
  diagnostic roles.

* Modernized the preferred `WdS.test()` interface. New code should use
  `group`, `adjustment`, and `data`, for example
  `WdS.test(dm, group = Enterotype, adjustment = ~ SeqTech, data = meta)`.
  Named legacy arguments `f`, `formula`, and `formula_data` still work through
  `...`, but now warn users to use the new names.

* Added tidy-friendly column lookup to `WdS.test()`. When `data` is supplied,
  `group` and `strata` can be vectors, unquoted column names, or length-one
  character column names.

* Unified `a.dist()` with `WdS.test()` using `adjustment` and `data`.
  Deprecated named aliases `formula` and `formula_data` are accepted through
  `...` with warnings; conflicting, duplicated, or unexpected arguments error.
  Omitted or `NULL` data uses the calling environment, and adjustment requires
  a one-sided formula. The first six positional arguments keep their meaning;
  legacy `formula_data` must now be named. Intercept handling, rank handling,
  eigenvalue rules, and the adjustment calculations are unchanged.

* Replaced remaining `mtcars`, random-matrix, and placeholder examples in the
  package documentation with `phyloseq::enterotype` microbiome examples.
  Examples now use Bray-Curtis distances, `Enterotype` as the tested factor,
  `SeqTech` as a categorical adjustment or stratification variable in separate
  examples, and an Age-complete subset when demonstrating continuous covariate
  adjustment.

* Polished help-file wording for low-level statistic and distance helpers,
  including clarifying that `Tw2()` and `WdS()` compute raw statistics while
  `Tw2.test()`, `WdS.test()`, and the generic test wrapper perform permutation
  testing.

* Corrected the documented return value of `dist.group.sigma2()` from a
  diagonal matrix to the named numeric vector returned by the function.

# WdStar 2.5.0 2026-08-26

* Added `dist.goodness.of.fit()` to compute a distance-based pseudo-coefficient
  of determination, \(R^2_{pseudo} = 1 - \sigma^2_{residual} /
  \sigma^2_{total}\), for users who compute original and residual distance
  matrices outside of `WdS.test()`. The helper now uses the more general
  `dm_residual` argument and validates that both inputs are `dist` objects with
  the same number of observations.

* `WdS.test()` now reports a separate `goodness.of.fit` component. This keeps
  `estimate` reserved for omega-squared (\(\omega^2\)) while reporting
  distance-based pseudo-\(R^2\) separately. The new `goodness` argument defaults
  to `"auto"`, which reports the tested factor's pseudo-\(R^2\) for unadjusted
  tests and the semi-partial pseudo-\(R^2\) for adjusted tests.

* Added `goodness.components` to `WdS.test()` results. This table stores the
  goodness-of-fit components computed along the way, including adjustment-only,
  full-model, semi-partial, and partial pseudo-\(R^2\) for adjusted tests.
  Users can request specific components or all available components with
  `goodness = "all"`.

* Added compact distance diagnostics to `a.dist()` and `WdS.test()`.
  `a.dist()` now attaches a `distance.diagnostics` attribute that records the
  eigendecomposition tolerance, raw positive/negative eigenvalue counts,
  eigenvalues zeroed by `tol`, and negative eigenvalues removed during distance
  reconstruction. `WdS.test()` collects these diagnostics in a separate
  `distance.diagnostics` output component for constructed residual distance
  matrices. Raw eigenvalue vectors are stored only when
  `keep.eigenvalues = TRUE`.

* Improved eigenvalue handling in `a.dist()`: eigenvalues smaller than `tol` and
  negative eigenvalues are now reported to the user and set to zero before the
  adjusted distance matrix is reconstructed.

* Improved `formula_data` handling for `a.dist()` and adjusted `WdS.test()`
  calls. Formula variables can be resolved from the caller environment by
  default, and `formula_data` now accepts environments, data frames, lists, and
  objects coercible to data frames, including `phyloseq::sample_data()` objects.

* Added `phyloseq` as a suggested package to support tests and examples involving
  `sample_data()` inputs.

* Updated documentation for `a.dist()`, `WdS.test()`, and
  `dist.goodness.of.fit()`, including accepted `formula_data` input types,
  parent-frame formula lookup, goodness-of-fit component definitions, distance
  diagnostics, and direct pseudo-\(R^2\) calculation from original and residual
  distance matrices.

* Added regression tests for parent-frame formula lookup,
  `phyloseq::sample_data()` compatibility, adjusted-test goodness-of-fit, direct
  `dist.goodness.of.fit()` use, `a.dist()` tolerance messages, and negative
  eigenvalue handling.

* Added GitHub issue templates for bug reports and feature requests.

# WdStar 2.3.0 2025-09-21

* Added package versioning.

* Documented installation from GitHub using `remotes::install_github()`.

* Added covariate-adjusted testing support through `WdS.test()`.

* Added omega-squared as an effect size estimate.

* Added between-group degrees of freedom to `WdS.test()` output.

* Updated `a.dist()` parameters to avoid requesting duplicate objects.

* Improved error handling and data type handling.

* Improved `a.dist()` and `WdS.test()` outputs.

* Revamped help files and updated examples for `a.dist()` and `WdS.test()`.
