# WdStar 2.4.0


* Added `dist.goodness.of.fit()` to compute Goodness of Fit as the coefficient
  of prediction, \(R^2 = 1 - \sigma^2_{residual} / \sigma^2_{total}\), for
  users who compute raw and adjusted distance matrices outside of `WdS.test()`.
  The helper validates that both inputs are `dist` objects with the same number
  of observations.

* `WdS.test()` now reports a separate `goodness.of.fit` component for adjusted
  tests. This keeps `estimate` reserved for omega squared (\(\omega^2\)) while
  reporting adjustment Goodness of Fit separately.

* Improved eigenvalue handling in `a.dist()`: eigenvalues smaller than `tol` and
  negative eigenvalues are now reported to the user and set to zero before the
  adjusted distance matrix is reconstructed.

* Added regression tests for adjusted-test Goodness of Fit, direct
  `dist.goodness.of.fit()` use, `a.dist()` tolerance messages, and negative
  eigenvalue handling.

* Added covariate-adjusted testing support for `WdS.test()` via `formula` and
  `formula_data`, allowing users to adjust the input distance matrix through
  `a.dist()` before computing the WdS test statistic.

* Improved `formula_data` handling for `a.dist()` and adjusted `WdS.test()`
  calls. Formula variables can be resolved from the caller environment by
  default, and `formula_data` now accepts environments, data frames, lists, and
  objects coercible to data frames, including `phyloseq::sample_data()` objects.

* Added `phyloseq` as a suggested package to support tests and examples involving
  `sample_data()` inputs.

* Updated `a.dist()` documentation to describe accepted `formula_data` input
  types, parent-frame lookup, and covariate-adjusted distance matrix behavior.

* Updated `WdS.test()` documentation and examples for adjusted tests, stratified
  tests, and the relationship between `WdS.test()` and `a.dist()`.

* Added test coverage for distance utilities, parent-frame formula lookup,
  adjusted `WdS.test()` calls, and `phyloseq::sample_data()` compatibility.

* Refreshed package documentation generated from roxygen comments, including
  exports for distance utility functions.

* Added GitHub issue templates for bug reports and feature requests.
