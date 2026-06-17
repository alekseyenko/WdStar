# WdStar 2.4.0 2026-06-17

* Added `dist.goodness.of.fit()` to compute goodness-of-fit as the coefficient
  of determination, \(R^2 = 1 - \sigma^2_{residual} /
  \sigma^2_{total}\), for users who compute raw and adjusted distance matrices
  outside of `WdS.test()`. The helper validates that both inputs are `dist`
  objects with the same number of observations.

* `WdS.test()` now reports a separate `goodness.of.fit` component for adjusted
  tests. This keeps `estimate` reserved for omega-squared (\(\omega^2\)) while
  reporting adjustment goodness-of-fit separately. For unadjusted tests,
  `goodness.of.fit` is `NULL`.

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
  parent-frame formula lookup, adjusted-test output, and direct Goodness of Fit
  calculation from raw and adjusted distance matrices.

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
