test_that("inherited eigenvalues coexist with new diagnostic summaries", {
  set.seed(410)
  dm <- dist(matrix(rnorm(48), 12))
  meta <- data.frame(g = factor(rep(1:3, each = 4)), z = rnorm(12))
  for (input_keeps in c(FALSE, TRUE)) for (output_keeps in c(FALSE, TRUE)) {
    adjusted <- suppressMessages(a.dist(dm, ~ z, meta, keep.eigenvalues = input_keeps))
    fit <- WdS.test(adjusted, g, data = meta, nrep = 0,
                    goodness = "all", keep.eigenvalues = output_keeps)
    d <- fit$distance.diagnostics
    expect_identical(d$component, c("input", "group.only"))
    expect_identical(d$model[1], "~ z")
    expect_identical(d$keep.eigenvalues, c(input_keeps, output_keeps))
    if (input_keeps || output_keeps) {
      expect_identical(lengths(d$eigenvalues), ifelse(d$keep.eigenvalues, 12L, 0L))
      if (input_keeps) expect_identical(d$eigenvalues[[1]],
        attr(adjusted, "distance.diagnostics")$eigenvalues[[1]])
    } else {
      expect_false("eigenvalues" %in% names(d))
    }
    # Optional diagnostic storage must not change numerical test results.
    baseline <- WdS.test(adjusted, g, data = meta, nrep = 0,
                         goodness = "all", distance.diagnostics = FALSE)
    expect_identical(fit$statistic, baseline$statistic)
    expect_identical(fit$estimate, baseline$estimate)
    expect_equal(fit$goodness.components, baseline$goodness.components)
  }
})
