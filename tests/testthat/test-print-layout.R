test_that("printed summaries follow the agreed layout in all model configurations", {
  set.seed(222)
  meta <- data.frame(Enterotype = factor(rep(1:3, each = 5)),
                     Age = rnorm(15), Subject = factor(rep(1:5, 3)))
  dm <- dist(matrix(rnorm(45), 15))
  attr(dm, "method") <- "bray"
  for (adjusted in c(FALSE, TRUE)) for (stratified in c(FALSE, TRUE)) {
    adjustment <- if (adjusted) ~ Age else NULL
    strata <- if (stratified) meta$Subject else NULL
    set.seed(223)
    res <- suppressMessages(WdS.test(dm, Enterotype, adjustment = adjustment,
                      strata = strata, data = meta, nrep = 3))
    before <- res
    output <- capture.output(returned <- print(res))
    positions <- vapply(c("Data:", "Parameters:", "WdS =", "P-value ",
                         "Effect size estimate", "Goodness-of-fit (pseudo-R\u00B2):"),
      function(prefix) which(startsWith(output, prefix))[1], integer(1))
    expect_false(anyNA(positions))
    expect_true(all(diff(positions) > 0))
    expect_equal(unname(positions[2] - positions[1]), 1L)
    expect_true(any(grepl("Enterotype (3 levels)", output, fixed = TRUE)))
    expect_true(any(grepl("Bray-Curtis distances", output, fixed = TRUE)))
    if (stratified) expect_true(any(grepl(
      "Stratification affects only the p-value", output, fixed = TRUE)))
    expect_identical(res, before)
    expect_identical(returned, before)
    expect_identical(res$data.info$group, "Enterotype")
    expect_true(all(names(res) == tolower(names(res))))
    # Same observed quantities as the unstratified model; display changes must
    # not affect the statistic, effect size, or distance-based components.
    plain <- suppressMessages(WdS.test(dm, meta$Enterotype, adjustment = adjustment,
                                      data = meta, nrep = 0))
    expect_identical(res$statistic, plain$statistic)
    expect_identical(res$estimate, plain$estimate)
    expect_equal(res$goodness.components$pseudo.R2, plain$goodness.components$pseudo.R2)
    expect_error(print(res, digits = 0), "digits")
  }
})

test_that("group labels are resolved and truncated only at print time", {
  set.seed(224)
  dm <- dist(matrix(rnorm(30), 10))
  long <- "Enterotype_label_longer_than_thirty_characters"
  meta <- setNames(data.frame(factor(rep(1:2, each = 5))), long)
  res <- WdS.test(dm, group = long, data = meta, nrep = 0)
  expect_identical(res$data.info$group, long)
  expect_true(grepl(long, res$data.name, fixed = TRUE))
  output <- capture.output(print(res))
  expect_true(any(grepl(paste0(substr(long, 1, 27), "..."), output, fixed = TRUE)))
  expect_false(any(grepl(long, output, fixed = TRUE)))
})

test_that("multiple goodness components preserve default and requested ordering", {
  set.seed(225)
  dm <- dist(matrix(rnorm(36), 12))
  meta <- data.frame(g = factor(rep(1:3, 4)), z = rnorm(12))
  for (selection in list("all", c("partial", "group.only"))) {
    res <- suppressMessages(WdS.test(dm, g, adjustment = ~ z, data = meta,
                                    nrep = 0, goodness = selection))
    output <- capture.output(print(res))
    labels <- WdStar:::.wdstest_goodness_labels(res)
    positions <- vapply(labels, function(label)
      which(startsWith(tolower(trimws(output)), tolower(substr(label, 1, 15))))[1], integer(1))
    expect_true(all(diff(positions) > 0))
    expect_true(any(output == "Goodness-of-fit (pseudo-R\u00B2):"))
    expect_false(any(grepl("Pseudo-R", output, fixed = TRUE)))
  }
})

test_that("parameter wrapping respects output width", {
  set.seed(226)
  dm <- dist(matrix(rnorm(36), 12))
  meta <- data.frame(g = factor(rep(1:3, 4)), z = rnorm(12))
  res <- suppressMessages(WdS.test(dm, g, adjustment = ~ z, data = meta, nrep = 0))
  old <- options(width = 60)
  on.exit(options(old))
  output <- capture.output(print(res))
  expect_true(all(nchar(output, type = "width") <= 60))
})
