test_that("WdS accepts prepared factors and agrees with WdS.test", {
  set.seed(92)
  dm <- dist(matrix(rnorm(48), nrow = 12))
  metadata <- data.frame(treatment = factor(rep(c("A", "B"), each = 6)))
  group_vector <- metadata$treatment
  expected <- WdS(dm, group_vector)
  expect_equal(WdS(dm, group = group_vector), expected)
  expect_warning(legacy <- WdS(dm, f = group_vector), "'f' is deprecated in WdS()",
                 fixed = TRUE)
  expect_equal(legacy, expected)
  # The permutation engine still calls the statistic positionally.
  result <- WdS.test(dm, treatment, data = metadata, nrep = 3, goodness = "none")
  expect_equal(unname(result$statistic), unname(expected))
})

test_that("WdS rejects ambiguous and invalid grouping arguments", {
  dm <- dist(matrix(seq_len(24), nrow = 12))
  group <- factor(rep(1:2, each = 6))
  expect_error(WdS(dm), "'group' must be provided", fixed = TRUE)
  expect_error(WdS(dm, group, f = group), "Use only one")
  expect_error(WdS(dm, f = group, f = group), "more than once")
  expect_error(WdS(dm, group, typo = TRUE), "Only")
  expect_error(WdS(dm, group, data = data.frame(group)), "Only")
  expect_error(WdS(dm, as.character(group)), "must be a factor")
  expect_error(WdS(dm, group[-1]), "must have length 12")
  expect_error(WdS(dm, rep("A", 12)), "at least two levels")
})

test_that("univariate distance statistics agree with standard Welch tests", {
  set.seed(86)
  y <- c(rnorm(6), rnorm(7, 1, 2), rnorm(8, 2, 3))
  group <- factor(rep(1:3, c(6, 7, 8)))
  expect_equal(unname(WdS(dist(y), group)),
               unname(stats::oneway.test(y ~ group)$statistic))
  keep <- group != "3"
  two_groups <- droplevels(group[keep])
  expect_equal(unname(Tw2(dist(y[keep]), two_groups)),
               unname(stats::t.test(y[keep] ~ two_groups)$statistic)^2)
})

test_that("post hoc comparisons retain correct group sizes and observed statistics", {
  set.seed(87)
  y <- rnorm(18)
  group <- factor(rep(c("A", "B", "C"), c(5, 6, 7)))
  dm <- dist(y)
  pairs <- Tw2.posthoc.tests(dm, group, nrep = 3)
  expect_equal(dim(pairs), c(3L, 7L))
  for (i in seq_len(nrow(pairs))) {
    keep <- group %in% pairs[i, c("Level1", "Level2")]
    expect_equal(as.numeric(pairs[i, c("N1", "N2")]),
                 as.numeric(table(droplevels(group[keep]))))
    expect_equal(as.numeric(pairs[i, "tw2.stat"]),
                 unname(Tw2(as.dist(as.matrix(dm)[keep, keep]), droplevels(group[keep]))))
  }
  one_vs_all <- Tw2.posthoc.1vsAll.tests(dm, group, nrep = 3)
  expect_equal(dim(one_vs_all), c(3L, 5L))
  expect_equal(as.numeric(one_vs_all[, "N2"]), c(5, 6, 7))
})
