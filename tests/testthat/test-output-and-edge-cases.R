test_that("printing keeps effect sizes on one line without changing results", {
  set.seed(83)
  dm <- dist(matrix(rnorm(48), 12))
  group <- factor(rep(1:2, each = 6))
  for (goodness in list("none", "auto", "all")) {
    x <- WdS.test(dm, group, nrep = 3, goodness = goodness)
    original <- x
    output <- capture.output(returned <- print(x))
    estimate_line <- output[grepl("Effect size estimate", output, fixed = TRUE)]
    expect_length(estimate_line, 1)
    expect_match(estimate_line, " = [-0-9.NA]")
    expect_false(any(grepl("=  =", output, fixed = TRUE)))
    expect_false(any(grepl("sample estimates:", output, fixed = TRUE)))
    expect_identical(returned, original)
    expect_identical(x, original)
    for (label in WdStar:::.wdstest_goodness_labels(x)) {
      expect_true(any(grepl(tolower(label), tolower(paste(trimws(output), collapse = " ")), fixed = TRUE)))
    }
  }
})

test_that("method labels reflect adjustment and stratification independently", {
  set.seed(84)
  dm <- dist(matrix(rnorm(48), 12))
  meta <- data.frame(group = factor(rep(1:2, 6)), covariate = rnorm(12),
                     block = factor(rep(1:3, each = 4)))
  for (adjusted in c(FALSE, TRUE)) {
    adjustment <- if (adjusted) ~ covariate else NULL
    plain <- suppressMessages(WdS.test(dm, group, adjustment = adjustment,
                                      data = meta, nrep = 3, goodness = "all"))
    blocked <- suppressMessages(WdS.test(dm, group, adjustment = adjustment,
                        strata = block, data = meta, nrep = 3, goodness = "all"))
    expect_match(blocked$method, if (adjusted) "Adjusted and stratified" else "Stratified distance")
    expect_identical(blocked$statistic, plain$statistic)
    expect_identical(blocked$estimate, plain$estimate)
    expect_equal(blocked$goodness.components, plain$goodness.components)
  }
})

test_that("restricted permutations preserve singleton numeric strata", {
  f <- c(5, 7, 9, 11)
  strata <- c("single", "pair", "pair", "other")
  statistic <- function(dm, f) {
    expect_equal(f[c(1, 4)], c(5, 11))
    expect_equal(sort(f[2:3]), c(7, 9))
    sum(f * seq_along(f))
  }
  result <- generic.distance.permutation.test(statistic, NULL, f, nrep = 9, strata = strata)
  expect_true(result$p.value > 0 && result$p.value <= 1)
  for (nrep in list(-1, 1.5, NA_real_, Inf, "9")) {
    expect_error(generic.distance.permutation.test(sum, NULL, f, nrep), "nrep")
  }
  expect_error(generic.distance.permutation.test(sum, NULL, f, strata = c(1, NA, 2, 2)), "strata")
})

test_that("adjustment rejects missing or mismatched metadata and supports no intercept", {
  set.seed(85)
  dm <- dist(matrix(rnorm(24), 8))
  meta <- data.frame(z = rnorm(8))
  expect_s3_class(suppressMessages(a.dist(dm, ~ z - 1, meta)), "dist")
  meta$z[1] <- NA
  expect_error(a.dist(dm, ~ z, meta), "missing")
  expect_error(a.dist(dm, ~ z, data.frame(z = 1:7)), "one row")
  expect_error(a.dist(dm, ~ z, meta, tol = Inf), "tol")
})

test_that("tapply permutations restore interleaved sample order and factor levels", {
  f <- factor(c("A", "B", "C", "A", "B", "C"),
              levels = c("C", "B", "A", "unused"))
  strata <- factor(c("second", "first", "single", "second", "first", "last"),
                   levels = c("first", "last", "second", "single", "unused"))
  statistic <- function(dm, permuted) {
    expect_identical(levels(permuted), levels(f))
    expect_identical(unname(as.character(permuted[c(3, 6)])), c("C", "C"))
    for (block in unique(strata)) {
      idx <- which(strata == block)
      expect_equal(sort(as.character(permuted[idx])), sort(as.character(f[idx])))
    }
    sum(as.integer(permuted) * seq_along(permuted))
  }
  set.seed(88)
  result <- generic.distance.permutation.test(statistic, NULL, f, nrep = 3, strata = strata)
  expect_true(is.finite(result$p.value))
})
