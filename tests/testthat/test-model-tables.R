test_that("component tables use actual model labels and the simplified schema", {
  set.seed(301)
  meta <- data.frame(Enterotype = factor(rep(1:3, each = 5)), Age = rnorm(15),
                     SeqTech = factor(rep(1:2, length.out = 15)))
  dm <- dist(matrix(rnorm(45), 15))
  fit <- function(adjustment = ~ SeqTech, goodness = "all") suppressMessages(
    WdS.test(dm, Enterotype, adjustment = adjustment, data = meta,
             nrep = 0, goodness = goodness, keep.eigenvalues = TRUE))
  x <- fit()
  expect_named(x$goodness.components, c("component", "model", "pseudo.R2", "displayed"))
  expect_identical(x$goodness.components$model, c("~ Enterotype", "~ SeqTech",
    "~ SeqTech + Enterotype", "(Enterotype | SeqTech) / total",
    "(Enterotype | SeqTech) / residual"))
  expect_identical(x$distance.diagnostics$component, c("group.only", "adjustment.only", "full"))
  expect_identical(x$distance.diagnostics$model, x$goodness.components$model[1:3])
  expect_false(any(c("role", "formula", "interpretation") %in% names(x$distance.diagnostics)))
  multiple <- fit(~ Age + SeqTech)
  expect_identical(multiple$goodness.components$model[4:5], c(
    "(Enterotype | (Age + SeqTech)) / total", "(Enterotype | (Age + SeqTech)) / residual"))
  none <- fit(goodness = "none")
  expect_null(none$goodness.components)
  expect_identical(none$distance.diagnostics$component, "adjustment.only")
  expect_identical(none$distance.diagnostics$model, "~ SeqTech")
  unadjusted <- fit(NULL)
  expect_identical(unadjusted$goodness.components$model,
    c("~ Enterotype", "~ Enterotype", "(Enterotype | none) / total",
      "(Enterotype | none) / residual"))
  expect_equal(unadjusted$goodness.components$pseudo.R2,
    rep(unadjusted$goodness.components$pseudo.R2[1], 4))
  # Counts distinguish tolerance zeroing from removal of remaining negatives.
  d <- x$distance.diagnostics
  expect_equal(d$n.negative.raw, d$n.negative.zeroed.by.tol + d$n.negative.removed)
  for (i in seq_len(nrow(d))) {
    raw <- d$eigenvalues[[i]]
    expect_equal(d$n.negative.removed[i], sum(raw < 0 & abs(raw) >= d$tol[i]))
  }
})

test_that("direct and inherited diagnostics retain model provenance", {
  set.seed(302)
  meta <- data.frame(Enterotype = factor(rep(1:3, each = 5)), Age = rnorm(15))
  dm <- dist(matrix(rnorm(45), 15))
  residual <- suppressMessages(a.dist(dm, ~ Age, data = meta))
  direct <- attr(residual, "distance.diagnostics")
  expect_identical(direct$component, "adjustment.only")
  expect_identical(direct$model, "~ Age")
  x <- WdS.test(residual, Enterotype, data = meta, goodness = "all", nrep = 0)
  expect_identical(x$distance.diagnostics$component, c("input", "group.only"))
  expect_identical(x$distance.diagnostics$model, c("~ Age", "~ Enterotype"))
  no_goodness <- WdS.test(residual, Enterotype, data = meta, goodness = "none", nrep = 0)
  expect_identical(no_goodness$distance.diagnostics$model, "~ Age")
})

test_that("model labels preserve expressions, nonsyntactic names and intercept conventions", {
  set.seed(303)
  meta <- data.frame(check.names = FALSE, "group name" = factor(rep(1:3, each = 5)), z = rnorm(15))
  dm <- dist(matrix(rnorm(45), 15))
  x <- suppressMessages(WdS.test(dm, "group name", adjustment = ~ z - 1,
                               data = meta, goodness = "all", nrep = 0))
  expect_identical(x$goodness.components$model[1], "~ `group name`")
  expect_identical(x$goodness.components$model[3], "~ z + `group name` - 1")
  expect_identical(x$goodness.components$model[4], "(`group name` | (z - 1)) / total")
  expression <- suppressMessages(WdS.test(dm, factor(meta[["group name"]]),
                                adjustment = ~ I(z^2), data = meta, goodness = "all", nrep = 0))
  expect_true(grepl('factor(meta[["group name"]])', expression$goodness.components$model[1], fixed=TRUE))
  expect_true(grepl("I(z^2)", expression$goodness.components$model[4], fixed=TRUE))
  names(meta)[1] <- "a+b"
  unusual <- WdS.test(dm, "a+b", data = meta, nrep = 0)
  expect_identical(unusual$goodness.components$model, "~ `a+b`")
})
