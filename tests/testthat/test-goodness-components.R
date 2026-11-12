test_that("component names, default order, and explicit selections are consistent", {
  set.seed(117)
  dm <- dist(matrix(rnorm(60), 15))
  meta <- data.frame(group = factor(rep(1:3, each = 5)), z = rnorm(15))
  fit <- function(goodness, adjustment = ~ z) suppressMessages(
    WdS.test(dm, group, data = meta, adjustment = adjustment,
             goodness = goodness, nrep = 0))
  expected <- c("group.only", "adjustment.only", "full", "semi.partial", "partial")
  all <- fit("all")
  expect_identical(all$goodness.components$component, expected)
  expect_identical(names(all$goodness.of.fit),
                   unname(WdStar:::.goodness_component_label(expected, TRUE)))
  selected <- fit(c("partial", "group.only"))
  expect_identical(names(selected$goodness.of.fit),
    unname(WdStar:::.goodness_component_label(c("partial", "group.only"), TRUE)))
  expect_identical(selected$goodness.components$component, expected)
  values <- setNames(all$goodness.components$pseudo.R2, expected)
  expect_equal(unname(values["semi.partial"]), unname(values["full"] - values["adjustment.only"]))
  expect_equal(unname(values["partial"]), unname(values["semi.partial"] / (1 - values["adjustment.only"])))
  expect_equal(unname(selected$goodness.of.fit), unname(values[c("partial", "group.only")]))
  expect_error(fit("adjustment"), "Unknown 'goodness'")
  expect_error(fit("adjustment.only", NULL), "requires an 'adjustment' formula")
  unadjusted <- fit("all", NULL)
  expect_identical(unadjusted$goodness.components$component, expected[-2])
  expect_equal(unname(unadjusted$goodness.of.fit), rep(unname(unadjusted$goodness.of.fit[1]), 4))
  expect_equal(fit("auto")$goodness.of.fit, all$goodness.of.fit[4])
  expect_identical(fit("auto")$goodness.components$component, expected[-1])
})
