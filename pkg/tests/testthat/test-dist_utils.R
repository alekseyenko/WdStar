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
