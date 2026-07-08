test_that("GRASP runs without errors", {
  set.seed(123)
  X <- matrix(rnorm(30 * 3), 30, 3)
  colnames(X) <- c("x1", "x2", "x3")
  beta <- c(1, 0.5, 0)
  y <- rbinom(30, 1, plogis(X %*% beta))
  data <- data.frame(y, X)

  groups <- list(g1 = c("x1", "x2"), g2 = c("x3"))

  expect_no_error({
    fit <- grasp(y ~ ., data = data, groups = groups,
                 n_iter = 100, burnin = 50, verbose = FALSE)
  })
})
