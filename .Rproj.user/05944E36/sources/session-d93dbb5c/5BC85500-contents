#' Bayesian Grouped Horseshoe (BGHS) for Logistic Regression
#'
#' Fits a Bayesian grouped logistic regression model using the BGHS prior.
#'
#' @param formula A model formula (e.g., \code{y ~ x1 + x2}).
#' @param data A data frame containing the variables.
#' @param groups A list of variable names or a vector of group IDs.
#' @param n_iter Total number of MCMC iterations.
#' @param burnin Number of burn-in iterations.
#' @param thin Thinning interval (default = 1).
#' @param init Starting values (optional).
#' @param verbose Print progress messages (default = TRUE).
#' @param print_every Print progress every N iterations (default = 100).
#' @param ... Additional arguments passed to \code{bghs_gibbs}.
#'
#' @return An object of class \code{bghs}
#'
#' @importFrom stats model.frame model.response model.matrix coef
#' @importFrom stats plogis quantile sd
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- bghs(y ~ ., data = my_data, groups = my_groups)
#' summary(fit)
#' }
bghs <- function(formula, data, groups, n_iter = 20000, burnin = 10000,
                 thin = 1, init = NULL, verbose = TRUE, print_every = 100, ...) {

  # 1. Extract response and design matrix
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, data)[, -1, drop = FALSE]

  n <- nrow(X)
  p <- ncol(X)

  # 2. Process groups
  if (is.list(groups)) {
    G <- length(groups)
    group_sizes <- sapply(groups, length)
    var_names <- colnames(X)
    all_names <- unlist(groups)
    if (!all(all_names %in% var_names)) {
      stop("Some variable names in groups not found in data")
    }
    X <- X[, all_names, drop = FALSE]
  } else if (is.numeric(groups) && length(groups) == p) {
    G <- length(unique(groups))
    group_sizes <- as.vector(table(groups))
    ord <- order(groups)
    X <- X[, ord, drop = FALSE]
  } else {
    stop("groups must be a list of variable names or a vector of group IDs")
  }

  # 3. Prepare data for Gibbs
  data_gibbs <- cbind(y, X)
  colnames(data_gibbs)[1] <- "y"

  # 4. Create initial values if not provided
  if (is.null(init)) {
    init <- create_bghs_init(n, p, G, group_sizes)
  }

  # 5. Run Gibbs sampler
  if (verbose) {
    cat("\nRunning BGHS Gibbs sampler...\n")
  }

  start_time <- Sys.time()

  posterior <- bghs_gibbs(
    data = data_gibbs,
    nsamples = n_iter,
    burnin = burnin,
    init = init,
    group_sizes = group_sizes,
    print_every = print_every,
    ...
  )

  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # 6. Extract beta columns
  beta_cols <- grep("^beta_", colnames(posterior))
  beta_names <- colnames(X)
  if (length(beta_cols) == length(beta_names)) {
    colnames(posterior)[beta_cols] <- beta_names
  }

  # 7. Compute summaries
  beta_summary <- t(apply(posterior[, beta_cols, drop = FALSE], 2,
                          function(x) c(mean = mean(x), sd = sd(x),
                                        quantile(x, c(0.025, 0.5, 0.975)))))

  colnames(beta_summary) <- c("Mean", "SD", "2.5%", "50%", "97.5%")

  # 8. Build output
  out <- list(
    posterior = posterior,
    beta = posterior[, beta_cols, drop = FALSE],
    coefficients = beta_summary,
    formula = formula,
    data = data,
    groups = groups,
    group_sizes = group_sizes,
    n_iter = n_iter,
    burnin = burnin,
    thin = thin,
    method = "BGHS",
    elapsed_time = elapsed_time
  )

  class(out) <- "bghs"
  return(out)
}

#' Create Initial Values for BGHS
#'
#' Internal function to create initial values for the BGHS sampler.
#'
#' @noRd
create_bghs_init <- function(n, p, G, group_sizes) {
  c(rep(0.5, n),          # w
    rep(0, p),            # beta
    rep(1, p),            # lambda_sq
    rep(1, p),            # t
    rep(1, G),            # delta_sq
    rep(1, G),            # c
    1,                    # tau_sq
    1)                    # v
}
