#' @importFrom stats coef model.matrix plogis quantile sd
#' @importFrom graphics plot
#' @importFrom grDevices dev.off pdf
NULL
#' Print Method for BGHS Objects
#'
#' @param x Object of class \code{bghs}.
#' @param ... Additional arguments.
#' @export
print.bghs <- function(x, ...) {
  cat("\n========================================\n")
  cat("Bayesian Grouped Horseshoe (BGHS) Model\n")
  cat("========================================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Iterations: %d\n", x$n_iter))
  cat(sprintf("Burn-in: %d\n", x$burnin))
  cat(sprintf("Posterior samples: %d\n", nrow(x$posterior)))
  cat(sprintf("Groups: %d\n", length(x$group_sizes)))
  cat(sprintf("Predictors: %d\n", nrow(x$coefficients)))
  cat(sprintf("Elapsed time: %.1f seconds\n", x$elapsed_time))
  cat("========================================\n")
}

#' Print Method for GRASP Objects
#'
#' @param x Object of class \code{grasp}.
#' @param ... Additional arguments.
#' @export
print.grasp <- function(x, ...) {
  cat("\n========================================\n")
  cat("GRASP Model\n")
  cat("========================================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Iterations: %d\n", x$n_iter))
  cat(sprintf("Burn-in: %d\n", x$burnin))
  cat(sprintf("Posterior samples: %d\n", nrow(x$posterior)))
  cat(sprintf("Groups: %d\n", length(x$group_sizes)))
  cat(sprintf("Predictors: %d\n", nrow(x$coefficients)))
  cat(sprintf("Elapsed time: %.1f seconds\n", x$elapsed_time))
  cat("========================================\n")
}

#' Print Method for GIGG Objects
#'
#' @param x Object of class \code{gigg}.
#' @param ... Additional arguments.
#' @export
print.gigg <- function(x, ...) {
  cat("\n========================================\n")
  cat("GIGG Model\n")
  cat("========================================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Iterations: %d\n", x$n_iter))
  cat(sprintf("Burn-in: %d\n", x$burnin))
  cat(sprintf("Posterior samples: %d\n", nrow(x$posterior)))
  cat(sprintf("Groups: %d\n", length(x$group_sizes)))
  cat(sprintf("Predictors: %d\n", nrow(x$coefficients)))
  cat(sprintf("Elapsed time: %.1f seconds\n", x$elapsed_time))
  cat("========================================\n")
}

#' Summary Method for BGHS Objects
#'
#' @param object Object of class \code{bghs}.
#' @param ... Additional arguments.
#' @export
summary.bghs <- function(object, ...) {
  cat("\n========================================\n")
  cat("Bayesian Grouped Horseshoe (BGHS) Summary\n")
  cat("========================================\n")
  cat(sprintf("Method: %s\n", object$method))
  cat(sprintf("Iterations: %d\n", object$n_iter))
  cat(sprintf("Burn-in: %d\n", object$burnin))
  cat(sprintf("Posterior samples: %d\n", nrow(object$posterior)))
  cat(sprintf("Groups: %d\n", length(object$group_sizes)))
  cat(sprintf("Predictors: %d\n", nrow(object$coefficients)))
  cat(sprintf("Elapsed time: %.1f seconds\n", object$elapsed_time))
  cat("\nCoefficients:\n")
  print(object$coefficients)
}

#' Summary Method for GRASP Objects
#'
#' @param object Object of class \code{grasp}.
#' @param ... Additional arguments.
#' @export
summary.grasp <- function(object, ...) {
  cat("\n========================================\n")
  cat("GRASP Summary\n")
  cat("========================================\n")
  cat(sprintf("Method: %s\n", object$method))
  cat(sprintf("Iterations: %d\n", object$n_iter))
  cat(sprintf("Burn-in: %d\n", object$burnin))
  cat(sprintf("Posterior samples: %d\n", nrow(object$posterior)))
  cat(sprintf("Groups: %d\n", length(object$group_sizes)))
  cat(sprintf("Predictors: %d\n", nrow(object$coefficients)))
  cat(sprintf("Elapsed time: %.1f seconds\n", object$elapsed_time))
  cat("\nMH Acceptance Rates:\n")
  print(object$acceptance)
  cat("\nCoefficients:\n")
  print(object$coefficients)
}

#' Summary Method for GIGG Objects
#'
#' @param object Object of class \code{gigg}.
#' @param ... Additional arguments.
#' @export
summary.gigg <- function(object, ...) {
  cat("\n========================================\n")
  cat("GIGG Summary\n")
  cat("========================================\n")
  cat(sprintf("Method: %s\n", object$method))
  cat(sprintf("Iterations: %d\n", object$n_iter))
  cat(sprintf("Burn-in: %d\n", object$burnin))
  cat(sprintf("Posterior samples: %d\n", nrow(object$posterior)))
  cat(sprintf("Groups: %d\n", length(object$group_sizes)))
  cat(sprintf("Predictors: %d\n", nrow(object$coefficients)))
  cat(sprintf("Elapsed time: %.1f seconds\n", object$elapsed_time))
  cat("\nMH Acceptance Rates:\n")
  print(object$acceptance)
  cat("\nCoefficients:\n")
  print(object$coefficients)
}

#' Coef Method for BGHS Objects
#'
#' @param object Object of class \code{bghs}.
#' @param ... Additional arguments.
#' @export
coef.bghs <- function(object, ...) {
  return(object$coefficients[, "Mean"])
}

#' Coef Method for GRASP Objects
#'
#' @param object Object of class \code{grasp}.
#' @param ... Additional arguments.
#' @export
coef.grasp <- function(object, ...) {
  return(object$coefficients[, "Mean"])
}

#' Coef Method for GIGG Objects
#'
#' @param object Object of class \code{gigg}.
#' @param ... Additional arguments.
#' @export
coef.gigg <- function(object, ...) {
  return(object$coefficients[, "Mean"])
}

#' Fitted Method for BGHS Objects
#'
#' @param object Object of class \code{bghs}.
#' @param ... Additional arguments.
#' @export
fitted.bghs <- function(object, ...) {
  X <- model.matrix(object$formula, object$data)[, -1, drop = FALSE]
  beta_mean <- coef(object)
  if (length(beta_mean) != ncol(X)) {
    var_names <- colnames(X)
    match_idx <- match(var_names, names(beta_mean))
    beta_mean <- beta_mean[match_idx]
  }
  eta <- as.vector(X %*% beta_mean)
  return(plogis(eta))
}

#' Fitted Method for GRASP Objects
#'
#' @param object Object of class \code{grasp}.
#' @param ... Additional arguments.
#' @export
fitted.grasp <- function(object, ...) {
  X <- model.matrix(object$formula, object$data)[, -1, drop = FALSE]
  beta_mean <- coef(object)
  if (length(beta_mean) != ncol(X)) {
    var_names <- colnames(X)
    match_idx <- match(var_names, names(beta_mean))
    beta_mean <- beta_mean[match_idx]
  }
  eta <- as.vector(X %*% beta_mean)
  return(plogis(eta))
}

#' Fitted Method for GIGG Objects
#'
#' @param object Object of class \code{gigg}.
#' @param ... Additional arguments.
#' @export
fitted.gigg <- function(object, ...) {
  X <- model.matrix(object$formula, object$data)[, -1, drop = FALSE]
  beta_mean <- coef(object)
  if (length(beta_mean) != ncol(X)) {
    var_names <- colnames(X)
    match_idx <- match(var_names, names(beta_mean))
    beta_mean <- beta_mean[match_idx]
  }
  eta <- as.vector(X %*% beta_mean)
  return(plogis(eta))
}

#' Predict Method for BGHS Objects
#'
#' @param object Object of class \code{bghs}.
#' @param newdata Data frame for predictions.
#' @param type "link" for log-odds or "response" for probabilities.
#' @param ... Additional arguments.
#' @export
predict.bghs <- function(object, newdata, type = c("link", "response"), ...) {
  type <- match.arg(type)
  X_new <- model.matrix(object$formula, newdata)[, -1, drop = FALSE]
  beta_mean <- coef(object)
  if (length(beta_mean) != ncol(X_new)) {
    var_names <- colnames(X_new)
    match_idx <- match(var_names, names(beta_mean))
    beta_mean <- beta_mean[match_idx]
  }
  eta <- as.vector(X_new %*% beta_mean)
  if (type == "link") return(eta)
  return(plogis(eta))
}

#' Predict Method for GRASP Objects
#'
#' @param object Object of class \code{grasp}.
#' @param newdata Data frame for predictions.
#' @param type "link" for log-odds or "response" for probabilities.
#' @param ... Additional arguments.
#' @export
predict.grasp <- function(object, newdata, type = c("link", "response"), ...) {
  type <- match.arg(type)
  X_new <- model.matrix(object$formula, newdata)[, -1, drop = FALSE]
  beta_mean <- coef(object)
  if (length(beta_mean) != ncol(X_new)) {
    var_names <- colnames(X_new)
    match_idx <- match(var_names, names(beta_mean))
    beta_mean <- beta_mean[match_idx]
  }
  eta <- as.vector(X_new %*% beta_mean)
  if (type == "link") return(eta)
  return(plogis(eta))
}

#' Predict Method for GIGG Objects
#'
#' @param object Object of class \code{gigg}.
#' @param newdata Data frame for predictions.
#' @param type "link" for log-odds or "response" for probabilities.
#' @param ... Additional arguments.
#' @export
predict.gigg <- function(object, newdata, type = c("link", "response"), ...) {
  type <- match.arg(type)
  X_new <- model.matrix(object$formula, newdata)[, -1, drop = FALSE]
  beta_mean <- coef(object)
  if (length(beta_mean) != ncol(X_new)) {
    var_names <- colnames(X_new)
    match_idx <- match(var_names, names(beta_mean))
    beta_mean <- beta_mean[match_idx]
  }
  eta <- as.vector(X_new %*% beta_mean)
  if (type == "link") return(eta)
  return(plogis(eta))
}

#' Plot Method for BGHS Objects
#'
#' @param x Object of class \code{bghs}.
#' @param type "interval", "trace", or "density".
#' @param ... Additional arguments.
#' @export
plot.bghs <- function(x, type = c("interval", "trace", "density"), ...) {
  type <- match.arg(type)

  if (type == "interval") {
    beta_cols <- grep("^beta_", colnames(x$posterior))
    beta_means <- colMeans(x$posterior[, beta_cols])
    beta_lower <- apply(x$posterior[, beta_cols], 2, quantile, 0.025)
    beta_upper <- apply(x$posterior[, beta_cols], 2, quantile, 0.975)

    plot_data <- data.frame(
      variable = 1:length(beta_means),
      mean = beta_means,
      lower = beta_lower,
      upper = beta_upper
    )

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = variable, y = mean)) +
        ggplot2::geom_point() +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::labs(title = paste(x$method, "Credible Intervals"),
                      x = "Coefficient Index", y = "Estimate") +
        ggplot2::theme_minimal()
      print(p)
    } else {
      cat("Install ggplot2 for plots.\n")
    }
  } else {
    cat("Trace and density plots not yet implemented.\n")
  }
}

#' Plot Method for GRASP Objects
#'
#' @param x Object of class \code{grasp}.
#' @param type "interval", "trace", or "density".
#' @param ... Additional arguments.
#' @export
plot.grasp <- function(x, type = c("interval", "trace", "density"), ...) {
  type <- match.arg(type)

  if (type == "interval") {
    beta_cols <- grep("^beta_", colnames(x$posterior))
    beta_means <- colMeans(x$posterior[, beta_cols])
    beta_lower <- apply(x$posterior[, beta_cols], 2, quantile, 0.025)
    beta_upper <- apply(x$posterior[, beta_cols], 2, quantile, 0.975)

    plot_data <- data.frame(
      variable = 1:length(beta_means),
      mean = beta_means,
      lower = beta_lower,
      upper = beta_upper
    )

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = variable, y = mean)) +
        ggplot2::geom_point() +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::labs(title = paste(x$method, "Credible Intervals"),
                      x = "Coefficient Index", y = "Estimate") +
        ggplot2::theme_minimal()
      print(p)
    } else {
      cat("Install ggplot2 for plots.\n")
    }
  } else {
    cat("Trace and density plots not yet implemented.\n")
  }
}

#' Plot Method for GIGG Objects
#'
#' @param x Object of class \code{gigg}.
#' @param type "interval", "trace", or "density".
#' @param ... Additional arguments.
#' @export
plot.gigg <- function(x, type = c("interval", "trace", "density"), ...) {
  type <- match.arg(type)

  if (type == "interval") {
    beta_cols <- grep("^beta_", colnames(x$posterior))
    beta_means <- colMeans(x$posterior[, beta_cols])
    beta_lower <- apply(x$posterior[, beta_cols], 2, quantile, 0.025)
    beta_upper <- apply(x$posterior[, beta_cols], 2, quantile, 0.975)

    plot_data <- data.frame(
      variable = 1:length(beta_means),
      mean = beta_means,
      lower = beta_lower,
      upper = beta_upper
    )

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = variable, y = mean)) +
        ggplot2::geom_point() +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::labs(title = paste(x$method, "Credible Intervals"),
                      x = "Coefficient Index", y = "Estimate") +
        ggplot2::theme_minimal()
      print(p)
    } else {
      cat("Install ggplot2 for plots.\n")
    }
  } else {
    cat("Trace and density plots not yet implemented.\n")
  }
}

#' As.mcmc Method for BGHS Objects
#'
#' @param x Object of class \code{bghs}.
#' @param ... Additional arguments.
#' @export
as.mcmc.bghs <- function(x, ...) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for this method")
  }
  return(coda::as.mcmc(x$posterior))
}

#' As.mcmc Method for GRASP Objects
#'
#' @param x Object of class \code{grasp}.
#' @param ... Additional arguments.
#' @export
as.mcmc.grasp <- function(x, ...) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for this method")
  }
  return(coda::as.mcmc(x$posterior))
}

#' As.mcmc Method for GIGG Objects
#'
#' @param x Object of class \code{gigg}.
#' @param ... Additional arguments.
#' @export
as.mcmc.gigg <- function(x, ...) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for this method")
  }
  return(coda::as.mcmc(x$posterior))
}
