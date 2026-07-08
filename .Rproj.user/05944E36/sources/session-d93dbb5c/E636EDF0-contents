# ====================================================================
# -------------------- Gibbs sampling for BGHS -----------------------
# ====================================================================

#' @noRd
bghs_gibbs = function(data, nsamples = 10000, burnin = 1000,
                      init, group_sizes, print_every = 100){

  n = nrow(data)
  p = ncol(data) - 1
  X = as.matrix(data[,-1])
  y = as.vector(data[,1])
  G = length(group_sizes)
  expected_len = n + 3*p + 2*G + 2

  # Few preliminary checks
  if(length(init) != expected_len) stop("init length does not match number of parameters.")
  if(burnin < 0) stop("burnin must be non-negative.")
  if(burnin >= nsamples) stop("burnin must be smaller nsamples.")
  if(sum(group_sizes) != p) stop("sum(group_sizes) must equal no. of covariates.")

  # Creating the sample dataframe
  colname = c(paste0("w_",1:n), paste0("beta_",1:p),
              unlist(Map(function(i,n) paste0("lam_sq_",i,1:n), 1:G, group_sizes)),
              unlist(Map(function(i,n) paste0("t_",i,1:n), 1:G, group_sizes)),
              paste0("delta_sq_",1:G),
              paste0("c_",1:G),
              "tau_sq","v"
  )
  sample = matrix(NA, nrow = nsamples, ncol = length(colname))
  colnames(sample) = colname

  # Tweaking initials
  current_w = init[1:n]
  current_beta = init[(n+1) :(n+p)]
  current_lambda_sq = init[(n+p+1) : (n+2*p)]
  current_t = init[(n+2*p+1) : (n+3*p)]
  current_delta_sq = init[(n+3*p+1) : (n+3*p+G)]
  current_c = init[(n+3*p+G+1) : (n+3*p+2*G)]
  current_tau_sq = init[n+3*p+2*G+1]
  current_v = init[n+3*p+2*G+2]

  # initial row
  sample[1,] = init

  # Few predefined things
  k = y - 0.5
  idx = split(1:length(current_beta), rep(1:G, times = group_sizes)) # for generating delta square

  # Time related (Initial)
  start_time = Sys.time()
  cat("\n===============================================================\n")
  cat("                Starting BGHS Gibbs Sampler\n")
  cat("===============================================================\n")
  cat("\n\n")
  cat("----------------------------------------\n")
  cat("Summary of experimental setting\n")
  cat("----------------------------------------\n\n")
  cat(sprintf("Observations    : %d\n", n))
  cat(sprintf("Covariates      : %d\n", p))
  cat(sprintf("Number of groups: %d\n", G))
  cat(sprintf("Group sizes     : %s\n",paste(group_sizes,collapse = ", ")))

  cat("\n----------------------------------------\n\n")

  for(i in 2:nsamples){

    # progress printing
    if(i %% print_every == 0){
      elapsed = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      avg_time = elapsed/(i-1)
      eta_remaining = avg_time * (nsamples - i)
      pct = round(100*i/nsamples,1)
      current_time = format(Sys.time(),"%H:%M:%S")

      cat(sprintf(
        "Time: %s => Samples: [%d/%d] | Completed: %.1f%% | elapsed: %.1fs | ETA: %.1fs\n",
        current_time, i, nsamples, pct, elapsed, eta_remaining
      ))
    }

    # for generating w
    eta = as.vector(X %*% current_beta)
    current_w = BayesLogit::rpg(num= n, h = 1, z = eta)

    # for generating beta
    XtOmegaX = crossprod(X, X*current_w)
    D_inv = diag(1 / (rep(current_delta_sq,times = group_sizes) * current_lambda_sq * current_tau_sq))
    A = XtOmegaX + D_inv
    R = chol(A)
    Sigma_beta = chol2inv(R)
    mu_beta = Sigma_beta %*% crossprod(X,k)
    current_beta = MASS::mvrnorm(n = 1, mu = as.vector(mu_beta),
                                 Sigma = Sigma_beta)

    # for generating lambda square
    current_lambda_sq = MCMCpack::rinvgamma(p, shape = 1, scale = (current_beta^2)/(2*current_tau_sq*rep(current_delta_sq, group_sizes)) + (1/current_t))

    # for generating t
    current_t = MCMCpack::rinvgamma(p, shape = 1, scale = 1 + (1/current_lambda_sq))

    # for generating delta square
    current_delta_sq = sapply(1:G,function(g){

      ind = idx[[g]]
      shape = (group_sizes[g] + 1)/2
      scale = (1/(2*current_tau_sq))*sum(current_beta[ind]^2 / current_lambda_sq[ind]) + 1/current_c[g]
      MCMCpack::rinvgamma(1, shape = shape, scale = scale)
    })

    # for generating c
    current_c = MCMCpack::rinvgamma(G, shape = 1, scale = 1 + (1/current_delta_sq))

    # for generating tau
    current_tau_sq = MCMCpack::rinvgamma(1,
                                         shape = (p+1)/2,
                                         scale =
                                           0.5*sum(current_beta^2/
                                                     (current_lambda_sq *
                                                        rep(current_delta_sq,
                                                            times = group_sizes))) +
                                           (1/current_v))

    # for generating v
    current_v = MCMCpack::rinvgamma(1,shape = 1, scale = 1 +
                                      (1/current_tau_sq))

    new_row = c(current_w, current_beta, current_lambda_sq, current_t, current_delta_sq, current_c, current_tau_sq, current_v)
    sample[i,] = new_row # adding new row
  }

  sample = as.data.frame(sample)
  post_idx = seq(burnin + 1, nsamples, by = 1)
  posterior_sample = sample[post_idx,]
  rownames(posterior_sample) =seq_len(nrow(posterior_sample))

  # Time related (End)
  elapsed =as.numeric(difftime(Sys.time(),start_time,units = "secs"))
  cat("\n----------------------------------------\n")
  cat("Run Summary\n")
  cat("----------------------------------------\n")
  cat(sprintf("Total samples : %d\n",nsamples))
  cat(sprintf("Burn-in       : %d\n",burnin))
  cat(sprintf("Posterior kept: %d\n",nrow(posterior_sample)))
  elapsed_min = floor(elapsed / 60)
  elapsed_sec = round(elapsed %% 60)
  cat(sprintf("Elapsed time  : %.1f seconds (%dm %02ds)\n",elapsed, elapsed_min, elapsed_sec))
  cat("----------------------------------------\n\n")
  return(posterior_sample)
}
