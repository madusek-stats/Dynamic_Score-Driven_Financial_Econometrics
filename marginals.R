# Dynamic Conditional Score (DCS) Burr XII marginal model
# Input observations must be chronological, finite and strictly positive.
# Optional inference requires the numDeriv and MASS packages.

# Construct a range-based volatility proxy from caller-supplied prices.
log_price_range <- function(high, low) {
  stopifnot(is.numeric(high), is.numeric(low),
            length(high) == length(low), length(high) > 0L,
            all(is.finite(high)), all(is.finite(low)),
            all(low > 0), all(high > low))
  log(high) - log(low)
}

# Burr XII conditional log-density.
burr_logpdf <- function(y, c, k, lambda) {
  z <- y / lambda
  log(c) + log(k) - log(lambda) +
    (c - 1) * log(z) - (k + 1) * log1p(z^c)
}

# Analytical scores with respect to f = log(lambda) and g = log(k).
score_f <- function(y, c, k, lambda) {
  zc <- (y / lambda)^c
  -c + (k + 1) * c * zc / (1 + zc)
}

score_g <- function(y, c, k, lambda) {
  zc <- (y / lambda)^c
  1 - k * log1p(zc)
}

# Diagonal Fisher information scaling as in the original implementation.
I_fisher_f <- function(c, k) c^2 * k / (k + 2)
I_fisher_g <- function() 1

# Burr XII CDF for standardised observations z = y / lambda.
burr_cdf <- function(z, c, k) {
  -expm1(-k * log1p(z^c))
}

# Score-driven filter: static shape c and dynamic scale lambda and tail k.
# theta = (logc, omega_f, alpha_f, atanh_beta_f,
#           omega_g, alpha_g, atanh_beta_g).
dcs_burr_filter <- function(theta, y) {
  stopifnot(is.numeric(theta), length(theta) == 7L,
            is.numeric(y), length(y) > 0L,
            all(is.finite(y)), all(y > 0))
  if (any(!is.finite(theta))) return(NULL)

  c <- exp(theta[1])
  omega_f <- theta[2]; alpha_f <- theta[3]; beta_f <- tanh(theta[4])
  omega_g <- theta[5]; alpha_g <- theta[6]; beta_g <- tanh(theta[7])
  if (!is.finite(c) || c <= 0 || abs(beta_f) >= 1 || abs(beta_g) >= 1) {
    return(NULL)
  }

  n <- length(y)
  f <- g <- lambda <- k <- ll <- numeric(n)

  # Initialise at the unconditional state means.
  f[1] <- omega_f / (1 - beta_f)
  g[1] <- omega_g / (1 - beta_g)
  lambda[1] <- exp(f[1])
  k[1] <- exp(g[1])

  for (t in seq_len(n)) {
    if (!is.finite(lambda[t]) || !is.finite(k[t]) ||
        lambda[t] <= 0 || k[t] <= 0) return(NULL)
    ll[t] <- burr_logpdf(y[t], c, k[t], lambda[t])
    if (!is.finite(ll[t])) return(NULL)

    if (t < n) {
      s_f <- score_f(y[t], c, k[t], lambda[t]) / I_fisher_f(c, k[t])
      s_g <- score_g(y[t], c, k[t], lambda[t]) / I_fisher_g()
      if (!all(is.finite(c(s_f, s_g)))) return(NULL)

      f[t + 1L] <- omega_f + beta_f * f[t] + alpha_f * s_f
      g[t + 1L] <- omega_g + beta_g * g[t] + alpha_g * s_g
      lambda[t + 1L] <- exp(f[t + 1L])
      k[t + 1L] <- exp(g[t + 1L])
    }
  }

  list(c = unname(c), lambda = lambda, k = k, ll = ll)
}

# Negative log-likelihood for numerical optimisation.
negloglik <- function(theta, y) {
  out <- dcs_burr_filter(theta, y)
  if (is.null(out) || !is.finite(sum(out$ll))) return(1e12)
  -sum(out$ll)
}

# Joint maximum likelihood estimation, optional inference and PIT generation.

fit_dcs_burr <- function(y, dates = NULL, asset = NULL, init = NULL,
                         inference = TRUE,
                         control = list(maxit = 5000, reltol = 1e-10)) {
  stopifnot(is.numeric(y), length(y) > 200L,
            all(is.finite(y)), all(y > 0),
            is.logical(inference), length(inference) == 1L, !is.na(inference))
  if (!is.null(dates)) {
    stopifnot(length(dates) == length(y), !anyNA(dates),
              !is.unsorted(dates, strictly = TRUE))
  }
  if (!is.null(asset)) stopifnot(is.character(asset), length(asset) == 1L)

  # Original parameter initialisation.
  if (is.null(init)) {
    init <- c(logc = log(2),
              omega_f = log(median(y)) * 0.05,
              alpha_f = 0.05,
              atanh_beta_f = atanh(0.95),
              omega_g = log(2) * 0.1,
              alpha_g = 0.05,
              atanh_beta_g = atanh(0.90))
  }
  stopifnot(is.numeric(init), length(init) == 7L, all(is.finite(init)))
  names(init) <- c("logc", "omega_f", "alpha_f", "atanh_beta_f",
                   "omega_g", "alpha_g", "atanh_beta_g")

  fit <- stats::optim(par = init, fn = negloglik, y = y,
                      method = "BFGS", control = control)
  if (fit$convergence != 0L) stop("BFGS optimisation did not converge.")
  flt <- dcs_burr_filter(fit$par, y)
  if (is.null(flt)) stop("Invalid fitted filter.")

  # Observed Hessian, covariance matrix and delta-method SE for shape c.
  H <- vcov <- se <- se_c <- NULL
  if (inference) {
    if (!requireNamespace("numDeriv", quietly = TRUE) ||
        !requireNamespace("MASS", quietly = TRUE)) {
      stop("Inference requires numDeriv and MASS; or set inference = FALSE.")
    }
    H <- numDeriv::hessian(func = negloglik, x = fit$par, y = y)
    H <- (H + t(H)) / 2
    dimnames(H) <- list(names(fit$par), names(fit$par))
    if (!all(is.finite(H))) stop("Non-finite likelihood Hessian.")
    if (min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      warning("Hessian is not positive definite; standard errors may be unreliable.")
    }
    vcov <- tryCatch(solve(H), error = function(e) {
      warning("Singular Hessian; using a generalised inverse.")
      MASS::ginv(H)
    })
    dimnames(vcov) <- dimnames(H)
    # Do not disguise negative variance estimates as zero standard errors.
    variances <- diag(vcov)
    se <- sqrt(ifelse(variances >= 0, variances, NA_real_))
    se_c <- unname(flt$c * se[1])
  }

  # Probability integral transforms for subsequent copula estimation.
  u <- burr_cdf(y / flt$lambda, flt$c, flt$k)
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
  pit <- data.frame(u = u, lambda = flt$lambda, k = flt$k)
  if (!is.null(dates)) pit <- data.frame(Date = dates, pit)

  # Retain the original fitted-object field names for downstream modules.
  list(asset = asset, dates = dates, y = y, c = flt$c,
       lambda_t = flt$lambda, k_t = flt$k, theta = fit$par,
       vcov = vcov, logLik = -fit$value, pit = pit,
       hessian = H, se = se, se_c = se_c,
       convergence = fit$convergence)
}
