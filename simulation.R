# Monte Carlo predictive returns: Gaussian copula on Burr volatility ranges.
# Columns are explicitly matched by name. Dates label TARGET observations.

burr_q <- function(u, c, k, lambda) {
  stopifnot(all(is.finite(u)), all(u >= 0 & u <= 1),
            all(is.finite(c) & c > 0), all(is.finite(k) & k > 0),
            all(is.finite(lambda) & lambda > 0))
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  lambda * exp(log(expm1(-log1p(-u) / k)) / c)
}

# Parkinson range-to-standard-deviation conversion.
parkinson_sigma <- function(range) range / sqrt(4 * log(2))

# Optional conditional mean uses only returns BEFORE the target row.
lagged_rolling_mean <- function(returns, window = 252L,
                                min_obs = max(30L, floor(0.6 * window))) {
  returns <- as.matrix(returns)
  stopifnot(is.numeric(returns), nrow(returns) > 0L, window >= 1,
            min_obs >= 1, min_obs <= window)
  out <- matrix(NA_real_, nrow(returns), ncol(returns), dimnames = dimnames(returns))
  for (i in seq_len(nrow(returns))[-1L]) {
    x <- returns[max(1L, i - window):(i - 1L), , drop = FALSE]
    for (j in seq_len(ncol(x))) {
      v <- x[is.finite(x[, j]), j]
      if (length(v) >= min_obs) out[i, j] <- mean(v)
    }
  }
  out
}

simulate_returns <- function(R, shape, tail, scale, mu, B = 20000L) {
  R <- as.matrix(R); d <- ncol(R)
  stopifnot(is.numeric(R), d == 3L, nrow(R) == d,
            all(is.finite(R)), isSymmetric(unname(R)),
            all(abs(diag(R) - 1) < 1e-8), B >= 2, B == as.integer(B))
  nm <- colnames(R)
  if (is.null(nm)) stop('Supply named correlation matrix columns.')
  stopifnot(identical(rownames(R), nm), !anyDuplicated(nm))
  align <- function(x) {
    stopifnot(is.numeric(x), length(x) == d, !anyDuplicated(names(x)),
              setequal(names(x), nm), all(is.finite(x)))
    x[nm]
  }
  shape <- align(shape); tail <- align(tail); scale <- align(scale); mu <- align(mu)
  Z <- matrix(stats::rnorm(B * d), ncol = d) %*% chol(R)
  U <- stats::pnorm(Z)
  ranges <- vapply(seq_len(d), function(j) burr_q(U[, j], shape[j], tail[j], scale[j]), numeric(B))
  eps <- matrix(stats::rnorm(B * d), ncol = d)
  r <- sweep(parkinson_sigma(ranges) * eps, 2, mu, '+')
  colnames(r) <- nm
  if (any(!is.finite(r))) stop('Non-finite simulated returns; inspect tail parameters.')
  r
}

prepare_simulation <- function(fits, copula_fit, returns, return_dates = NULL,
                               window = 252L, min_obs = max(30L, floor(0.6 * window))) {
  if (!is.null(copula_fit$nu)) stop("This simulator requires a Gaussian copula fit.")
  nm <- copula_fit$assets; R <- copula_fit$R_t
  returns <- as.matrix(returns)
  stopifnot(length(nm) == 3L, !anyDuplicated(nm), setequal(names(fits), nm),
            setequal(colnames(returns), nm), !anyDuplicated(colnames(returns)),
            identical(dimnames(R)[[1]], nm), identical(dimnames(R)[[2]], nm))
  fits <- fits[nm]; returns <- returns[, nm, drop = FALSE]
  n <- dim(R)[3]; dates <- copula_fit$dates
  mu_all <- lagged_rolling_mean(returns, window, min_obs)
  if (is.null(dates)) {
    stopifnot(is.null(return_dates), all(vapply(fits, function(f) is.null(f$dates), logical(1))),
              nrow(returns) == n, all(vapply(fits, function(f) length(f$lambda_t) == n, logical(1))))
    idx <- rep(list(seq_len(n)), 3L); ir <- seq_len(n)
  } else {
    check_dates <- function(x, n) {
      stopifnot(length(x) == n, !anyNA(x), !is.unsorted(x, strictly = TRUE))
      as.character(x)
    }
    keys <- check_dates(dates, n)
    ir <- match(keys, check_dates(return_dates, nrow(returns)))
    idx <- lapply(fits, function(f) match(keys, check_dates(f$dates, length(f$lambda_t))))
    if (anyNA(ir) || any(vapply(idx, anyNA, logical(1)))) stop('Missing target dates.')
  }
  extract <- function(field) {
    out <- vapply(seq_along(fits), function(j) fits[[j]][[field]][idx[[j]]], numeric(n))
    colnames(out) <- nm; out
  }
  list(dates = dates, assets = nm, R_t = R,
       shape = vapply(fits, function(f) f$c, numeric(1)),
       tail = extract('k_t'), scale = extract('lambda_t'), mu = mu_all[ir, , drop = FALSE],
       returns = returns[ir, , drop = FALSE])
}

simulate_returns_t <- function(state, t_idx, B = 20000L) {
  stopifnot(length(t_idx) == 1L, t_idx %in% seq_len(dim(state$R_t)[3]))
  simulate_returns(state$R_t[, , t_idx], state$shape, state$tail[t_idx, ],
                   state$scale[t_idx, ], state$mu[t_idx, ], B)
}
