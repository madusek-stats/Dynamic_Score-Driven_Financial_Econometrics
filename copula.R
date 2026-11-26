# Dynamic score-driven Gaussian copula (three series).
# Requires Matrix only if a correlation matrix needs projection.
# Input columns must be aligned chronologically; their order is preserved.

# Align PITs returned by fit_dcs_burr() in marginals.R.
# A named list controls the asset order. Dates are joined by intersection;
# without dates, the caller must ensure observations correspond by row.

align_marginal_pits <- function(fits) {
  stopifnot(is.list(fits), length(fits) == 3L)
  nm <- names(fits)
  if (is.null(nm)) nm <- paste0("series", seq_along(fits))
  stopifnot(!anyNA(nm), all(nzchar(nm)), !anyDuplicated(nm))
  pits <- lapply(fits, function(f) {
    u <- f$pit$u
    stopifnot(is.numeric(u), length(u) > 0L,
              all(is.finite(u)), all(u >= 0 & u <= 1))
    u
  })
  dates <- lapply(fits, function(f) f$dates)
  dated <- !vapply(dates, is.null, logical(1))
  if (any(dated) && !all(dated)) stop("Supply dates for all or none of the fits.")
  if (all(dated)) {
    for (j in seq_along(dates)) {
      dates[[j]] <- as.Date(dates[[j]])
      stopifnot(length(dates[[j]]) == length(pits[[j]]),
                !anyNA(dates[[j]]), !is.unsorted(dates[[j]], strictly = TRUE))
    }
    keys <- lapply(dates, as.character)
    common <- sort(Reduce(intersect, keys))
    if (!length(common)) stop("No common dates among marginal fits.")
    U <- do.call(cbind, lapply(seq_along(pits), function(j) {
      pits[[j]][match(common, keys[[j]])]
    }))
    dates_out <- as.Date(common)
  } else {
    stopifnot(length(unique(lengths(pits))) == 1L)
    U <- do.call(cbind, pits)
    dates_out <- NULL
  }
  colnames(U) <- nm
  list(U = U, dates = dates_out)
}

# Correlations follow pair order (1,2), (1,3), (2,3).
make_R <- function(rho) {
  stopifnot(is.numeric(rho), length(rho) == 3L, all(is.finite(rho)),
            all(abs(rho) <= 1))
  R <- diag(3L)
  R[1, 2] <- R[2, 1] <- rho[1]
  R[1, 3] <- R[3, 1] <- rho[2]
  R[2, 3] <- R[3, 2] <- rho[3]
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Install Matrix for positive-definite correlation projection.")
    }
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
  }
  R
}


score_rho_gaussian <- function(x, R) {
  Ri <- solve(R)
  G <- 0.5 * (Ri %*% tcrossprod(x) %*% Ri - Ri)
  c(G[1, 2], G[1, 3], G[2, 3])
}

# par = (omega[1:3], alpha[1:3], atanh_beta[1:3]).
# X consists of Gaussian scores qnorm(U), in the caller's column order.
# The original sample-correlation initialisation uses all supplied rows.
dcs_gauss_copula_filter <- function(par, X, rho_init = NULL) {
  stopifnot(is.numeric(par), length(par) == 9L,
            is.matrix(X), is.numeric(X), ncol(X) == 3L,
            nrow(X) >= 2L, all(is.finite(X)))
  if (any(!is.finite(par))) return(NULL)
  omega <- unname(par[1:3])
  alpha <- unname(par[4:6])
  beta <- tanh(unname(par[7:9]))
  n <- nrow(X)
  h <- rho <- matrix(0, n, 3L)
  ll <- numeric(n)
  if (is.null(rho_init)) {
    R0 <- stats::cor(X)
    rho_init <- c(R0[1, 2], R0[1, 3], R0[2, 3])
  }
  stopifnot(is.numeric(rho_init), length(rho_init) == 3L,
            all(is.finite(rho_init)), all(abs(rho_init) <= 1))
  h[1, ] <- atanh(pmin(pmax(rho_init, -0.95), 0.95))
  for (t in seq_len(n)) {
    rho[t, ] <- tanh(h[t, ])
    R <- make_R(rho[t, ])
    x <- X[t, ]
    Ri <- solve(R)
    # Original Gaussian likelihood kernel (constant in parameters omitted).
    ll[t] <- -0.5 * (as.numeric(determinant(R, logarithm = TRUE)$modulus) +
                      drop(crossprod(x, Ri %*% x)))
    if (!is.finite(ll[t])) return(NULL)
    if (t < n) {
      s <- score_rho_gaussian(x, R) * (1 - rho[t, ]^2)
      if (any(!is.finite(s))) return(NULL)
      h[t + 1L, ] <- pmax(pmin(omega + beta * h[t, ] + alpha * s, 5), -5)
    }
  }
  list(h = h, rho = rho, ll = ll, beta = beta)
}

negloglik_gauss_copula <- function(par, X, rho_init = NULL) {
  out <- dcs_gauss_copula_filter(par, X, rho_init)
  if (is.null(out) || !is.finite(sum(out$ll))) return(1e12)
  -sum(out$ll)
}

# Fit by BFGS. U is an n x 3 matrix of aligned marginal PITs.
fit_gaussian_copula <- function(U, dates = NULL, init = NULL,
                                rho_init = NULL, eps = 1e-6,
                                control = list(maxit = 5000, reltol = 1e-10)) {
  U <- as.matrix(U)
  stopifnot(is.numeric(U), ncol(U) == 3L, nrow(U) > 200L,
            all(is.finite(U)), all(U >= 0 & U <= 1),
            is.numeric(eps), length(eps) == 1L, is.finite(eps),
            eps > 0, eps < 0.5)
  nm <- colnames(U)
  if (is.null(nm)) nm <- paste0("series", 1:3)
  stopifnot(!anyNA(nm), all(nzchar(nm)), !anyDuplicated(nm))
  colnames(U) <- nm
  if (!is.null(dates)) {
    stopifnot(length(dates) == nrow(U), !anyNA(dates),
              !is.unsorted(dates, strictly = TRUE))
  }
  U[] <- pmin(pmax(U, eps), 1 - eps)
  X <- U
  X[] <- stats::qnorm(U)
  stopifnot(all(apply(X, 2, stats::sd) > 0))
  if (is.null(rho_init)) {
    R0 <- stats::cor(X)
    rho_init <- c(R0[1, 2], R0[1, 3], R0[2, 3])
  }
  if (is.null(init)) init <- c(rep(0, 3), rep(0.05, 3), rep(atanh(0.95), 3))
  stopifnot(is.numeric(init), length(init) == 9L, all(is.finite(init)))
  pairs <- c(paste(nm[1], nm[2], sep = "__"),
             paste(nm[1], nm[3], sep = "__"),
             paste(nm[2], nm[3], sep = "__"))
  names(init) <- c(paste0("omega_", pairs), paste0("alpha_", pairs),
                   paste0("atanh_beta_", pairs))
  fit <- stats::optim(init, negloglik_gauss_copula, X = X, rho_init = rho_init,
                      method = "BFGS", control = control)
  if (fit$convergence != 0L) stop("BFGS optimisation did not converge.")
  flt <- dcs_gauss_copula_filter(fit$par, X, rho_init)
  if (is.null(flt)) stop("Invalid fitted copula filter.")
  n <- nrow(U)
  R_t <- array(NA_real_, c(3L, 3L, n), dimnames = list(nm, nm, NULL))
  rho_t <- flt$rho
  colnames(rho_t) <- colnames(flt$rho) <- colnames(flt$h) <- pairs
  for (t in seq_len(n)) {
    R <- make_R(flt$rho[t, ])
    R_t[, , t] <- R
    rho_t[t, ] <- c(R[1, 2], R[1, 3], R[2, 3])
  }
  # The actual copula log-density adds 0.5*x*x to the original kernel.
  # This constant does not alter the fitted parameters for fixed U.
  copula_ll <- flt$ll + 0.5 * rowSums(X^2)
  list(dates = dates, assets = nm, U = U, rho_t = rho_t, R_t = R_t,
       rho_raw = flt$rho, h = flt$h, params = fit$par, beta = flt$beta,
       rho_init = rho_init, logLik = sum(copula_ll), logLik_kernel = -fit$value,
       log_density = copula_ll, convergence = fit$convergence, eps = eps)
}

# Sequential conditional PITs in the supplied column order.
rosenblatt_gaussian_tv <- function(U, R_t, eps = 1e-6) {
  U <- as.matrix(U)
  stopifnot(is.numeric(U), ncol(U) == 3L, nrow(U) > 0L,
            all(is.finite(U)), all(U >= 0 & U <= 1),
            identical(dim(R_t), c(3L, 3L, nrow(U))),
            all(is.finite(R_t)), length(eps) == 1L, is.finite(eps),
            eps > 0, eps < 0.5)
  nm <- colnames(U)
  rn <- dimnames(R_t)[[1L]]
  cn <- dimnames(R_t)[[2L]]
  if (!is.null(nm) && (!is.null(rn) || !is.null(cn))) {
    stopifnot(identical(nm, rn), identical(nm, cn))
  }
  if (is.null(nm)) nm <- paste0("series", 1:3)
  U[] <- pmin(pmax(U, eps), 1 - eps)
  Z <- U
  Z[] <- stats::qnorm(U)
  V <- U
  for (t in seq_len(nrow(U))) {
    R <- R_t[, , t]
    z <- Z[t, ]
    r12 <- R[1, 2]
    V[t, 2] <- stats::pnorm((z[2] - r12 * z[1]) /
                             sqrt(max(1 - r12^2, 1e-12)))
    R12 <- R[1:2, 1:2, drop = FALSE]
    r3 <- R[3, 1:2, drop = FALSE]
    mu3 <- drop(r3 %*% solve(R12, z[1:2]))
    var3 <- 1 - drop(r3 %*% solve(R12, t(r3)))
    V[t, 3] <- stats::pnorm((z[3] - mu3) / sqrt(max(var3, 1e-12)))
  }
  colnames(V) <- c(paste0("v_", nm[1]),
                   paste0("v_", nm[2], "_given_", nm[1]),
                   paste0("v_", nm[3], "_given_", nm[1], "_", nm[2]))
  V[] <- pmin(pmax(V, eps), 1 - eps)
  V
}

# Uniformity and serial-dependence diagnostics. 
rosenblatt_diagnostics <- function(V, lag = 20L) {
  V <- as.matrix(V)
  stopifnot(is.numeric(V), ncol(V) == 3L, all(is.finite(V)),
            all(V > 0 & V < 1), length(lag) == 1L, is.finite(lag),
            lag == as.integer(lag), lag >= 1L, nrow(V) > lag + 1L)
  Z <- V
  Z[] <- stats::qnorm(V)
  nm <- colnames(V)
  if (is.null(nm)) nm <- paste0("v", 1:3)
  data.frame(Series = nm,
             KS_p_uniform = apply(V, 2, function(x) stats::ks.test(x, "punif")$p.value),
             LjungBox_p_z = apply(Z, 2, function(x) stats::Box.test(x, lag = lag, type = "Ljung-Box")$p.value),
             LjungBox_p_z2 = apply(Z^2, 2, function(x) stats::Box.test(x, lag = lag, type = "Ljung-Box")$p.value),
             row.names = NULL)
}
