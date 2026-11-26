# Dynamic Gaussian versus Student-t copula model selection (three series).
# Load copula.R first; Matrix is required if correlation projection is needed.
# No input/output operations: supply an aligned matrix of marginal PITs.

# Full Student-t copula log-density at t-quantiles x. Q may be supplied
# by the filter to avoid repeating matrix operations.
student_t_copula_logpdf <- function(x, R, nu, Q = NULL) {
  stopifnot(is.numeric(x), length(x) == 3L, all(is.finite(x)),
            is.matrix(R), identical(dim(R), c(3L, 3L)),
            length(nu) == 1L, is.finite(nu), nu > 2)
  d <- length(x)
  C <- chol(R)
  logdet <- 2 * sum(log(diag(C)))
  if (is.null(Q)) Q <- sum(forwardsolve(t(C), x)^2)
  joint <- lgamma((nu + d) / 2) - lgamma(nu / 2) -
    (d / 2) * log(nu * pi) - 0.5 * logdet -
    ((nu + d) / 2) * log1p(Q / nu)
  joint - sum(stats::dt(x, df = nu, log = TRUE))
}


score_rho_student_t <- function(x, R, nu, Ri = NULL, Q = NULL) {
  if (is.null(Ri)) Ri <- chol2inv(chol(R))
  if (is.null(Q)) Q <- drop(crossprod(x, Ri %*% x))
  w <- (nu + length(x)) / (nu + Q)
  G <- 0.5 * (w * (Ri %*% tcrossprod(x) %*% Ri) - Ri)
  c(G[1, 2], G[1, 3], G[2, 3])
}

# Internal recursion on already transformed PITs. Parameter order is
# omega[1:3], alpha[1:3], atanh_beta[1:3], log(nu - 2).
# Use common Gaussian-score initial correlations for both model families.
.student_t_filter <- function(par, U, rho_init, keep_R = FALSE) {
  if (any(!is.finite(par))) return(NULL)
  nu <- 2 + exp(par[10L])
  if (!is.finite(nu) || nu <= 2) return(NULL)
  omega <- unname(par[1:3])
  alpha <- unname(par[4:6])
  beta <- tanh(unname(par[7:9]))
  X <- U
  X[] <- stats::qt(U, df = nu)
  if (any(!is.finite(X))) return(NULL)
  n <- nrow(U)
  h <- rho <- rho_raw <- matrix(0, n, 3L)
  ll <- numeric(n)
  R_t <- if (keep_R) array(NA_real_, c(3L, 3L, n),
                           dimnames = list(colnames(U), colnames(U), NULL)) else NULL
  h[1L, ] <- atanh(pmin(pmax(rho_init, -0.95), 0.95))
  d <- ncol(U)
  # Terms depending on nu must remain in the estimation objective.
  constant <- lgamma((nu + d) / 2) - lgamma(nu / 2) -
    (d / 2) * log(nu * pi)
  marginal_ll <- rowSums(matrix(stats::dt(X, df = nu, log = TRUE), nrow = n))
  for (t in seq_len(n)) {
    rho_raw[t, ] <- tanh(h[t, ])
    R <- make_R(rho_raw[t, ])
    C <- chol(R)
    Ri <- chol2inv(C)
    x <- X[t, ]
    Q <- sum(forwardsolve(t(C), x)^2)
    ll[t] <- constant - sum(log(diag(C))) -
      ((nu + d) / 2) * log1p(Q / nu) - marginal_ll[t]
    if (!is.finite(ll[t])) return(NULL)
    rho[t, ] <- c(R[1, 2], R[1, 3], R[2, 3])
    if (keep_R) R_t[, , t] <- R
    if (t < n) {
      s <- score_rho_student_t(x, R, nu, Ri, Q) * (1 - rho_raw[t, ]^2)
      if (any(!is.finite(s))) return(NULL)
      h[t + 1L, ] <- pmax(pmin(omega + beta * h[t, ] + alpha * s, 5), -5)
    }
  }
  list(h = h, rho = rho, rho_raw = rho_raw, R_t = R_t,
       ll = ll, beta = beta, nu = nu)
}

# Validate PIT inputs once rather than at every likelihood evaluation.
.prepare_student_t_inputs <- function(U, dates, rho_init, eps) {
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
  stopifnot(is.numeric(rho_init), length(rho_init) == 3L,
            all(is.finite(rho_init)), all(abs(rho_init) <= 1))
  list(U = U, dates = dates, assets = nm, rho_init = rho_init)
}

# Fit the dynamic Student-t copula by constrained maximum likelihood.
# Finite nu bounds improve stability close to the Gaussian limit.
# Default bounds are configurable; a boundary estimate is flagged.
fit_student_t_copula <- function(U, dates = NULL, init = NULL,
                                 rho_init = NULL, eps = 1e-6,
                                 nu_bounds = c(2.01, 500),
                                 control = list(maxit = 5000, factr = 1e7)) {
  if (!exists("make_R", mode = "function")) stop("Load copula.R first.")
  z <- .prepare_student_t_inputs(U, dates, rho_init, eps)
  U <- z$U
  rho_init <- z$rho_init
  nm <- z$assets
  stopifnot(is.numeric(nu_bounds), length(nu_bounds) == 2L,
            all(is.finite(nu_bounds)), nu_bounds[1] > 2,
            nu_bounds[2] > nu_bounds[1])
  if (is.null(init)) {
    nu0 <- max(nu_bounds[1], min(10, nu_bounds[2]))
    init <- c(rep(0, 3), rep(0.05, 3), rep(atanh(0.95), 3), log(nu0 - 2))
  }
  lower <- c(rep(-Inf, 9), log(nu_bounds[1] - 2))
  upper <- c(rep(Inf, 9), log(nu_bounds[2] - 2))
  stopifnot(is.numeric(init), length(init) == 10L, all(is.finite(init)),
            all(init >= lower & init <= upper))
  pairs <- c(paste(nm[1], nm[2], sep = "__"),
             paste(nm[1], nm[3], sep = "__"),
             paste(nm[2], nm[3], sep = "__"))
  names(init) <- c(paste0("omega_", pairs), paste0("alpha_", pairs),
                   paste0("atanh_beta_", pairs), "log_nu_minus2")
  objective <- function(par) {
    flt <- tryCatch(.student_t_filter(par, U, rho_init), error = function(e) NULL)
    if (is.null(flt) || !is.finite(sum(flt$ll))) return(1e12)
    -sum(flt$ll)
  }
  fit <- stats::optim(init, objective, method = "L-BFGS-B",
                      lower = lower, upper = upper, control = control)
  if (fit$convergence != 0L) stop("Student-t optimisation did not converge: ", fit$message)
  flt <- .student_t_filter(fit$par, U, rho_init, keep_R = TRUE)
  if (is.null(flt)) stop("Invalid fitted Student-t filter.")
  at_boundary <- any(abs(fit$par[10L] - c(lower[10L], upper[10L])) < 1e-4)
  if (at_boundary) warning("Degrees of freedom reached an estimation bound; interpret AIC cautiously.")
  colnames(flt$rho) <- colnames(flt$rho_raw) <- colnames(flt$h) <- pairs
  list(dates = dates, assets = nm, U = U, rho_t = flt$rho,
       R_t = flt$R_t, rho_raw = flt$rho_raw, h = flt$h,
       params = fit$par, beta = flt$beta, rho_init = rho_init,
       nu = flt$nu, logLik = sum(flt$ll), log_density = flt$ll,
       convergence = fit$convergence, nu_bounds = nu_bounds,
       nu_at_boundary = at_boundary, eps = eps)
}

# Summary statistics for the actual (post-projection) correlations.
correlation_summary <- function(rho_t) {
  stopifnot(is.matrix(rho_t), is.numeric(rho_t), ncol(rho_t) == 3L,
            nrow(rho_t) > 1L, all(is.finite(rho_t)))
  out <- t(apply(rho_t, 2, function(x) {
    c(Mean = mean(x), SD = stats::sd(x), Min = min(x),
      Q05 = stats::quantile(x, 0.05, names = FALSE),
      Q95 = stats::quantile(x, 0.95, names = FALSE), Max = max(x))
  }))
  out
}

# Fit both families to exactly the same PIT observations and compare AIC.
compare_dynamic_copulas <- function(U, dates = NULL, rho_init = NULL,
                                    eps = 1e-6, nu_bounds = c(2.01, 500),
                                    gaussian_init = NULL, student_init = NULL,
                                    gaussian_control = list(maxit = 5000, reltol = 1e-10),
                                    student_control = list(maxit = 5000, factr = 1e7)) {
  if (!exists("fit_gaussian_copula", mode = "function")) stop("Load copula.R first.")
  g <- fit_gaussian_copula(U, dates = dates, init = gaussian_init,
                           rho_init = rho_init, eps = eps, control = gaussian_control)
  t_fit <- fit_student_t_copula(g$U, dates = g$dates, init = student_init,
                                rho_init = g$rho_init, eps = eps,
                                nu_bounds = nu_bounds, control = student_control)
  k <- c(length(g$params), length(t_fit$params))
  ll <- c(g$logLik, t_fit$logLik)
  comparison <- data.frame(Model = c("Dynamic Gaussian", "Dynamic Student-t"),
                           logLik = ll, k = k, AIC = -2 * ll + 2 * k)
  comparison <- comparison[order(comparison$AIC), , drop = FALSE]
  rownames(comparison) <- NULL
  comparison$delta_AIC <- comparison$AIC - min(comparison$AIC)
  list(comparison = comparison, gaussian = g, student_t = t_fit,
       correlation_summary = list(gaussian = correlation_summary(g$rho_t),
                                  student_t = correlation_summary(t_fit$rho_t)))
}
