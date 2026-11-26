# Distributional and score-residual diagnostics, with no plots or I/O.
# Source marginals.R and copula.R for score functions and fitted objects.
# Optional Cramer-von Mises and Anderson-Darling tests require goftest.

# Conditional Gaussian AR(1) likelihood ratio on matched observations.
berkowitz_lr <- function(z) {
  stopifnot(is.numeric(z), all(is.finite(z)))
  n <- length(z)
  if (n < 50L) return(list(LR = NA_real_, p = NA_real_))
  z1 <- z[-1L]; z0 <- z[-n]
  fit <- stats::lm(z1 ~ z0)
  phi <- unname(stats::coef(fit)[2]); intercept <- unname(stats::coef(fit)[1])
  sigma <- sqrt(mean(stats::residuals(fit)^2))
  if (!is.finite(phi) || !is.finite(sigma) || sigma <= 0) return(list(LR = NA_real_, p = NA_real_))
  ll0 <- sum(stats::dnorm(z1, log = TRUE))
  ll1 <- sum(stats::dnorm(stats::residuals(fit), sd = sigma, log = TRUE))
  LR <- max(0, 2 * (ll1 - ll0))
  list(LR = LR, p = stats::pchisq(LR, df = 3, lower.tail = FALSE), phi = phi,
       mu = if (abs(phi) < 1) intercept / (1 - phi) else NA_real_, sigma = sigma)
}

score_acf_diagnostics <- function(x, lag = 20L) {
  stopifnot(is.numeric(x), all(is.finite(x)), length(x) > lag + 1L, lag >= 1)
  if (stats::sd(x) == 0) return(list(n = length(x), mean = mean(x), acf = NULL, ljung_box = NULL))
  list(n = length(x), mean = mean(x), acf = stats::acf(x, lag.max = lag, plot = FALSE),
       ljung_box = stats::Box.test(x, lag = lag, type = 'Ljung-Box'))
}

# P-values are conventional diagnostic approximations
pit_diagnostics <- function(u, alpha_tail = 0.05, lag = 20L, eps = 1e-6) {
  stopifnot(is.numeric(u), all(is.finite(u)), all(u >= 0 & u <= 1),
            length(u) > lag + 1L, alpha_tail > 0, alpha_tail < 0.5, eps > 0, eps < 0.5)
  z <- stats::qnorm(pmin(pmax(u, eps), 1 - eps)); n <- length(u)
  extra <- if (requireNamespace('goftest', quietly = TRUE)) {
    list(cvm = goftest::cvm.test(u, null = 'punif'), ad = goftest::ad.test(u, null = 'punif'))
  } else NULL
  list(n = n, ks = stats::ks.test(u, 'punif'), extra = extra,
       lower_tail = stats::binom.test(sum(u < alpha_tail), n, p = alpha_tail),
       upper_tail = stats::binom.test(sum(u > 1 - alpha_tail), n, p = alpha_tail),
       normal_scores = score_acf_diagnostics(z, lag), squared_scores = score_acf_diagnostics(z^2, lag),
       berkowitz = berkowitz_lr(z), qq = data.frame(theoretical = stats::ppoints(n), empirical = sort(u)))
}

marginal_score_diagnostics <- function(fit, lag = 20L) {
  y <- fit$y; c <- fit$c; k <- fit$k_t; lambda <- fit$lambda_t
  stopifnot(length(y) == length(k), length(k) == length(lambda))
  scale <- score_f(y, c, k, lambda) / I_fisher_f(c, k)
  tail <- score_g(y, c, k, lambda) / I_fisher_g()
  list(dates = fit$dates, scores = cbind(scale = scale, tail = tail),
       scale = score_acf_diagnostics(scale, lag), tail = score_acf_diagnostics(tail, lag),
       pit = pit_diagnostics(fit$pit$u, lag = lag))
}

# Match the half-gradient / chain-rule convention in copula.R. Diagnostics
# describe the actual filter score, not the derivative of the nearPD projection.
copula_score_diagnostics <- function(fit, lag = 20L) {
  U <- fit$U; R <- fit$R_t; raw <- fit$rho_raw
  stopifnot(is.matrix(U), ncol(U) == 3L, identical(dim(R), c(3L, 3L, nrow(U))),
            identical(dim(raw), c(nrow(U), 3L)))
  X <- stats::qnorm(U)
  S <- t(vapply(seq_len(nrow(U)), function(i) score_rho_gaussian(X[i, ], R[, , i]) *
                      (1 - raw[i, ]^2), numeric(3)))
  colnames(S) <- colnames(raw)
  tests <- lapply(seq_len(3L), function(j) score_acf_diagnostics(S[, j], lag))
  names(tests) <- colnames(S)
  list(dates = fit$dates, scores = S, tests = tests)
}

# Rosenblatt tests for every permutation, using the same fitted Gaussian model.
rosenblatt_all_orders <- function(fit, lag = 20L) {
  orders <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  out <- lapply(orders, function(j) {
    U <- fit$U[, j, drop = FALSE]; R <- fit$R_t[j, j, , drop = FALSE]
    V <- rosenblatt_gaussian_tv(U, R)
    list(order = colnames(U), transforms = V, diagnostics = rosenblatt_diagnostics(V, lag))
  })
  names(out) <- vapply(orders, function(j) paste(fit$assets[j], collapse = '__'), character(1))
  out
}
