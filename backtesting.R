# Upper-tail loss convention: L = -return; alpha is exceedance probability.

hac_mean_test <- function(x, lag = 21L, min_obs = 200L) {
  stopifnot(is.numeric(x), all(is.finite(x)), lag >= 0, min_obs >= 2)
  n <- length(x)
  if (n < min_obs) return(list(n = n, mean = if (n) mean(x) else NA_real_, t = NA_real_, p = NA_real_))
  if (!requireNamespace('sandwich', quietly = TRUE) || !requireNamespace('lmtest', quietly = TRUE))
    stop('Install sandwich and lmtest for HAC inference.')
  if (stats::sd(x) == 0) return(list(n = n, mean = mean(x), t = NA_real_, p = NA_real_))
  fit <- stats::lm(x ~ 1)
  ct <- lmtest::coeftest(fit, vcov. = sandwich::NeweyWest(fit, lag = min(lag, n - 1L),
                                                         prewhite = FALSE, adjust = TRUE))
  list(n = n, mean = mean(x), t = unname(ct[1, 3]), p = unname(ct[1, 4]))
}

check_loss_forecasts <- function(L, q, e, alpha) {
  stopifnot(is.numeric(L), is.numeric(q), is.numeric(e), length(L) == length(q),
            length(q) == length(e), length(alpha) == 1L, alpha > 0, alpha < 1,
            all(is.finite(L)), all(is.finite(q)), all(is.finite(e)), all(e >= q))
}


AS_joint_test_loss <- function(L, VaR, ES, alpha = 0.05, nw_lag = 21L, min_obs = 200L) {
  check_loss_forecasts(L, VaR, ES, alpha)
  Z <- pmax(L - VaR, 0) / alpha - (ES - VaR)
  out <- hac_mean_test(Z, nw_lag, min_obs)
  out$meanZ <- out$mean; out
}


FZ0 <- function(L, q, e, alpha = 0.05) {
  check_loss_forecasts(L, q, e, alpha)
  if (any(e <= 0)) stop('FZ0 requires strictly positive loss ES.')
  pmax(L - q, 0) / (alpha * e) + q / e + log(e) - 1
}

dm_hac <- function(d, lag = 21L, min_obs = 200L) hac_mean_test(d, lag, min_obs)

# Standard comparative forecast test: both forecasts predict the SAME loss
# series. Negative mean score difference favours model A.
compare_risk_forecasts <- function(L, q_a, e_a, q_b, e_b, alpha = 0.05,
                                   nw_lag = 21L, min_obs = 200L) {
  dm_hac(FZ0(L, q_a, e_a, alpha) - FZ0(L, q_b, e_b, alpha), nw_lag, min_obs)
}

backtest_allocation <- function(returns, forecasts, dates = NULL,
                                nw_lag = 21L, min_obs = 200L, freq = 252) {
  returns <- as.matrix(returns); nm <- forecasts$assets
  stopifnot(identical(forecasts$timing, 'target_date'), setequal(colnames(returns), nm),
            !anyDuplicated(colnames(returns)), is.numeric(returns), all(is.finite(returns)))
  returns <- returns[, nm, drop = FALSE]
  if (!is.null(forecasts$dates)) {
    stopifnot(!is.null(dates), identical(as.character(dates), as.character(forecasts$dates)))
  }
  n <- nrow(returns); ports <- names(forecasts$weights)
  stopifnot(nrow(forecasts$VaR) == n, nrow(forecasts$ES) == n,
            identical(colnames(forecasts$VaR), ports), identical(colnames(forecasts$ES), ports))
  realised <- vapply(forecasts$weights, function(w) {
    stopifnot(identical(dim(w), dim(returns)), identical(colnames(w), nm))
    rowSums(w * returns)
  }, numeric(n))
  realised <- matrix(realised, nrow = n, ncol = length(ports), dimnames = list(NULL, ports))
  keep <- which(rowSums(!is.finite(realised) | !is.finite(forecasts$VaR) | !is.finite(forecasts$ES)) == 0)
  if (!length(keep)) stop('No complete evaluation rows.')
  if (length(keep) > 1L && any(diff(keep) != 1L)) stop('Internal gaps in evaluation rows.')
  r <- realised[keep, , drop = FALSE]; L <- -r
  q <- forecasts$VaR[keep, , drop = FALSE]; e <- forecasts$ES[keep, , drop = FALSE]
  alpha <- forecasts$alpha
  calibration <- lapply(ports, function(p) AS_joint_test_loss(L[, p], q[, p], e[, p], alpha, nw_lag, min_obs))
  names(calibration) <- ports
  scores <- matrix(NA_real_, nrow(r), length(ports), dimnames = list(NULL, ports))
  for (p in ports) if (all(e[, p] > 0)) scores[, p] <- FZ0(L[, p], q[, p], e[, p], alpha)
  # Scores of different portfolio outcomes are descriptive
  score_summary <- vapply(ports, function(p) if (all(is.finite(scores[, p]))) mean(scores[, p]) else NA_real_, numeric(1))
  list(dates = if (is.null(dates)) NULL else dates[keep], rows = keep,
       returns = r, losses = L, VaR = q, ES = e, calibration = calibration,
       exceedance_rate = colMeans(L > q), mean_forecast_ES = colMeans(e),
       scores = scores, descriptive_mean_scores = score_summary,
       performance = compare_portfolio_risk(r, freq = freq, alpha = alpha))
}
