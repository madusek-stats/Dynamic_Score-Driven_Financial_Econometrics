# Return-based risk, wealth and recovery metrics. Inputs are log returns.
# Missing returns are rejected rather than silently treated as zero.
check_returns <- function(r) {
  stopifnot(is.numeric(r), length(r) > 0L, all(is.finite(r)))
  as.numeric(r)
}

VaR_ES_loss <- function(L, alpha = 0.05) {
  L <- check_returns(L)
  stopifnot(length(alpha) == 1L, alpha > 0, alpha < 1)
  q <- unname(stats::quantile(L, 1 - alpha, type = 7))
  list(VaR = q, ES = mean(L[L >= q]))
}
VaR_realised <- function(r, alpha = 0.05) VaR_ES_loss(-r, alpha)$VaR
ES_realised <- function(r, alpha = 0.05) VaR_ES_loss(-r, alpha)$ES

make_wealth <- function(r) exp(cumsum(check_returns(r)))
# Include initial wealth of one in the high-water mark.
make_drawdown <- function(r) {
  w <- make_wealth(r)
  w / cummax(c(1, w))[-1L] - 1
}
max_drawdown <- function(r) min(make_drawdown(r))
prop_time_underwater <- function(r) mean(make_drawdown(r) < 0)
CDaR_realised <- function(r, alpha = 0.10) VaR_ES_loss(-make_drawdown(r), alpha)$ES
omega_ratio <- function(r, threshold = 0) {
  r <- check_returns(r)
  den <- mean(pmax(threshold - r, 0))
  if (den <= 1e-12) return(NA_real_)
  mean(pmax(r - threshold, 0)) / den
}

# Episode lengths are observation counts, not calendar days.
# Unrecovered episodes are right-censored, not assigned fictitious recovery dates.
drawdown_episodes <- function(r, dates = NULL) {
  dd <- make_drawdown(r); n <- length(dd)
  if (!is.null(dates)) stopifnot(length(dates) == n, !anyNA(dates), !is.unsorted(dates, strictly = TRUE))
  runs <- rle(dd < 0); ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  out <- data.frame(start = integer(), trough = integer(), end = integer(),
                    depth = numeric(), duration = integer(), recovery = integer(), censored = logical())
  for (i in which(runs$values)) {
    s <- starts[i]; e <- ends[i]; tr <- s - 1L + which.min(dd[s:e]); recovered <- e < n
    out[nrow(out) + 1L, ] <- list(s, tr, if (recovered) e + 1L else NA_integer_,
                                 dd[tr], e - s + 1L, if (recovered) e + 1L - tr else NA_integer_, !recovered)
  }
  if (!is.null(dates)) {
    out$start_date <- dates[out$start]; out$trough_date <- dates[out$trough]; out$end_date <- dates[out$end]
  }
  out
}

# Penalise unrecovered episodes by full window length.
# Explicit configurable convention for optimisation only.
recovery_penalty <- function(r, censor_penalty = length(r)) {
  ep <- drawdown_episodes(r)
  if (!nrow(ep)) return(0)
  mean(ifelse(ep$censored, censor_penalty, ep$recovery))
}

portfolio_risk_summary <- function(r, freq = 252, alpha = 0.05, alpha_cdar = 0.10,
                                   threshold = 0) {
  r <- check_returns(r); ep <- drawdown_episodes(r)
  ve <- VaR_ES_loss(-r, alpha); vol <- stats::sd(r) * sqrt(freq)
  avg <- function(x) if (length(x) && any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
  c(n = length(r), ann_log_return = mean(r) * freq, ann_vol = vol,
    sharpe_zero_rf = if (is.finite(vol) && vol > 0) mean(r) * freq / vol else NA_real_,
    VaR = ve$VaR, ES = ve$ES, max_drawdown = max_drawdown(r), CDaR = CDaR_realised(r, alpha_cdar),
    omega = omega_ratio(r, threshold), time_underwater = prop_time_underwater(r),
    avg_duration = avg(ep$duration), avg_recovery = avg(ep$recovery),
    censored_episodes = sum(ep$censored), terminal_wealth = tail(make_wealth(r), 1L))
}

# User-supplied event windows, not embedded empirical dates.
crisis_window_analysis <- function(returns, dates, windows) {
  returns <- as.matrix(returns)
  stopifnot(nrow(returns) == length(dates), !anyNA(dates), !is.unsorted(dates, strictly = TRUE),
            all(c('start', 'end') %in% names(windows)), !is.null(colnames(returns)))
  rows <- list()
  for (i in seq_len(nrow(windows))) {
    idx <- which(dates >= windows$start[i] & dates <= windows$end[i])
    if (!length(idx)) next
    for (j in seq_len(ncol(returns))) {
      r <- check_returns(returns[idx, j])
      rows[[length(rows) + 1L]] <- data.frame(window = i, portfolio = colnames(returns)[j],
          n = length(r), cumulative_return = expm1(sum(r)), max_drawdown = max_drawdown(r))
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

# In-memory comparison
compare_portfolio_risk <- function(returns, ...) {
  returns <- as.matrix(returns)
  stopifnot(!is.null(colnames(returns)), ncol(returns) > 0L)
  out <- do.call(rbind, lapply(seq_len(ncol(returns)), function(j) portfolio_risk_summary(returns[, j], ...)))
  rownames(out) <- colnames(returns); out
}
