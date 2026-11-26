
candidate_grid <- function(w_min, w_max, step = 0.01) {
  stopifnot(length(w_min) == 3L, length(w_max) == 3L,
            all(is.finite(c(w_min, w_max))), all(w_min >= 0), all(w_max <= 1),
            all(w_min <= w_max), sum(w_min) <= 1 + 1e-10, sum(w_max) >= 1 - 1e-10,
            step > 0, step <= 1)
  a <- expand.grid(w1 = seq(w_min[1], w_max[1], by = step),
                   w2 = seq(w_min[2], w_max[2], by = step))
  W <- cbind(as.matrix(a), w3 = 1 - a$w1 - a$w2)
  keep <- W[, 3] >= w_min[3] - 1e-10 & W[, 3] <= w_max[3] + 1e-10
  W <- W[keep, , drop = FALSE]; W[abs(W) < 1e-12] <- 0
  if (!nrow(W)) stop('No grid-feasible portfolios; change step or constraints.')
  colnames(W) <- names(w_min); W
}

validate_weights <- function(w, assets) {
  stopifnot(is.numeric(w), length(w) == length(assets), !anyDuplicated(names(w)),
            setequal(names(w), assets), all(is.finite(w)), all(w >= 0), abs(sum(w) - 1) < 1e-8)
  w[assets]
}

optim_es_3asset <- function(r_sim, alpha = 0.05, w_min, w_max, w_bench,
                            step = 0.01, ret_tol = 0, fallback_relax = TRUE) {
  r_sim <- as.matrix(r_sim); nm <- colnames(r_sim)
  stopifnot(ncol(r_sim) == 3L, !is.null(nm), !anyDuplicated(nm), all(is.finite(r_sim)),
            setequal(names(w_min), nm), setequal(names(w_max), nm))
  w_bench <- validate_weights(w_bench, nm)
  W <- candidate_grid(w_min[nm], w_max[nm], step)
  mu <- drop(W %*% colMeans(r_sim)); mu_bench <- sum(w_bench * colMeans(r_sim))
  idx <- which(mu >= mu_bench + ret_tol - 1e-12); relaxed <- FALSE
  if (!length(idx) && fallback_relax && ret_tol > -1e-4) {
    ret_tol <- -1e-4; idx <- which(mu >= mu_bench + ret_tol - 1e-12); relaxed <- TRUE
  }
  if (!length(idx)) stop('No portfolio meets the return floor.')
  risk <- vapply(idx, function(i) unlist(VaR_ES_loss(-drop(r_sim %*% W[i, ]), alpha)), numeric(2))
  j <- which.min(risk[2, ]); i <- idx[j]
  list(w = W[i, ], VaR = risk[1, j], ES = risk[2, j], mu = mu[i],
       mu_bench = mu_bench, return_tolerance = ret_tol, relaxed = relaxed)
}

# Every strategy supplies list(min = named vector, max = named vector).
# Exclude an asset by setting its min and max to zero.
# Target-date weights and forecasts are aligned with same-row realised returns.
run_mc_allocation <- function(state, constraints, w_bench, alpha = 0.05,
                              B = 20000L, rebalance_every = 21L, step = 0.01,
                              refresh_daily = TRUE, seed = NULL, fallback_relax = TRUE) {
  nm <- state$assets; n <- nrow(state$mu)
  w_bench <- validate_weights(w_bench, nm)
  stopifnot(length(constraints) > 0L, !is.null(names(constraints)),
            !anyDuplicated(names(constraints)), !('benchmark' %in% names(constraints)),
            rebalance_every >= 1, rebalance_every == as.integer(rebalance_every))
  if (!is.null(seed)) set.seed(seed)
  ports <- c(names(constraints), 'benchmark')
  W <- setNames(lapply(ports, function(p) matrix(NA_real_, n, 3L, dimnames = list(NULL, nm))), ports)
  q <- e <- m <- matrix(NA_real_, n, length(ports), dimnames = list(NULL, ports))
  relaxed <- matrix(FALSE, n, length(ports), dimnames = list(NULL, ports))
  current <- NULL; last_rebalance <- -Inf
  for (t in seq_len(n)) {
    if (any(!is.finite(state$mu[t, ]))) {
      current <- NULL # Restart rather than carrying stale risk forecasts across gaps.
      next
    }
    rebalance <- is.null(current) || t - last_rebalance >= rebalance_every
    if (rebalance || refresh_daily) {
      sim <- simulate_returns_t(state, t, B)
      if (rebalance) {
        sols <- lapply(constraints, function(z) optim_es_3asset(sim, alpha, z$min, z$max,
                          w_bench, step, fallback_relax = fallback_relax))
        current <- c(lapply(sols, function(s) s$w), list(benchmark = w_bench))
        relaxed[t, names(constraints)] <- vapply(sols, function(s) s$relaxed, logical(1))
        last_rebalance <- t
      }
      for (p in ports) {
        rp <- drop(sim %*% current[[p]]); ve <- VaR_ES_loss(-rp, alpha)
        q[t, p] <- ve$VaR; e[t, p] <- ve$ES; m[t, p] <- mean(rp)
      }
    } else {
      q[t, ] <- q[t-1L, ]; e[t, ] <- e[t-1L, ]; m[t, ] <- m[t-1L, ]
    }
    for (p in ports) W[[p]][t, ] <- current[[p]]
  }
  list(dates = state$dates, assets = nm, weights = W, VaR = q, ES = e, mu = m,
       alpha = alpha, relaxed = relaxed, timing = 'target_date', refresh_daily = refresh_daily)
}

# Historical return/recovery strategy, adapted from Portfolio_Rule.R.
# Strict constraints: do not silently relax ES or return floors.
select_portfolios <- function(R_window, W_grid, w_bench, alpha = 0.05, alpha_cdar = 0.10,
                              ES_tolerance_to_benchmark = 0.75, return_tie_tol = 0.00002) {
  nm <- colnames(R_window); w_bench <- validate_weights(w_bench, nm)
  stopifnot(identical(colnames(W_grid), nm), all(is.finite(R_window)),
            ES_tolerance_to_benchmark >= 0, ES_tolerance_to_benchmark <= 1, return_tie_tol >= 0)
  stats <- t(vapply(seq_len(nrow(W_grid)), function(i) {
    r <- drop(R_window %*% W_grid[i, ])
    c(mu = mean(r), ES = ES_realised(r, alpha), CDaR = CDaR_realised(r, alpha_cdar),
      Recovery = recovery_penalty(r), Omega = omega_ratio(r))
  }, numeric(5)))
  mu_b <- mean(drop(R_window %*% w_bench)); es_b <- ES_realised(drop(R_window %*% w_bench), alpha)
  ix <- which(stats[, 'mu'] >= mu_b - 1e-12)
  if (!length(ix)) stop('No candidate meets benchmark return floor.')
  es_i <- ix[which.min(stats[ix, 'ES'])]; es_min <- stats[es_i, 'ES']
  threshold <- es_min + ES_tolerance_to_benchmark * (es_b - es_min)
  feasible <- ix[stats[ix, 'ES'] <= threshold + 1e-12]
  if (!length(feasible)) stop('No ES-feasible portfolio; benchmark ES is below grid minimum.')
  near <- feasible[stats[feasible, 'mu'] >= max(stats[feasible, 'mu']) - return_tie_tol]
  final <- near[order(stats[near, 'Recovery'], stats[near, 'CDaR'])][1]
  rc <- ix[which.max(stats[ix, 'mu'] / pmax(stats[ix, 'CDaR'], 1e-12))]
  om <- ix[is.finite(stats[ix, 'Omega'])]
  if (!length(om)) stop('Omega ratio undefined for all feasible candidates.')
  om <- om[which.max(stats[om, 'Omega'])]
  list(weights = list(ES_only = W_grid[es_i, ], ReturnCDaR = W_grid[rc, ],
                      Omega = W_grid[om, ], FinalRule = W_grid[final, ]),
       diagnostics = c(mu_bench = mu_b, ES_bench = es_b, ES_min = es_min,
                       ES_threshold = threshold, n_feasible = length(feasible), n_tie = length(near)),
       candidates = stats)
}

# Rolling allocation uses rows strictly before each target date; no double lag.
run_return_recovery <- function(returns, w_min, w_max, w_bench, dates = NULL,
                                lookback = 252L, rebalance_every = 21L, step = 0.01, ...) {
  returns <- as.matrix(returns); nm <- colnames(returns); n <- nrow(returns)
  stopifnot(length(nm) == 3L, all(is.finite(returns)), n > lookback,
            lookback >= 2, rebalance_every >= 1,
            setequal(names(w_min), nm), setequal(names(w_max), nm))
  if (!is.null(dates)) stopifnot(length(dates) == n, !anyNA(dates), !is.unsorted(dates, strictly = TRUE))
  w_bench <- validate_weights(w_bench, nm); grid <- candidate_grid(w_min[nm], w_max[nm], step)
  ports <- c('ES_only', 'ReturnCDaR', 'Omega', 'FinalRule', 'benchmark')
  W <- setNames(lapply(ports, function(p) matrix(NA_real_, n, 3L, dimnames = list(NULL, nm))), ports)
  diagnostics <- list()
  for (t in seq.int(lookback + 1L, n)) {
    if ((t - lookback - 1L) %% rebalance_every == 0L) {
      fit <- select_portfolios(returns[(t-lookback):(t-1L), , drop = FALSE], grid, w_bench, ...)
      current <- c(fit$weights, list(benchmark = w_bench)); diagnostics[[as.character(t)]] <- fit$diagnostics
    }
    for (p in ports) W[[p]][t, ] <- current[[p]]
  }
  realised <- sapply(W, function(w) rowSums(w * returns))
  list(dates = dates, assets = nm, weights = W, returns = realised, diagnostics = diagnostics,
       timing = 'target_date')
}
