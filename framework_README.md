# Dynamic score-driven financial econometrics

Data-independent R functions adapted from a dissertation research implementation. No original market data, fitted results, private paths or automatic file I/O are included. The code is a research prototype, not a validated trading system.

## Modules and loading order

Load the files in this order in an R session (using the location of your downloaded files):

1. `marginals.R`
2. `copula.R`
3. `model_selection.R`
4. `risk_metrics.R`
5. `simulation.R`
6. `optimisation.R`
7. `backtesting.R`
8. `diagnostics.R`

Dependencies: Matrix (correlation projection), numDeriv and MASS (optional marginal inference), sandwich and lmtest (HAC tests), goftest (optional distributional tests). Base R supplies other functions. No dependencies are installed automatically.

## Interfaces

Supply three chronological, positive volatility-range series and a named, aligned matrix of log returns. Supply dates for all series, or align by row yourself. Each named marginal fit is returned by `fit_dcs_burr()`. Use `align_marginal_pits()` and `fit_gaussian_copula()` for the Gaussian model, or `compare_dynamic_copulas()` for model selection.

- `prepare_simulation(fits, copula_fit, returns, return_dates)` matches columns by name and dates explicitly; it calculates conditional means from prior returns only. `simulate_returns_t(state, t_idx)` generates predictive scenarios for a target row. Call `set.seed()` for reproducibility.
- `run_mc_allocation(state, constraints, w_bench, seed = ...)` accepts a named list of strategies, each with `min` and `max` named three-element weight vectors; `w_bench` is a named fully invested benchmark. Exclude an asset by giving it zero lower and upper bounds. The result carries target-date weights and VaR/ES. `refresh_daily = FALSE` carries forecasts between rebalances as in the dissertation. The default refreshes risk forecasts daily.
- `backtest_allocation(state$returns, forecasts, dates = state$dates)` evaluates aligned forecasts and weights against realised outcomes. Initial incomplete rows are dropped, not replaced with zeros.
- `run_return_recovery()` evaluates rolling historical allocations and return/recovery objectives, independently of the copula simulation. It uses only previous rows for decisions. `compare_portfolio_risk()` and `crisis_window_analysis()` summarise user-supplied portfolio log returns.
- `marginal_score_diagnostics()`, `copula_score_diagnostics()` and `rosenblatt_all_orders()` evaluate fitted-model diagnostics.

## Methodological notes and validation status

- Earlier marginal and copula fits use the whole supplied sample. The new allocation functions align filtered states and weights to their target observations, but this alone does **not** remove look-ahead bias from estimated parameters or sample-based initialisation. For honest out-of-sample forecasting, estimate exclusively on information available before each target (fixed training, rolling or expanding estimation), and construct forward states. That full re-estimation pipeline is not implemented here.
- Gaussian-copula simulation preserves the original dependence of positive volatility ranges and independent Gaussian return innovations. It does not impose the copula correlations on signed returns. Student-t is retained for model comparison, not used by the simulator.
- The model-selection module uses a full Student-t copula likelihood. The backtesting module uses a corrected upper-tail FZ0 loss and a HAC test of the joint calibration moment; it is not an exact replication of the original Acerbi–Szekely test. Formal Diebold–Mariano comparisons use forecasts of the same realised loss series, rather than different portfolios.
- The return/recovery module enforces its stated return and ES constraints rather than silently widening them. The Monte Carlo optimiser retains the original optional return-floor relaxation and reports it. Drawdowns include initial wealth, and unrecovered episodes are marked censored; the original full-window censoring penalty is retained only for optimisation.
- Weighted asset log returns are used as in the original research (an approximation to exact portfolio log returns); transaction costs, turnover costs, market impact and execution are not modelled. Empirical quantile-tail averaging follows the original convention, including ties.
- Original dissertation results are not reproduced or claimed by these modules. Code has been statically reviewed, not executed or unit tested.


