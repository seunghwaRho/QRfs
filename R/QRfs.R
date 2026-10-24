# QRfs: fixed-smoothing (fixed-b and fixed-K) HAR inference for quantile regression.
# Package source, built from library/qrfs_functions_reduced.R of the review repository by
# agentClaude/scripts/2026-10-19/build_package_claude.py: identical except that the Adrian data
# is the dataset 'adrian' (data/adrian.rda) rather than a function.
# The paper: Chen, K., Galvao, A. F., Rho, S., Vogelsang, T. J. & Yoon, J.,
#   "Fixed-smoothing Uniform Inference for Quantile Regression". Section, equation
#   and algorithm numbers below refer to it. qrhac is the conventional counterpart of the stacked
#   test: the small-b kernel HAC of Galvao & Yoon (2024) or the large-K OS HAR, chi-square critical value.
# Dependencies, always called with their namespace: quantreg (rq and its summary), Matrix (bdiag),
# mvtnorm -- or pbivnorm in its place -- for the bivariate normal probabilities of .phi2. As a
# package: Imports quantreg, Matrix, mvtnorm; Suggests pbivnorm. No code runs at the top level.

# QR scores z_t(tau) = x_t (tau - 1{u_t(tau) <= 0}) (paper, Section 2). Residuals that are zero only up to
# rounding, of order 1e-16 and of arbitrary sign, would be classified by that noise under a bare <= 0, so a
# residual within 1e-8 * sd(residuals) of zero counts as non-positive.
compute_score_matrix <- function(x_mat, resid_vec, tau) {
  x_mat * (tau - (resid_vec <= 1e-8 * stats::sd(resid_vec)))
}

# Fit the quantile grid; return model components, inverse Hessians and cached scores.
# arg_hinv_method: quantreg 'ker' (Powell's kernel sandwich, the paper's eq. (4); default) or 'nid'
#   (the Hendricks-Koenker sparsity sandwich).
# arg_rq_method: the solver of quantreg::rq -- 'br' (Barrodale-Roberts simplex, quantreg's default),
#   'fn' (Frisch-Newton interior point, the solver of Adrian et al.'s own code). The
#   nid Hessian's refits inside summary.rq use the same solver.
fit_qr_grid <- function(arg_formula, arg_data, arg_tau_grid,
                         arg_hinv_method = "ker", arg_rq_method = "br") {
  if (!is.character(arg_hinv_method) || length(arg_hinv_method) != 1L ||
      is.na(arg_hinv_method) || !(arg_hinv_method %in% c("nid", "ker"))) {
    stop("arg_hinv_method must be 'nid' or 'ker'.", call. = FALSE)
  }
  arg_rq_method <- match.arg(arg_rq_method, c("br", "fn"))

  fit <- quantreg::rq(formula = arg_formula, data = arg_data,
                      tau = arg_tau_grid, method = arg_rq_method, model = TRUE)

  model_frame <- stats::model.frame(fit)
  x_mat <- stats::model.matrix(stats::terms(fit), data = model_frame,
                                contrasts.arg = fit$contrasts)
  coef_mat <- as.matrix(stats::coef(fit))

  summaries <- summary(fit, se = arg_hinv_method, covariance = TRUE)
  if (length(fit$tau) == 1L) summaries <- list(summaries)
  Hinv_list <- lapply(summaries, function(s) s$Hinv)

  residuals <- as.matrix(fit$residuals)
  m <- length(fit$tau); p <- ncol(x_mat)
  score_mat <- matrix(0, nrow(x_mat), m * p)
  for (j in seq_len(m)) {
    ii <- ((j - 1L) * p + 1L):(j * p)
    score_mat[, ii] <- compute_score_matrix(x_mat, residuals[, j], fit$tau[j])
  }

  list(
    fit = fit,
    Hinv_list = Hinv_list,
    coef_mat = coef_mat,
    residuals = residuals,
    x_mat = x_mat,
    score_mat = score_mat,
    nobs = nrow(model_frame),
    p = p,
    m = m,
    arg_hinv_method = arg_hinv_method
  )
}

# Build the T x K DCT-II cosine basis phi_k(r) = sqrt(2) cos(pi k (r - 1/(2T))) at r = t/T, as in
# Section 5 of the paper. It is orthonormal in sample for K <= T - 1 and aliased above it.
make_os_basis <- function(arg_T, arg_K) {
  T_n <- as.integer(arg_T); K <- as.integer(arg_K)
  if (length(T_n) != 1L || is.na(T_n) || T_n < 2L ||
      length(K) != 1L || is.na(K) || K < 1L)
    stop("make_os_basis requires scalar integers T >= 2 and K >= 1.", call. = FALSE)
  if (K > T_n - 1L)
    stop(sprintf(paste0("make_os_basis: K=%d exceeds the orthonormality limit %d for T=%d -- ",
      "the cosine basis is not orthonormal above it (aliasing/degeneracy)."), K, T_n - 1L, T_n), call. = FALSE)
  sqrt(2) * cos(pi * outer((1:T_n) - 0.5, 1:K) / T_n)
}

# Wald statistic for a deviation and covariance matrix.
# A variance vector instead gives separate scalar Wald statistics.
compute_wald <- function(arg_deviation, arg_covariance) {
  if (is.null(dim(arg_covariance))) {
    if (any(!is.finite(arg_covariance)) || any(arg_covariance == 0))
      stop("The simulated Wald covariance is non-finite or singular.", call. = FALSE)
    return(as.numeric(arg_deviation)^2 / arg_covariance)
  }
  covariance <- (arg_covariance + t(arg_covariance)) / 2
  as.numeric(crossprod(arg_deviation, solve(covariance, arg_deviation)))
}

# Sum-scale score covariance: Bartlett with numeric arg_M, or OS with arg_K.
# arg_phi_mat optionally supplies a precomputed OS basis.
.compute_Jhat_from_scores <- function(arg_scores, arg_har_method, arg_M,
                                      arg_K, arg_phi_mat) {
  if (arg_har_method == "ker-har") {
    if (arg_M == 1L) return(crossprod(arg_scores))
    n <- nrow(arg_scores)
    # Moving sums, including partial windows at both ends, give the Bartlett weights.
    cumulative_scores <- rbind(0, apply(arg_scores, 2L, cumsum))
    r <- seq_len(n + arg_M - 1L)
    window_sums <- cumulative_scores[pmin(r, n) + 1L, , drop = FALSE] -
      cumulative_scores[pmax(r - arg_M, 0L) + 1L, , drop = FALSE]
    return(crossprod(window_sums) / arg_M)
  }
  if (arg_har_method == "os-har") {
    Phi <- arg_phi_mat
    if (is.null(Phi)) Phi <- make_os_basis(nrow(arg_scores), arg_K)
    projections <- crossprod(arg_scores, Phi)
    return(tcrossprod(projections) / ncol(Phi))
  }
  stop("arg_har_method must be 'ker-har' or 'os-har'.", call. = FALSE)
}

# Stacked score covariance, ordered by coefficient within quantile. Sum scale: J = T * Lambda-hat.
compute_Jhat_stacked <- function(arg_qr_fit, arg_har_method, arg_M, arg_K) {
  .compute_Jhat_from_scores(arg_qr_fit$score_mat, arg_har_method, arg_M,
                            arg_K, NULL)
}

# Select integer K: numeric arg_type, 'simple', 'test_optimal',
# 'mse_optimal_var', or 'mse_optimal_ar'. Test-optimal selection needs stacked arg_R.
# arg_K_min is the floor the automatic rules apply: q + 4 for a fixed-K test (the F(q, K - q + 1) limit of
# equation (19) has a finite variance only then), 5 for a shape test, m*p + 4 for the nuisance covariance,
# and q for qrhac, whose chi-square limit needs only an invertible estimate (paper, footnote 2).
# arg_R, arg_alpha and arg_K_min have no default; only 'test_optimal' reads arg_R and arg_alpha.
compute_K <- function(arg_qr_fit, arg_type, arg_R, arg_alpha, arg_K_min,
                      arg_threshold_method = "theoretical",
                      arg_var_fit = "gaussian",
                      arg_var_stationary = "andrews_monahan",
                      arg_kappa = NULL) {
  T_n   <- arg_qr_fit$nobs
  K_max <- T_n - 1L
  if (is.numeric(arg_type)) {
    if (arg_K_min > K_max)
      warning(sprintf("Required minimum K (%d) exceeds the maximum number of basis functions (%d).",
                      arg_K_min, K_max))
    return(as.integer(arg_type))
  }
  K_min <- arg_K_min
  if (K_min > K_max)
    stop(sprintf("Required minimum K (%d) exceeds the maximum number of basis functions (%d).",
                 K_min, K_max), call. = FALSE)
  switch(arg_type,
    "simple" = as.integer(min(K_max, max(K_min, floor(2 * T_n^(1/3))))),
    "test_optimal" = compute_test_optimal_K(
      arg_qr_fit, arg_R, arg_alpha = arg_alpha,
      arg_kappa = arg_kappa, arg_K_min = arg_K_min,
      arg_threshold_method = arg_threshold_method, arg_var_fit = arg_var_fit, arg_var_stationary = arg_var_stationary),
    "mse_optimal_var" = compute_mse_optimal_K_var(
      arg_qr_fit, arg_K_min = arg_K_min,
      arg_threshold_method = arg_threshold_method, arg_var_fit = arg_var_fit, arg_var_stationary = arg_var_stationary),
    "mse_optimal_ar" = compute_mse_optimal_K_ar(
      arg_qr_fit, arg_K_min = arg_K_min,
      arg_threshold_method = arg_threshold_method, arg_var_fit = arg_var_fit,
      arg_var_stationary = arg_var_stationary),
    stop("arg_type must be numeric or one of 'test_optimal', 'mse_optimal_var', ",
         "'mse_optimal_ar', 'simple'.", call. = FALSE)
  )
}

# Bivariate standard-normal CDF at (z1, z2), with correlation rho.
.phi2 <- function(z1, z2, rho) {
  if (rho >= 1) return(stats::pnorm(min(z1, z2)))
  if (rho <= -1) return(max(0, stats::pnorm(z1) + stats::pnorm(z2) - 1))
  if (requireNamespace("mvtnorm", quietly = TRUE)) {
    as.numeric(mvtnorm::pmvnorm(upper = c(z1, z2),
                                corr = matrix(c(1, rho, rho, 1), 2, 2)))
  } else if (requireNamespace("pbivnorm", quietly = TRUE)) {
    as.numeric(pbivnorm::pbivnorm(z1, z2, rho))
  } else {
    stop("Neither mvtnorm nor pbivnorm is available.", call. = FALSE)
  }
}

# Bivariate-normal approximation to the quantile-score transition matrix.
# The paper's alternative estimator of A (Section 4.5), adapting Galvao & Yoon (2024).
# arg_threshold_method: 'theoretical' normal or 'empirical' standardized-residual quantiles.
.estimate_A_r <- function(arg_qr_fit, arg_threshold_method) {
  m <- arg_qr_fit$m
  list_tau <- arg_qr_fit$fit$tau
  E <- as.matrix(arg_qr_fit$residuals)
  nobs <- arg_qr_fit$nobs
  rho_hat <- stats::cor(E[-1L, , drop = FALSE], E[-nobs, , drop = FALSE])

  q_tau <- if (identical(arg_threshold_method, "theoretical")) {
    stats::qnorm(list_tau)
  } else {
    E_std <- scale(E)
    vapply(seq_len(m), function(j) {
      as.numeric(stats::quantile(E_std[, j], probs = list_tau[j],
                                 names = FALSE, type = 7))
    }, numeric(1L))
  }

  # Lag-one score covariance under the bivariate-normal approximation.
  E_lag1 <- matrix(0, m, m)
  for (i in 1:m) for (j in 1:m) {
    E_lag1[i, j] <- .phi2(q_tau[i], q_tau[j], rho_hat[i, j]) -
      list_tau[i] * list_tau[j]
  }

  # Contemporaneous quantile-score covariance.
  E_lag0 <- outer(list_tau, list_tau, pmin) - outer(list_tau, list_tau, "*")

  E_lag1 %*% solve(E_lag0)
}

# Estimate the regressor VAR(1) transition matrix A_x.
.estimate_A_x <- function(arg_qr_fit) {
  X <- arg_qr_fit$x_mat
  T_n <- nrow(X)
  X_t   <- X[-1, , drop = FALSE]
  X_lag <- X[-T_n, , drop = FALSE]
  # QR least squares avoids forming the lagged-regressor cross-product matrix.
  t(stats::lm.fit(X_lag, X_t, singular.ok = FALSE)$coefficients)
}

# The stationarity bound of the working models: arg_var_stationary as a number -- 'andrews_monahan'
# 0.97 (Andrews & Monahan 1992), 'sul_phillips_choi' 1 - 1/sqrt(T) (Sul, Phillips & Choi 2005), a number
# in (0, 1) as given -- or NA for 'unrestricted'. The VAR's spectral radius and each AR(1) coefficient
# are held within it; unrestricted, a fit that is not stationary stops instead.
.qrfs_stationary_bound <- function(arg_var_stationary, arg_nobs) {
  if (identical(arg_var_stationary, "unrestricted")) return(NA_real_)
  bound <- if (identical(arg_var_stationary, "andrews_monahan")) 0.97 else
    if (identical(arg_var_stationary, "sul_phillips_choi")) 1 - 1 / sqrt(arg_nobs) else arg_var_stationary
  if (!is.numeric(bound) || length(bound) != 1L || !is.finite(bound) || bound <= 0 || bound >= 1)
    stop('arg_var_stationary must be "andrews_monahan", "sul_phillips_choi", ',
         '"unrestricted", or a finite scalar in (0, 1).', call. = FALSE)
  bound
}

# Estimate the VAR(1) working model of the stacked scores: A and the innovation covariance Sigma_eps.
# arg_var_fit 'gaussian': A = A_r %x% A_x, Galvao & Yoon's separation of the indicator's and the
# regressors' dependence; 'ols': the VAR(1) fitted to the scores by least squares (Sun 2014).
# arg_var_stationary caps the spectral radius of A (see .qrfs_stationary_bound).
.estimate_stacked_var <- function(arg_qr_fit,
                                  arg_threshold_method,
                                  arg_var_stationary,
                                  arg_var_fit) {
  T_n <- arg_qr_fit$nobs
  s_stack <- arg_qr_fit$score_mat

  if (identical(arg_var_fit, "ols")) {
    # The VAR(1) of the mp-dimensional score by least squares, which needs T - 1 > m p observations.
    mp <- ncol(s_stack)
    if (T_n - 1L <= mp)
      stop("arg_var_fit = \"ols\" fits a VAR(1) with ", mp, " x ", mp, " coefficients and needs ",
           "T - 1 > m * p; here T = ", T_n, " and m * p = ", mp, ". Use \"gaussian\".", call. = FALSE)
    lagged <- s_stack[-T_n, , drop = FALSE]
    # QR least squares avoids forming the lagged-score cross-product matrix.
    A_raw <- t(stats::lm.fit(lagged, s_stack[-1L, , drop = FALSE],
                            singular.ok = FALSE)$coefficients)
  } else {
    A_r <- .estimate_A_r(arg_qr_fit, arg_threshold_method = arg_threshold_method)
    A_x <- .estimate_A_x(arg_qr_fit)
    A_raw <- kronecker(A_r, A_x)
  }

  bound <- .qrfs_stationary_bound(arg_var_stationary, T_n)
  radius_raw <- max(Mod(eigen(A_raw, only.values = TRUE)$values))
  # Unrestricted, the fit must be stationary on its own; a bound rescales A only when it binds.
  if (is.na(bound) && radius_raw >= 1)
    stop(sprintf("The fitted VAR(1) of the scores is not stationary (spectral radius %.3f) and ", radius_raw),
         "arg_var_stationary = \"unrestricted\" imposes no bound; choose \"andrews_monahan\" (0.97), ",
         "\"sul_phillips_choi\" (1 - 1/sqrt(T)), or a number in (0, 1).", call. = FALSE)
  scale <- if (!is.na(bound) && radius_raw > bound) bound / radius_raw else 1
  A <- scale * A_raw

  eps <- s_stack[-1L, , drop = FALSE] - s_stack[-T_n, , drop = FALSE] %*% t(A)
  Sigma_eps <- crossprod(eps) / (T_n - 1L)
  list(A = A, Sigma_eps = Sigma_eps)
}

# VAR long-run covariance and OS bias; arg_bias_constant scales the bias.
# The bias polynomial is that of the paper, Section 4.5, from Hwang & Valdes (2025).
.compute_os_moments <- function(arg_A, arg_Sigma_eps, arg_bias_constant) {
  A <- arg_A
  E <- arg_Sigma_eps
  mp <- nrow(A)
  inv_IA <- solve(diag(mp) - A)
  Omega_hat <- inv_IA %*% E %*% t(inv_IA)
  Omega_hat <- (Omega_hat + t(Omega_hat)) / 2

  At   <- t(A)
  A2   <- A %*% A
  At2  <- At %*% At
  inner <- A %*% E +
           A2 %*% E %*% At +
           A2 %*% E -
           6 * A %*% E %*% At +
           E %*% At2 +
           A %*% E %*% At2 +
           E %*% At
  inv_IA3  <- inv_IA %*% inv_IA %*% inv_IA
  # Second-order lag moment: sum_h h^2 Gamma(h).
  S2_hat <- inv_IA3 %*% inner %*% t(inv_IA3)
  S2_hat <- (S2_hat + t(S2_hat)) / 2
  B_tilde <- -arg_bias_constant * S2_hat
  list(Omega_hat = Omega_hat, B_tilde = B_tilde)
}

# VAR plug-in MSE-optimal K, with the caller's minimum arg_K_min and a basis-specific cap.
# MSE criterion of Hwang & Valdes (2025, Section 6) with an identity weight.
compute_mse_optimal_K_var <- function(arg_qr_fit, arg_K_min,
                                      arg_threshold_method,
                                      arg_var_stationary,
                                      arg_var_fit) {
  T_n <- arg_qr_fit$nobs
  K_max <- T_n - 1L
  bias_constant <- pi^2 / 6   # Hwang & Valdes (2025), DCT-II
  var_approx <- .estimate_stacked_var(arg_qr_fit,
                               arg_threshold_method = arg_threshold_method, arg_var_fit = arg_var_fit,
                               arg_var_stationary = arg_var_stationary)
  os_moments <- .compute_os_moments(var_approx$A, var_approx$Sigma_eps,
                                     bias_constant)
  Omega_hat <- os_moments$Omega_hat
  B_tilde <- os_moments$B_tilde

  # AMSE = V / K + squared bias * (K / T)^4.
  V_scalar <- sum(diag(Omega_hat))^2 + sum(Omega_hat^2)
  B_scalar <- sum(B_tilde^2)

  K_min <- arg_K_min
  K_tilde <- (V_scalar / (4 * B_scalar))^(1/5) * T_n^(4/5)
  as.integer(min(max(K_min, ceiling(K_tilde)), K_max))
}

# AR(1) coefficient and innovation variance of each stacked score column, for the AR-based MSE rules
# (Galvao & Yoon 2024, Section 5). arg_var_fit 'gaussian': the indicator's coefficient from the
# bivariate normal at the threshold with the residuals' lag-one correlation -- GY's eq. (2) with
# Phi^-1(tau) ('theoretical'), or their alternative estimator with the standardised residuals'
# tau-quantile ('empirical') -- times the regressor's own AR(1) coefficient (1 for the intercept);
# 'ols': GY's natural estimator, the least-squares AR(1) coefficient of the score itself.
# Each coefficient is held within the stationarity bound of arg_var_stationary, as the VAR's radius is.
.estimate_score_ar <- function(arg_qr_fit, arg_threshold_method, arg_var_fit, arg_var_stationary) {
  nobs <- arg_qr_fit$nobs
  Z <- arg_qr_fit$score_mat
  if (identical(arg_var_fit, "ols")) {
    phi <- vapply(seq_len(ncol(Z)), function(k) {
      z1 <- Z[-1L, k] - mean(Z[-1L, k]); z0 <- Z[-nobs, k] - mean(Z[-nobs, k])
      denom <- sum(z0^2)
      sum(z1 * z0) / denom
    }, numeric(1L))
  } else {
    taus <- arg_qr_fit$fit$tau
    X <- arg_qr_fit$x_mat
    E <- as.matrix(arg_qr_fit$residuals)
    p <- ncol(X); m <- length(taus)
    # Model-matrix assignment 0 denotes the intercept.
    non_intercept <- which(attr(X, "assign") != 0L)
    rho_x <- rep(1, p)
    rho_x[non_intercept] <- vapply(non_intercept, function(j) stats::cor(X[-1L, j], X[-nobs, j]),
                                   numeric(1L))
    rho_e <- vapply(seq_len(m), function(j) stats::cor(E[-1L, j], E[-nobs, j]), numeric(1L))
    q_tau <- if (identical(arg_threshold_method, "theoretical")) {
      stats::qnorm(taus)
    } else {
      E_std <- scale(E)
      vapply(seq_len(m), function(j) {
        as.numeric(stats::quantile(E_std[, j], probs = taus[j], names = FALSE, type = 7))
      }, numeric(1L))
    }
    p11 <- mapply(function(q, r) .phi2(q, q, r), q_tau, rho_e)
    phi_star <- (p11 - taus^2) / (taus * (1 - taus))
    # Stack coefficients within each quantile.
    phi <- as.vector(outer(rho_x, phi_star))
  }
  bound <- .qrfs_stationary_bound(arg_var_stationary, nobs)
  if (is.na(bound) && any(abs(phi) >= 1))
    stop(sprintf("An AR(1) fit of the scores is not stationary (largest |coefficient| %.3f) and ", max(abs(phi))),
         "arg_var_stationary = \"unrestricted\" imposes no bound; choose \"andrews_monahan\" (0.97), ",
         "\"sul_phillips_choi\" (1 - 1/sqrt(T)), or a number in (0, 1).", call. = FALSE)
  if (!is.na(bound)) phi <- pmin(pmax(phi, -bound), bound)
  list(phi = phi, sigma2 = apply(Z, 2L, stats::var) * (1 - phi^2))
}

# Componentwise AR(1) MSE-optimal K, with the caller's minimum arg_K_min and a basis-specific cap.
# MSE criterion of Hwang & Valdes (2025, Section 6) with an identity weight.
compute_mse_optimal_K_ar <- function(arg_qr_fit, arg_K_min,
                                     arg_threshold_method,
                                     arg_var_fit,
                                     arg_var_stationary) {
  T_n <- arg_qr_fit$nobs
  K_max <- T_n - 1L
  bias_constant <- pi^2 / 6   # Hwang & Valdes (2025), DCT-II
  ar1 <- .estimate_score_ar(arg_qr_fit, arg_threshold_method, arg_var_fit, arg_var_stationary)
  phi <- ar1$phi; s2 <- ar1$sigma2
  omega <- s2 / (1 - phi)^2
  bias <- -2 * bias_constant * s2 * phi / (1 - phi)^4

  # AMSE = V / K + squared bias * (K / T)^4.
  V_scalar <- sum(omega)^2 + sum(omega^2)
  B_scalar <- sum(bias^2)

  K_min <- arg_K_min
  K_tilde <- (V_scalar / (4 * B_scalar))^(1/5) * T_n^(4/5)
  as.integer(min(max(K_min, ceiling(K_tilde)), K_max))
}

# MSE-optimal Bartlett M from the AR(1) coefficients of the stacked score columns.
# Galvao & Yoon (2024, Corollary 5.1); 1.1447 is the Bartlett constant of Andrews (1991).
compute_mse_optimal_M_ar <- function(arg_qr_fit, arg_threshold_method,
                                     arg_var_fit,
                                     arg_var_stationary) {
  nobs <- arg_qr_fit$nobs
  ar1 <- .estimate_score_ar(arg_qr_fit, arg_threshold_method, arg_var_fit, arg_var_stationary)
  phi_vec <- ar1$phi; sigma2_vec <- ar1$sigma2

  denom <- sum(sigma2_vec^2 / (1 - phi_vec)^4)
  numer <- sum(4 * phi_vec^2 * sigma2_vec^2 / ((1 - phi_vec)^6 * (1 + phi_vec)^2))
  theta_hat <- numer / denom
  S_opt <- 1.1447 * (theta_hat * nobs)^(1/3)

  min(max(as.integer(round(S_opt)), 1L), as.integer(nobs))
}

# Stationary covariance \Gamma_0 = Var(Z_t) from \Gamma_0 - A \Gamma_0 A' = \Sigma_\epsilon for a VAR(1) score process.
# Required to evaluate the Bartlett spectral bias moment \sum |h| \Gamma(h) in kernel bandwidth selection.
.solve_discrete_lyapunov <- function(A, Q) {
  n <- nrow(A)
  schur <- Matrix::Schur(A)
  U <- schur$Q
  S <- schur$T
  D <- crossprod(U, Q %*% U)
  blocks <- list()
  last <- n
  while (last > 0L) {
    first <- if (last > 1L && S[last, last - 1L] != 0) last - 1L else last
    blocks[[length(blocks) + 1L]] <- first:last
    last <- first - 1L
  }
  # Solve the quasi-triangular blocks from the bottom right.
  Y <- SY <- matrix(0, n, n)
  for (J in blocks) {
    after_j <- if (max(J) < n) seq.int(max(J) + 1L, n) else integer()
    C <- D[, J, drop = FALSE]
    if (length(after_j))
      C <- C + SY[, after_j, drop = FALSE] %*% t(S[J, after_j, drop = FALSE])
    S_jj <- S[J, J, drop = FALSE]
    for (I in blocks) {
      after_i <- if (max(I) < n) seq.int(max(I) + 1L, n) else integer()
      rhs <- C[I, , drop = FALSE]
      if (length(after_i))
        rhs <- rhs + S[I, after_i, drop = FALSE] %*%
          Y[after_i, J, drop = FALSE] %*% t(S_jj)
      block_system <- diag(length(I) * length(J)) - kronecker(S_jj, S[I, I, drop = FALSE])
      Y[I, J] <- matrix(solve(block_system, c(rhs)), length(I), length(J))
    }
    SY[, J] <- S %*% Y[, J, drop = FALSE]
  }
  U %*% Y %*% t(U)
}

# VAR plug-in MSE-optimal Bartlett M, returned as an integer in [1, nobs].
# Andrews (1991) MSE-optimal Bartlett rule with a VAR(1) plug-in.
compute_mse_optimal_M_var <- function(arg_qr_fit,
                                      arg_threshold_method,
                                      arg_var_stationary,
                                      arg_var_fit) {
  nobs <- arg_qr_fit$nobs
  var_approx <- .estimate_stacked_var(arg_qr_fit,
                               arg_threshold_method = arg_threshold_method, arg_var_fit = arg_var_fit,
                               arg_var_stationary = arg_var_stationary)
  mp <- arg_qr_fit$m * arg_qr_fit$p
  A <- var_approx$A; Sigma_eps <- var_approx$Sigma_eps
  inv_IA <- solve(diag(mp) - A)
  Omega_hat <- inv_IA %*% Sigma_eps %*% t(inv_IA)
  Omega_hat <- (Omega_hat + t(Omega_hat)) / 2

  Gamma0 <- .solve_discrete_lyapunov(A, Sigma_eps)
  Gamma0 <- (Gamma0 + t(Gamma0)) / 2
  inv_IA2 <- inv_IA %*% inv_IA
  # Bartlett bias moment: sum_h |h| Gamma(h).
  Lam1 <- A %*% inv_IA2 %*% Gamma0 + Gamma0 %*% t(A) %*% t(inv_IA2)
  Lam1 <- (Lam1 + t(Lam1)) / 2

  V_scalar <- sum(diag(Omega_hat))^2 + sum(Omega_hat^2)
  B_scalar <- sum(Lam1^2)
  M_tilde <- (3 * B_scalar / V_scalar)^(1/3) * nobs^(1/3)
  as.integer(max(1, min(round(M_tilde), nobs)))
}

# Test-optimal K
# The paper, Section 4.5, adapting Hwang & Valdes (2025); kappa as in Hwang & Valdes.
# arg_K_min is the floor, set by the caller independently of nrow(arg_R).
# arg_kappa = user selected size-inflation factor. if not specified (NULL) 1.05 for T >= 500, otherwise 1.10.
compute_test_optimal_K <- function(arg_qr_fit, arg_R,
                                   arg_alpha,
                                   arg_kappa,
                                   arg_K_min,
                                   arg_threshold_method,
                                   arg_var_stationary,
                                   arg_var_fit) {
  T_n <- arg_qr_fit$nobs
  K_max <- T_n - 1L
  bias_constant <- pi^2 / 6   # Hwang & Valdes (2025), DCT-II
  var_approx <- .estimate_stacked_var(arg_qr_fit,
                               arg_threshold_method = arg_threshold_method, arg_var_fit = arg_var_fit,
                               arg_var_stationary = arg_var_stationary)
  os_moments <- .compute_os_moments(var_approx$A, var_approx$Sigma_eps,
                                     bias_constant)
  Omega_hat <- os_moments$Omega_hat
  B_tilde <- os_moments$B_tilde

  Hinv_block <- as.matrix(Matrix::bdiag(arg_qr_fit$Hinv_list))

  if (is.null(arg_kappa)) arg_kappa <- if (T_n >= 500) 1.05 else 1.10

  arg_R <- as.matrix(arg_R)
  d_a <- nrow(arg_R)

  # Average bias in the tested directions.
  R <- arg_R
  RD <- R %*% Hinv_block
  num <- RD %*% B_tilde %*% t(RD)
  den <- RD %*% Omega_hat %*% t(RD)
  num <- (num + t(num)) / 2
  den <- (den + t(den)) / 2
  B_bar_a <- sum(diag(solve(den, num))) / d_a

  K_min <- arg_K_min
  Q_chi2 <- stats::qchisq(1 - arg_alpha, df = d_a)
  if (B_bar_a > 0) {
    delta2 <- stats::uniroot(
      function(d2) stats::pchisq(Q_chi2, df = d_a, ncp = d2,
                          lower.tail = FALSE) - 0.75,
      lower = 1e-6, upper = 200, extendInt = "upX"
    )$root
    g_da_d2 <- stats::dchisq(Q_chi2, df = d_a,     ncp = delta2)
    g_da2_d <- stats::dchisq(Q_chi2, df = d_a + 2, ncp = delta2)
    K_tilde <- (delta2 * g_da2_d / (4 * g_da_d2 * B_bar_a))^(1/3) * T_n^(2/3)
  } else {
    # abs() avoids negative zero when the bias is zero.
    K_tilde <- ((arg_kappa - 1) * arg_alpha /
                (stats::dchisq(Q_chi2, df = d_a) * Q_chi2 * abs(B_bar_a)))^(1/2) * T_n
  }

  # Bound before integer conversion; zero bias selects K_max.
  as.integer(min(max(K_min, ceiling(K_tilde)), K_max))
}

# Test-optimal Bartlett M
# The paper, Section 4.5, adapting Sun (2014).
# arg_kappa = user selected size-inflation factor. if not specified (NULL) 1.05 for T >= 500, otherwise 1.10.
compute_test_optimal_M <- function(arg_qr_fit, arg_R,
                                   arg_alpha,
                                   arg_kappa,
                                   arg_threshold_method,
                                   arg_var_stationary,
                                   arg_var_fit) {
  T_n <- arg_qr_fit$nobs
  if (is.null(arg_kappa)) arg_kappa <- if (T_n >= 500) 1.05 else 1.10
  arg_R <- as.matrix(arg_R)
  d_a <- nrow(arg_R)

  var_approx <- .estimate_stacked_var(arg_qr_fit, arg_threshold_method,
                               arg_var_stationary = arg_var_stationary,
                               arg_var_fit = arg_var_fit)
  mp <- arg_qr_fit$m * arg_qr_fit$p
  A <- var_approx$A
  Sigma_eps <- var_approx$Sigma_eps

  inv_IA    <- solve(diag(mp) - A)
  Omega_hat <- inv_IA %*% Sigma_eps %*% t(inv_IA)
  Omega_hat <- (Omega_hat + t(Omega_hat)) / 2
  Gamma0 <- .solve_discrete_lyapunov(A, Sigma_eps)
  Gamma0 <- (Gamma0 + t(Gamma0)) / 2
  inv_IA2 <- inv_IA %*% inv_IA
  L <- A %*% inv_IA2 %*% Gamma0
  B_tilde <- -(L + t(L))

  Hinv_block <- as.matrix(Matrix::bdiag(arg_qr_fit$Hinv_list))
  RD  <- arg_R %*% Hinv_block
  num <- RD %*% B_tilde   %*% t(RD)
  den <- RD %*% Omega_hat %*% t(RD)
  num <- (num + t(num)) / 2
  den <- (den + t(den)) / 2
  B_bar_a <- sum(diag(solve(den, num))) / d_a

  Q_chi2 <- stats::qchisq(1 - arg_alpha, df = d_a)
  if (B_bar_a > 0) {
    delta2 <- stats::uniroot(
      function(d2) stats::pchisq(Q_chi2, df = d_a, ncp = d2,
                          lower.tail = FALSE) - 0.75,
      lower = 1e-6, upper = 200, extendInt = "upX"
    )$root
    g_p  <- stats::dchisq(Q_chi2, df = d_a,     ncp = delta2)
    g_p2 <- stats::dchisq(Q_chi2, df = d_a + 2, ncp = delta2)
    M_tilde <- sqrt(3 * g_p * B_bar_a / (delta2 * g_p2)) * sqrt(T_n)
  } else {
    M_tilde <- stats::dchisq(Q_chi2, df = d_a) * Q_chi2 * (-B_bar_a) /
      ((arg_kappa - 1) * arg_alpha)
  }

  # Bound before integer conversion; zero bias selects M = 1.
  as.integer(max(1, min(floor(M_tilde), T_n)))
}

# The stack-Wald statistic of the paper, Section 4.2, eq. (16).
compute_stacked_wald <- function(arg_qr_fit, arg_R, arg_null_val,
                                 arg_har_method,
                                 arg_M, arg_K) {
  R <- as.matrix(arg_R)
  Hinv_block <- as.matrix(Matrix::bdiag(arg_qr_fit$Hinv_list))
  RH <- R %*% Hinv_block
  # Estimate covariance directly for the tested coefficient combinations.
  tested_scores <- arg_qr_fit$score_mat %*% t(RH)
  V_hat <- .compute_Jhat_from_scores(tested_scores, arg_har_method, arg_M, arg_K, NULL)
  deviation <- R %*% as.vector(arg_qr_fit$coef_mat) - as.numeric(arg_null_val)
  compute_wald(deviation, V_hat)
}


# The sup statistic of the paper, Section 4.1.
compute_sup_wald <- function(arg_qr_fit, arg_R, arg_null_val,
                             arg_har_method,
                             arg_M, arg_K) {
  R <- as.matrix(arg_R)
  m <- arg_qr_fit$m; p <- arg_qr_fit$p
  # Reuse one OS basis across quantiles.
  Phi <- if (arg_har_method == "os-har")
    make_os_basis(arg_qr_fit$nobs, arg_K)
  null <- matrix(as.numeric(arg_null_val), ncol = 1)
  wald_tau <- numeric(m)
  # Standard errors of the tested combinations, for the uniform band of a one-restriction test.
  se_tau_mat <- matrix(NA_real_, nrow(R), m)
  for (j in seq_len(m)) {
    ii <- ((j - 1L) * p + 1L):(j * p)
    RH <- R %*% arg_qr_fit$Hinv_list[[j]]
    tested_scores <- arg_qr_fit$score_mat[, ii, drop = FALSE] %*% t(RH)
    V_hat <- .compute_Jhat_from_scores(tested_scores, arg_har_method, arg_M, arg_K, Phi)
    deviation <- R %*% arg_qr_fit$coef_mat[, j, drop = FALSE] - null
    wald_tau[j] <- compute_wald(deviation, V_hat)
    se_tau_mat[, j] <- sqrt(diag(V_hat))
  }
  list(sup_wald = max(wald_tau), wald_tau = wald_tau,
       se_tau_mat = se_tau_mat)
}

# Bartlett lagged cross-products, interpolated when b falls between grid points.
# arg_diagonal returns marginal integrals only.
.bartlett_shift_crossprod <- function(arg_path, arg_b, arg_nr, arg_diagonal) {
  path <- as.matrix(arg_path)
  lag_grid <- as.numeric(arg_nr) * arg_b
  lag_lo <- floor(lag_grid)
  frac <- lag_grid - lag_lo
  tol <- 32 * .Machine$double.eps * max(1, abs(lag_grid))
  if (frac <= tol) frac <- 0
  if (1 - frac <= tol) { lag_lo <- lag_lo + 1L; frac <- 0 }

  if (frac == 0) {
    n_use <- as.integer(arg_nr) + 1L - as.integer(lag_lo)
    idx_lo <- seq_len(n_use)
    shifted <- path[idx_lo + as.integer(lag_lo), , drop = FALSE]
  } else {
    n_use <- as.integer(arg_nr) - as.integer(lag_lo)
    idx_lo <- seq_len(n_use)
    shifted <- (1 - frac) * path[idx_lo + as.integer(lag_lo), , drop = FALSE] +
               frac * path[idx_lo + as.integer(lag_lo) + 1L, , drop = FALSE]
    # The last interval has width (1 - frac) / arg_nr.
    shifted[n_use, ] <- (1 - frac) * shifted[n_use, ]
  }
  if (arg_diagonal) colSums(path[idx_lo, , drop = FALSE] * shifted) else
    crossprod(path[idx_lo, , drop = FALSE], shifted)
}

# Bartlett fixed-b covariance functional of a tied-down path; b = M/T.
# The P functional of Kiefer & Vogelsang (2005); the paper, Theorem 3.2.
# arg_nr is the grid interval count, arg_dr the step; arg_diagonal returns only variances.
compute_P_functional <- function(arg_g_tilde, arg_b, arg_dr, arg_nr,
                                 arg_diagonal = FALSE) {
  base <- if (arg_diagonal) colSums(arg_g_tilde^2) else crossprod(arg_g_tilde)
  shifted <- .bartlett_shift_crossprod(arg_g_tilde, arg_b, arg_nr, arg_diagonal)
  if (arg_diagonal) (2 * arg_dr / arg_b) * (base - shifted) else
    as.matrix((2 * arg_dr / arg_b) * base - (arg_dr / arg_b) * (shifted + t(shifted)))
}

# Retain leading positive covariance eigenvalues and return their Gaussian factor.
# The truncated eigen-decomposition of the paper, Algorithm 1, step (i).
# arg_min_rank is the minimum retained dimension; arg_coverage is the variance fraction.
truncate_lambda_eigendecomp <- function(arg_lambda_hat, arg_min_rank, arg_coverage) {
  eigs_all <- eigen(arg_lambda_hat, symmetric = TRUE)
  pos_idx <- which(eigs_all$values > 0)

  if (length(pos_idx) < arg_min_rank) {
    stop(sprintf(
      "arg_lambda_hat does not meet the rank condition. Only %d positive eigenvalues are available; requested at least %d.",
      length(pos_idx), arg_min_rank
    ), call. = FALSE)
  }

  pos_evals <- eigs_all$values[pos_idx]
  total_pos <- sum(pos_evals)
  cum_prop <- cumsum(pos_evals) / total_pos
  num_eigs <- which(cum_prop >= arg_coverage)[1]
  num_eigs <- min(length(pos_idx), max(num_eigs, arg_min_rank))

  eigs <- pos_evals[1:num_eigs]
  eigvecs_mat <- eigs_all$vectors[, pos_idx[1:num_eigs], drop = FALSE]
  # Fix signs by the largest absolute component for stable seeded draws.
  signs <- apply(eigvecs_mat, 2L, function(v) sign(v[which.max(abs(v))]))
  eigvecs_mat <- sweep(eigvecs_mat, 2L, signs, "*")

  list(
    sqrt_eigs_eigvecs = sqrt(eigs) * t(eigvecs_mat),
    num_eigs = num_eigs
  )
}

# Brownian path with arg_nr intervals of length arg_dr, including the initial zero row.
# arg_eig_trunc supplies the Gaussian factor, optionally projected onto restrictions.
simulate_gaussian_path <- function(arg_eig_trunc, arg_nr, arg_dr) {
  num_eigs <- arg_eig_trunc$num_eigs
  sqrt_eigs_eigvecs <- arg_eig_trunc$sqrt_eigs_eigvecs

  bm <- matrix(0, nrow = arg_nr + 1, ncol = num_eigs)
  d_bm <- matrix(stats::rnorm(arg_nr * num_eigs), nrow = arg_nr, ncol = num_eigs)
  bm[2:(arg_nr + 1), ] <- apply(d_bm * sqrt(arg_dr), 2, cumsum)

  bm %*% sqrt_eigs_eigvecs
}

# Select integer Bartlett M: numeric arg_type, 'simple', 'test_optimal',
# 'mse_optimal_var', or 'mse_optimal_ar'. Test-optimal selection needs stacked arg_R.
# arg_R and arg_alpha have no default; only 'test_optimal' reads them.
compute_M <- function(arg_qr_fit, arg_type, arg_R, arg_alpha,
                      arg_threshold_method = "theoretical", arg_var_stationary = "andrews_monahan",
                      arg_kappa = NULL, arg_var_fit = "gaussian") {
  nobs <- arg_qr_fit$nobs
  if (is.numeric(arg_type)) return(as.integer(arg_type))
  switch(arg_type,
    "simple" = as.integer(min(nobs, max(1L, floor(2 * nobs^(1/3))))),
    "test_optimal" = compute_test_optimal_M(
      arg_qr_fit, arg_R, arg_alpha = arg_alpha, arg_kappa = arg_kappa,
      arg_threshold_method = arg_threshold_method, arg_var_fit = arg_var_fit, arg_var_stationary = arg_var_stationary),
    "mse_optimal_var" = compute_mse_optimal_M_var(
      arg_qr_fit, arg_threshold_method = arg_threshold_method, arg_var_fit = arg_var_fit, arg_var_stationary = arg_var_stationary),
    "mse_optimal_ar" = compute_mse_optimal_M_ar(
      arg_qr_fit, arg_threshold_method = arg_threshold_method, arg_var_fit = arg_var_fit,
      arg_var_stationary = arg_var_stationary),
    stop("arg_type must be numeric or one of 'test_optimal', 'mse_optimal_var', ",
         "'mse_optimal_ar', 'simple'.", call. = FALSE)
  )
}

# Select the test and nuisance bandwidth M or basis count K from validated choices.
# arg_K_min is the floor of the test K; the nuisance K is floored at m*p + 4. NULL skips nuisance selection.
# The nuisance rules read no restriction, so arg_R = NULL is passed for them.
.qrfs_select_M_K <- function(arg_qr_fit, arg_har_method, arg_test_choice,
                              arg_nuisance_choice, arg_R, arg_K_min, arg_alpha,
                              arg_kappa, arg_threshold_method, arg_var_fit,
                              arg_var_stationary) {
  M <- K <- M_lambda <- K_lambda <- NA_integer_
  if (arg_har_method == "ker-har") {
    M <- compute_M(arg_qr_fit, arg_test_choice, arg_R = arg_R, arg_alpha = arg_alpha,
      arg_kappa = arg_kappa, arg_threshold_method = arg_threshold_method,
      arg_var_fit = arg_var_fit, arg_var_stationary = arg_var_stationary)
    if (!is.null(arg_nuisance_choice))
      M_lambda <- compute_M(arg_qr_fit, arg_nuisance_choice, arg_R = NULL, arg_alpha = arg_alpha,
        arg_kappa = arg_kappa, arg_threshold_method = arg_threshold_method,
        arg_var_fit = arg_var_fit, arg_var_stationary = arg_var_stationary)
  } else {
    K <- compute_K(arg_qr_fit, arg_test_choice, arg_R = arg_R, arg_alpha = arg_alpha,
      arg_K_min = arg_K_min, arg_kappa = arg_kappa,
      arg_threshold_method = arg_threshold_method, arg_var_fit = arg_var_fit,
      arg_var_stationary = arg_var_stationary)
    if (!is.null(arg_nuisance_choice))
      K_lambda <- compute_K(arg_qr_fit, arg_nuisance_choice, arg_R = NULL, arg_alpha = arg_alpha,
        arg_K_min = arg_qr_fit$m * arg_qr_fit$p + 4L, arg_kappa = arg_kappa,
        arg_threshold_method = arg_threshold_method, arg_var_fit = arg_var_fit,
        arg_var_stationary = arg_var_stationary)
  }
  list(M = M, K = K, M_lambda = M_lambda, K_lambda = K_lambda)
}

# Nonpivotal sup CV and p-value, preserving dependence across quantiles: the paper's Algorithm 1,
# with the fixed-b functional (Theorem 3.2) for the kernel and the fixed-K limit (Theorem 3.3) for the
# series. arg_R applies at each quantile; arg_M / arg_K is the statistic's smoothing, arg_M_lambda /
# arg_K_lambda the nuisance's. The Gaussian factor of the truncated nuisance covariance is projected
# onto the tested combinations once; every draw is then an endpoint with its covariance.
direct_cv_sup <- function(arg_qr_fit, arg_R, arg_har_method, arg_M, arg_K,
                          arg_M_lambda, arg_K_lambda, arg_stat, arg_n_cv, arg_n_r,
                          arg_coverage, arg_alpha) {
  kernel <- arg_har_method == "ker-har"
  m <- arg_qr_fit$m; p <- arg_qr_fit$p; q <- nrow(arg_R); nobs <- arg_qr_fit$nobs
  J_hat <- compute_Jhat_stacked(arg_qr_fit, arg_har_method, arg_M_lambda, arg_K_lambda)
  eig_trunc <- truncate_lambda_eigendecomp(J_hat / nobs, q, arg_coverage)
  factor <- eig_trunc$sqrt_eigs_eigvecs
  eig_trunc$sqrt_eigs_eigvecs <- do.call(cbind, lapply(seq_len(m), function(j) {
    idx <- (j - 1L) * p + seq_len(p)
    factor[, idx, drop = FALSE] %*% t(arg_R %*% arg_qr_fit$Hinv_list[[j]])
  }))
  # With one restriction every quantile's Wald term is a ratio, so the m variances suffice; with more,
  # each quantile's q x q covariance is inverted.
  blocks <- if (q > 1L) lapply(seq_len(m), function(j) (j - 1L) * q + seq_len(q))
  if (kernel) {
    b <- arg_M / nobs
    dr <- 1 / arg_n_r
    r_grid <- seq(0, 1, length.out = arg_n_r + 1L)
  }
  draws <- numeric(arg_n_cv)
  for (sim in seq_len(arg_n_cv)) {
    if (kernel) {  # Endpoint and Brownian bridge, the P functional of the bridge as covariance.
      path <- simulate_gaussian_path(eig_trunc, arg_n_r, dr)
      endpoint <- path[arg_n_r + 1L, ]
      bridge <- path - r_grid %o% endpoint
      walds <- if (q == 1L) {
        compute_wald(endpoint, compute_P_functional(bridge, b, dr, arg_n_r, arg_diagonal = TRUE))
      } else {
        vapply(blocks, function(idx) compute_wald(endpoint[idx],
          compute_P_functional(bridge[, idx, drop = FALSE], b, dr, arg_n_r,
                               arg_diagonal = FALSE)), numeric(1))
      }
    } else {  # Independent Gaussian endpoint and K projections.
      endpoint <- as.numeric(stats::rnorm(eig_trunc$num_eigs) %*% eig_trunc$sqrt_eigs_eigvecs)
      projections <- matrix(stats::rnorm(arg_K * eig_trunc$num_eigs), nrow = arg_K) %*%
        eig_trunc$sqrt_eigs_eigvecs
      walds <- if (q == 1L) {
        compute_wald(endpoint, colSums(projections^2) / arg_K)
      } else {
        vapply(blocks, function(idx) compute_wald(endpoint[idx],
          crossprod(projections[, idx, drop = FALSE]) / arg_K), numeric(1))
      }
    }
    draws[sim] <- max(walds)
  }
  list(cv = as.numeric(stats::quantile(draws, 1 - arg_alpha, type = 1)),
       p_value = mean(draws >= arg_stat))
}

# Pivotal stacked OS CV and p-value from the scaled F(q, K - q + 1) distribution.
# The paper, Section 4.2, eq. (19); Sun (2013).
# q is the restriction count, K the basis count, and arg_stat the observed Wald statistic.
compute_fixed_K_cv <- function(q, K, alpha, arg_stat) {
  df2 <- K - q + 1L
  scale <- q * K / df2
  cv <- scale * stats::qf(1 - alpha, df1 = q, df2 = df2)
  names(cv) <- paste0("alpha=", alpha)

  list(cv = cv,
       p_value = stats::pf(arg_stat / scale, df1 = q, df2 = df2, lower.tail = FALSE))
}

# Convert linear hypothesis strings to R, null_val and labels.
# Coefficients use name[tau]; a stacked hypothesis names a quantile for every coefficient.
# Sup hypotheses apply separately at each quantile and cannot specify tau.
parse_hypothesis <- function(arg_hypothesis, arg_qr_fit, arg_type) {
  coef_names <- rownames(arg_qr_fit$coef_mat)
  taus <- arg_qr_fit$fit$tau
  m <- arg_qr_fit$m
  p <- arg_qr_fit$p
  n_par <- if (arg_type == "sup") p else m * p

  # Position in stacked theta (or in beta(tau) for a sup test).
  position <- function(name, tau) {
    k <- match(name, coef_names)
    if (is.na(k))
      stop("Unknown coefficient '", name, "'. The model has: ",
           paste(coef_names, collapse = ", "), ".", call. = FALSE)
    if (arg_type == "sup") {
      if (!is.null(tau)) stop("A sup hypothesis cannot name a quantile.", call. = FALSE)
      return(k)
    }
    if (is.null(tau))
      stop("A stacked hypothesis names the quantiles it is about, as in ",
           "\"NFCI[0.1] = NFCI[0.9]\".", call. = FALSE)
    j <- which(abs(taus - tau) < 1e-8)
    if (length(j) != 1L) stop("Quantile ", tau, " was not fitted.", call. = FALSE)
    (j - 1L) * p + k
  }

  # Substitute coefficient positions in theta.
  to_theta <- function(e) {
    if (is.numeric(e)) return(e)
    if (is.symbol(e)) return(bquote(.theta[.(position(as.character(e), NULL))]))
    op <- as.character(e[[1]])
    if (op == "[")
      return(bquote(.theta[.(position(as.character(e[[2]]), eval(e[[3]], baseenv())))]))
    if (!op %in% c("+", "-", "*", "/", "("))
      stop("A hypothesis may use only + - * / and parentheses.", call. = FALSE)
    e[-1] <- lapply(as.list(e)[-1], to_theta)
    e
  }

  # Expand a = b = c into a = b and b = c.
  sides <- lapply(strsplit(arg_hypothesis, "==|="), trimws)
  if (any(lengths(sides) < 2L)) stop("A hypothesis needs an '='.", call. = FALSE)
  pairs <- unlist(lapply(sides, function(s) paste(s[-length(s)], "=", s[-1])))

  R <- NULL
  null_val <- NULL
  label <- NULL
  for (h in pairs) {
    side <- strsplit(h, " = ", fixed = TRUE)[[1]]
    e <- call("-", str2lang(side[1]), call("(", str2lang(side[2])))
    e_j <- to_theta(e)
    f <- function(theta) eval(e_j, list(.theta = theta), baseenv())
    f_0 <- f(numeric(n_par))
    row <- vapply(seq_len(n_par),
                  function(k) f(replace(numeric(n_par), k, 1)) - f_0, numeric(1))
    if (abs(f(rep(2, n_par)) - f_0 - 2 * sum(row)) > 1e-10)
      stop("The hypothesis is not linear in the coefficients: ", h, call. = FALSE)
    R <- rbind(R, row, deparse.level = 0)
    null_val <- c(null_val, -f_0)
    label <- c(label, h)
  }
  list(R = R, null_val = null_val, label = label)
}

# Validate alpha, the smoothing rules or values and the truncation coverage once per call.
# NULL arg_nuisance skips the nuisance and coverage checks, which a stacked test has no use for;
# arg_q is the minimum test K.
.qrfs_check_inputs <- function(arg_har_method, arg_alpha, arg_smoothing, arg_nuisance,
                               arg_coverage, arg_nobs, arg_q, arg_kappa) {
  kernel <- arg_har_method == "ker-har"
  if (!is.numeric(arg_alpha) || length(arg_alpha) != 1L || !is.finite(arg_alpha) ||
      arg_alpha <= 0 || arg_alpha >= 1)
    stop("arg_alpha must be one finite number in (0, 1).", call. = FALSE)
  if (!is.null(arg_kappa) &&
      (!is.numeric(arg_kappa) || length(arg_kappa) != 1L || !is.finite(arg_kappa) ||
       arg_kappa <= 1))
    stop("arg_kappa must be NULL or one finite number greater than 1.", call. = FALSE)
  rules <- c("simple", "test_optimal", "mse_optimal_var", "mse_optimal_ar")
  if (is.numeric(arg_smoothing)) {
    if (length(arg_smoothing) != 1L || !is.finite(arg_smoothing) || arg_smoothing != floor(arg_smoothing))
      stop("Test smoothing must be one finite integer or a supported rule.", call. = FALSE)
    if (kernel && (arg_smoothing < 1 || arg_smoothing > arg_nobs))
      stop(sprintf("Test M = %g: it must be between 1 and nobs = %d.", arg_smoothing, arg_nobs), call. = FALSE)
    if (!kernel && (arg_smoothing < arg_q || arg_smoothing > arg_nobs))
      stop(sprintf("Test K = %g: it must be between the number of restrictions (%d) and nobs = %d.",
                   arg_smoothing, arg_q, arg_nobs), call. = FALSE)
  } else if (!is.character(arg_smoothing) || length(arg_smoothing) != 1L || !arg_smoothing %in% rules) {
    stop("Unknown smoothing rule for the test.", call. = FALSE)
  }
  if (!is.null(arg_nuisance)) {
    if (is.numeric(arg_nuisance)) {
      if (length(arg_nuisance) != 1L || !is.finite(arg_nuisance) ||
          arg_nuisance != floor(arg_nuisance) || arg_nuisance < 1)
        stop("Nuisance smoothing must be one positive integer or an MSE/simple rule.", call. = FALSE)
      if (arg_nuisance > arg_nobs)
        stop(sprintf("Nuisance %s = %g: it must be at most nobs = %d.",
                     if (kernel) "M" else "K", arg_nuisance, arg_nobs), call. = FALSE)
    } else if (!is.character(arg_nuisance) || length(arg_nuisance) != 1L ||
               !arg_nuisance %in% c("simple", "mse_optimal_var", "mse_optimal_ar")) {
      stop("Nuisance smoothing must be numeric, 'simple', 'mse_optimal_var', or 'mse_optimal_ar'.", call. = FALSE)
    }
    if (!is.numeric(arg_coverage) || length(arg_coverage) != 1L ||
        !is.finite(arg_coverage) || arg_coverage <= 0 || arg_coverage > .99)
      stop("arg_coverage must be one finite fraction in (0, 0.99].", call. = FALSE)
  }
  invisible(NULL)
}

# Extracts the unique quantiles named in a hypothesis (e.g., "x[0.1] = x[0.9]" -> c(0.1, 0.9)).
.qrfs_taus_named <- function(arg_hypothesis) {
  if (!is.character(arg_hypothesis)) return(numeric())
  found <- numeric()
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    if (identical(e[[1]], as.name("[")) && length(e) == 3L && is.numeric(e[[3]]))
      found <<- c(found, e[[3]])
    for (z in as.list(e)[-1]) walk(z)
  }
  for (h in arg_hypothesis)
    for (side in trimws(strsplit(h, "==|=")[[1]])) walk(str2lang(side))
  sort(unique(found))
}

# Resolves the quantile grid and the label of the range the null covers, for the uniform sup test
# and the shape tests. arg_tau_grid is NULL for 0.1 0.3 0.5 0.7 0.9 over (0, 1), the quantiles
# themselves, or list(range = c(lo, hi), n = points) for a generated grid in that range.
# The list keeps range and n together so that neither can be given without the other, and it is what
# lets the null be stated over the range asked about rather than over the span of the points chosen.
# A range reaching 0 or 1 excludes that endpoint. Explicit grids are sorted and duplicates removed.
# A single quantile is refused: that is a pointwise test, qrfstest_stacked's, and two doors must not
# answer it with two slightly different critical values.
.qrfs_tau_grid <- function(arg_tau_grid) {
  if (is.null(arg_tau_grid))
    return(list(grid = c(0.1, 0.3, 0.5, 0.7, 0.9), range_label = "(0, 1)"))

  if (!is.list(arg_tau_grid)) {
    if (!is.numeric(arg_tau_grid) || !length(arg_tau_grid) || any(!is.finite(arg_tau_grid)) ||
        any(arg_tau_grid <= 0) || any(arg_tau_grid >= 1))
      stop("arg_tau_grid must contain finite quantiles strictly inside (0, 1); ",
           "use list(range = c(0, 0.25), n = 5) to specify a range.", call. = FALSE)
    grid <- sort(unique(arg_tau_grid))
    if (length(grid) < 2L)
      stop("A uniform or shape test needs more than one quantile; to test at a single quantile ",
           "use the stacked test, e.g. qrfstest_stacked(..., \"NFCI[0.5] = 0\").", call. = FALSE)
    return(list(grid = grid, range_label = sprintf("[%g, %g]", grid[1L], grid[length(grid)])))
  }

  if (!setequal(names(arg_tau_grid), c("range", "n")))
    stop("A list arg_tau_grid must have exactly the elements 'range' and 'n'.", call. = FALSE)
  tau_range <- arg_tau_grid$range
  n_tau <- arg_tau_grid$n
  if (!is.numeric(tau_range) || length(tau_range) != 2L || any(!is.finite(tau_range)) ||
      tau_range[1] >= tau_range[2] || tau_range[1] < 0 || tau_range[2] > 1)
    stop("arg_tau_grid$range must be two increasing numbers within [0, 1].", call. = FALSE)
  if (!is.numeric(n_tau) || length(n_tau) != 1L || !is.finite(n_tau) ||
      n_tau < 2 || n_tau != floor(n_tau))
    stop("arg_tau_grid$n must be one integer of at least 2.", call. = FALSE)

  # An endpoint at 0 or 1 is dropped, so the points stay interior to the range.
  open_low <- tau_range[1] <= 0
  open_high <- tau_range[2] >= 1
  grid <- seq(tau_range[1], tau_range[2], length.out = n_tau + open_low + open_high)
  if (open_low) grid <- grid[-1]
  if (open_high) grid <- grid[-length(grid)]
  # Preserve the first and last quantiles; round only the middle ones.
  grid[-c(1L, n_tau)] <- round(grid[-c(1L, n_tau)], 2)
  if (any(diff(grid) <= 0))
    stop("Rounding the middle quantiles to two decimals does not give an increasing grid; ",
         "use fewer points or give the quantiles directly.", call. = FALSE)
  list(grid = grid,
       range_label = sprintf("%s%g, %g%s", if (open_low) "(" else "[",
                             tau_range[1], tau_range[2], if (open_high) ")" else "]"))
}

# Joint test of restrictions that name their quantiles: the stacked Wald statistic of equation (16)
qrfstest_stacked <- function(arg_formula, arg_data, arg_hypothesis, arg_har_method,
                             arg_M = "test_optimal", arg_K = "test_optimal",
                             arg_alpha = 0.05, arg_kappa = NULL,
                             arg_threshold_method = "theoretical",
                             arg_var_fit = "gaussian",
                             arg_var_stationary = "andrews_monahan",
                             arg_hinv_method = "ker", arg_rq_method = "br",
                             arg_seed = NULL) {
  if (!is.null(arg_seed)) set.seed(arg_seed)
  arg_har_method <- match.arg(arg_har_method, c("ker-har", "os-har"))
  arg_threshold_method <- match.arg(arg_threshold_method, c("theoretical", "empirical"))
  arg_var_fit <- match.arg(arg_var_fit, c("gaussian", "ols"))
  tau_grid <- .qrfs_taus_named(arg_hypothesis)
  if (!length(tau_grid))
    stop("A stacked hypothesis names the quantiles it is about, as in \"NFCI[0.1] = NFCI[0.9]\".", call. = FALSE)
  if (any(tau_grid <= 0 | tau_grid >= 1))
    stop("Quantiles in a stacked hypothesis must lie strictly inside (0, 1).", call. = FALSE)
  qr_fit <- fit_qr_grid(arg_formula, arg_data, tau_grid, arg_hinv_method, arg_rq_method)
  hypothesis <- parse_hypothesis(arg_hypothesis, qr_fit, "stacked")
  R <- hypothesis$R
  q <- nrow(R)
  if (qr(R)$rank != q) stop("The restrictions are redundant: R must have independent rows.", call. = FALSE)
  null <- as.numeric(hypothesis$null_val)
  if (any(!is.finite(null))) stop("null_val must be finite.", call. = FALSE)

  smoothing <- if (arg_har_method == "ker-har") arg_M else arg_K
  # No nuisance rule or truncation coverage to check.
  .qrfs_check_inputs(arg_har_method, arg_alpha, smoothing, NULL, NULL,
                     qr_fit$nobs, q, arg_kappa)
  # K floor q + 4: the F(q, K - q + 1) limit of equation (19) has a finite variance only then (Section 4.5).
  selected <- .qrfs_select_M_K(qr_fit, arg_har_method, smoothing, NULL, R, q + 4L,
    arg_alpha, arg_kappa, arg_threshold_method, arg_var_fit, arg_var_stationary)
  M <- selected$M; K <- selected$K

  stat <- compute_stacked_wald(qr_fit, R, null, arg_har_method, arg_M = M, arg_K = K)
  # Critical-value simulation: draws and grid intervals per Brownian path. Fixed
  # here rather than offered as arguments; they set its numerical accuracy, not the test.
  n_cv <- 10000L; n_r <- 1000L
  cv_result <- if (arg_har_method == "os-har")
    compute_fixed_K_cv(q, K, arg_alpha, arg_stat = stat) else
    compute_pivotal_cv_stacked(q, M / qr_fit$nobs, arg_alpha, arg_n_cv = n_cv,
                               arg_n_r = n_r, arg_stat = stat)
  cv <- cv_result$cv; names(cv) <- paste0("alpha=", arg_alpha)

  out <- list(
    hypothesis = hypothesis$label, coefficients = qr_fit$coef_mat, stat = stat, cv = cv,
    p_value = cv_result$p_value, reject = as.integer(stat > cv),
    M_effective = M, K_effective = K, M_lambda = NA_integer_, K_lambda = NA_integer_,
    q = q, type = "stacked", har_method = arg_har_method, cv_method = "pivotal",
    alpha = arg_alpha, n_cv = if (arg_har_method == "ker-har") n_cv else NA,
    nobs = qr_fit$nobs, tau = qr_fit$fit$tau, smoothing_rule = smoothing, tau_range = NULL,
    R = R, null_val = null)
  class(out) <- "qrfstest"
  out
}

# The same joint test with the conventional critical value: the small-b kernel HAC of Galvao & Yoon
# (2024) or the large-K OS HAR. Under traditional asymptotics, b -> 0 and K -> infinity as T -> infinity,
# the limits rely on consistency of the variance estimator and do not depend on the kernel, the
# bandwidth or K (paper, Section 2): the stacked Wald statistic of equation (16) is chi-square with
# q degrees of freedom. The default rule for M or K is therefore MSE-optimal; "test_optimal" or a
# number may be given instead. Nothing is simulated, so there is no arg_seed.
qrhac <- function(arg_formula, arg_data, arg_hypothesis, arg_har_method,
                  arg_M = "mse_optimal_ar", arg_K = "mse_optimal_var",
                  arg_alpha = 0.05, arg_kappa = NULL,
                  arg_threshold_method = "theoretical",
                  arg_var_fit = "gaussian",
                  arg_var_stationary = "andrews_monahan",
                  arg_hinv_method = "ker", arg_rq_method = "br") {
  arg_har_method <- match.arg(arg_har_method, c("ker-har", "os-har"))
  arg_threshold_method <- match.arg(arg_threshold_method, c("theoretical", "empirical"))
  arg_var_fit <- match.arg(arg_var_fit, c("gaussian", "ols"))
  tau_grid <- .qrfs_taus_named(arg_hypothesis)
  if (!length(tau_grid))
    stop("A stacked hypothesis names the quantiles it is about, as in \"NFCI[0.1] = NFCI[0.9]\".", call. = FALSE)
  if (any(tau_grid <= 0 | tau_grid >= 1))
    stop("Quantiles in a stacked hypothesis must lie strictly inside (0, 1).", call. = FALSE)
  qr_fit <- fit_qr_grid(arg_formula, arg_data, tau_grid, arg_hinv_method, arg_rq_method)
  hypothesis <- parse_hypothesis(arg_hypothesis, qr_fit, "stacked")
  R <- hypothesis$R
  q <- nrow(R)
  if (qr(R)$rank != q) stop("The restrictions are redundant: R must have independent rows.", call. = FALSE)
  null <- as.numeric(hypothesis$null_val)
  if (any(!is.finite(null))) stop("null_val must be finite.", call. = FALSE)

  smoothing <- if (arg_har_method == "ker-har") arg_M else arg_K
  # No nuisance rule or truncation coverage to check.
  .qrfs_check_inputs(arg_har_method, arg_alpha, smoothing, NULL, NULL,
                     qr_fit$nobs, q, arg_kappa)
  # K floor q, not q + 4: the chi-square limit needs only an invertible estimate (paper, footnote 2).
  selected <- .qrfs_select_M_K(qr_fit, arg_har_method, smoothing, NULL, R, q,
    arg_alpha, arg_kappa, arg_threshold_method, arg_var_fit, arg_var_stationary)
  M <- selected$M; K <- selected$K

  stat <- compute_stacked_wald(qr_fit, R, null, arg_har_method, arg_M = M, arg_K = K)
  # The chi-square limit with q degrees of freedom: the small-b limit of equation (18) and the
  # large-K limit of equation (19).
  cv <- stats::qchisq(1 - arg_alpha, q); names(cv) <- paste0("alpha=", arg_alpha)
  p_value <- stats::pchisq(stat, q, lower.tail = FALSE)

  out <- list(
    hypothesis = hypothesis$label, coefficients = qr_fit$coef_mat, stat = stat, cv = cv,
    p_value = p_value, reject = as.integer(stat > cv),
    M_effective = M, K_effective = K, M_lambda = NA_integer_, K_lambda = NA_integer_,
    q = q, type = "stacked", har_method = arg_har_method, cv_method = "chisq",
    alpha = arg_alpha, n_cv = NA,
    nobs = qr_fit$nobs, tau = qr_fit$fit$tau, smoothing_rule = smoothing, tau_range = NULL,
    R = R, null_val = null)
  class(out) <- "qrfstest"
  out
}

# Uniform test over a range of quantiles
# hypothesis may not name a quantile; the grid discretizes the range. arg_tau_grid is NULL for
# 0.1 0.3 0.5 0.7 0.9 over (0, 1), the quantiles themselves, or list(range =, n =).
# A single restriction also gets the uniform confidence band of Section 4.3 (Algorithm 2), returned
# as band and printed; arg_band_plot = TRUE draws it as well (see plot.qrfstest).
qrfstest_sup <- function(arg_formula, arg_data, arg_hypothesis, arg_har_method,
                         arg_tau_grid = NULL,
                         arg_M = "test_optimal", arg_K = "test_optimal",
                         arg_M_lambda = "mse_optimal_ar",
                         arg_K_lambda = "mse_optimal_var",
                         arg_alpha = 0.05, arg_kappa = NULL,
                         arg_threshold_method = "theoretical",
                         arg_var_fit = "gaussian",
                         arg_var_stationary = "andrews_monahan",
                         arg_coverage = 0.99,
                         arg_hinv_method = "ker", arg_rq_method = "br",
                         arg_band_plot = FALSE, arg_seed = NULL) {
  if (!is.null(arg_seed)) set.seed(arg_seed)
  arg_har_method <- match.arg(arg_har_method, c("ker-har", "os-har"))
  arg_threshold_method <- match.arg(arg_threshold_method, c("theoretical", "empirical"))
  arg_var_fit <- match.arg(arg_var_fit, c("gaussian", "ols"))
  if (!is.character(arg_hypothesis))
    stop("arg_hypothesis must be one or more strings, e.g. \"NFCI = 0\".", call. = FALSE)
  tau_spec <- .qrfs_tau_grid(arg_tau_grid)
  qr_fit <- fit_qr_grid(arg_formula, arg_data, tau_spec$grid, arg_hinv_method, arg_rq_method)
  hypothesis <- parse_hypothesis(arg_hypothesis, qr_fit, "sup")
  R <- hypothesis$R
  q <- nrow(R)
  if (qr(R)$rank != q) stop("The restrictions are redundant: R must have independent rows.", call. = FALSE)
  null <- as.numeric(hypothesis$null_val)
  if (any(!is.finite(null))) stop("null_val must be finite.", call. = FALSE)

  smoothing <- if (arg_har_method == "ker-har") arg_M else arg_K
  nuisance <- if (arg_har_method == "ker-har") arg_M_lambda else arg_K_lambda
  .qrfs_check_inputs(arg_har_method, arg_alpha, smoothing, nuisance, arg_coverage,
                     qr_fit$nobs, q, arg_kappa)
  # The paper gives no test-optimal rule for a maximum, so its stacked rule is used, evaluated at
  # I_m (x) R0. The K floor stays q0 + 4 (arg_K_min = q + 4 below), the rank of each per-quantile
  # statistic. Only that rule reads the restriction.
  R_final <- if (smoothing == "test_optimal") kronecker(diag(qr_fit$m), R) else R
  selected <- .qrfs_select_M_K(qr_fit, arg_har_method, smoothing, nuisance, R_final, q + 4L,
    arg_alpha, arg_kappa, arg_threshold_method, arg_var_fit, arg_var_stationary)
  M <- selected$M; K <- selected$K
  M_lambda <- selected$M_lambda; K_lambda <- selected$K_lambda

  result <- compute_sup_wald(qr_fit, R, null, arg_har_method, arg_M = M, arg_K = K)
  stat <- result$sup_wald
  # Critical-value simulation: draws and grid intervals per Brownian path. Fixed
  # here rather than offered as arguments; they set its numerical accuracy, not the test.
  n_cv <- 10000L; n_r <- 1000L
  cv_result <- direct_cv_sup(qr_fit, R, arg_har_method, M, K, M_lambda, K_lambda,
                               arg_stat = stat, arg_n_cv = n_cv, arg_n_r = n_r,
                               arg_coverage = arg_coverage, arg_alpha = arg_alpha)
  cv <- cv_result$cv; names(cv) <- paste0("alpha=", arg_alpha)

  # The uniform confidence band of Section 4.3 for a single restriction: the critical value of the
  # sup-Wald statistic is the square of the sup-|t| one, so the band is R0 b(tau) +- sqrt(cv) s.e.(tau)
  # and the null value leaves it somewhere exactly when the test rejects.
  band <- NULL
  if (q == 1L) {
    estimate <- as.numeric(R %*% qr_fit$coef_mat)
    half_width <- sqrt(as.numeric(cv)) * as.numeric(result$se_tau_mat)
    band <- data.frame(tau = qr_fit$fit$tau, estimate = estimate,
                       lower = estimate - half_width, upper = estimate + half_width)
  }

  out <- list(
    hypothesis = hypothesis$label, coefficients = qr_fit$coef_mat, stat = stat, cv = cv,
    p_value = cv_result$p_value, reject = as.integer(stat > cv),
    M_effective = M, K_effective = K, M_lambda = M_lambda, K_lambda = K_lambda,
    q = q, type = "sup", har_method = arg_har_method, cv_method = "nonpivotal-sim",
    alpha = arg_alpha, n_cv = n_cv, nobs = qr_fit$nobs, tau = qr_fit$fit$tau,
    smoothing_rule = smoothing, tau_range = tau_spec$range_label, R = R, null_val = null,
    wald_tau = result$wald_tau, band = band)
  class(out) <- "qrfstest"
  if (arg_band_plot) {
    if (is.null(band)) warning("A band is drawn only for a single restriction; nothing is plotted.")
    else plot.qrfstest(out)
  }
  out
}

# Print test results and the selected bandwidth fraction b or basis count K.
print.qrfstest <- function(x, ...) {
  kernel <- x$har_method == "ker-har"
  # A qrhac result carries the chi-square calibration and is headed by its method.
  hac <- identical(x$cv_method, "chisq")

  # A Wald statistic for the stacked and sup tests, a t statistic for the shape tests.
  title <- switch(x$type,
    stacked = "Stacked Wald test",
    sup = "Uniform sup Wald test",
    homogeneity = "Homogeneity test",
    monotonicity_adjacent = "Monotonicity test, adjacent differences",
    monotonicity_allpair = "Monotonicity test, all pairwise differences")
  stat_name <- if (x$type %in% c("stacked", "sup")) "W" else "t"

  # With no simulated exceedances, print p < 1 / n_cv to indicate simulation resolution.
  # Small positive p-values use scientific notation.
  p_text <- if (!is.na(x$n_cv) && x$p_value == 0) {
    sprintf("p-val < %s", formatC(1 / x$n_cv, format = "fg"))
  } else if (x$p_value < 1e-4) {
    sprintf("p-val = %.1e", x$p_value)
  } else {
    sprintf("p-val = %.4f", x$p_value)
  }
  rule <- strrep("-", 74)

  cat(sprintf("\n%s, %s\n",
              if (hac) { if (kernel) "Small-b Bartlett kernel HAC" else "Large-K orthonormal series HAR" }
              else sprintf("HAR %s", if (kernel) "Bartlett kernel" else "orthonormal series"), title))
  # A uniform null holds over a range and is evaluated on the grid printed on the next line;
  # a stacked null names its quantiles in the hypothesis itself.
  cat("H0:",
      paste(x$hypothesis, collapse = ";  "),
      if (identical(x$type, "sup")) sprintf("for all tau in %s", x$tau_range),
      if (!identical(x$type, "stacked"))
        sprintf("\ntau grid { %s }", paste(format(x$tau), collapse = ", ")), "\n")
  cat(rule, "\n", sep = "")
  cat(sprintf("%s-stat = %.4f,  %s,  T = %d\n", stat_name, x$stat, p_text, x$nobs))
  cat(sprintf("CV at the %g%% significance level = %.4f\n", 100 * x$alpha, x$cv))
  cat(rule, "\n", sep = "")
  # A number the user supplied is not data dependent. The rule is not named: the sup test-optimal
  # value is the stacked rule reused, not a rule derived for a maximum.
  origin <- if (is.numeric(x$smoothing_rule)) "user provided" else "data dependent"
  cat(if (kernel) sprintf("%s b = %.3f\n", origin, x$M_effective / x$nobs) else
      sprintf("%s K = %d\n", origin, x$K_effective))
  if (!is.null(x$t_ratio))
    cat(sprintf("t-stat of each %s difference:",
                if (identical(x$type, "monotonicity_allpair")) "pairwise" else "adjacent"),
        paste(sprintf("%.2f", x$t_ratio), collapse = "  "), "\n")
  # The terms of the sup statistic, so the quantile attaining the maximum is visible; with a single
  # restriction they come with the uniform confidence band, whose reference is the null value.
  if (!is.null(x$band)) {
    cat(sprintf("%g%% uniform confidence band\n",
                100 * (1 - x$alpha)))
    print(round(cbind(x$band, W = x$wald_tau), 3), row.names = FALSE)
  } else if (!is.null(x$wald_tau)) {
    cat("W at each tau:", paste(sprintf("%.2f", x$wald_tau), collapse = "  "), "\n")
  }
  invisible(x)
}

# Summary of a test result, in the usual form: summary() returns an object of class
# "summary.qrfstest" and its print method shows what print.qrfstest leaves out -- the estimates,
# the calibration of the critical value, and the nuisance smoothing of a sup or shape test.
summary.qrfstest <- function(object, ...) {
  structure(list(result = object), class = "summary.qrfstest")
}

print.summary.qrfstest <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  r <- x$result
  kernel <- r$har_method == "ker-har"
  hac <- identical(r$cv_method, "chisq")
  title <- switch(r$type,
    stacked = "Stacked Wald test",
    sup = "Uniform sup Wald test",
    homogeneity = "Homogeneity test",
    monotonicity_adjacent = "Monotonicity test, adjacent differences",
    monotonicity_allpair = "Monotonicity test, all pairwise differences")
  cat(sprintf("\n%s, %s\n",
              if (hac) { if (kernel) "Small-b Bartlett kernel HAC" else "Large-K orthonormal series HAR" }
              else sprintf("HAR %s", if (kernel) "Bartlett kernel" else "orthonormal series"), title))
  cat(sprintf("T = %d, tau grid { %s }\n", r$nobs, paste(format(r$tau), collapse = ", ")))
  cat("\nH0:", paste(r$hypothesis, collapse = ";  "),
      if (identical(r$type, "sup")) sprintf("for all tau in %s", r$tau_range), "\n")

  cat("\nCoefficients:\n")
  if (is.matrix(r$coefficients)) print(r$coefficients, digits = digits)
  else print(stats::setNames(r$coefficients, paste0("tau=", format(r$tau))), digits = digits)
  if (!is.null(r$band)) {
    cat(sprintf("\n%g%% uniform confidence band (Section 4.3), with W at each tau:\n",
                100 * (1 - r$alpha)))
    print(cbind(r$band, W = r$wald_tau), digits = digits, row.names = FALSE)
  } else if (!is.null(r$wald_tau)) {
    cat("\nW at each tau:\n")
    print(stats::setNames(r$wald_tau, paste0("tau=", format(r$tau))), digits = digits)
  }
  if (!is.null(r$t_ratio)) {
    cat(sprintf("\n%s differences and their t-statistics:\n",
                if (identical(r$type, "monotonicity_allpair")) "Pairwise" else "Adjacent"))
    print(rbind(difference = r$difference, t = r$t_ratio), digits = digits)
  }

  # The calibration, named as the paper names it.
  calibration <- if (hac) {
    sprintf("chi-square limit with %d degrees of freedom, %s asymptotics", r$q,
            if (kernel) "small-b" else "large-K")
  } else if (identical(r$type, "stacked")) {
    if (kernel) sprintf("fixed-b limit, eq. (18), simulated with %d draws", r$n_cv)
    else sprintf("fixed-K limit, eq. (19): scaled F(%d, %d)", r$q, r$K_effective - r$q + 1L)
  } else sprintf("Algorithm 1, simulated with %d draws", r$n_cv)
  p_text <- if (!is.na(r$n_cv) && r$p_value == 0) {
    sprintf("p-value < %s", formatC(1 / r$n_cv, format = "fg"))
  } else sprintf("p-value = %s", format(r$p_value, digits = digits))
  cat("\nTest:\n")
  cat(sprintf("  %s = %s, %s\n", if (r$type %in% c("stacked", "sup")) "W" else "t",
              format(r$stat, digits = digits), p_text))
  cat(sprintf("  critical value at the %g%% level = %s\n",
              100 * r$alpha, format(as.numeric(r$cv), digits = digits)))
  cat("  calibration:", calibration, "\n")

  cat("\nSmoothing:\n")
  origin <- if (is.numeric(r$smoothing_rule)) "user provided" else r$smoothing_rule
  if (kernel) cat(sprintf("  test M = %d (b = %.3f), %s\n", r$M_effective,
                          r$M_effective / r$nobs, origin))
  else cat(sprintf("  test K = %d, %s\n", r$K_effective, origin))
  if (kernel && !is.na(r$M_lambda)) cat(sprintf("  nuisance M = %d\n", r$M_lambda))
  if (!kernel && !is.na(r$K_lambda)) cat(sprintf("  nuisance K = %d\n", r$K_lambda))
  invisible(x)
}

# Draw the uniform confidence band of a sup test with a single restriction: the estimated
# combination at each quantile with the band of Section 4.3, and the null value as a dashed line.
# arg_result is the object qrfstest_sup returns, which carries the hypothesis the band is for.
# Several restrictions have no band: the dual of the joint statistic is an ellipse at each
# quantile, and the paper's Remark 4.2 rejects the projected rectangle as conservative.
uniform_CI_band <- function(arg_result, main = NULL, xlab = expression(tau), ylab = NULL,
                      ylim = NULL, ...) {
  x <- arg_result
  if (!inherits(x, "qrfstest")) stop("arg_result must be the object a qrfs test returns.", call. = FALSE)
  if (is.null(x$band)) stop("Only a sup test with a single restriction has a band to plot.", call. = FALSE)
  band <- x$band
  if (is.null(ylab)) {
    # R identifies the plotted combination, including coefficients on either side of '='.
    selected <- x$R[1L, ] != 0
    weights <- x$R[1L, selected]
    terms <- rownames(x$coefficients)[selected]
    terms <- ifelse(abs(weights) == 1, terms,
                    paste0(format(abs(weights), trim = TRUE), " * ", terms))
    ylab <- sub("^\\+ ", "", paste0(ifelse(weights < 0, "- ", "+ "), terms, collapse = " "))
  }
  if (is.null(ylim)) ylim <- range(band$lower, band$upper, x$null_val)
  if (is.null(main))
    main <- sprintf("%g%% uniform confidence band, HAR %s", 100 * (1 - x$alpha),
                    if (x$har_method == "ker-har") "Bartlett kernel" else "orthonormal series")
  plot(band$tau, band$estimate, type = "n", ylim = ylim,
       xlab = xlab, ylab = ylab, main = main, ...)
  graphics::polygon(c(band$tau, rev(band$tau)), c(band$lower, rev(band$upper)), col = "grey85", border = NA)
  graphics::lines(band$tau, band$estimate, lwd = 2)
  graphics::points(band$tau, band$estimate, pch = 16)
  graphics::abline(h = x$null_val, lty = 2)
  invisible(band)
}

# The plot method keeps the generic's arguments plot(x, y, ...); y has no meaning for one series.
plot.qrfstest <- function(x, y, ...) {
  uniform_CI_band(x, ...)
  invisible(x)
}

# The estimates a result carries: the p x m coefficient matrix of a Wald test, or the tested
# coefficient at each quantile for a shape test.
coef.qrfstest <- function(object, ...) object$coefficients

# The simulated critical value and p-value of a shape statistic: Algorithm 1 of the paper applied to
# the tested differences, whose projection arg_RH carries. The transformation of a draw must be the
# one the observed statistic used, so the caller states it: |t| two-sided, or arg_alternative_sign times t
# one-sided, with arg_alternative_sign +1 when a rise supports the alternative and -1 when a fall does.
.qrfs_shape_cv <- function(arg_qr_fit, arg_RH, arg_har_method, arg_M, arg_K,
                           arg_M_lambda, arg_K_lambda, arg_stat, arg_two_sided,
                           arg_alternative_sign, arg_n_cv, arg_n_r, arg_coverage, arg_alpha) {
  kernel <- arg_har_method == "ker-har"
  nobs <- arg_qr_fit$nobs
  J_nuisance <- compute_Jhat_stacked(arg_qr_fit, arg_har_method, arg_M_lambda, arg_K_lambda)
  eig_trunc <- truncate_lambda_eigendecomp(J_nuisance / nobs, 1L, arg_coverage)
  # Project the Gaussian factors onto the tested differences once.
  eig_trunc$sqrt_eigs_eigvecs <- eig_trunc$sqrt_eigs_eigvecs %*% t(arg_RH)
  # The Brownian grid belongs to the kernel path; the series draws are finite-dimensional.
  if (kernel) r_grid <- seq(0, 1, length.out = arg_n_r + 1L)
  draws <- numeric(arg_n_cv)
  for (sim in seq_len(arg_n_cv)) {
    if (kernel) {  # Endpoint and Brownian bridge.
      path <- simulate_gaussian_path(eig_trunc, arg_n_r, 1 / arg_n_r)
      endpoint <- path[arg_n_r + 1L, ]
      bridge <- path - r_grid %o% endpoint
      variance <- compute_P_functional(bridge, arg_M / nobs, 1 / arg_n_r, arg_n_r,
                                       arg_diagonal = TRUE)
    } else {  # Independent Gaussian endpoint and K projections.
      endpoint <- as.numeric(stats::rnorm(eig_trunc$num_eigs) %*% eig_trunc$sqrt_eigs_eigvecs)
      normals <- matrix(stats::rnorm(arg_K * eig_trunc$num_eigs), nrow = arg_K)
      variance <- colSums((normals %*% eig_trunc$sqrt_eigs_eigvecs)^2) / arg_K
    }
    draw <- endpoint / sqrt(variance)
    draws[sim] <- max(if (arg_two_sided) abs(draw) else arg_alternative_sign * draw)
  }
  cv <- as.numeric(stats::quantile(draws, 1 - arg_alpha, type = 1))
  names(cv) <- paste0("alpha=", arg_alpha)
  list(cv = cv, p_value = mean(draws >= arg_stat))
}
# Shape tests of one coefficient over the quantile grid, on the t scale: each contrast is a single
# restriction, so its t-ratio is reported rather than its square. The paper, Section 4.4,
# Proposition 4.2 and Algorithm 3.
# arg_test: "homogeneity" -- the coefficient is constant; t_H = max |t_i| over the adjacent
#   differences, two-sided. "monotonicity_adjacent" -- one-sided over the same differences, as in
#   the paper's simulations and Section 4.5 (R_M = R_H). "monotonicity_allpair" -- one-sided over
#   every pair of quantiles. With the default "decreasing" alternative, the statistic
#   is -t_M, the negative of the paper's minimum in Section 4.4; large values reject.
# arg_alternative names what a rejection of a monotonicity test would support, as in stats::t.test:
#   "increasing" tests the non-increasing null and rejects when some step rises, "decreasing" the
#   reverse. The homogeneity test is two-sided and ignores it, with a warning if it is given.
# arg_tau_grid takes the same three forms as in qrfstest_sup.
qrfstest_shape <- function(arg_formula, arg_data, arg_variable, arg_har_method,
                           arg_test,
                           arg_alternative = c("decreasing", "increasing"),
                           arg_tau_grid = NULL,
                           arg_M = "test_optimal", arg_K = "test_optimal",
                           arg_M_lambda = "mse_optimal_ar",
                           arg_K_lambda = "mse_optimal_var",
                           arg_alpha = 0.05, arg_kappa = NULL,
                           arg_threshold_method = "theoretical",
                           arg_var_fit = "gaussian",
                           arg_var_stationary = "andrews_monahan",
                           arg_coverage = 0.99,
                           arg_hinv_method = "ker", arg_rq_method = "br",
                           arg_seed = NULL) {
  if (!is.null(arg_seed)) set.seed(arg_seed)
  if (missing(arg_test))
    stop("arg_test must be \"homogeneity\", \"monotonicity_adjacent\" or \"monotonicity_allpair\".", call. = FALSE)
  arg_test <- match.arg(arg_test, c("homogeneity", "monotonicity_adjacent", "monotonicity_allpair"))
  arg_har_method <- match.arg(arg_har_method, c("ker-har", "os-har"))
  arg_threshold_method <- match.arg(arg_threshold_method, c("theoretical", "empirical"))
  arg_var_fit <- match.arg(arg_var_fit, c("gaussian", "ols"))
  kernel <- arg_har_method == "ker-har"
  # Two-sided for homogeneity: under the null the coefficient is constant, and a departure in either
  # direction counts, so the statistic is max |t_i| -- nonnegative, and large when any pair differs.
  # One-sided for monotonicity: only a step in the direction of the alternative counts -- a rise
  # under "increasing", a fall under "decreasing" -- and the sign turns that direction into a
  # positive term, so the statistic is nonpositive when no step supports the alternative.
  two_sided <- arg_test == "homogeneity"
  if (two_sided && !missing(arg_alternative))
    warning("The homogeneity test is always two-sided; arg_alternative is ignored.")
  arg_alternative <- match.arg(arg_alternative)
  alternative_sign <- if (arg_alternative == "increasing") 1 else -1

  tau_spec <- .qrfs_tau_grid(arg_tau_grid)
  arg_tau_grid <- tau_spec$grid
  range_label <- tau_spec$range_label
  qr_fit <- fit_qr_grid(arg_formula, arg_data, arg_tau_grid, arg_hinv_method, arg_rq_method)
  taus <- qr_fit$fit$tau
  m <- qr_fit$m
  p <- qr_fit$p
  nobs <- qr_fit$nobs
  k <- match(arg_variable, rownames(qr_fit$coef_mat))
  if (is.na(k))
    stop("Unknown coefficient '", arg_variable, "'. The model has: ",
         paste(rownames(qr_fit$coef_mat), collapse = ", "), ".", call. = FALSE)
  # Each difference is one restriction, of rank 1.
  .qrfs_check_inputs(arg_har_method, arg_alpha, if (kernel) arg_M else arg_K,
                     if (kernel) arg_M_lambda else arg_K_lambda, arg_coverage,
                     nobs, 1L, arg_kappa)

  # The quantile pairs compared, lower index first: adjacent quantiles, or every pair.
  pairs <- if (arg_test == "monotonicity_allpair") utils::combn(m, 2L) else
    rbind(seq_len(m - 1L), seq_len(m - 1L) + 1L)
  R <- matrix(0, ncol(pairs), m * p)
  for (i in seq_len(ncol(pairs))) {
    R[i, (pairs[1L, i] - 1L) * p + k] <- -1
    R[i, (pairs[2L, i] - 1L) * p + k] <- 1
  }
  RH <- R %*% as.matrix(Matrix::bdiag(qr_fit$Hinv_list))

  # The M / K rules are evaluated at the adjacent differences for every test (Section 4.5,
  # R_M = R_H): the all-pair differences are their linear combinations and add no restriction.
  R_adjacent <- R[pairs[2L, ] - pairs[1L, ] == 1L, , drop = FALSE]
  selected <- .qrfs_select_M_K(qr_fit, arg_har_method,
    if (kernel) arg_M else arg_K, if (kernel) arg_M_lambda else arg_K_lambda,
    R_adjacent, 5L, arg_alpha, arg_kappa, arg_threshold_method, arg_var_fit, arg_var_stationary)
  M <- selected$M; K <- selected$K
  M_lambda <- selected$M_lambda; K_lambda <- selected$K_lambda

  tested_scores <- qr_fit$score_mat %*% t(RH)
  V_hat <- .compute_Jhat_from_scores(tested_scores, arg_har_method, M, K, NULL)
  difference <- as.numeric(R %*% as.vector(qr_fit$coef_mat))
  t_ratio <- difference / sqrt(diag(V_hat))
  names(difference) <- names(t_ratio) <-
    paste0(format(taus)[pairs[2L, ]], "-", format(taus)[pairs[1L, ]])
  stat <- max(if (two_sided) abs(t_ratio) else alternative_sign * t_ratio)

  # Critical-value simulation: draws and grid intervals per Brownian path. Fixed
  # here rather than offered as arguments; they set its numerical accuracy, not the test.
  n_cv <- 10000L; n_r <- 1000L
  cv_result <- .qrfs_shape_cv(qr_fit, RH, arg_har_method, M, K, M_lambda, K_lambda,
                                arg_stat = stat, arg_two_sided = two_sided,
                                arg_alternative_sign = alternative_sign,
                                arg_n_cv = n_cv, arg_n_r = n_r,
                                arg_coverage = arg_coverage, arg_alpha = arg_alpha)
  cv <- cv_result$cv

  out <- list(
    hypothesis = sprintf("%s %s in tau over %s", arg_variable,
      if (two_sided) "is constant" else if (arg_alternative == "decreasing")
        "is non-decreasing" else "is non-increasing", range_label),
    coefficients = qr_fit$coef_mat[k, ], stat = stat, cv = cv,
    p_value = cv_result$p_value,
    reject = as.integer(stat > cv),
    M_effective = M, K_effective = K, M_lambda = M_lambda, K_lambda = K_lambda,
    smoothing_rule = if (kernel) arg_M else arg_K, tau_range = range_label,
    q = nrow(R),
    type = arg_test,
    har_method = arg_har_method, alpha = arg_alpha, n_cv = n_cv, nobs = nobs,
    tau = taus, difference = difference, t_ratio = t_ratio)
  class(out) <- "qrfstest"
  out
}

# Pivotal Bartlett CV and optional p-value for stacked Wald inference.
# The paper, Section 4.2, eq. (18): the fixed-b limit of Kiefer & Vogelsang (2005).
# q is the restriction count; b = M/T; arg_stat is the observed Wald statistic.
# arg_n_cv draws use arg_n_r grid intervals.
compute_pivotal_cv_stacked <- function(q, b, alpha,
                                        arg_n_cv = 10000, arg_n_r=1000,
                                        arg_stat = NULL) {
  nr1 <- arg_n_r + 1
  r_grid <- seq(0, 1, length.out=nr1)
  dr <- 1 / arg_n_r

  vals <- numeric(arg_n_cv)
  for (sim in seq_len(arg_n_cv)) {
    BM <- matrix(0, nrow=nr1, ncol=q)
    Z_inc <- matrix(stats::rnorm(arg_n_r * q), nrow=arg_n_r, ncol=q)
    BM[2:nr1, ] <- apply(Z_inc * sqrt(dr), 2, cumsum)

    W1 <- BM[nr1, ]
    W_tilde <- BM - r_grid %o% W1

    P_mat <- compute_P_functional(W_tilde, b, dr, arg_n_r, arg_diagonal = FALSE)

    vals[sim] <- compute_wald(W1, P_mat)
  }

  # type = 1, the inverse empirical CDF, as in every other simulated critical value here.
  cv_out <- as.numeric(stats::quantile(vals, 1 - alpha, type = 1))
  names(cv_out) <- paste0("alpha=", alpha)
  list(cv = cv_out,
       p_value = if (is.null(arg_stat)) NA_real_ else
         mean(vals >= arg_stat))
}

# Run and print Adrian examples; return results invisibly.
# arg_test: 'stacked', 'sup', or 'shape'; arg_har_method picks HAR.
# A non-NULL arg_seed is passed to every command as its own arg_seed, so each printed command
# reproduces its result by itself; NULL leaves the commands unseeded (the current random-number stream).
# arg_var_fit is passed to every command: "gaussian" uses Galvao & Yoon's approximation;
# "ols" fits AR(1) or VAR(1) directly to the scores, depending on the smoothing rule.
qrfs_example <- function(arg_test = c("stacked", "sup", "shape"),
                         arg_har_method = "ker-har", arg_var_fit = "gaussian",
                         arg_seed = NULL) {
  arg_test <- match.arg(arg_test)
  har <- arg_har_method
  vf <- arg_var_fit
  # The Adrian, Boyarchenko and Giannone (2019) data shipped with the package (data/adrian.rda).
  env <- new.env(); utils::data("adrian", package = "QRfs", envir = env); data <- env$adrian
  commands <- switch(arg_test,
    stacked = list(
      bquote(qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, data,
                              "NFCI[0.1] = NFCI[0.5] = NFCI[0.9]", .(har), arg_var_fit = .(vf))),
      bquote(qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, data,
                              "ggdp[0.1] = ggdp[0.5] = ggdp[0.9]", .(har), arg_var_fit = .(vf))),
      bquote(qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, data,
                              c("NFCI[0.1] = NFCI[0.5] = NFCI[0.9]",
                                "ggdp[0.1] = ggdp[0.5] = ggdp[0.9]"), .(har), arg_var_fit = .(vf)))),
    sup = list(
      bquote(qrfstest_sup(ggdp_h1 ~ ggdp + NFCI, data, "ggdp = 0", .(har), arg_var_fit = .(vf))),
      bquote(qrfstest_sup(ggdp_h1 ~ ggdp + NFCI, data, "NFCI = 0", .(har), arg_var_fit = .(vf),
                          arg_tau_grid = c(0.05, 0.1, 0.15, 0.2, 0.25)))),
    shape = list(
      bquote(qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, data, "ggdp", .(har), arg_var_fit = .(vf),
                            arg_test = "homogeneity")),
      bquote(qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, data, "NFCI", .(har), arg_var_fit = .(vf),
                            arg_test = "monotonicity_adjacent")),
      bquote(qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, data, "NFCI", .(har), arg_var_fit = .(vf),
                            arg_test = "monotonicity_adjacent", arg_alternative = "increasing")),
      bquote(qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, data, "NFCI", .(har), arg_var_fit = .(vf),
                            arg_test = "monotonicity_allpair", arg_alternative = "increasing"))))

  if (!is.null(arg_seed))
    commands <- lapply(commands, function(command) { command$arg_seed <- arg_seed; command })

  cat("Adrian, Boyarchenko and Giannone (2019): GDP growth in quarter t + 1 on GDP growth\n",
      "and the NFCI in quarter t, 1973Q1-2015Q3, T = ", nrow(data), ".\n",
      "> data <- adrian\n", sep = "")
  results <- vector("list", length(commands))
  for (i in seq_along(commands)) {
    cat("\n> ", paste(deparse(commands[[i]], width.cutoff = 70), collapse = "\n  "), "\n", sep = "")
    results[[i]] <- eval(commands[[i]])
    print(results[[i]])
    # Print estimates only when they change.
    estimates <- round(results[[i]]$coefficients, 3)
    if (i == 1 || !identical(estimates, shown)) {
      cat("estimates:\n")
      print(estimates)
    }
    shown <- estimates
  }
  invisible(results)
}
