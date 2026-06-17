#' Pure-R ridge + permutation kernel
#'
#' Solves \eqn{\beta = (X'X + \lambda I)^{-1} X'Y} via Cholesky, then
#' estimates null via Y-row permutations. Precomputes the projection
#' \eqn{T = (X'X + \lambda I)^{-1} X'} once so each permutation is a
#' single p-by-m dgemm instead of two triangular solves plus a crossprod.
#'
#' @param X Numeric matrix, n x p, column-scaled.
#' @param Y Numeric matrix, n x m, column-scaled.
#' @param lambda Ridge penalty.
#' @param nrand Number of permutations.
#' @return list(beta, se, zscore, pvalue) — each a p x m matrix.
#' @keywords internal
.ridge_pureR <- function(X, Y, lambda, nrand) {
  n <- nrow(Y); p <- ncol(X); m <- ncol(Y)

  A <- crossprod(X) + lambda * diag(p)
  R <- chol(A)
  Rt <- t(R)
  T_proj <- backsolve(R, forwardsolve(Rt, t(X)))   # p x n
  beta <- T_proj %*% Y

  perm_table <- .gsl_mt19937_perm_table(n, nrand)

  aver <- matrix(0, p, m)
  aver_sq <- matrix(0, p, m)
  pvalue_count <- matrix(0, p, m)

  for (i in seq_len(nrand)) {
    perm <- perm_table[i, ] + 1L
    beta_rand <- T_proj %*% Y[perm, , drop = FALSE]
    aver <- aver + beta_rand
    aver_sq <- aver_sq + beta_rand^2
    pvalue_count <- pvalue_count + (abs(beta_rand) >= abs(beta))
  }

  aver <- aver / nrand
  aver_sq <- aver_sq / nrand
  se <- sqrt(aver_sq - aver^2)
  zscore <- (beta - aver) / se
  pvalue <- (pvalue_count + 1) / (nrand + 1)

  dn <- list(colnames(X), colnames(Y))
  dimnames(beta)   <- dn
  dimnames(se)     <- dn
  dimnames(zscore) <- dn
  dimnames(pvalue) <- dn

  list(beta = beta, se = se, zscore = zscore, pvalue = pvalue)
}

#' Resolve backend based on installed FlashReg and GPU availability
#'
#' Legacy backend names ("gpu", "cpu-fast", "cpu-pure") are still
#' accepted and remap onto the corresponding FlashReg backends so
#' downstream code keeps working without changes.
#' @keywords internal
.resolve_backend <- function(backend = c("auto", "gpu", "cpu-fast", "cpu-pure")) {
  backend <- match.arg(backend)
  if (backend != "auto") return(backend)

  if (requireNamespace("FlashReg", quietly = TRUE)) {
    gpu_ok <- tryCatch(FlashReg::cuda_available(), error = function(e) FALSE)
    if (isTRUE(gpu_ok)) return("gpu")
    return("cpu-fast")
  }
  "cpu-pure"
}

# Map SecAct's legacy backend name onto FlashReg's backend name.
.secact_to_flashreg_backend <- function(chosen) {
  switch(chosen,
         gpu        = "cuda_native",
         `cpu-fast` = "omp",
         `cpu-pure` = "pure_r",
         stop("internal: unknown SecAct backend label '", chosen, "'"))
}

#' Dispatch ridge+permutation call to FlashReg or the in-house pure-R loop
#'
#' Picks GPU (\code{FlashReg::ridge(backend="cuda_native")}) > CPU-fast
#' (\code{FlashReg::ridge(backend="omp")}) > CPU-pure (this package's
#' \code{.ridge_pureR}). FlashReg replaces the historical pair of
#' optional packages (RidgeFast for CPU, RidgeCuda for GPU); legacy
#' backend names are kept as aliases for backward compatibility.
#'
#' With \code{rng_method="mt19937"} and \code{ncores=1}, FlashReg's
#' backends are bit-identical to the pure-R loop.
#'
#' @keywords internal
.ridge_dispatch <- function(X, Y, lambda, nrand,
                            ncores = 1L,
                            rng_method = "mt19937",
                            backend = c("auto", "gpu", "cpu-fast", "cpu-pure")) {
  rng_method <- match.arg(tolower(rng_method), c("mt19937", "srand"))
  chosen <- .resolve_backend(backend)
  if (isTRUE(getOption("SecAct.verbose", FALSE))) {
    message("[SecAct] ridge backend: ", chosen, " (rng_method=", rng_method,
            ", ncores=", ncores, ")")
  }

  if (chosen %in% c("gpu", "cpu-fast")) {
    if (!requireNamespace("FlashReg", quietly = TRUE)) {
      stop("backend='", chosen, "' requires the FlashReg package. ",
           "Install it from https://github.com/data2intelligence/FlashReg ",
           "or set backend='cpu-pure' to use the in-house pure-R loop.")
    }
    FlashReg::ridge(X = X, Y = Y, lambda = lambda, nrand = nrand,
                    backend = .secact_to_flashreg_backend(chosen),
                    ncores = ncores, rng_method = rng_method)
  } else {
    if (rng_method != "mt19937") {
      stop("rng_method='", rng_method, "' requires FlashReg; ",
           "pure-R supports only 'mt19937'.")
    }
    .ridge_pureR(X, Y, lambda, nrand)
  }
}
