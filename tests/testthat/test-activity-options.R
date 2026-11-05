sig_fixture <- function(n_ref = 200, n_in = 60, p = 6, m = 4, seed = 7L) {
  set.seed(seed)
  X <- matrix(rnorm(n_ref * p), n_ref, p,
              dimnames = list(paste0("g", seq_len(n_ref)), paste0("sig", seq_len(p))))
  # sig1 and sig2 agree on the input genes only
  X[seq_len(n_in), "sig2"] <- X[seq_len(n_in), "sig1"] + rnorm(n_in, sd = 0.05)
  f <- tempfile(fileext = ".tsv")
  write.table(X, f, sep = "\t", quote = FALSE)
  Y <- matrix(rnorm(n_in * m), n_in, m,
              dimnames = list(paste0("g", seq_len(n_in)), paste0("s", seq_len(m))))
  list(file = f, X = X, Y = Y)
}

run_act <- function(d, Y = d$Y, ...) {
  SecAct.activity.inference(Y, is.differential = TRUE, sigMatrix = d$file,
                            nrand = 30, backend = "cpu-pure", ...)
}

test_that("lambda defaults to NULL so per-signature defaults apply", {
  for (f in list(SecAct.activity.inference, SecAct.activity.inference.ST,
                 SecAct.activity.inference.scRNAseq, SecAct.CCC.scRNAseq)) {
    expect_null(formals(f)$lambda)
  }
  d <- sig_fixture()
  expect_equal(SecAct:::load_sig_matrix(d$file, NULL)$lambda, 5e5)
  expect_equal(run_act(d, is.group.sig = FALSE), run_act(d, is.group.sig = FALSE, lambda = 5e5))
  expect_error(SecAct:::load_sig_matrix(d$file, "a"), "single non-negative number")
})

test_that("signatures are grouped on the overlapping genes", {
  d <- sig_fixture()
  res <- run_act(d, lambda = 100)
  expect_identical(res$beta["sig1", ], res$beta["sig2", ])
  expect_false(identical(res$beta["sig1", ], res$beta["sig3", ]))
})

test_that("lambda_scale = 'overlap' scales lambda by (n_overlap-1)/(n_ref-1)", {
  d <- sig_fixture()
  lam <- 1000
  res <- run_act(d, is.group.sig = FALSE, lambda = lam, lambda_scale = "overlap")
  X <- scale(d$X[rownames(d$Y), ]); Y <- scale(d$Y)
  ref <- SecAct:::.ridge_dispatch(X, Y, lam * 59 / 199, 30, backend = "cpu-pure")
  expect_equal(res$beta, ref$beta)
  expect_equal(res$zscore, ref$zscore)
  expect_equal(run_act(d, is.group.sig = FALSE, lambda = lam),
               run_act(d, is.group.sig = FALSE, lambda = lam, lambda_scale = "none"))
  expect_error(run_act(d, lambda_scale = "bogus"))
})

test_that("pvalue = 'normal' uses 2*pnorm(-|z|); permutation p is floored", {
  d <- sig_fixture()
  perm <- run_act(d, is.group.sig = FALSE, lambda = 100)
  norm <- run_act(d, is.group.sig = FALSE, lambda = 100, pvalue = "normal")
  expect_equal(norm$zscore, perm$zscore)
  expect_equal(norm$pvalue, 2 * pnorm(-abs(perm$zscore)))
  expect_gte(min(perm$pvalue), 1 / 31)
})

test_that("NA genes are dropped with a warning", {
  d <- sig_fixture()
  Yna <- d$Y; Yna[c(3, 9), 2] <- NA
  expect_warning(res <- run_act(d, Y = Yna, is.group.sig = FALSE, lambda = 100),
                 "2 of 60 overlapping genes have NA")
  expect_false(anyNA(res$zscore))
  expect_equal(res, run_act(d, Y = d$Y[-c(3, 9), ], is.group.sig = FALSE, lambda = 100))
})

test_that("duplicated gene symbols are deduplicated with a warning", {
  d <- sig_fixture()
  Yd <- rbind(d$Y, d$Y[5, , drop = FALSE] + 10)
  expect_warning(res <- run_act(d, Y = Yd, is.group.sig = FALSE, lambda = 100),
                 "1 duplicated gene symbol")
  Ymax <- Yd[-5, ]  # the kept row stays in its original position
  expect_equal(res, run_act(d, Y = Ymax, is.group.sig = FALSE, lambda = 100))
})

test_that("zero-variance samples warn", {
  d <- sig_fixture()
  Yz <- d$Y; Yz[, 3] <- 1
  expect_warning(run_act(d, Y = Yz, is.group.sig = FALSE, lambda = 100), "zero variance")
})
