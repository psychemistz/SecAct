## SecAct (development version)

### Output-changing fixes

* **`sigMatrix = "CytoSig"` now uses its documented lambda (10000).**
  `SecAct.activity.inference()`, `.ST()`, `.scRNAseq()` and
  `SecAct.CCC.scRNAseq()` defaulted to `lambda = 5e+05`, so the
  per-signature default in `load_sig_matrix()` never applied and CytoSig
  ran with 5e+05. The default is now `lambda = NULL`: 10000 for CytoSig,
  5e+05 for SecAct, SecAct-XXX and custom signature files. To reproduce
  earlier CytoSig outputs pass `lambda = 5e+05`. SecAct outputs are
  unchanged.
* **Signature grouping (`is.group.sig = TRUE`) is computed on the genes
  shared by the input and the signature matrix**, not on the full
  signature gene set. Groups, and therefore activities, can change for
  inputs that do not cover the full panel.
* Inputs with NA: overlapping genes with NA in any sample are dropped
  from all samples, with a warning. Previously every signature of an
  affected sample was NA.
* Duplicated gene symbols in `inputProfile` / `inputProfile_control` are
  deduplicated as in the ST/scRNA-seq paths (keep the row with the
  largest sum), with a warning. Previously the first row was used
  silently.

### New options (defaults keep earlier behavior)

* `lambda_scale = c("none", "overlap")`: `"overlap"` multiplies lambda by
  `(n_overlap - 1) / (n_ref - 1)`, `n_ref` being the number of genes in
  the signature matrix, so shrinkage does not depend on panel size.
* `pvalue = c("permutation", "normal")`: `"normal"` returns
  `2 * pnorm(-abs(zscore))`. The permutation p-value is floored at
  `1 / (nrand + 1)`; with `nrand = 1000` every |z| above about 3.3 ties
  near 0.001, so FDR on it cannot rank strong hits.

### Other

* A zero-variance sample (over the overlapping genes) now warns; its
  outputs are NaN as before.
* A custom signature file with `lambda = NULL` now uses 5e+05 instead of
  failing with an obscure error; a non-numeric lambda is rejected.
* Documented that the permutation null SD is the population SD (divisor
  `nrand`), as in SecAct 1.0.0; unchanged.

## SecAct 1.1.0

*2026-07-28*

* Add a new tutorial for analyzing Visium HD data.

## SecAct 1.0.0

*2024-08-03*

* Release SecAct package.
