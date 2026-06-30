# Branch Output Format

This note summarizes the output schemas for the current MethylRegion branch
implementations. It was checked against the source files under `R/` and with
small example calls from `tests/api.R`.

## Wrapper Dispatch

| Entry point | Arguments | Implementation branch |
| --- | --- | --- |
| `mr_bi()` | `data.type = "independent"`, `cov.mod = NULL` | `dmr_case_control()` |
| `mr_bi()` | `data.type = "independent"`, `cov.mod` non-`NULL` | `dmr_case_control_cov()` |
| `mr_bi()` | `data.type = "paired"`, `cov.mod = NULL` | `dmr_paired()` |
| `mr_bi()` | `data.type = "paired"`, `cov.mod` non-`NULL` | `dmr_longitudinal()` |
| `mr_bi()` | `data.type = "longitudinal"` | `dmr_longitudinal()` |
| `mr_continuous()` | `data.type = "independent"` | `amr_continuous()` |
| `mr_continuous()` | `data.type = "longitudinal"` | `amr_longitudinal()` |

The dispatch depends only on `data.type` and whether `cov.mod` is `NULL`. It does
not depend on the balance of 0/1 labels in `y`. If the first column of `y`
contains more than two groups, `mr_bi()` runs each group pair through the same
binary branch after subsetting samples and recoding that pair to 0 and 1.

## Output Columns By Branch

| Branch | Columns, in return order |
| --- | --- |
| `dmr_case_control()` | `chr`, `start`, `end`, `N.CpGs`, `ks_stat`, `mean_diff`, `p_value`, `methX`, `methY`, `e_value`, `FDR`, `e_adjust` |
| `dmr_case_control_cov()` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_lm_group`, `p_value`, `methX`, `methY`, `FDR` |
| `dmr_paired()` | `chr`, `start`, `end`, `N.CpGs`, `ks_stat`, `mean_diff`, `p_value`, `methX`, `methY`, `FDR` |
| `dmr_longitudinal()` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_lmm`, `p_value`, `methX`, `methY`, `FDR` |
| `amr_continuous()` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_lm`, `p_value`, `methX`, `methY`, `FDR` |
| `amr_longitudinal()` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_meth`, `p_value`, `methX`, `methY`, `FDR` |

Only `dmr_case_control()` returns e-value columns (`e_value`, `e_adjust`).
The other branches do not currently return `e_value`, `e_adjust`, or
`e_bh_significant`.
Multi-group `mr_bi()` results prepend `comparison`, `groupA`, and `groupB` to
the corresponding branch schema, with adjusted statistics computed separately
inside each pairwise comparison. `attr(result, "comparisons")` records all
attempted pair labels, including pairs with no reported regions.

## Column Meaning

| Column | Meaning |
| --- | --- |
| `comparison` | Pairwise group label such as `A_vs_B` for multi-group `mr_bi()` output. |
| `groupA` | First group in a multi-group pairwise `mr_bi()` comparison, recoded to 0. |
| `groupB` | Second group in a multi-group pairwise `mr_bi()` comparison, recoded to 1. |
| `chr` | Chromosome label copied from `input_dat$chr`. |
| `start` | Region start coordinate, computed as first CpG `pos - 1`. |
| `end` | Region end coordinate, computed as last CpG `pos`. |
| `N.CpGs` | Number of CpGs in the reported region. |
| `ks_stat` | 2D-KS statistic, used by non-covariate binary DMR branches. |
| `mean_diff` | Mean methylation difference, usually test/case minus control. |
| `cor_est` | Correlation estimate used by regression/model-based branches. |
| `coef_lm_group` | Group coefficient from the covariate-adjusted case-control model. |
| `coef_lmm` | Phenotype coefficient from the longitudinal binary model. |
| `coef_lm` | Phenotype coefficient from the independent continuous model. |
| `coef_meth` | Methylation coefficient from the longitudinal continuous model. |
| `p_value` | Branch-level region p value. |
| `methX` | Mean methylation value across region methylation entries. |
| `methY` | Mean phenotype/group value from `y`. |
| `e_value` | Evidence-value statistic from `dmr_case_control()`. |
| `FDR` | Benjamini-Hochberg adjusted `p_value`. |
| `e_adjust` | Benjamini-Hochberg style adjustment of `e_value`. |

## Empty Result Behavior

When no candidate region is found, outputs are not guaranteed to preserve the
full branch schema. The public wrapper `.methylregion_result()` converts `NULL`
or non-data-frame empty results to `data.frame()`, but a branch may also return
a 0-row data frame after partially adding post-processing columns.

Observed with the small `tests/api.R`-style data and `controlist = list(mincpgs
= 2, trend = 0)`:

| Branch | Observed empty output shape |
| --- | --- |
| `dmr_case_control()` | 0 rows with `FDR`, `e_adjust` only |
| `dmr_case_control_cov()` | 0 rows and 0 columns |
| `amr_continuous()` | 0 rows and 0 columns |
| `amr_longitudinal()` | 0 rows and 0 columns |

Callers that need stable schemas should normalize empty outputs explicitly.

## Source Locations Checked

| Branch | Source lines for final column naming |
| --- | --- |
| `dmr_case_control()` | `R/impl-dmr_case_control.R`, final `colnames(...)`, `FDR`, `e_adjust`, and `return(...[, -4])` block |
| `dmr_case_control_cov()` | `R/impl-dmr_case_control_cov.R`, final `colnames(...)`, `FDR`, and `return(...[, -4])` block |
| `dmr_paired()` | `R/impl-dmr_paired.R`, final `colnames(...)`, `FDR`, and `return(...[, -4])` block |
| `dmr_longitudinal()` | `R/impl-dmr_longitudinal.R`, final `colnames(...)`, `FDR`, and `return(...[, -4])` block |
| `amr_continuous()` | `R/impl-amr_continuous.R`, final `colnames(...)`, `FDR`, and `return(...[, -4])` block |
| `amr_longitudinal()` | `R/impl-amr_longitudinal.R`, final `colnames(...)`, `FDR`, and `return(...[, -4])` block |
