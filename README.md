# MethylRegion

**MethylRegion: A unified framework for differential methylation regions across designs, traits, and covariates**

MethylRegion consolidates the AMRfinder development branches into one R package
with a small, design-oriented API. The recommended entry points are
`mr_bi()` for binary phenotypes and `mr_continuous()` for continuous
phenotypes. The original branch-specific functions remain exported for direct
use and backward compatibility.

## Functions

### Main entry points

| Function | Phenotype | `data` values | Covariate dispatch |
| --- | --- | --- | --- |
| `mr_bi()` | Binary | `"independent"`, `"paired"`, `"longitudinal"` | Uses `is.null(cov.mod)` to choose no-covariate versus covariate-aware branches where implemented. |
| `mr_continuous()` | Continuous | `"independent"`, `"longitudinal"` | Passes `cov.mod` through; it may be `NULL` or a covariate data frame. |

For `mr_bi()`:

| `data` | `cov.mod` | Computing function |
| --- | --- | --- |
| `"independent"` | `NULL` | `dmr_case_control()` |
| `"independent"` | non-`NULL` | `dmr_case_control_cov()` |
| `"paired"` | `NULL` | `dmr_paired()` |
| `"paired"` | non-`NULL` | Not implemented yet; throws an explicit error. |
| `"longitudinal"` | `NULL` or non-`NULL` | Not implemented yet; throws an explicit error. |

For `mr_continuous()`:

| `data` | `cov.mod` | Computing function |
| --- | --- | --- |
| `"independent"` | `NULL` or non-`NULL` | `amr_continuous()` |
| `"longitudinal"` | `NULL` or non-`NULL` | `amr_longitudinal()` |

### Branch-specific functions

| Function | Design | Trait | Covariates | Source branch |
| --- | --- | --- | --- | --- |
| `dmr_case_control()` | Independent case-control | Binary | No | `dmr.no.cov.2dks-mwu` |
| `dmr_case_control_cov()` | Independent case-control | Binary | Yes | `dmr.with.cov.lm` |
| `dmr_paired()` | Pre/post treatment or paired case-control | Binary | No | `dmr.no.cov.2dks-paired-wilcox` |
| `amr_continuous()` | Independent EWAS / phenotype association | Continuous | Yes/No | `main` |
| `amr_longitudinal()` | Longitudinal or repeated-measure AMR | Continuous | Yes/No | `dmr.with.cov.lm.phenotype-model` |

`amr_longitudinal()` currently supports continuous phenotypes only. It does not
support discrete or binary phenotypes.

## Installation

```r
remotes::install_github("godot4k/MethylRegion")
```

For local development:

```sh
R CMD INSTALL MethylRegion
```

## Input

The methylation matrix is supplied as a data frame with genomic columns followed
by one methylation column per sample:

```r
input_dat <- data.frame(
  chr = rep("chr1", 6),
  pos = seq(100, 600, by = 100),
  sample1 = c(0.1, 0.2, 0.2, 0.7, 0.8, 0.8),
  sample2 = c(0.2, 0.2, 0.3, 0.6, 0.7, 0.7),
  sample3 = c(0.8, 0.7, 0.7, 0.2, 0.2, 0.1),
  sample4 = c(0.7, 0.8, 0.8, 0.1, 0.2, 0.2)
)
```

## Examples

Independent binary DMR analysis without covariates:

```r
library(MethylRegion)

y <- data.frame(group = c(0, 0, 1, 1))

dmr <- mr_bi(
  input_dat,
  y,
  data = "independent",
  controlist = list(mincpgs = 2, trend = 0)
)
```

Independent binary DMR analysis with covariates:

```r
covariates <- data.frame(age = c(43, 51, 39, 58))

dmr_cov <- mr_bi(
  input_dat,
  y,
  data = "independent",
  cov.mod = covariates,
  controlist = list(mincpgs = 2, trend = 0)
)
```

Paired pre/post DMR analysis:

```r
y_paired <- data.frame(
  state = c(0, 1, 0, 1),
  pair_id = c("p1", "p1", "p2", "p2")
)

dmr_prepost <- mr_bi(
  input_dat,
  y_paired,
  data = "paired",
  controlist = list(mincpgs = 2, trend = 0)
)
```

Binary paired designs with `cov.mod` and binary longitudinal designs are not
implemented yet:

```r
try(mr_bi(input_dat, y_paired, data = "paired", cov.mod = covariates))
try(mr_bi(input_dat, y, data = "longitudinal"))
```

Continuous phenotype AMR analysis:

```r
y_continuous <- data.frame(phenotype = c(1.2, 2.0, 4.2, 5.1))

amr <- mr_continuous(
  input_dat,
  y_continuous,
  data = "independent",
  controlist = list(mincpgs = 2, trend = 0)
)
```

Longitudinal continuous-phenotype AMR analysis:

```r
y_long <- data.frame(
  phenotype = c(1.2, 1.8, 3.9, 4.5),
  subject_id = c("s1", "s1", "s2", "s2"),
  visit = c("pre", "post", "pre", "post")
)

amr_long <- mr_continuous(
  input_dat,
  y_long,
  data = "longitudinal",
  controlist = list(mincpgs = 2, trend = 0)
)
```

The lower-level branch-specific functions can still be called directly:

```r
dmr_case_control(input_dat, y, controlist = list(mincpgs = 2, trend = 0))
dmr_case_control_cov(input_dat, y, cov.mod = covariates)
dmr_paired(input_dat, y_paired)
amr_continuous(input_dat, y_continuous)
amr_longitudinal(input_dat, y_long)
```

## Output notes

`dmr_case_control()` and therefore `mr_bi(data = "independent", cov.mod = NULL)`
retain the current e-value output columns. Other current branches do not compute
or report `e_value`, `e_adjust`, or `e_bh_significant`.

## Development notes

The main entry points provide a stable user-facing API. The branch-specific
functions intentionally preserve their own statistics and output columns so the
original prototype behavior remains available inside one package namespace.
