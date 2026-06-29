# MethylRegion

**MethylRegion: A unified framework for differential methylation regions across designs, traits, and covariates**

The recommended entry points are `mr_bi()` for binary phenotypes and
`mr_continuous()` for continuous phenotypes.

## Installation

```r
remotes::install_github("godot4k/MethylRegion")
```

For local development:

```sh
R CMD INSTALL MethylRegion
```

## Main Entry Points

| Function | Phenotype | `data.type` values | Covariate dispatch |
| --- | --- | --- | --- |
| `mr_bi(input_dat, y, data.type, cov.mod)` | Binary | `"independent"`, `"paired"`, `"longitudinal"` | Uses `is.null(cov.mod)` to choose no-covariate versus covariate-aware branches where implemented. |
| `mr_continuous(input_dat, y, data.type, cov.mod)` | Continuous | `"independent"`, `"longitudinal"` | Passes `cov.mod` through; it may be `NULL` or a covariate data frame. |

## Input File Requirements

### Methylation Matrix: `input_dat`

| Column | Type | Description |
| --- | --- | --- |
| `chr` | Character | Chromosome label for each CpG site. |
| `pos` | Integer or numeric | Genomic coordinate for each CpG site. Rows should be ordered by chromosome and position. |
| Sample columns | Numeric | One methylation column per sample. Values should be methylation proportions or beta values, typically in `[0, 1]`. |

Sample columns must be in the same order as rows in `y` and `cov.mod`.

### Phenotype Table: `y`

| Entry point | `data.type` | Requirements |
| --- | --- | --- |
| `mr_bi()` | `"independent"` | One column: binary group. Controls should be coded `0`; cases/tests should be coded `1`. |
| `mr_bi()` | `"paired"` | Two columns: binary state and optional `pair_id`. Binary state should be coded `0` and `1`. Include `pair_id` for paired analysis. |
| `mr_bi()` | `"longitudinal"` | Two columns: binary phenotype and a subject/family column. Binary phenotype should be coded `0` and `1`. A family-like identifier such as `family`, `family_id`, or `subject_id` enables random-intercept modeling. |
| `mr_continuous()` | `"independent"` | One column: continuous phenotype. Numeric phenotype values, one row per sample. |
| `mr_continuous()` | `"longitudinal"` | Two columns: continuous phenotype and a subject column. Use numeric phenotype values. A subject-like identifier such as `subject_id` enables longitudinal modeling. |

### Covariate Table: `cov.mod`

| Argument  | Description |
| --- | --- |
| `cov.mod = NULL` | Runs the no-covariate branch where available. |
| `cov.mod` non-`NULL` |  One row per sample. Columns are sample-level covariates such as age, gender, or batch. |

## Output Columns Description

| Entry point | Controls | Columns in results |
| --- | --- | --- |
| `mr_bi(...)` | `data.type = "independent"`, `cov.mod = NULL` | `chr`, `start`, `end`, `N.CpGs`, `ks_stat`, `mean_diff`, `p_value`, `methX`, `methY`, `e_value`, `FDR`, `e_adjust` |
| `mr_bi(...)` | `data.type = "paired"`, `cov.mod = NULL` | `chr`, `start`, `end`, `N.CpGs`, `ks_stat`, `mean_diff`, `p_value`, `methX`, `methY`, `FDR` |
| `mr_bi(...)` | `data.type = "independent"`, `cov.mod` non-`NULL` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_lm_group`, `p_value`, `methX`, `methY`, `FDR` |
| `mr_bi(...)` | `data.type = "paired"`, `cov.mod` non-`NULL`; or `data.type = "longitudinal"` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_lmm`, `p_value`, `methX`, `methY`, `FDR` |
| `mr_continuous(...)` | `data.type = "independent"` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_lm`, `p_value`, `methX`, `methY`, `FDR` |
| `mr_continuous(...)` | `data.type = "longitudinal"` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_meth`, `p_value`, `methX`, `methY`, `FDR` |

> Only `mr_bi(..., data.type = "independent", cov.mod = NULL)` returns e-value
columns (`e_value`, `e_adjust`). The other entry-point configurations do not
currently return `e_value` or `e_adjust`.

### Column Meaning

| Column | Meaning |
| --- | --- |
| `chr` | Chromosome. |
| `start` | The start position of a methylation region. |
| `end` | The end position of a methylation region. |
| `N.CpGs` | The number of CpG sites. |
| `ks_stat` | The estimated 2D-KS statistic. |
| `mean_diff` | The mean difference of methylation between two groups. |
| `cor_est` | The estimated correlation coefficient. |
| `coef_lm_group` | The group coefficient derived from covariate-adjusted regression analysis. |
| `coef_lmm` | The phenotype coefficient derived from longitudinal mixed-model analysis. |
| `coef_lm` | The methylation coefficient derived from linear regression analysis. |
| `coef_meth` | The methylation coefficient derived from longitudinal continuous-phenotype analysis. |
| `p_value` | The methylation P value derived from the branch-specific statistical test. |
| `methX` | Mean DNA methylation level within each region across all samples. |
| `methY` | Mean phenotype value within each region across all samples. |
| `e_value` | E-value against the null hypothesis. |
| `FDR` | False discovery rate-adjusted P value using the BH approach. |
| `e_adjust` | Adjusted E-value using the BH approach. |

## Demo in Application Scenarios

### Binary Phenotype, Independent Samples

For differentially methylated region (DMRs) analysis without covariates in case/control data, the following code performs a calculation with a simple simulation setting.

```{R}
library(MethylRegion)
bulk.independent.data <- readRDS("./TestData/bulk.sub.txt.20.Rds")
y <- data.frame(y = rbinom(20, 1, 0.5))
cov.mod <- NULL
nfo <- mr_bi(bulk.independent.data, y, data.type = "independent", cov.mod = cov.mod)
```

The result is as follows:

| chr | start | end | N.CpGs | ks_stat | mean_diff | p_value | methX | methY | e_value | FDR | e_adjust |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr21 | 9437431 | 9437472 | 7 | 0.16450216450216454 | -0.013689889164760882 | 0.8900256987555473 | 0.5643629904121975 | 0.55 | 4.6311027629749155 | 0.9065689273569515 | 2.320139782575192 |
| chr21 | 9647715 | 9648024 | 16 | 0.13825757575757575 | 0.03219346291614744 | 0.07073041521948592 | 0.5795481425014936 | 0.55 | 2.6608306521707155 | 0.1732588079343753 | 1.5364256511221384 |
| chr21 | 9704362 | 9704508 | 10 | 0.13686868686868686 | 0.024799507226123674 | 0.1740722678630579 | 0.6685971431430848 | 0.55 | 2.6132289157665136 | 0.2944185271264066 | 1.5355104213080608 |
| chr21 | 9825466 | 9825663 | 26 | 0.1668609168609169 | 0.035602086446362646 | 3.8590475161541835e-4 | 0.28645697337726633 | 0.55 | 157.63873603884426 | 0.011772106935516786 | 35.67007895769469 |
| chr21 | 9825669 | 9825796 | 20 | 0.13636363636363638 | 0.044960734639017896 | 0.006798019357858196 | 0.3521189578815905 | 0.55 | 37.4086016046296 | 0.05321878011580416 | 10.512636217359415 |
| chr21 | 9825827 | 9825874 | 10 | 0.12424242424242424 | -0.0127511118284514 | 0.3006328633307783 | 0.36266710358403925 | 0.55 | 2.312619148092159 | 0.4268052049359236 | 1.4016335732545946 |

A detailed description of the output columns is provided in the `Column Meaning` section above.

### Continuous Phenotype, Independent Samples

For associated methylated regions (AMRs) analysis between phenotype and CpGs  without covariates in independent samples, the following code performs a calculation with a simple simulation setting.

```{R}
library(MethylRegion)
bulk.independent.data <- readRDS("./TestData/bulk.sub.txt.20.Rds")
y <- data.frame(y = rpois(20, 50))
cov.mod <- NULL
nfo <- mr_continuous(bulk.independent.data, y, data.type = "independent", cov.mod = cov.mod)
```

The result is as follows:

| chr | start | end | N.CpGs | cor_est | coef_lm | p_value | methX | methY | FDR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr21 | 9704320 | 9704392 | 10 | 0.36120935927417863 | 56.125576189511904 | 0.11764498673810235 | 0.6664195806128757 | 50.35 | 0.31305500886652343 |
| chr21 | 9709073 | 9709187 | 12 | 0.2608101311142022 | 50.825517836520014 | 0.2667248331457196 | 0.6393514542017502 | 50.35 | 0.42986720681737334 |
| chr21 | 9825466 | 9825716 | 35 | -0.4254471130382529 | -45.35237648892144 | 0.06146380906501976 | 0.32055122295901234 | 50.35 | 0.283943371021882 |
| chr21 | 9825717 | 9825748 | 5 | -0.21321161332966002 | -8.368715083755447 | 0.05757684713137738 | 0.28239092437323654 | 50.35 | 0.283943371021882 |
| chr21 | 9825756 | 9825793 | 5 | -0.4323876130982823 | -44.66298790221596 | 0.056907389515607094 | 0.31143659569301807 | 50.35 | 0.283943371021882 |
| chr21 | 9825795 | 9825871 | 15 | -0.3893348273742662 | -36.74142259281716 | 0.0897434344126131 | 0.32038108918764824 | 50.35 | 0.302336581767808 |

A detailed description of the output columns is also provided in the `Column Meaning` section above.
