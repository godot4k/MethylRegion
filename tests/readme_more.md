
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

| Entry point | Controls / dispatch | Columns, in return order |
| --- | --- | --- |
| `mr_bi(...)` | `data.type = "independent"`, `cov.mod = NULL` | `chr`, `start`, `end`, `N.CpGs`, `ks_stat`, `mean_diff`, `p_value`, `methX`, `methY`, `e_value`, `FDR`, `e_adjust` |
| `mr_bi(...)` | `data.type = "paired"`, `cov.mod = NULL` | `chr`, `start`, `end`, `N.CpGs`, `ks_stat`, `mean_diff`, `p_value`, `methX`, `methY`, `FDR` |
| `mr_bi(...)` | `data.type = "independent"`, `cov.mod` non-`NULL` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_lm_group`, `p_value`, `methX`, `methY`, `FDR` |
| `mr_bi(...)` | `data.type = "paired"`, `cov.mod` non-`NULL`; or `data.type = "longitudinal"` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_lmm`, `p_value`, `methX`, `methY`, `FDR` |
| `mr_continuous(...)` | `data.type = "independent"` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_lm`, `p_value`, `methX`, `methY`, `FDR` |
| `mr_continuous(...)` | `data.type = "longitudinal"` | `chr`, `start`, `end`, `N.CpGs`, `cor_est`, `coef_meth`, `p_value`, `methX`, `methY`, `FDR` |

Only `mr_bi(..., data.type = "independent", cov.mod = NULL)` returns e-value columns (`e_value`, `e_adjust`).
The other entry-point configurations do not currently return `e_value` or `e_adjust`.
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
| `e_value` | E value against null hypothesis. |
| `FDR` | False discovery rate adjusted P value using BH approach. |
| `e_adjust` | Adjusted E value using BH approach. |



## Demo

With `MethylRegion::mr_contious`, one could perform the following calculation:
### Demo with binary phenotype

```{R}
library(MethylRegion)
bulk.independent.data <- readRDS("./TestData/bulk.sub.txt.20.Rds")
y <- data.frame(y = rpois(20, 50))
cov.mod <- NULL
nfo <- mr_continuous(bulk.independent.data, y, data.type = "independent", cov.mod = cov.mod)
```

The results is as following:

| chr   | start   | end     | N.CpGs | ks_stat             | mean_diff             | p_value               | methX               | methY | e_value            | FDR                  | e_adjust           |
| ----- | ------- | ------- | ------ | ------------------- | --------------------- | --------------------- | ------------------- | ----- | ------------------ | -------------------- | ------------------ |
| chr21 | 9437431 | 9437472 | 7      | 0.16450216450216454 | -0.013689889164760882 | 0.8900256987555473    | 0.5643629904121975  | 0.55  | 4.6311027629749155 | 0.9065689273569515   | 2.320139782575192  |
| chr21 | 9647715 | 9648024 | 16     | 0.13825757575757575 | 0.03219346291614744   | 0.07073041521948592   | 0.5795481425014936  | 0.55  | 2.6608306521707155 | 0.1732588079343753   | 1.5364256511221384 |
| chr21 | 9704362 | 9704508 | 10     | 0.13686868686868686 | 0.024799507226123674  | 0.1740722678630579    | 0.6685971431430848  | 0.55  | 2.6132289157665136 | 0.2944185271264066   | 1.5355104213080608 |
| chr21 | 9825466 | 9825663 | 26     | 0.1668609168609169  | 0.035602086446362646  | 3.8590475161541835e-4 | 0.28645697337726633 | 0.55  | 157.63873603884426 | 0.011772106935516786 | 35.67007895769469  |
| chr21 | 9825669 | 9825796 | 20     | 0.13636363636363638 | 0.044960734639017896  | 0.006798019357858196  | 0.3521189578815905  | 0.55  | 37.4086016046296   | 0.05321878011580416  | 10.512636217359415 |
| chr21 | 9825827 | 9825874 | 10     | 0.12424242424242424 | -0.0127511118284514   | 0.3006328633307783    | 0.36266710358403925 | 0.55  | 2.312619148092159  | 0.4268052049359236   | 1.4016335732545946 |



### Demo with countious phenotype
```{R}
library(MethylRegion)
bulk.independent.data <- readRDS("./TestData/bulk.sub.txt.20.Rds")
y <- data.frame(y = rpois(20, 50))
cov.mod <- NULL
nfo <- mr_continuous(bulk.independent.data, y, data.type = "independent", cov.mod = cov.mod)
```

The Result is as following:

| chr   | start   | end     | N.CpGs | cor_est              | coef_lm            | p_value              | methX               | methY | FDR                 |
| ----- | ------- | ------- | ------ | -------------------- | ------------------ | -------------------- | ------------------- | ----- | ------------------- |
| chr21 | 9704320 | 9704392 | 10     | 0.36120935927417863  | 56.125576189511904 | 0.11764498673810235  | 0.6664195806128757  | 50.35 | 0.31305500886652343 |
| chr21 | 9709073 | 9709187 | 12     | 0.2608101311142022   | 50.825517836520014 | 0.2667248331457196   | 0.6393514542017502  | 50.35 | 0.42986720681737334 |
| chr21 | 9825466 | 9825716 | 35     | -0.4254471130382529  | -45.35237648892144 | 0.06146380906501976  | 0.32055122295901234 | 50.35 | 0.283943371021882   |
| chr21 | 9825717 | 9825748 | 5      | -0.21321161332966002 | -8.368715083755447 | 0.05757684713137738  | 0.28239092437323654 | 50.35 | 0.283943371021882   |
| chr21 | 9825756 | 9825793 | 5      | -0.4323876130982823  | -44.66298790221596 | 0.056907389515607094 | 0.31143659569301807 | 50.35 | 0.283943371021882   |
| chr21 | 9825795 | 9825871 | 15     | -0.3893348273742662  | -36.74142259281716 | 0.0897434344126131   | 0.32038108918764824 | 50.35 | 0.302336581767808   |
