.methylregion_default_control <- function() {
  list(
    maxdist = 300,
    method = "pearson",
    maxseg = -1,
    mincpgs = 5,
    threads = 1,
    mode = 1,
    mtc = 1,
    name = "sample",
    trend = 0.6,
    minNo = -1,
    minFactor = 0.8,
    valley = 0.7,
    minMethDist = 0.1,
    randomseed = 26061981
  )
}

.methylregion_control <- function(controlist) {
  if (is.null(controlist)) {
    controlist <- list()
  }
  if (!is.list(controlist)) {
    stop("controlist must be a list.", call. = FALSE)
  }
  defaults <- .methylregion_default_control()
  for (name in names(controlist)) {
    defaults[[name]] <- controlist[[name]]
  }
  defaults
}

.methylregion_input <- function(input_dat, intput_dat = NULL) {
  if (missing(input_dat)) {
    if (is.null(intput_dat)) {
      stop("input_dat is required.", call. = FALSE)
    }
    input_dat <- intput_dat
  }
  input_dat <- as.data.frame(input_dat)
  required <- c("chr", "pos")
  missing_cols <- setdiff(required, names(input_dat))
  if (length(missing_cols) > 0) {
    stop(
      "input_dat must contain columns: ",
      paste(required, collapse = ", "),
      call. = FALSE
    )
  }
  input_dat
}

.methylregion_result <- function(result) {
  if (is.null(result) || (!is.data.frame(result) && length(result) == 0)) {
    return(data.frame())
  }
  result
}

.methylregion_y_frame <- function(y) {
  y_df <- as.data.frame(y)
  if (ncol(y_df) < 1) {
    stop("y must contain at least one column.", call. = FALSE)
  }
  y_df
}

.methylregion_group_levels <- function(y) {
  group <- as.character(y[[1]])
  unique(group[!is.na(group)])
}

.methylregion_bind_rows_fill <- function(rows) {
  rows <- rows[vapply(rows, is.data.frame, logical(1))]
  if (length(rows) == 0) {
    return(data.frame())
  }
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  if (length(cols) == 0) {
    return(data.frame())
  }
  rows <- lapply(rows, function(x) {
    missing_cols <- setdiff(cols, names(x))
    for (col in missing_cols) {
      x[[col]] <- rep(NA, nrow(x))
    }
    x[, cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

.methylregion_annotate_pair_result <- function(result, group_a, group_b) {
  result <- .methylregion_result(result)
  comparison <- paste(group_a, group_b, sep = "_vs_")
  result$comparison <- rep(comparison, nrow(result))
  result$groupA <- rep(group_a, nrow(result))
  result$groupB <- rep(group_b, nrow(result))
  leading <- c("comparison", "groupA", "groupB")
  result[, c(leading, setdiff(names(result), leading)), drop = FALSE]
}

.methylregion_mr_bi_binary <- function(input_dat, y, data.type, cov.mod, controlist) {
  if (data.type == "independent") {
    if (is.null(cov.mod)) {
      return(dmr_case_control(input_dat, y, controlist = controlist))
    }
    return(dmr_case_control_cov(input_dat, y, cov.mod = cov.mod, controlist = controlist))
  }

  if (data.type == "paired") {
    if (!is.null(cov.mod)) {
      return(dmr_longitudinal(input_dat, y, cov.mod = cov.mod, controlist = controlist))
    }
    return(dmr_paired(input_dat, y, controlist = controlist))
  }

  dmr_longitudinal(input_dat, y, cov.mod = cov.mod, controlist = controlist)
}

.methylregion_mr_bi_pairwise <- function(input_dat, y, data.type, cov.mod, controlist) {
  y_df <- .methylregion_y_frame(y)
  group <- as.character(y_df[[1]])
  levels <- .methylregion_group_levels(y_df)
  if (length(levels) < 2) {
    stop("At least two non-missing groups are required in y.", call. = FALSE)
  }
  if (ncol(input_dat) - 2 != nrow(y_df)) {
    stop("The number of methylation sample columns must match the number of rows in y.", call. = FALSE)
  }
  if (!is.null(cov.mod) && nrow(as.data.frame(cov.mod)) != nrow(y_df)) {
    stop("cov.mod must have one row per sample.", call. = FALSE)
  }

  pairs <- combn(levels, 2, simplify = FALSE)
  comparison_labels <- vapply(pairs, paste, character(1), collapse = "_vs_")
  results <- lapply(pairs, function(pair) {
    keep <- which(group %in% pair)
    pair_input <- input_dat[, c(1, 2, keep + 2), drop = FALSE]
    pair_y <- y_df[keep, , drop = FALSE]
    pair_y[[1]] <- ifelse(as.character(pair_y[[1]]) == pair[1], 0, 1)
    pair_cov <- NULL
    if (!is.null(cov.mod)) {
      pair_cov <- as.data.frame(cov.mod)[keep, , drop = FALSE]
    }
    result <- .methylregion_mr_bi_binary(pair_input, pair_y, data.type, pair_cov, controlist)
    .methylregion_annotate_pair_result(result, pair[1], pair[2])
  })
  result <- .methylregion_bind_rows_fill(results)
  attr(result, "comparisons") <- comparison_labels
  result
}

#' Detect methylation regions for a binary phenotype
#'
#' `mr_bi()` is the main entry point for binary phenotype designs. If the first
#' column of `y` contains more than two groups, `mr_bi()` runs all pairwise
#' group comparisons and returns one combined result with `comparison`,
#' `groupA`, and `groupB` columns. The `data.type` argument selects the study
#' structure, while `cov.mod = NULL` versus a non-`NULL` covariate data frame
#' determines whether a covariate-adjusted branch is used when that branch
#' exists.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Data frame, matrix, or vector whose first column is the binary group
#'   or state indicator coded as 0 and 1. If the first column contains more than
#'   two distinct non-missing values, all pairwise comparisons are run after
#'   recoding each pair to 0 and 1.
#' @param data.type Study structure. One of `"independent"`, `"paired"`, or
#'   `"longitudinal"`.
#' @param cov.mod Optional data frame of sample-level covariates. For
#'   `data.type = "independent"`, `NULL` dispatches to `dmr_case_control()` and
#'   non-`NULL` dispatches to `dmr_case_control_cov()`. For binary paired
#'   designs with covariates and binary longitudinal designs, `mr_bi()`
#'   dispatches to `dmr_longitudinal()`.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions. Multi-group
#'   `mr_bi()` calls include `comparison`, `groupA`, and `groupB` columns.
mr_bi <- function(input_dat, y, data.type = c("independent", "paired", "longitudinal"),
                  cov.mod = NULL, controlist = list(), intput_dat = NULL) {
  data.type <- match.arg(data.type)
  y_df <- .methylregion_y_frame(y)
  if (length(.methylregion_group_levels(y_df)) > 2) {
    input_dat <- .methylregion_input(input_dat, intput_dat)
    return(.methylregion_mr_bi_pairwise(input_dat, y_df, data.type, cov.mod, controlist))
  }

  if (data.type == "independent") {
    if (is.null(cov.mod)) {
      return(dmr_case_control(input_dat, y, controlist = controlist, intput_dat = intput_dat))
    }
    return(dmr_case_control_cov(input_dat, y, cov.mod = cov.mod, controlist = controlist, intput_dat = intput_dat))
  }

  if (data.type == "paired") {
    if (!is.null(cov.mod)) {
      return(dmr_longitudinal(input_dat, y, cov.mod = cov.mod, controlist = controlist, intput_dat = intput_dat))
    }
    return(dmr_paired(input_dat, y, controlist = controlist, intput_dat = intput_dat))
  }

  dmr_longitudinal(input_dat, y, cov.mod = cov.mod, controlist = controlist, intput_dat = intput_dat)
}

#' Detect methylation regions for a continuous phenotype
#'
#' `mr_continuous()` is the main entry point for continuous phenotype designs.
#' The `data.type` argument selects independent versus longitudinal/repeated-measure
#' analyses. `cov.mod` is passed through and may be `NULL`.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Data frame, matrix, or vector whose first column is a numeric
#'   continuous phenotype.
#' @param data.type Study structure. One of `"independent"` or `"longitudinal"`.
#' @param cov.mod Optional data frame of sample-level covariates.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
mr_continuous <- function(input_dat, y, data.type = c("independent", "longitudinal"),
                          cov.mod = NULL, controlist = list(), intput_dat = NULL) {
  data.type <- match.arg(data.type)

  if (data.type == "independent") {
    return(amr_continuous(input_dat, y, cov.mod = cov.mod, controlist = controlist, intput_dat = intput_dat))
  }

  amr_longitudinal(input_dat, y, cov.mod = cov.mod, controlist = controlist, intput_dat = intput_dat)
}

#' Detect DMRs for a binary case-control design without covariates
#'
#' `dmr_case_control()` compares two independent sample groups coded as 0 and 1.
#' It uses the 2D-KS/Mann-Whitney branch implementation.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Data frame, matrix, or vector whose first column is the binary group
#'   indicator. Controls must be coded 0 and cases/tests must be coded 1.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
dmr_case_control <- function(input_dat, y, controlist = list(), intput_dat = NULL) {
  input_dat <- .methylregion_input(input_dat, intput_dat)
  .methylregion_result(
    .dmr_case_control_impl(input_dat, y, cov.mod = NULL, controlist = .methylregion_control(controlist))
  )
}

#' Detect DMRs for a binary case-control design with covariates
#'
#' `dmr_case_control_cov()` tests the group effect while adjusting for
#' sample-level covariates using the covariate linear-model branch.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Binary group variable; controls are coded 0 and cases/tests are coded 1.
#' @param cov.mod Optional data frame of sample-level covariates.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
dmr_case_control_cov <- function(input_dat, y, cov.mod = NULL, controlist = list(), intput_dat = NULL) {
  input_dat <- .methylregion_input(input_dat, intput_dat)
  .methylregion_result(
    .dmr_case_control_cov_impl(input_dat, y, cov.mod = cov.mod, controlist = .methylregion_control(controlist))
  )
}

#' Detect paired pre/post treatment DMRs
#'
#' `dmr_paired()` compares paired binary states coded as 0 and 1. If `y`
#' contains a `pair_id` column, the region p value is calculated with a paired
#' Wilcoxon test; otherwise it falls back to the unpaired rank-sum test.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Data frame whose first column is the binary pre/post state. Include
#'   `pair_id` for paired testing.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
dmr_paired <- function(input_dat, y, controlist = list(), intput_dat = NULL) {
  input_dat <- .methylregion_input(input_dat, intput_dat)
  .methylregion_result(
    .dmr_paired_impl(input_dat, y, cov.mod = NULL, controlist = .methylregion_control(controlist))
  )
}

#' Detect longitudinal DMRs for a binary phenotype
#'
#' `dmr_longitudinal()` tests a binary phenotype effect with a linear mixed
#' model of the form `meth ~ phenotype + (1 | family)` when a family-like random
#' effect column is available. Additional non-random columns in `y` or `cov.mod`
#' are included as fixed effects, matching the model-building behavior of
#' `amr_longitudinal()`. If mixed-model fitting is not available, the function
#' falls back to a fixed-effect linear model.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Data frame whose first column is the binary phenotype coded 0 and 1.
#'   Include a family-like column such as `family` or `family_id` for the random
#'   intercept.
#' @param cov.mod Optional data frame of additional sample-level covariates.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
dmr_longitudinal <- function(input_dat, y, cov.mod = NULL, controlist = list(), intput_dat = NULL) {
  input_dat <- .methylregion_input(input_dat, intput_dat)
  .methylregion_result(
    .dmr_longitudinal_impl(input_dat, y, cov.mod = cov.mod, controlist = .methylregion_control(controlist))
  )
}

#' Detect AMRs associated with a continuous phenotype
#'
#' `amr_continuous()` identifies methylation regions associated with a numeric
#' phenotype, optionally adjusting for sample-level covariates.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Numeric continuous phenotype in the first column.
#' @param cov.mod Optional data frame of sample-level covariates.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
amr_continuous <- function(input_dat, y, cov.mod = NULL, controlist = list(), intput_dat = NULL) {
  input_dat <- .methylregion_input(input_dat, intput_dat)
  .methylregion_result(
    .amr_continuous_impl(input_dat, y, cov.mod = cov.mod, controlist = .methylregion_control(controlist))
  )
}

#' Detect longitudinal AMRs associated with a continuous phenotype
#'
#' `amr_longitudinal()` identifies regions associated with a numeric continuous
#' phenotype while allowing adjustment variables and repeated-measure identifiers.
#' If a subject/id column is available and `lmerTest` is installed, the branch
#' implementation attempts a random-intercept mixed model; otherwise it falls
#' back to a fixed-effect linear model. Discrete phenotypes are not supported.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Data frame whose first column is a numeric continuous phenotype.
#'   Additional columns are adjustment variables or repeated-measure identifiers.
#' @param cov.mod Optional data frame of additional sample-level covariates.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
amr_longitudinal <- function(input_dat, y, cov.mod = NULL, controlist = list(), intput_dat = NULL) {
  input_dat <- .methylregion_input(input_dat, intput_dat)
  .methylregion_result(
    .amr_longitudinal_impl(input_dat, y, cov.mod = cov.mod, controlist = .methylregion_control(controlist))
  )
}
