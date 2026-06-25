\name{dmr_case_control}
\alias{dmr_case_control}
\alias{dmr_case_control_cov}
\alias{dmr_paired}
\alias{amr_continuous}
\alias{amr_longitudinal}
\title{MethylRegion study-design functions}
\description{
MethylRegion provides separate entry points for methylation-region analyses
across binary DMR designs and continuous-phenotype AMR designs.
}
\usage{
dmr_case_control(input_dat, y, controlist = list(), intput_dat = NULL)

dmr_case_control_cov(input_dat, y, cov.mod = NULL, controlist = list(),
  intput_dat = NULL)

dmr_paired(input_dat, y, controlist = list(), intput_dat = NULL)

amr_continuous(input_dat, y, cov.mod = NULL, controlist = list(),
  intput_dat = NULL)

amr_longitudinal(input_dat, y, cov.mod = NULL, controlist = list(),
  intput_dat = NULL)
}
\arguments{
\item{input_dat}{A data frame with columns \code{chr}, \code{pos}, and one
methylation column per sample.}
\item{y}{Study-design response data. For DMR functions, the first column is a
binary group or state coded 0 and 1. For AMR functions, the first column is a
numeric continuous phenotype.}
\item{controlist}{Optional list of segmentation controls. Missing entries are
filled from MethylRegion defaults.}
\item{intput_dat}{Deprecated misspelling retained for compatibility with the
original AMRfinder prototypes.}
\item{cov.mod}{Optional data frame of sample-level covariates.}
}
\details{
\code{dmr_case_control()} handles independent binary case-control comparisons
without covariates using the 2D-KS/Mann-Whitney branch.

\code{dmr_case_control_cov()} handles binary case-control comparisons with
sample-level covariates using the covariate linear-model branch.

\code{dmr_paired()} handles pre/post or paired binary designs. If \code{y}
contains \code{pair_id}, the region p value is calculated with a paired
Wilcoxon test.

\code{amr_continuous()} handles continuous phenotype-associated methylation
regions.

\code{amr_longitudinal()} handles longitudinal or repeated-measure AMRs for
continuous phenotypes only. If an id-like column is present and \pkg{lmerTest}
is installed, the implementation attempts a random-intercept mixed model;
otherwise it uses a fixed-effect linear model.
}
\value{
A data frame of candidate methylation regions. Output columns vary by design
because each function preserves the corresponding branch-specific statistics.
}
\examples{
dat <- data.frame(
  chr = rep("chr1", 6),
  pos = seq(100, 600, by = 100),
  sample1 = c(0.1, 0.2, 0.2, 0.7, 0.8, 0.8),
  sample2 = c(0.2, 0.2, 0.3, 0.6, 0.7, 0.7),
  sample3 = c(0.8, 0.7, 0.7, 0.2, 0.2, 0.1),
  sample4 = c(0.7, 0.8, 0.8, 0.1, 0.2, 0.2)
)
y_bin <- data.frame(group = c(0, 0, 1, 1))
res <- dmr_case_control(dat, y_bin, controlist = list(mincpgs = 2, trend = 0))

y_cont <- data.frame(phenotype = c(1, 2, 4, 5))
res2 <- amr_continuous(dat, y_cont, controlist = list(mincpgs = 2, trend = 0))
}
\keyword{models}
\keyword{htest}
