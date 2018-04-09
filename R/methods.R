#' Testing for neutrality on cancer sequencing data
#'
#' \code{neutralitytest} returns a neutralitytest object which contains the result of
#' various test statistics to test for neutrality as described in Williams et al. Nature Genetics 2018.
#'
#' @param VAF Vector of variant allele frequencies (VAFs) from a deep sequencing experiment,
#' numbers should be between 0 and 1
#' @param fmin Minimum VAF of integration range, default is 0.12
#' @param fmax Maximum VAF of integration range, default is 0.24
#' @return neutralitytest object which contains test statistics which tests
#' if the sequencing data is consistent a neutral evolutionary model.
#' Test statistics are area between theoretical and empirical curves, kolmogorov distance, mean distance and R^2 statistics
#' from linear model fit. Also returns an estimate of the mutation rate per tumour tumour doubling, the raw VAFs and
#' cumulative distribution
#' @examples
#' neutralitytest(runif(100))
#' neutralitytest(VAFselection, fmin = 0.1, fmax = 0.25)
#' @export
neutralitytest <- function(VAF, fmin = 0.12, fmax = 0.24) {

  cumulativefrequency <- analyse_vaf(VAF, fmin, fmax)

  A <- areametric(cumulativefrequency, fmin, fmax)

  Dk <- kolmogorovdist(data.frame(VAF = VAF), fmin, fmax)

  meanD <- meandist(cumulativefrequency)

  fit <- lsqfit(cumulativefrequency, fmin, fmax)

  out <- list(mutation.rate = fit$mu,
              model = fit,
              rsq = data.frame(metric = fit$rsq, pval = fit$pval),
              area = A,
              Dk = Dk,
              meanDist = meanD,
              cumulativefrequency = cumulativefrequency,
              VAF = VAF)

  class(out) <- "neutralitytest"

  invisible(out)

}

#' @export
print.neutralitytest <- function(x, ...){
  cat("Summary of neutrality metrics:","\n\n")

  cat("Area:\n ","value = ", x$area$metric, ", p-value = ", x$area$pval,"\n")
  cat("Kolmogorov Distance:\n ","value = ", x$Dk$metric, ", p-value = ", x$Dk$pval,"\n")
  cat("Mean distance:\n ","value = ", x$meanD$metric, ", p-value = ", x$meanD$pval,"\n")
  cat("R^2:\n ","value = ", x$rsq$metric, ", p-value = ", x$rsq$pval,"\n\n")

  cat("Effective mutation rate = ",x$mutation.rate)
}

#' @export
summary.neutralitytest <- function(object, ...){
  cat("Summary of neutrality metrics:","\n\n")

  cat("Area:\n ","value = ", object$area$metric, ", p-value = ", object$area$pval,"\n")
  cat("Kolmogorov Distance:\n ","value = ", object$Dk$metric, ", p-value = ", object$Dk$pval,"\n")
  cat("Mean distance:\n ","value = ", object$meanD$metric, ", p-value = ", object$meanD$pval,"\n")
  cat("R^2:\n ","value = ", object$rsq$metric, ", p-value = ", object$rsq$pval,"\n\n")

  cat("Effective mutation rate = ", object$mutation.rate)
}
