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
print.neutralitytest <- function(object){
  print("Neutrality test object")
}

#' @export
summary.neutralitytest <- function(object){
  cat("Summary of neutrality metrics:","\n\n")

  cat("Area:\n ","value = ", object$area$metric, ", p-value = ", object$area$pval,"\n")
  cat("Kolmogorov Distance:\n ","value = ", object$Dk$metric, ", p-value = ", object$Dk$pval,"\n")
  cat("Mean distance:\n ","value = ", object$meanD$metric, ", p-value = ", object$meanD$pval,"\n")
  cat("R^2:\n ","value = ", object$rsq$metric, ", p-value = ", object$rsq$pval,"\n\n")

  cat("Effective mutation rate = ",object$mutation.rate)
}
