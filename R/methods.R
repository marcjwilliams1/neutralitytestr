neutralitytest <- function(dfVAF, fmin = 0.12, fmax = 0.24, removepeak = FALSE, lims = c(0.2,0.8), N = 1000) {
  
  if (removepeak == FALSE) {
    cumulativefrequency <- analyse_vaf(dfVAF$VAF, fmin, fmax)
    vafhist <- hist(dfVAF$VAF, breaks = seq(0,1.0,0.01), plot = FALSE)
    dfhist <- data.frame(VAF = vafhist$mids, freq = vafhist$counts)
    post <- NULL
    codaout <- NULL
  }
  else{
    df1 <- dfVAF[dfVAF$VAF > lims[1] & dfVAF$VAF < lims[2],]
    data <- df1[c("depth", "counts", "VAF")]
    plims <- matrix(nrow=2,ncol=2, data=c(lims[1],lims[2], 0.0,0.2), byrow=T)
    mcmcout <- run.mcmc(data, initp = c((lims[2] - lims[1])/2, 0.0001), sc = c(1.0,1.0), nsamps = N,
                     p.log=function(x){ LLbetabin(data, x) + log( prior(x,plims)) } )
    post <- mcmcout[[1]]
    codaout <- mcmcout[[2]]
    dfhist <- removepeak(dfVAF, lims, 100, post, median(dfVAF$depth))
    dftemp <- removepeak(dfVAF, lims, 1000, post, median(dfVAF$depth))
    cumulativefrequency <- cumulativefreqpeak(dftemp, fmin, fmax)
  }

  A <- areametric(cumulativefrequency, fmin, fmax)

  Dk <- kolmogorovdist(cumulativefrequency, fmin, fmax)

  meanD <- meandist(cumulativefrequency)

  fit <- lsqfit(cumulativefrequency, fmin, fmax)

  out <- list(mutation.rate = fit$mu,
              rsq = data.frame(metric = fit$rsq, pval = fit$pval),
              area = A,
              Dk = Dk,
              meanDist = meanD,
              cumulativefrequency = cumulativefrequency,
              histogram = dfhist,
              clonalposterior = post,
              codaout = codaout)

  class(out) <- "neutralitytest"

  invisible(out)

}

print.neutralitytest <- function(object){
  print("Neutrality test object")
}

summary.neutralitytest <- function(object){
  cat("Summary of neutrality metrics:","\n\n")

  cat("Area:\n ","value = ", object$area$metric, ", p-value = ", object$area$pval,"\n")
  cat("Kolmogorov Distance:\n ","value = ", object$Dk$metric, ", p-value = ", object$Dk$pval,"\n")
  cat("Mean distance:\n ","value = ", object$meanD$metric, ", p-value = ", object$meanD$pval,"\n")
  cat("R^2:\n ","value = ", object$rsq$metric, ", p-value = ", object$rsq$pval,"\n\n")

  cat("Effective mutation rate = ",object$mutation.rate)
}

plot.neutralitytest <- function(object){

  p1 <- vaf_histogram(object)
  p2 <- lsq_cumulative_plot(object)
  p3 <- norm_cumulative_plot(object)

  layout <- matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE)

  p <- multiplot(p1, p2, p3, layout = layout)

  return(p)

}
