kolmogorovdist <- function(dfVAF, fmin = 0.12, fmax = 0.24){

  names(dfVAF) <- "VAF"

  dfVAF <- dplyr::filter(dfVAF,VAF > fmin, VAF < fmax)

  n = length(dfVAF$VAF)

  cdfs <- 1 - theoreticalcurve(sort(dfVAF$VAF),fmin , fmax)

  dp <- max((1:n) / n - cdfs)
  dn <- - min((0:(n-1)) / n - cdfs)
  d <- max(c(dn, dp))

  return(data.frame(metric = d, pval = metricp$pval[findInterval(d,metricp$Dkmetric)]))
}

meandist <- function(cumulativefrequency) {

  d = mean(abs(cumulativefrequency$tM_f - cumulativefrequency$nM_f))

  return(data.frame(metric = d, pval = metricp$pval[findInterval(d, metricp$meanDmetric)]))

}

areametric <- function(cumulativefrequency, fmin = 0.12, fmax = 0.24){

  # Calculate area of theoretical curve
  theoryA <- pracma::trapz(cumulativefrequency$inv_f, cumulativefrequency$tM_f)

  # Calculate area of emprical curve
  dataA <- pracma::trapz(cumulativefrequency$inv_f, cumulativefrequency$nM_f)

  # Take absolute difference between the two
  A <- abs(theoryA - dataA)

  # Normalize so that metric is invariant to chosen limits
  A <- A / (1 / fmin - 1 / fmax)

  return(data.frame(metric = A, pval = metricp$pval[findInterval(A, metricp$areametric)]))

}

lsqfit <- function(cumulativefrequency, fmin = 0.12, fmax = 0.24){

  # Least squares fit
  mylm <- lm( cumulativefrequency$M_f ~ cumulativefrequency$inv_f + 0 )

  # Extract R^2 valyes
  rsq <- summary(mylm)$r.squared

  # Extract coefficient from fit
  mu <- coefficients(mylm)[1]

  return(data.frame(rsq = rsq, mu = mu, pval = metricp$pval[findInterval(1 - rsq, metricp$invrsqmetric)]))

}
