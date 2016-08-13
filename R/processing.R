theoreticalcurve <- function(f, fmin = 0.12, fmax = 0.24){

  (1 / f - 1 / fmax) / (1 / fmin - 1 / fmax)

}

analyse_vaf <- function(VAF, fmin = 0.12, fmax = 0.24) {

  # Make vector of steps between fmin and fmax
  steps <- seq(fmax,fmin,-0.001)

  # Init cumulative frequency data.frame:
  cumulativefrequency <- data.frame( M_f=sapply(steps,FUN = function(x) sum(VAF>x)), f=steps )

  # Scale data and calculate the inverse:
  cumulativefrequency$M_f   <- cumulativefrequency$M_f - cumulativefrequency$M_f[1]
  cumulativefrequency$inv_f <- ( 1 / cumulativefrequency$f - 1 / fmax )

  # Add normalized M(f)
  cumulativefrequency$nM_f <- cumulativefrequency$M_f / max(cumulativefrequency$M_f)

  # Add theoretical M(f) prediction
  cumulativefrequency$tM_f <- theoreticalcurve(cumulativefrequency$f, fmin, fmax)

  return( cumulativefrequency )
}
