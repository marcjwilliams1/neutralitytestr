# neutralitytestr

### WARNING: It is recommended that you use [*mobster*](https://github.com/caravagnalab/mobster) for this type of analysis moving forward

[![Travis-CI Build Status](https://travis-ci.org/marcjwilliams1/neutralitytestr.svg?branch=master)](https://travis-ci.org/marcjwilliams1/neutralitytestr)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/neutralitytestr)](https://cran.r-project.org/package=neutralitytestr)
![](http://cranlogs.r-pkg.org/badges/grand-total/neutralitytestr)

This is an R package to analyse Variant Allele Frequency (VAF) distributions as reported from high throughput cancer sequencing. It reports 4 summary statistics and associated p-values based on a neutral model of tumour evolution, as well as functions to plot the VAF histogram and model fits.

## Getting Started
You can download the package from [CRAN](https://cran.r-project.org/package=neutralitytestr) in the usual way.
``` r
install.packages(neutralitytestr)
library(neutralitytestr)
```


To download the latest development version you'll need to use the ```devtools``` package:
```R
devtools::install_github("marcjwilliams1/neutralitytestr")
```

## Analysis
The package comes with some preloaded test data, generated from a simulation of tumour growth. These test data sets are called ```VAFselection``` and ```VAFneutral```.

The basic functionality of the ```neutralitytestr``` package is achieved by creating a ```neutralitytest``` object. The ```neutralitytest``` object contains a range of metrics to test for neutrality, and makes plotting histograms and cumulative distributions to visualize the output easy. The ```neutralitytest``` function takes a vector of VAFs and an upper and lower limit for the frequency range over which we wish to test whether the data is consistent with a neutral model, and then calculates all 4 metrics.
```R
out <- neutralitytest(VAFneutral, fmin = 0.1, fmax = 0.25)
```

The neutralitytest object can be summarised using the ```summary(out)``` command.
```
Summary of neutrality metrics:

Area:
  value =  0.03067203 , p-value =  0.679
Kolmogorov Distance:
  value =  0.09036603 , p-value =  0.641
Mean distance:
  value =  0.04131414 , p-value =  0.595
R^2:
  value =  0.98993 , p-value =  0.404

Effective mutation rate =  216.7985
```


The following commands will plot a VAF histogram, a least squares model fit and the normalized distributions. For more information see the vignette.

```R
vaf_histogram(out)
lsq_plot(out)
normalized_plot(out)
```

We can also input the read depth, cellularity, overdispersion rho and ploidy and let the package calculate an appropriate upper integration limit by considering the expected standard deviation of the clonal peak. Using this on the VAFselection data we would do the following:

```R
out <- neutralitytest(VAFselection, read_depth = 100.0, cellularity = 0.8, rho = 0.0, ploidy = 2)
plot_all(out) #this will plot all 3 of the above plots and combine into 1 figure.
```

For a more detailed introduction to the package see this [vignette]( https://CRAN.R-project.org/package=neutralitytestr/vignettes/neutraltytestr.html).

### Notes
Note that the p-values should be interpreted with care and are meant to serve as an approximation to guide the interpretation of the test statistics. These p-values were generated empirically from a simulated cohort of cancers with known ground truth and are derived from the same data that generated the ROC curves in supplementary figure 3 from the paper. This cohort of simulated tumours were "sequenced" to 100X and thus if a sample you are analysing is sequenced to much higher or lower depth the p-values may no longer be valid. We have also developed a Bayesian alternative to identifying neutral and non-neutral tumours, this is available [here](https://github.com/marcjwilliams1/SubClonalSelection.jl). Note that this Bayesian method is much more computationally expensive and can take upwards of 10 hours per sample.
