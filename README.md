# neutralitytestr

[![Travis-CI Build Status](https://travis-ci.org/marcjwilliams1/neutralitytestr.svg?branch=master)](https://travis-ci.org/marcjwilliams1/neutralitytestr)

This is an R package to analyse Variant Allele Frequency (VAF) distributions as reported from high throughput cancer sequencing. It reports 4 summary statistics and associated p-values based on a neutral model of tumour evolution, as well as functions to plot the VAF histogram and model fits.

## Getting Started
To download the package you will need the ```devtools``` package. Once this is installed you can then download the ```neutralitytestr``` package with the command:
```R
devtools::install_github("marcjwilliams1/neutralitytestr")
```
In your R session you can then start using the package with the usual command
```R
library(neutralitytestr)
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

We can do the same with the VAFselection data:

```R
out <- neutralitytest(VAFselection, fmin = 0.1, fmax = 0.25)
plot_all(out) #this will plot all 3 of the above plots and combine into 1 figure.
```


### Notes
Note that the p-values should be interpreted with care and are meant to serve as an approximation to guide the interpretation of the test statistics. These p-values were generated empirically from a simulated cohort of cancers with known ground truth and are derived from the same data that generated the ROC curves in supplementary figure 3 from the paper. This cohort of simulated tumours were "sequenced" to 100X and thus if a sample you are analysing is sequenced to much higher or lower depth the p-values may no longer be valid. We have also developed a Bayesian alternative to identifying neutral and non-neutral tumours, this is available [here](https://github.com/marcjwilliams1/SubClonalSelection.jl). Note that this Bayesian method is much more computationally expensive and can take upwards of 10 hours per sample.
