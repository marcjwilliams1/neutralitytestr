# neutralitytestr

This is an R package to analyse Variant Allele Frequency (VAF) distributions as reported from high throughput cancer sequencing. It reports 4 summary statistics and associated p-values based on a neutral model of tumour evolution, as well as functions to plot the VAF histogram and model fits.

## Getting Started
To download the package you will need the ```devtools``` package. Once this is installed you can then download the ```neutralitytestr``` package with the command:
```
devtools::install_github("marcjwilliams1/neutralitytestr")
```
In your R session you can then start using the package with the usual command
```
library(neutralitytestr)
```

## Analysis
The package comes with some preloaded test data, generated from a simulation of tumour growth. These test data sets are called ```VAFselection``` and ```VAFneutral```. We'll demonstrate the functionality of the package with the ``VAFselection`` data.

The basic functionality of the ```neutralitytestr``` package is achieved by creating a ```neutralitytest``` object. The ```neutralitytest``` object contains a range of metrics to test for neutrality, and makes plotting histograms and cumulative distributions to visualize the output easy. The ```neutralitytest``` function takes a vector of VAFs and an upper and lower limit for the frequency range over which we wish to test whether the data is consistent with a neutral model, and then calculates all 4 metrics.
```
out <- neutralitytest(VAFselection, fmin = 0.05, fmax = 0.4)
```

The neutralitytest object can be summarised using the ```summary(out)``` command.

The following commands will plot a VAF histogram, a least squares model fit and the normalized distributions. For more information see the vignette in ```vignettes/```.
```
vaf_histogram(out)
lsq_plot(out)
normalized_plot(out)
```
