# neutralitytestr

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
The package comes with some preloaded test data, generated from a simulation of tumour growth. These test data sets are called ```VAFselection``` and ```VAFneutral```. We'll demonstrate the functionality of the package with the ``VAFselection`` data.

The basic functionality of the ```neutralitytestr``` package is achieved by creating a ```neutralitytest``` object. The ```neutralitytest``` object contains a range of metrics to test for neutrality, and makes plotting histograms and cumulative distributions to visualize the output easy. The ```neutralitytest``` function takes a vector of VAFs and an upper and lower limit for the frequency range over which we wish to test whether the data is consistent with a neutral model, and then calculates all 4 metrics.
```R
out <- neutralitytest(VAFneutral, fmin = 0.05, fmax = 0.4)
```

The neutralitytest object can be summarised using the ```summary(out)``` command.
```
Summary of neutrality metrics:

Area:
  value =  0.001196022 , p-value =  0.988
Kolmogorov Distance:
  value =  0.04412539 , p-value =  0.953
Mean distance:
  value =  0.01950672 , p-value =  0.906
R^2:
  value =  0.9941321 , p-value =  0.585

Effective mutation rate =  33.84608
```

The following commands will plot a VAF histogram, a least squares model fit and the normalized distributions. For more information see the vignette in ```vignettes/```.
```R
vaf_histogram(out)
```
![plot](/figure/unnamed-chunk-7-1.png)
```R
lsq_plot(out)
```
![plot](/figure/unnamed-chunk-8-1.png)
```R
normalized_plot(out)
```
![plot](/figure/unnamed-chunk-9-1.png)
