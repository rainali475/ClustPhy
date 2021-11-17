
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ClustPhy

<!-- badges: start -->
<!-- badges: end -->

## Description

`ClustPhy` is an R package for clustering phylogenetic trees (using PAM
or EM clustering), comparing different clusterings (using gap
statistics), and visualizing the clusters (in a phylogenetic tree or in
a 2D biplot).

## Installation

To install the latest version of the package:

``` r
require("devtools")
devtools::install_github("rainali475/ClustPhy", build_vignettes = TRUE)
library("TestingPackage")
```

To run the Shiny app:

``` r
Under construction
```

## Overview

``` r
ls("package:TestingPackage")
data(package = "TestingPackage")
```

`ClustPhy` contains 7 functions… Refer to package vignettes for more
details.

``` r
browseVignettes("ClustPhy")
```

An overview of the package is illustrated below.

## Contributions

The author of this package is Yuzi (Raina) Li. …

## References

…

## Acknowledgement

This package was developed as part of an assessment for 2021 BCB410H:
Applied Bioinformatics, University of Toronto, Toronto, CANADA.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ClustPhy)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
