
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RCAmle

<!-- badges: start -->

[![R-CMD-check](https://github.com/jrvanderdoes/RCAmlr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jrvanderdoes/RCAmlr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This is a package examining change point in RCA(1) model. This package
is related to the paper ‘The maximally selected likelihood ratio test in
random coefficient models’. It contains all the data and code (written
in R) to reproduce the simulations, tables, and figures in the paper.
The code used in the paper is directly saved in paperCode vignette.

This was produced in R 4.2.0 and all necessary external libraries are
loaded in the libraries section. Each function has comments describing
the purpose of the function, the inputs, and the outputs.

The data used in the analysis are all located in the related data
folder–UK covid data and US housing index data.

## Installation

You can install the released version of RCAmle from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jrvanderdoes/RCAmle@*release")
```
or target a specific release with:
``` r
# install.packages("devtools")
devtools::install_github("jrvanderdoes/RCAmle@v1.0.0")
```

You can also install the development version of RCAmle from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jrvanderdoes/RCAmle")
```
