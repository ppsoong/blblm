# blblm R package

<!-- badges: start -->
<!-- badges: end -->

## Bag of Little Bootstraps for Linear Models

### Overview

The `blblm` package is used for statistical analysis of linear models. The functions in the package return the coefficients, sigma values, predicted values, and the confidence interval. The user can also enable parallelization and dynamic scheduling to increase performance on larger datasets. Additionally, some functions have been rewritten in C++ using Rcpp and RcppArmadillo to increase performance. The functions included in the package are based off linear regression and least squares techniques.

Please check `blblm-vignette.Rmd` for more information about the package and its functions.


### The "Bag of Little Bootstraps Algorithm

The "Bag of Little Bootstraps" algorithm first splits to reduce the original data into m number of subsamples. Then the algorithm takes B number of bootstrap samples for each subsample. After bootstrapping, the algorithm computes the statistics for each subsample and uses reduce to return an average for each statistic.

You can install the package by running the following command in your R Studio console:

```r
devtools::install_github("ppsoong/blblm")
```

### Examples

``` r
library(blblm)
library(parallel)

# Fitting a model
fit = blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)

# Fitting a model with parallelization
cl = makeCluster(2)
fit = par_blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, cl)
stopCluster(cl)

# Fitting a model with parallelization and load balancing
cl = makeCluster(2)
fit = parLB_blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, cl)
stopCluster(cl)

coef(fit)
#> (Intercept)          wt          hp       wt:hp 
#> 48.88428523 -7.88702986 -0.11576659  0.02600976

confint(fit, c("wt", "hp"))
#>           2.5%       97.5%
#> wt -10.7902240 -5.61586271
#> hp  -0.1960903 -0.07049867

sigma(fit)
#> [1] 1.838911
sigma(fit, confidence = TRUE)
#>    sigma      lwr      upr 
#> 1.838911 1.350269 2.276347

predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)))
#>        1        2 
#> 21.55538 18.80785

predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE)
#>        fit      lwr      upr
#> 1 21.55538 20.02457 22.48764
#> 2 18.80785 17.50654 19.71772

```
