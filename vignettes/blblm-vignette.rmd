---
title: "Bag of Little Bootstraps for Linear Models"
author: "Patrick Soong"
date: "6/9/2020"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Bag of Little Bootstraps}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(blblm)
```

## blblm package

### Topics
 * Overview
 * "Bag of Little Bootstraps" Algorithm
 * Main Functions
 * Examples
 * References
 * Changelog
 
 
### Overview

The `blblm` package is used for statistical analysis of linear models. The functions in the package return the coefficients, sigma values, predicted values, and the confidence interval. The user can also enable parallelization and dynamic scheduling to increase performance on larger datasets. The functions included in the package are based off linear regression and least squares techniques.

Please view changelog for full list of changes from the original package.


### The "Bag of Little Bootstraps Algorithm

The "Bag of Little Bootstraps" algorithm first splits the original data into subsamples. Then the algorithm takes bootstrap samples for each subsample. After bootstrapping, the algorithm computes the statistics for each subsample and uses reduce to return an average for each statistic.


### Functions

#### 1. par_blblm

```r
par_blblm = function(formula,
                      data,
                      m = 10,
                      B = 5000,
                      cl) {
  data_list = split_data(data, m)
  estimates =
    parLapply(cl, data_list, function(formula, data, n, B) {
      lm_each_subsample(
        formula = formula,
        data = data,
        n = nrow(data),
        B = B
      )
    },
    formula = formula, n = nrow(data), B = B)
  results = list(estimates = estimates, formula = formula)
  class(results) = "blblm"
  invisible(results)
}

```
**Input**

  * `formula`: A linear regression formula
  * `data`: A chosen dataframe
  * `m`: Number of subsamples
  * `B`: Number of bootstrap samples
  * `cl`: Number of cores for use in parallelization (Made using `cl = makeCluster()`)

**Returns**

Coefficients from a fitted linear model
  
**Description**

This function provides the coefficients for a fitted linear model. This function performs the same function as `blblm()`, but uses `parLapply()` to allow the user to utilize parallelization. One must first create a cluster use `makeCluster()` from R's `parallel` package.

#### 2. parLB_blblm

```r
parLB_blblm = function(formula,
                        data,
                        m = 10,
                        B = 5000,
                        cl) {
  data_list = split_data(data, m)
  estimates =
    parLapplyLB(cl, data_list, function(formula, data, n, B) {
      lm_each_subsample(
        formula = formula,
        data = data,
        n = nrow(data),
        B = B
      )
    },
    formula = formula, n = nrow(data), B = B)
  results = list(estimates = estimates, formula = formula)
  class(results) = "blblm"
  invisible(results)
}

```
**Input**

  * `formula`: A linear regression formula
  * `data`: A chosen dataframe
  * `m`: Number of subsamples
  * `B`: Number of bootstrap samples
  * `cl`: Number of cores for use in parallelization (Made using `cl = makeCluster()`)

**Returns**

Coefficients from a fitted linear model
  
**Description**

This function provides the coefficients for a fitted linear model. This function performs the same function as `par_blblm()`, but uses `parLapplyLB()` to allow the user to utilize parallelization and load balancing/dynamic scheduling. One must first create a cluster use `makeCluster()` from R's `parallel` package.

#### 3. LmC

``` r
lmC = function(formula, data, freqs) {
  environment(formula) = environment()
  X = model.matrix(reformulate(attr(terms(formula), "term.labels")), data)
  y = as.matrix(data[, all.vars(formula)[1]])
  w = freqs
  rw = sqrt(w)
  rw = as.vector(w)

  X_t = rw * X
  y_t = rw * y

  fit = fast_lm(X_t, y_t) # Calls C++ version

  list(
    formula = formula,
    coef = blbcoef(fit, formula),
    sigma = blbsigma(fit),
    stderr = fit$stderr
  )
}

```
**Input**

  * `formula`: A linear regression formula
  * `data`: A chosen dataframe
  * `freq`: Weights computed from `lm_each_boot()`

**Returns**

A list containing the following items:
  
  * `formula`: A linear regression formula
  * `coef`: The linear model's coefficients
  * `sigma`: The linear model's sigma value
  * `stderr`: The linear model's standard error

**Description**

A modified version of `lm1()` that calls a `fast_lm()` function written in C++. This function preprocesses the data for use in `fast_lm()` and applies the weights calculated from `lm_each_boot()` to both model matrix (X) and response vector (y) before passing them to `fast_lm()`. `lmC()` returns the computed statistics for the given data after calling `fast_lm()`.
  
#### 4. fast_lm (C++)

``` cpp
List fast_lm(const arma::mat& X, const arma::colvec& y) {
  int n = X.n_rows, k = X.n_cols;

  arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
  arma::colvec res  = y - X*coef;           // residuals

  // std.errors of coefficients
  double s2 =
    std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);

  arma::colvec std_err =
    arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));

  return List::create(Named("coefficients") = coef,
                      Named("stderr")       = std_err,
                      Named("residuals")  = res,
                      Named("df.residual")  = n - k);
}

```
**Input**

  * `const arma::vec &y`: Response vector
  * `const arma::mat &X`: Model matrix

**Returns**

A list containing the following items:
  
  * `coefficeints`: Linear model's coefficients
  * `stderr`: Linear model's standard error estimates
  * `residuals`: Linear model's residual values
  * `df.residuals`: Degrees of freedom of the residuals

**Description**

A linear model function implemented in C++ for improved performance. This function returns the statistics used in all other functions. Benchmarks are included in the relevant section. This function is modeled off of `fastLM()` from the "Fast Linear Models with Armaillo" page from the _RCPP Gallery_.


### Examples

``` r
library(blblm)
library(parallel)

# Fitting a model
fit = blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)

# Fitting a model with parallelization
cl = makeCluster(2)
fit = blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, cl)
stopCluster(cl)

# Fitting a model with parallelization and load balancing
cl = makeCluster(2)
fit = blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, cl)
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


### References

 1. *original `blblm` package* Randy Lai (STA 141C Spring 2020).
    * <https://github.com/ucdavis-sta141c-sq-2020/blblm>
 2. *Fast Linear Models with Armadillo.* Dirk Eddelbuettel, Dec 19, 2012.
    * <https://gallery.rcpp.org/articles/fast-linear-model-with-armadillo/>
 3. *FastLM.R* Dirk Eddelbuettel, Aug 19, 2017.
    * <https://github.com/RcppCore/RcppArmadillo/blob/master/R/fastLm.R>
 4. *FastLM.cpp* Dirk Eddelbuettel, Aug 19, 2017
    * <https://github.com/RcppCore/RcppArmadillo/blob/master/src/fastLm.cpp>
 5. *Linear Algebra with RcppArmadillo* Jonathan Olmsted (Q-APS), May 30, 2014.
    * <https://scholar.princeton.edu/sites/default/files/q-aps/files/slides_day4_am.pdf>
 6. *R: lm() result differs when using `weights` argument and when using manually reweighted data* Magean, Sep 4, 2016
    * <https://stackoverflow.com/questions/39315867/r-lm-result-differs-when-using-weights-argument-and-when-using-manually-rew>
 
  
### Changelog

#### Release 1.0.0

  * Added `par_blblm()` function
    * Applied parallelization to base `blblm()` function
    
  * Added `parLB_blblm()` function
    * Applied parallelization and load balancing to base `blblm()` function
    
  * Changed `lm1()` function to `lmC()` function
    * Modified `lm1()` function to work with `fast_lm()` C++ function
    * Function preproccesses data for use in `fast_lm()` function
    
  * Added `fasm_lm()` function
    * `fast_lm()` is a C++ implementation of the `lm()` function optimized for performance
    
  * Modified `blbsigma()`
    * Modified `blbsigma()` to compute sigma value from `fast_lm()` results
    
  * Modified `blbcoef()`
    * Modified `blbcoef()` to extract coefficients from `fast_lm()` results
    
  * Added additional tests
    * Added the following tests: `test-blblm.R`, `test-confint.R`, `test-par_blblm.R`, `test-sigma.R`
    
  * Added vignette for documentation purposes
  
  * Updated readme to include basic details about the package
    
  
  
