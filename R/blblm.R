#' @import purrr
#' @import stats
#' @import parallel
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#'
#' @aliases NULL
#' @details
#' Bag of Little Bootstraps for Linear Models
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c(".", "fit"))


#' @title Bag of Little Bootstraps Linear Model (single core)
#' @description Computes coefficeints for a linear model using blb (single core)
#'
#' @param formula Linear regression formula
#' @param data Chosen dataset
#' @param m Number of times to split data for bootstrapping
#' @param B Number of bootstrap samples
#' @return Model coefficients
#' @export
#'
blblm = function(formula, data, m = 10, B = 5000) {
  data_list = split_data(data, m)
  estimates = map(data_list,
                   ~ lm_each_subsample(
                     formula = formula,
                     data = .,
                     n = nrow(data),
                     B = B
                   ))
  res = list(estimates = estimates, formula = formula)
  class(res) = "blblm"
  invisible(res)
}


#' @title Bag of Little Bootstraps Linear Model (parallelization)
#' @description Computes coefficeints for a linear model using blb (multi core)
#'
#' @param formula Linear regression formula
#' @param data Chosen dataset
#' @param m Number of times to split data for bootstrapping
#' @param B umber of bootstrap samples
#' @param cl umber of cores to use (Made using `cl = makeCluster()`)
#' @return Model coefficients
#' @export
#'
par_blblm = function(formula,
                      data,
                      m = 10,
                      B = 5000,
                      cl = NULL) {
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


#' @title Bag of Little Bootstraps Linear Model (parallelization and load balancing/dynamic scheduling)
#' @description Computes coefficeints using blb and dynamic scheduling (multi core)
#'
#' @param formula Linear regression formula
#' @param data Chosen dataset
#' @param m umber of times to split data for bootstrapping
#' @param B Number of bootstrap samples
#' @param cl Number of cores to use (Made using `cl = makeCluster()`)
#' @return Model coefficients
#' @export
#'
parLB_blblm = function(formula,
                        data,
                        m = 10,
                        B = 5000,
                        cl = NULL) {
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


#' @title Split Data
#' @description Splits data into m parts of approximated equal sizes
#'
#' @param data Chosen dataset
#' @param m Number of times to split data for bootstrapping
#'
split_data = function(data, m) {
  idx = sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' @title Linear Model Subsample (C++)
#' @description Computes estimates
#'
#' @param formula Linear regression formula
#' @param data Chosen dataset
#' @param n Number of rows in the dataset
#' @param B Number of bootstrap samples
#'
lm_each_subsample = function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' @title Linear Model for Bootstrap
#' @description Computes regression estimates for a blb dataset
#'
#' @param formula Linear regression formula
#' @param data Chosen dataset
#' @param n Number of rows in the dataset
#'
lm_each_boot = function(formula, data, n) {
  freqs = rmultinom(1, n, rep(1, nrow(data)))
  lmC(formula, data, freqs)
}


#' @title Linear Model (in C++)
#' @description Estimate the regression estimates using a lm function written in C++
#'
#' @param formula Linear regression formula
#' @param data Chosen dataset
#' @param freqs Weights computed from lm_each_boot
#'
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


#' @title Sigma (C++)
#' @description Computes sigma value from fit with lmC
#'
#' @param object Fitted linear regression model
#' @return Sigma value
#'
blbsigma = function(object) {
  sqrt(sum(object$residuals ^2) / object$df.residual)
}


#' @title Coefficients (C++)
#' @description Computes coefficients from fit with lmC
#'
#' @param object Fitted linear regression model
#' @param formula Linear regression formula
#' @return Model coefficients
#' @export
#'
blbcoef = function(object, formula) {
  cf = coef(object)
  parm = attr(terms(formula), "term.labels") # Extracts labels
  names(cf) = c("intercept", parm) # Labels
  cf
}


#' @title Print blblm
#' @description Prints formula used in blblm
#'
#' @param x Fitted linear regression model
#' @param  ... Additional parameters
#' @export
#' @method print blblm
#'
print.blblm = function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @title Sigma blblm
#' @description Computes sigma for blblm
#'
#' @param object Fitted linear regression model
#' @param confidence Enables confidence interval
#' @param level Confidence level
#' @param  ... Additional parameters
#' @export
#' @method sigma blblm
#'
sigma.blblm = function(object,
                       confidence = FALSE,
                       level = 0.95, ...) {
  est = object$estimates
  sigma = mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha = 1 - 0.95
    limits = est %>%
      map_mean( ~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
    }
  else {
    return(sigma)
  }
}


#' @title Coefficients blblm
#' @description Extracts coefficients from the blblm
#'
#' @param object Fitted linear regression model
#' @param  ... Additional parameters
#' @export
#' @method coef blblm
#'
coef.blblm = function(object, ...) {
  est = object$estimates
  means = map_mean(est, ~ map_rbind(., "coef") %>%
                     colMeans())
  parm = attr(terms(object$formula), "term.labels") # Extract parameters
  labs = c("(intercept)", parm) # Add intercept to extracted parameters
  cf = (means)
  names(cf) = labs # Assign labels
  (cf)
}


#' @title Confidence Interval blblm
#' @description Computes confidence interval for a given confidence level
#'
#' @param object Fitted linear regression model
#' @param parm Model parameters
#' @param level Confidence level
#' @param  ... Additional parameters
#' @export
#' @method confint blblm
#'
confint.blblm = function(object,
                         parm = NULL,
                         level = 0.95,
                         ...) {
  if (is.null(parm)) {
    parm = attr(terms(object$formula), "term.labels")
  }

  alpha = 1 - level
  est = object$estimates
  out = map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>%
               quantile(c(alpha / 2, 1 - alpha / 2)))
  })

  if (is.vector(out)) {
    out = as.matrix(t(out))
  }
  dimnames(out)[[1]] = parm
  out
}


#' @title Predict
#' @description Computes a predicted value based on new data
#'
#' @param object Fitted linear regression model
#' @param new_data Dataset
#' @param confidence boolean: Enables confidence interval
#' @param level float: Confidence level
#' @param  ... Additional parameters
#' @return Predicted values; confidence interval (if TRUE)
#' @export
#' @method predict blblm
#'
predict.blblm = function(object,
                         new_data,
                         confidence = FALSE,
                         level = 0.95,
                         ...) {
  est = object$estimates
  X = model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)

  if (confidence) {
    map_mean(est,
             ~ map_cbind(., ~ X %*% .$coef[1, ]) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  }
  else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef[1, ]) %>%
               rowMeans())
  }
}


mean_lwr_upr = function(x, level = 0.95) {
  alpha = 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>%
      set_names(c("lwr", "upr")))
}


map_mean = function(.x, .f, ...) {
  (map(.x, .f, ...) %>%
     reduce(`+`)) / length(.x)
}


map_cbind = function(.x, .f, ...) {
  map(.x, .f, ...) %>%
    reduce(cbind)
}


map_rbind = function(.x, .f, ...) {
  map(.x, .f, ...) %>%
    reduce(rbind)
}

