% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{predict.blblm}
\alias{predict.blblm}
\title{Predict}
\usage{
\method{predict}{blblm}(object, new_data, confidence = FALSE, level = 0.95, ...)
}
\arguments{
\item{object}{Fitted linear regression model}

\item{new_data}{Dataset}

\item{confidence}{boolean: Enables confidence interval}

\item{level}{float: Confidence level}

\item{...}{Additional parameters}
}
\value{
Predicted values; confidence interval (if TRUE)
}
\description{
Computes a predicted value based on new data
}
