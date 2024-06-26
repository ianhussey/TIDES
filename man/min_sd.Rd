% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/min_sd.R
\name{min_sd}
\alias{min_sd}
\title{Calculate minimum Standard Deviation for a truncated and binned variable}
\usage{
min_sd(mean, n)
}
\arguments{
\item{mean}{numeric variable representing the reported mean.}

\item{n}{numeric variable representing the reported sample size.}
}
\value{
A numeric variable representing the minimum SD.
}
\description{
Explanation to be added
}
\details{
Assumptions to be checked:
does this assume the mean represents (a) a variable where the possible scores are integers as opposed to something with smaller granularity, and relatedly (b) the mean of a mean mean score rather than the mean of a sum score or single item
}
