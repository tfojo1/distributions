% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/master_class.R
\name{generate.random.marginal.samples}
\alias{generate.random.marginal.samples}
\title{Generate marginal random samples for each variable in a distribution}
\usage{
generate.random.marginal.samples(dist, n)
}
\arguments{
\item{dist}{An object of class Distribution or a subclass}

\item{n}{The number of random samples to be generated}
}
\value{
If n==1, a numeric vector with one value for each variable in the distribution. If there is only one variable in the distribution, a numeric vector with n values. If n>1 and there is more than one variable in the distribution, a matrix with n rows, where each row represents one random sample, and each column represents a variable in the distribution
}
\description{
Generates random samples for each variable in the distribution, but does NOT preserve the correlations across variables within samples. For some distribution implementations, this can be substantially faster than generating true, multivariate random samples, so this method is preferred if random samples are only needed for one variable at a time
}
\details{
Note to developers: this function should not be overridden. Override \code{\link{do.generate.random.marginal.samples}} instead
}
