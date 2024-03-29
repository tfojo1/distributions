% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/canonical_mixture_constructors.R
\name{Logitnormal.Mixture}
\alias{Logitnormal.Mixture}
\title{Create a Mixture Univariate Logit-Normal Distributions}
\usage{
Logitnormal.Mixture(
  meanslogit,
  sdslogit,
  weights = 1,
  lower = 0,
  upper = 1,
  var.name = NULL
)
}
\arguments{
\item{meanslogit, sdslogit}{The means and standard deviations (on the logit scale) of the components of the distribution}

\item{weights}{The weights applied to each component of the mixture. Weights need not be normalized to sum to 1. If a scalar value is given, all components of the mixture will have the same weight}

\item{lower, upper}{The lower and the upper bounds (on the original scale, NOT the logit scale), if this is a bounded distribution (the same bounds apply for all components of the mixture)}
}
\description{
Create a Mixture Univariate Logit-Normal Distributions
}
\seealso{
Other Canonical Mixture Distribution Constructors: 
\code{\link{Lognormal.Mixture}()},
\code{\link{Normal.Mixture}()},
\code{\link{Transformed.Normal.Mixture}()}

Other Distribution Constructors: 
\code{\link{Autoregressive.Multivariate.Normal.Distribution}()},
\code{\link{Bernoulli.Distribution}()},
\code{\link{Beta.Distribution}()},
\code{\link{Binomial.Distribution}()},
\code{\link{Canonical.Mixture.Distribution}()},
\code{\link{Compound.Symmetry.Multivariate.Normal.Distribution}()},
\code{\link{Constant.Distribution}()},
\code{\link{Discrete.Set.Distribution}()},
\code{\link{Empiric.Distribution}()},
\code{\link{Logitnormal.Distribution}()},
\code{\link{Logituniform.Distribution}()},
\code{\link{Lognormal.Distribution}()},
\code{\link{Lognormal.Mixture}()},
\code{\link{Loguniform.Distribution}()},
\code{\link{Multivariate.Correlated.Uniform.Distribution}()},
\code{\link{Multivariate.Logitnormal.Distribution}()},
\code{\link{Multivariate.Lognormal.Distribution}()},
\code{\link{Multivariate.Normal.Distribution}()},
\code{\link{Normal.Distribution}()},
\code{\link{Normal.Mixture}()},
\code{\link{Smoothed.Empiric.Distribution}()},
\code{\link{Transformed.Multivariate.Normal.Distribution}()},
\code{\link{Transformed.Normal.Distribution}()},
\code{\link{Transformed.Normal.Mixture}()},
\code{\link{Uniform.Distribution}()},
\code{\link{Univariate.Canonical.Distribution}()},
\code{\link{join.distributions}()}
}
\concept{Canonical Mixture Distribution Constructors}
\concept{Distribution Constructors}
