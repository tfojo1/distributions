#'@title Create a Mixture of Univariate Normal Distributions
#'
#'@inheritParams Normal.Distribtion
#'@param means,sds The means and standard deviations of each component of the mixture
#'@param weights The weights applied to each component of the mixture. Weights need not be normalized to sum to 1. If a scalar value is given, all components of the mixture will have the same weight
#'@param lower,upper The lower and the upper bounds, if this is a bounded distribution (the same bounds apply for all components of the mixture)
#'
#'@family Canonical Mixture Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Normal.Mixture <- function(means=0, sds=1, weights=1, lower=-Inf, upper=Inf, var.name=NULL)
{
    if (any(sds<0))
        stop("sds cannot be negative")

    if (lower==-Inf && upper==Inf)
    {
        mean.values = means
        variance.values = sds^2
    }
    else
    {
        alpha = (lower - means) / sds
        beta = (upper - means) / sds

        Z = pnorm(beta) - pnorm(alpha)
        mean.values = means + (dnorm(alpha) - dnorm(beta)) * sds / Z
        variance.values = sds^2 *
            ( 1 + (alpha*dnorm(alpha) - beta*dnorm(beta)) / Z -
                  ((dnorm(alpha) - dnorm(beta)) / Z)^2 )
    }

    Canonical.Mixture.Distribution(name=get.canonical.mixture.name('normal',
                                                                          n.components=max(length(means), length(sds)),
                                                                          lower=lower, upper=upper,
                                                                          min.lower=-Inf, max.upper=Inf),
                                          weights=weights,
                                          dist.name='norm',
                                          parameters=list(mean=means, sd=sds),
                                          var.name=var.name,
                                   support = Continuous.Support(lower, upper),
                                          mean.values=mean.values,
                                          variance.values=variance.values)
}

#'@title Create a Mixture of Univariate Log-Normal Distributions
#'
#'@param meanslog,sdslog The means and standard deviations (on the log scale) of the components of distribution
#'@param lower,upper The lower and the upper bounds (on the original scale, NOT the log scale), if this is a bounded distribution (the same bounds apply for all components of the mixture)
#'@inheritParams Normal.Mixture
#'
#'@family Canonical Mixture Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Lognormal.Mixture <- function(meanslog=0, sdslog=1, weights=1, lower=0, upper=Inf, var.name=NULL)
{
    if (lower < 0)
        stop("A log-normal distribution cannot have a lower bound less than zero")

    if (lower==0 && upper==Inf)
    {
        mean.values = exp(meanslog + sdslog^2/2)
        variance.values = (exp(sdslog^2)-1) * mean.values^2
    }
    else
        mean.values = variance.values = as.numeric(NA)

    if (any(sdslog<0))
        stop("sdslog cannot be negative")

    Canonical.Mixture.Distribution(name=get.canonical.mixture.name('log-normal',
                                                                             n.components=max(length(meanslog), length(sdslog)),
                                                                             lower=lower, upper=upper,
                                                                             min.lower=0, max.upper=Inf),
                                          weights=weights,
                                          dist.name='norm',
                                          parameters=list(mean=meanslog, sd=sdslog),
                                          var.name=var.name,
                                   support = Continuous.Support(lower,upper),
                                          transformation='log',
                                          mean.value=mean.values,
                                          variance.value=variance.values)
}

#'@title Create a Mixture Univariate Logit-Normal Distributions
#'
#'@param meanslogit,sdslogit The means and standard deviations (on the logit scale) of the components of the distribution
#'@param lower,upper The lower and the upper bounds (on the original scale, NOT the logit scale), if this is a bounded distribution (the same bounds apply for all components of the mixture)
#'@inheritParams Normal.Mixture
#'
#'@family Canonical Mixture Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Logitnormal.Mixture <- function(meanslogit, sdslogit, weights=1, lower=0, upper=1, var.name=NULL)
{
    if (lower < 0)
        stop("A logit-normal distribution cannot have a lower bound less than zero")
    if (upper > 1)
        stop("A logit-normal distribution cannot have an upper bound greater than one")

    if (any(sdslogit<0))
        stop("sdslogit cannot be negative")

    Canonical.Mixture.Distribution(get.canonical.mixture.name('logit-normal',
                                                                        n.components=max(length(meanslogit), length(sdslogit)),
                                                                        lower=lower, upper=upper,
                                                                        min.lower=0, max.upper=1),
                                          weights=weights,
                                          dist.name='norm',
                                          parameters=list(mean=meanslogit, sd=sdslogit),
                                          var.name=var.name,
                                   support = Continuous.Support(lower, upper),
                                          transformation='logit',
                                          mean.values=as.numeric(NA),
                                          variance.values=as.numeric(NA))
}

#'@title Create a Mixture of Transformed Normal Distributions
#'
#'@description Creates a mixture of distributions representing random variables that are normal after transformation (eg, mixture of normals, log-normals, logit-normals, etc.)
#'
#'@inheritParams Normal.Mixture
#'@param means,sds The means and standard deviations of each component of the mixture (after transformation has been applied)
#'@param lower,upper The lower and the upper bounds, if this is a bounded distribution. The bounds are on the ORIGINAL scale (before transformation). The same bounds apply for all components of the mixture
#'
#'@family Canonical Mixture Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Transformed.Normal.Mixture <- function(means=0, sds=1, weights=1, transformation=NULL,
                                       lower=-Inf, upper=Inf, var.name=NULL)
{
    if (any(sds<0))
        stop("sds cannot be negative")

    if (is.null(transformation))
        transformation.name = 'identity'
    else if (is(transformation, 'transformation'))
        transformation.name = transformation@name
    else if (is(transformation, 'character'))
        transformation.name = transformation
    else
        stop("'transformation' must be either an object of class 'transformation' or the name of a pre-defined transformation")

    if (transformation.name=='identity')
        Normal.Mixture(means=means, sds=sds, weights=weights, lower=lower, upper=upper, var.name=var.name)
    else if (transformation.name=='log')
        Lognormal.Mixture(meanslog=means, sdslog=sds, weights=weights, lower=lower, upper=upper, var.name=var.name)
    else if (transformation.name=='logit')
        Logitnormal.Mixture(meanslogit=means, sdslogit=sds, weights=weights, lower=lower, upper=upper, var.name=var.name)
    else
        Canonical.Mixture.Distribution(name=get.canonical.mixture.name(paste0('normal (after',
                                                                                    transformation.name,
                                                                                    ' transformation)'),
                                                                              n.components=max(length(means), length(sds)),
                                                                              lower=lower, upper=upper,
                                                                              min.lower=NA, max.upper=NA),
                                              weights=weights,
                                              dist.name='norm',
                                              parameters=list(mean=means, sd=sds),
                                              var.name=var.name,
                                       transformation=transformation,
                                       support = Continuous.Support(lower, upper),
                                              mean.values=as.numeric(NA),
                                              variance.values=as.numeric(NA))
}
