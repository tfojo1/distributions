
#'@title Create a Univariate Normal Distribution Object
#'
#'@param mean,sd The mean and standard deviation of the distribution
#'@param var.name The name of the single variable in the distribution, or NULL if no name is specified
#'@param lower,upper The lower and the upper bounds, if this is a bounded distribution
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Normal.Distribution <- function(mean=0, sd=1, lower=-Inf, upper=Inf, var.name=NULL)
{
    if ((!is(sd, 'numeric') && !is(sd, 'integer')) ||
        length(sd)!=1 || is.na(sd) || sd<0)
        stop("sd must be a scalar, non-negative, non-NA numeric value")

    if (lower==-Inf && upper==Inf)
    {
        mean.value = mean
        variance.value = sd^2
    }
    else
    {
        alpha = (lower - mean) / sd
        beta = (upper - mean) / sd

        Z = pnorm(beta) - pnorm(alpha)
        mean.value = mean + (dnorm(alpha) - dnorm(beta)) * sd / Z
        variance.value = sd^2 *
            ( 1 + (alpha*dnorm(alpha) - beta*dnorm(beta)) / Z -
                  ((dnorm(alpha) - dnorm(beta)) / Z)^2 )
    }

    Univariate.Canonical.Distribution(name=get.canonical.distribution.name('normal', mean, sd, lower=lower, upper=upper),
                                             dist.name='norm',
                                             parameters=list(mean=mean, sd=sd),
                                             var.name=var.name,
                                             support=Continuous.Support(lower, upper),
                                             mean.value=mean.value,
                                             variance.value=variance.value)
}

#'@title Create a Univariate Log-Normal Distribution Object
#'
#'@param meanlog,sdlog The mean and standard deviation (on the log scale) of the distribution
#'@param lower,upper The lower and the upper bounds (on the original scale, NOT the log scale), if this is a bounded distribution
#'@inheritParams Normal.Distribution
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Lognormal.Distribution <- function(meanlog=0, sdlog=1, lower=0, upper=Inf, var.name=NULL)
{
    if (lower < 0)
        stop("A log-normal distribution cannot have a lower bound less than zero")

    if (lower==0 && upper==Inf)
    {
        mean.value = exp(meanlog + sdlog^2/2)
        variance.value = (exp(sdlog^2)-1) * mean.value^2
    }
    else
        mean.value = variance.value = as.numeric(NA)

    if ((!is(sdlog, 'numeric') && !is(sdlog, 'integer')) ||
        length(sdlog)!=1 || is.na(sdlog) || sdlog<0)
        stop("sdlog must be a scalar, non-negative, non-NA numeric value")

    Univariate.Canonical.Distribution(name=get.canonical.distribution.name('log-normal', meanlog, sdlog, lower=lower, upper=upper, min.lower=0),
                                             dist.name='norm',
                                             parameters=list(mean=meanlog, sd=sdlog),
                                             var.name=var.name,
                                             support = Continuous.Support(lower, upper),
                                             transformation='log',
                                             mean.value=mean.value,
                                             variance.value=variance.value)
}

#'@title Create a Univariate Logit-Normal Distribution Object
#'
#'@param meanlogit,sdlogit The mean and standard deviation (on the logit scale) of the distribution
#'@param lower,upper The lower and the upper bounds (on the original scale, NOT the logit scale), if this is a bounded distribution
#'@inheritParams Normal.Distribution
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Logitnormal.Distribution <- function(meanlogit, sdlogit, lower=0, upper=1, var.name=NULL)
{
    if (lower < 0)
        stop("A logit-normal distribution cannot have a lower bound less than zero")
    if (upper > 1)
        stop("A logit-normal distribution cannot have an upper bound greater than one")

    if ((!is(sdlogit, 'numeric') && !is(sdlogit, 'integer')) ||
        length(sdlogit)!=1 || is.na(sdlogit) || sdlogit<0)
        stop("sdlogit must be a scalar, non-negative, non-NA numeric value")

    Univariate.Canonical.Distribution(name=get.canonical.distribution.name('logit-normal', meanlogit, sdlogit, lower=lower, upper=upper, min.lower=0, max.upper=1),
                                             dist.name='norm',
                                             parameters=list(mean=meanlogit, sd=sdlogit),
                                             var.name=var.name,
                                             support = Continuous.Support(lower, upper),
                                             transformation='logit',
                                             mean.value=as.numeric(NA),
                                             variance.value=as.numeric(NA))
}

#'@title Create a Distribution of a Random Variable that is Normal after an arbitrary transformation
#'
#'@inheritParams Normal.Distribution
#'@param mean,sd The mean and standard deviation of the random variable (after transformation has been applied)
#'@param lower,upper The lower and the upper bounds, if this is a bounded distribution. The bounds are on the ORIGINAL scale (before transformation).
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Transformed.Normal.Distribution <- function(mean=0, sd=1, transformation=NULL,
                                       lower=-Inf, upper=Inf, var.name=NULL)
{
    if ((!is(sd, 'numeric') && !is(sd, 'integer')) ||
        length(sd)!=1 || is.na(sd) || sd<0)
        stop("sd must be a scalar, non-negative, non-NA numeric value")

    if (is.null(transformation))
        transformation.name = 'identity'
    else if (is(transformation, 'transformation'))
        transformation.name = transformation@name
    else if (is(transformation, 'character'))
        transformation.name = transformation
    else
        stop("'transformation' must be either an object of class 'transformation' or the name of a pre-defined transformation")

    if (transformation.name=='identity')
        Normal.Distribution(mean=mean, sd=sd, lower=lower, upper=upper, var.name=var.name)
    else if (transformation.name=='log')
        Lognormal.Distribution(meanlog=mean, sdlog=sd, lower=lower, upper=upper, var.name=var.name)
    else if (transformation.name=='logit')
        Logitnormal.Distribution(meanlogit=mean, sdlogit=sd, weights=weights, lower=lower, upper=upper, var.name=var.name)
    else
        Univariate.Canonical.Distribution(name=get.canonical.distribution.name(paste0('normal (after',
                                                                                      transformation.name,
                                                                                      ' transformation)'),
                                                                               mean, sd, lower=lower, upper=upper, min.lower=NA, max.upper=NA),
                                          dist.name='norm',
                                          parameters=list(mean=mean, sd=sd),
                                          var.name=var.name,
                                          transformation=transformation,
                                          support = Continuous.Support(lower, upper),
                                          mean.value=as.numeric(NA),
                                          variance.value=as.numeric(NA))
}


#'@title Create a Univariate Beta Distribution Object
#'
#'@param alpha,beta The parameters to the beta distribution (AKA shape1, shape2 respectively in dbeta)
#'@param lower,upper The lower and the upper bounds
#'@inheritParams Normal.Distribution
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Beta.Distribution <- function(alpha=1, beta=1, lower=0, upper=1, var.name=NULL)
{
    if (lower < 0)
        stop("A beta distribution cannot have a lower bound less than zero")
    if (upper > 1)
        stop("A beta distribution cannot have an upper bound greater than one")

    if (alpha<=0 || beta<=0)
        stop("alpha and beta must be greater than 0")

    if (lower==0 && upper==1)
    {
        mean.value = alpha / (alpha+beta)
        variance.value = alpha * beta / (alpha+beta)^2 / (alpha+beta+1)
    }
    else
        mean.value = variance.value = as.numeric(NA)

    Univariate.Canonical.Distribution(name=get.canonical.distribution.name('beta', alpha, beta, lower=lower, upper=upper, min.lower=0, max.upper=1),
                                             dist.name='beta',
                                             parameters=list(shape1=alpha, shape2=beta),
                                             var.name=var.name,
                                      support = Continuous.Support(lower, upper),
                                             mean.value=mean.value,
                                             variance.value=variance.value)
}

#'@title Create a Uniform Distribution Object
#'
#'@param min,max The bounds of the distribution. If either is infinite, will create an improper uniform distribution (which returns 1 for all densities, but cannot evaluate cdfs, quantiles, or generate random samples)
#'@inheritParams Normal.Distribution
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Uniform.Distribution <- function(min=0, max=1, var.name=NULL)
{
    if (is.infinite(min) || is.infinite(max))
    {
        Univariate.Canonical.Distribution(name=paste0("Improper Uniform(", min, ", ", max, ")"),
                                                 dist.name='NA',
                                                 density.function.name='improper.uniform.density',
                                                 parameters=list(),
                                                 var.name=var.name,
                                          support = Continuous.Support(min, max),
                                                 is.improper=T)
    }
    else
    {
        Univariate.Canonical.Distribution(name=get.canonical.distribution.name('uniform', min, max, lower=min, upper=max, min.lower=NA, max.lower=NA),
                                                 dist.name='unif',
                                                 parameters=list(min=min, max=max),
                                                 var.name=var.name,
                                                 support = Continuous.Support(min, max),
                                                 mean.value=(max+min)/2,
                                                 variance.value=(max-min)^2/12)
    }
}

#'@title Create a Log-Uniform Distribution Object
#'
#'@param min,max The bounds of the distribution (on the original scale, NOT the log scale). If min is 0 or max is Inf, returns an improper log-uniform distribution (which returns 1/x for all densities, but cannot evaluate cdfs, quantiles, or generate random samples)
#'@inheritParams Normal.Distribution
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Loguniform.Distribution <- function(min=0, max=1, var.name=NULL)
{
    if (min < 0)
        stop("A log-uniform distribution cannot have a min less than zero")

    if (is.infinite(log(min)) || is.infinite(log(max)))
    {
        Univariate.Canonical.Distribution(name=paste0("Improper Log-Uniform(", min, ", ", max, ")"),
                                                 dist.name='NA',
                                                 density.function.name='improper.uniform.density',
                                                 parameters=list(),
                                                 var.name=var.name,
                                          support = Continuous.Support(min, max),
                                                 is.improper=T,
                                                 transformation='log')
    }
    else
    {
        Univariate.Canonical.Distribution(name=get.canonical.distribution.name('uniform', min, max, lower=min, upper=max, min.lower=NA, max.lower=NA),
                                                 dist.name='unif',
                                                 parameters=list(min=min, max=max),
                                                 var.name=var.name,
                                          support = Continuous.Support(min, max),
                                                 transformation='log')
    }
}

#'@title Create a Logit-Uniform Distribution Object
#'
#'@param min,max The bounds of the distribution (on the original scale, NOT the logit scale). If min is 0 or max is 1, returns an improper logit-uniform distribution (which returns 1/x/(1-x) for all densities, but cannot evaluate cdfs, quantiles, or generate random samples)
#'@inheritParams Normal.Distribution
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Logituniform.Distribution <- function(min=0, max=1, var.name=NULL)
{
    if (min < 0)
        stop("A logit-uniform distribution cannot have a min less than zero")
    if (max > 1)
        stop("A logit-uniform distribution cannot have a max greater than one")

    if (min==0 || max==1)
    {
        Univariate.Canonical.Distribution(name=paste0("Improper Logit-Uniform(", min, ", ", max, ")"),
                                                 dist.name='NA',
                                                 density.function.name='improper.uniform.density',
                                                 parameters=list(),
                                                 var.name=var.name,
                                          support = Continuous.Support(min, max),
                                                 is.improper=T,
                                                 transformation='logit')
    }
    else
    {
        Univariate.Canonical.Distribution(name=get.canonical.distribution.name('uniform', min, max, lower=min, upper=max, min.lower=NA, max.lower=NA),
                                                 dist.name='unif',
                                                 parameters=list(min=min, max=max),
                                                 var.name=var.name,
                                          support = Continuous.Support(min, max),
                                                 transformation='logit')
    }
}


#'@title Create a Univariate Binomial Distribution Object
#'
#'@param n,p The size and probability parameters to the distribution (AKA size, prob in dbinom)
#'@param lower,upper The lower and the upper bounds
#'@inheritParams Normal.Distribution
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Binomial.Distribution <- function(n=10, p=0.5, lower=0, upper=n, var.name=NULL)
{
    if (lower < 0)
        stop("A binomial distribution cannot have a lower bound less than zero")
    if (upper > n)
        stop("A binomial distribution cannot have a upper bound greater than n")

    if (p<0 || p>1)
        stop("p must be between 0 and 1")
    if (n<0)
        stop("n cannot be negative")

    if (lower==0 && upper==n)
    {
        mean.value = n*p
        variance.value = n*p*(1-p)
    }
    else
        mean.value = variance.value = as.numeric(NA)

    Univariate.Canonical.Distribution(name=get.canonical.distribution.name('binomial', n, p, lower=lower, upper=upper, min.lower=0, max.upper=n),
                                             dist.name='binom',
                                             parameters=list(size=n, prob=p),
                                             var.name=var.name,
                                      support = Integer.Support(lower, upper),
                                             mean.value=mean.value,
                                             variance.value=variance.value)
}

#'@title Create a Univariate Bernoulli Distribution Object
#'
#'@param p The probability of success in the distribution
#'@inheritParams Normal.Distribution
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@export
Bernoulli.Distribution <- function(p=0.5, var.name=NULL)
{
    if (p<0 || p>1)
        stop("p must be between 0 and 1")

    Univariate.Canonical.Distribution(name=paste0('bernoulli(', p, ')'),
                                             dist.name='binom',
                                             parameters=list(size=1, prob=p),
                                             var.name=var.name,
                                             support = Integer.Support(0,1),
                                             mean.value=p,
                                             variance.value=p*(1-p))
}
