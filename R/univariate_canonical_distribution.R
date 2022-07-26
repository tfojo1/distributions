
##--------------------------------------##
##-- CLASS DEFINITION and CONSTRUCTOR --##
##--------------------------------------##

#'@title The Univariate_Canonical_Distribution class
#'
#'@description A class that represents a univariate canonical distribution, such as a normal or beta distribution
#'
#'@slot name A short, descriptive name for the distribution
#'@slot density.function.name,cdf.function.name,quantile.function.name Names of the (elsewhere defined) functions that calculate density, cdf, and quantiles (eg 'dnorm', 'pnorm', and 'qnorm)
#'@slot parameters A list of parameters that would be passed to density.function, cdf.function, and quantile.function. For example, for a normal(0,1), parameters=list(mean=0, sd=1)
#'@slot p.in.bounds,log.p.in.bounds The cumulative probability and of the log cumulative probability between the lower and upper bounds
#'@slot p.lower.bound The value of the cdf evaluated at the lower bound
#'@slot transformation An object of class \link{transformation}, defining the transformation applied to the variable before applying the distribution (eg, )
#'@slot mean.value,variance If known, the mean and variance of the distribution. May be NA
#'
#'@seealso The family of functions that create specific distributions, such as \link{Normal.Distribution}
#'
#'@name Univariate_Canonical_Distribution
#'@rdname Univariate_Canonical_Distribution
#'@aliases Univariate_Canonical_Distribution-class
#'@exportClass Univariate_Canonical_Distribution
#'@export
setClass('Univariate_Canonical_Distribution',
         contains='Distribution',
         representation=list(name='character',
                             density.function.name='character',
                             cdf.function.name='character',
                             quantile.function.name='character',
                             parameters='list',
                             p.in.bounds='numeric',
                             log.p.in.bounds='numeric',
                             p.lower.bound='numeric',
                             transformation='transformation',
                             mean='numeric',
                             variance='numeric'))


#'@title Create a distribution object based off of a canonical, univariate distribution
#'
#'@description
#'
#'@param name A short, descriptive name for the distribution
#'@param dist.name The root name of the distribution as used by d/r/p functions in R. For example, to create a normal distribution, dist.name='norm' (like dnorm, rnorm, pnorm). To create a beta distribution, dist.name='beta' (dbeta, rbeta, qbeta). This argument is only necessary if density.function, cdf.function, and quantile.function are not passed explicitly
#'@param parameters A list of parameters that would be passed to d/r/p functions in R. For example, for a normal(0,1), parameters=list(mean=0, sd=1)
#'@param var.name The name of the variable in this distribution
#'@param support A \link{support} object that
#'@param density.function,cdf.function,quantile.function The functions that calculate the density, cdf, and quantile respectively. These should take the standard parameters of such functions in R. density.function should take 'x', distribution-specific parameters, and 'log'. cdf.function should take 'q', distribution-specific parameters, 'lower.tail' and 'log.p'. quantile.function should take 'p', distribution-specific parameters, and 'lower.tail' and 'log.p'
#'@param transformation A transformation object, if the named distribution operates on a transformation of the random variable (can also pass NULL for no transformation or the name of a predefined transformation - see \code{\link{get.defined.transformation}}). For example, to create a log-normal distribution, use dist.name='norm' with transformation='log
#'@param is.improper Logicals indicating whether the distribution is improper/discrete
#'@param mean.value,variance If known, these values will be returned for calls to get.means, get.sds, get.covariance.matrix. If passed NA, random sampling will be used to estimate distribution parameters if requested
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@seealso \code{\link{get.defined.transformation}}, \code{\link{create.transformation}}
#'
#'@export
Univariate.Canonical.Distribution <- function(name,
                                              dist.name,
                                              parameters,
                                              var.name=NULL,
                                              support=Continuous.Support(),
                                              density.function.name=paste0('d', dist.name),
                                              cdf.function.name=paste0('p', dist.name),
                                              quantile.function.name=paste0('q', dist.name),
                                              transformation=NULL,
                                              is.improper=F,
                                              mean.value=as.numeric(NA),
                                              variance.value=as.numeric(NA))
{
    new('Univariate_Canonical_Distribution',
        name=name,
        parameters=parameters,
        var.name=var.name,
        support=support,
        density.function.name=density.function.name,
        cdf.function.name=cdf.function.name,
        quantile.function.name=quantile.function.name,
        transformation=transformation,
        is.improper=is.improper,
        mean.value=mean.value,
        variance.value=variance.value)
}

setMethod('initialize',
          signature(.Object='Univariate_Canonical_Distribution'),
def=function(.Object,
             name,
             parameters,
             var.name=NULL,
             support=Continuous.Support(),
             density.function.name,
             cdf.function.name,
             quantile.function.name,
             transformation=NULL,
             is.improper=F,
             mean.value=as.numeric(NA),
             variance.value=as.numeric(NA))
{
    .Object = callNextMethod(.Object,
                             var.names=var.name,
                             n.var=1,
                             support=support,
                             is.improper=is.improper)


    #-- Set up transformation functions --#

    if (is.null(transformation))
        .Object@transformation = get.defined.transformation('identity')
    else if (is(transformation, 'character'))
        .Object@transformation = get.defined.transformation(transformation, throw.error.if.no.match=T)
    else if (is(transformation, 'transformation'))
        .Object@transformation = transformation
    else
        stop("transformation must be either a transformation object, a character name of a predefined transformation, or NULL")


    #-- Check parameters --#
    if (!is(parameters, 'list'))
        stop("'parameters' must be a list")

    #-- Set Up Bound Probabilities --#
    if (!.Object@support@n.var==1)
        stop("The support for a Univariate_Canonical_Distribution must be univariate")

    if (!.Object@support@is.contiguous)
        stop("To create a Canonical_Univariate_Distribution, the support must be over a contiguous range")

    bounds = get.support.bounds(support)
    if (bounds[1]==bounds[2])
        stop("The lower bounds of support cannot equal the upper bound")

    transformed.lower.bound = .Object@transformation@transform(bounds[1])
    if (is.na(transformed.lower.bound))
        stop(paste0("The transformation function ('",
                    .Object@transformation@name,
                    "') produces NA when evaluated at the lower bound (",
                    bounds[1], ")"))
    transformed.upper.bound = .Object@transformation@transform(bounds[2])
    if (is.na(transformed.upper.bound))
        stop(paste0("The transformation function ('",
                    .Object@transformation@name,
                    "') produces NA when evaluated at the upper bound (",
                    bounds[2], ")"))

    if (!is.improper)
        p.bounds = do.call(cdf.function.name,
                                args=c(list(q=c(transformed.lower.bound, transformed.upper.bound)), parameters))
    else
        p.bounds = c(0,1)

    .Object@p.in.bounds = p.bounds[2] - p.bounds[1]
    .Object@log.p.in.bounds = log(.Object@p.in.bounds)
    .Object@p.lower.bound = p.bounds[1]

    #-- Set Other Member Variables --#

    .Object@name = name
    .Object@density.function.name = density.function.name
    .Object@cdf.function.name = cdf.function.name
    .Object@quantile.function.name = quantile.function.name
    .Object@parameters = parameters
    .Object@mean = mean.value
    .Object@variance = variance.value

    #-- Return --#
    .Object
})

##------------------------------##
##-- HELPERS for CONSTRUCTORS --##
##------------------------------##

improper.uniform.density <- function(x, ...)
{
    rep(1, length(x))
}

get.canonical.distribution.name <- function(name,
                                ...,
                                lower,
                                upper,
                                min.lower=-Inf,
                                max.upper=Inf)
{
    parameters = list(...)

    rv = paste0(name, "(",
                paste0(unlist(parameters), collapse=', '),
                ")")

    if ((!is.na(min.lower) && lower != min.lower) ||
        (!is.na(max.upper) && upper != max.upper))
        rv = paste0(rv, ", bounded on [",
                    lower, ", ", upper, "]")

    rv
}

##----------------------------##
##-- METHODS IMPLEMENTATION --##
##----------------------------##


setMethod('do.calculate.density',
          signature(dist='Univariate_Canonical_Distribution'),
def=function(dist, x, log=F, n.sim=1000)
{
    x = as.numeric(x)

    in.bounds = is.supported(dist@support, x)

    orig.x = x
    x[in.bounds] = dist@transformation@transform(x[in.bounds])

    if (log)
        rv = rep(-Inf, length(x))
    else
        rv = rep(0, length(x))


    if (log)
    {
        log.derivative = dist@transformation@log.abs.transformation.derivative(orig.x[in.bounds])
        log.derivative[log.derivative==Inf] = 0

        rv[in.bounds] =  do.call(dist@density.function.name,
                                 args=c(list(x=x[in.bounds], log=log), dist@parameters)) +
                            log.derivative - dist@log.p.in.bounds
    }
    else
    {
        derivative = abs(dist@transformation@transformation.derivative(orig.x[in.bounds]))
        derivative[derivative==Inf] = 0

        rv[in.bounds] = do.call(dist@density.function.name,
                                args=c(list(x=x[in.bounds], log=log), dist@parameters)) *
                            derivative / dist@p.in.bounds
    }

    rv
})


setMethod('do.calculate.cdf',
          signature(dist='Univariate_Canonical_Distribution'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000)
{
    q = as.numeric(q)
    if (dist@is.improper)
    {
        return(rep(NA, length(q)))
    }

    bounds = get.support.bounds(dist@support)
    below.lower = q <= bounds[1]
    above.upper = q >= bounds[2]
    in.bounds = is.supported(dist@support, q)

    q[in.bounds] = dist@transformation@transform(q[in.bounds])

    rv = rep(0, length(q))
    rv[above.upper] = 1

    if (!lower.tail)
        rv = 1-rv

    if (log.p)
        rv = log(rv)

    rv[in.bounds] = (do.call(dist@cdf.function.name,
                            args=c(list(q=q[in.bounds], lower.tail=lower.tail, log.p=log.p), dist@parameters)) -
                         dist@p.lower.bound)

    if (log.p)
        rv[in.bounds] = rv[in.bounds] - dist@log.p.in.bounds
    else
        rv[in.bounds] = rv[in.bounds] / dist@p.in.bounds

    rv
})


setMethod('do.get.quantiles',
          signature(dist='Univariate_Canonical_Distribution'),
def=function(dist, p, lower.tail=T, log.p=F, n.sim=1000)
{
    p = as.numeric(p)
    if (dist@is.improper)
    {
        return(rep(NA, length(p)))
    }

    if (log.p)
        p = exp(p)
    if (!lower.tail)
        p = 1-p

    p = dist@p.lower.bound + p * dist@p.in.bounds

    rv = do.call(dist@quantile.function.name,
            args=c(list(p=p, lower.tail=T, log.p=F), dist@parameters))

    rv = dist@transformation@reverse.transform(rv)

    rv
})

setMethod('do.generate.random.samples',
          signature(dist='Univariate_Canonical_Distribution'),
def=function(dist, n)
{

    if (dist@is.improper)
    {
        return(rep(NA, n))
    }

    p = runif(n)
    do.get.quantiles(dist, p=p)
})



setMethod('do.get.covariance.matrix',
          signature(dist='Univariate_Canonical_Distribution'),
def=function(dist, n.sim=1000)
{
    if (any(is.na(dist@variance)))
        callNextMethod(dist, n.sim=n.sim)
    else
        dist@variance
})

setMethod('do.get.means',
          signature(dist='Univariate_Canonical_Distribution'),
def=function(dist, n.sim=200)
{
    if (any(is.na(dist@mean)))
        callNextMethod(dist, n.sim=n.sim)
    else
        dist@mean
})

setMethod('get.description',
          signature(object='Univariate_Canonical_Distribution'),
          def=function(object){
                object@name
          })

