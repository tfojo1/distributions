
##--------------------------------------##
##-- CLASS DEFINITION and CONSTRUCTOR --##
##--------------------------------------##

#'@title The Canonical_Mixture_Distribution class
#'
#'@description A class that represents a univariate mixture of different instances of a canonical distribution (the same distribution, but with different parameters)
#'
#'@slot name A short, descriptive name for the distribution
#'@slot n.components The number of component distributions in this mixture
#'@slot weights The weights for the component distributions in the mixture (sum(weights)=1)
#'@slot density.function.name,cdf.function.name,quantile.function.name Names of the (elsewhere defined) functions that calculate density, cdf, and quantiles (eg 'dnorm', 'pnorm', and 'qnorm)
#'@slot parameters A list of vectors; each element of the list represents a parameter to the distribution, and the values of the vector represent the parameter values for each component of the mixture
#'@slot p.in.bounds,log.p.in.bounds A vector cumulative probabilities and of the log cumulative probabilities between the lower and upper bounds for each component of the mixture
#'@slot p.lower.bound A vector of values of the cdf evaluated at the lower bound for each components of the distribution
#'@slot transformation An object of class \link{transformation}, defining the transformation applied to the variable before applying the distribution (eg, )
#'@slot means,variances If known, the means and variances of each component of the distribution. May be NA
#'
#'@seealso The family of functions that create specific mixture distributions, such as \link{normal.mixture}
#'
#'@name Canonical_Mixture_Distribution
#'@rdname Canonical_Mixture_Distribution
#'@aliases Canonical_Mixture_Distribution-class
#'@exportClass Canonical_Mixture_Distribution
#'@export
setClass('Canonical_Mixture_Distribution',
         contains='Distribution',
         representation=list(name='character',
                             n.components='integer',
                             weights='numeric',
                             log.weights='numeric',
                             density.function.name='character',
                             cdf.function.name='character',
                             quantile.function.name='character',
                             parameters='list',
                             p.in.bounds='numeric',
                             log.p.in.bounds='numeric',
                             p.lower.bound='numeric',
                             transformation='transformation',
                             means='numeric',
                             variances='numeric'))


#'@title Create a mixture distribution object based off of a canonical, univariate distribution
#'
#'@description
#'
#'@param name A short, descriptive name for the distribution
#'@param weights A weight for each component distribution in the mixture. The weights need not sum to 1. If a scalar is passed, all components are given the same weight
#'@param dist.name The root name of the distribution as used by d/r/p functions in R. For example, to create a normal distribution, dist.name='norm' (like dnorm, rnorm, pnorm). To create a beta distribution, dist.name='beta' (dbeta, rbeta, qbeta). This argument is only necessary if density.function, cdf.function, and quantile.function are not passed explicitly
#'@param parameters A list of vectors, where each element has parameter values that would be passed to d/r/p functions in R. For example, for a mixture of normal(0,1) and normal(1,1) distributions, parameters=list(mean=c(0,1), sd=1) or parameters=list(mean=c(0,1), sd=c(1,1))
#'@param var.name The name of the variable in this distribution
#'@param lower.bound,upper.bound If this is a bounded distribution, the lower and upper bounds
#'@param density.function,cdf.function,quantile.function The functions that calculate the density, cdf, and quantile respectively. These should take the standard parameters of such functions in R. density.function should take 'x', distribution-specific parameters, and 'log'. cdf.function should take 'q', distribution-specific parameters, 'lower.tail' and 'log.p'. quantile.function should take 'p', distribution-specific parameters, and 'lower.tail' and 'log.p'
#'@param transformation A transformation object, if the named distribution operates on a transformation of the random variable (can also pass NULL for no transformation or the name of a predefined transformation - see \code{\link{get.defined.transformation}}). For example, to create a log-normal distribution, use dist.name='norm' with transformation='log
#'@param is.improper,is.discrete Logicals indicating whether the distribution is improper/discrete
#'@param mean.values,variance.values If known, vectors of values for the means and variances of each component of the distribution
#'
#'@family Univariate Canonical Distribution Constructors
#'@family Distribution Constructors
#'
#'@seealso \code{\link{get.defined.transformation}}, \code{\link{create.transformation}}
#'
#'@export
Canonical.Mixture.Distribution <- function(name,
                                           dist.name,
                                           parameters,
                                           weights=1,
                                           var.name=NULL,
                                           support = Continuous.Support(),
                                           density.function.name=paste0('d', dist.name),
                                           cdf.function.name=paste0('p', dist.name),
                                           quantile.function.name=paste0('q', dist.name),
                                           transformation=NULL,
                                           is.improper=F,
                                           mean.values=as.numeric(NA),
                                           variance.values=as.numeric(NA))
{
    new('Canonical_Mixture_Distribution',
        name=name,
        parameters=parameters,
        var.name=var.name,
        support=support,
        density.function.name=density.function.name,
        cdf.function.name=cdf.function.name,
        quantile.function.name=quantile.function.name,
        transformation=transformation,
        is.improper=is.improper,
        mean.values=mean.values,
        variance.values=variance.values)
}
setMethod('initialize',
          signature(.Object='Canonical_Mixture_Distribution'),
def = function(.Object,
               name,
               parameters,
               weights=1,
               var.name=NULL,
               support = Continuous.Support(),
               density.function.name,
               cdf.function.name,
               quantile.function.name,
               transformation=NULL,
               is.improper=F,
               mean.values=as.numeric(NA),
               variance.values=as.numeric(NA))
{
    .Object = callNextMethod(.Object,
                             var.names=var.name,
                             n.var=1,
                             support=support,
                             is.improper=is.improper)


    #-- Check parameters and set up n.components --#
    n.components = 1

    if (!is(parameters, 'list'))
        stop("'parameters' must be a list")

    for (i in 1:length(parameters))
    {
        elem = parameters[[i]]
        if (length(elem) == 0)
            stop("The elements of 'parameters' cannot have length 0")
        n.components = max(n.components, length(elem))
        if (length(elem) != n.components && length(elem) != 1)
            stop("All the elements of 'parameters' must have the same length as the other elements or be of length 1")
    }

    .Object@n.components = as.integer(n.components)

    #-- Set up weights --#
    if (length(weights)==1)
        weights = rep(1, n.components)
    if (length(weights) != n.components)
        stop(paste0("There must be one value of weights for each component of the mixture (", n.components, ")"))

    .Object@weights = weights / sum(weights)
    .Object@log.weights = log(.Object@weights)

    #-- Set up transformation functions --#

    if (is.null(transformation))
        .Object@transformation = get.defined.transformation('identity')
    else if (is(transformation, 'character'))
        .Object@transformation = get.defined.transformation(transformation, throw.error.if.no.match=T)
    else if (!is(transformation, 'transformation'))
        stop("transformation must be either a transformation object, a character name of a predefined transformation, or NULL")

    #-- Set Up Bound Probabilities --#

    if (!.Object@support@n.var==1)
        stop("The support for a Canonical_Mixture_Distribution must be univariate")

    if (!.Object@support@is.contiguous)
        stop("To create a Canonical_Univariate_Distribution, the support must be over a contiguous range")

    bounds = get.support.bounds(.Object@support)
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
        p.lower.bound = sum(.Object@weights * do.call(cdf.function.name,
                                args=c(list(q=transformed.lower.bound), parameters)))
    else
        p.lower.bound = rep(0, .Object@n.components)

    if (!is.improper)
        p.upper.bound = sum(.Object@weights * do.call(cdf.function.name,
                                args=c(list(q=transformed.upper.bound), parameters)))
    else
        p.upper.bound = rep(1, .Object@n.components)

    .Object@p.in.bounds = p.upper.bound - p.lower.bound
    .Object@log.p.in.bounds = log(.Object@p.in.bounds)
    .Object@p.lower.bound = p.lower.bound

    #-- Set Other Member Variables --#

    .Object@name = name
    .Object@density.function.name = density.function.name
    .Object@cdf.function.name = cdf.function.name
    .Object@quantile.function.name = quantile.function.name
    .Object@parameters = parameters
    .Object@means = mean.values
    .Object@variances = variance.values

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

get.canonical.mixture.name <- function(component.name,
                                       n.components,
                                       lower,
                                       upper,
                                       min.lower=-Inf,
                                       max.upper=Inf)
{
    rv = paste0("A mixture of ", n.components, " ",
                component.name, " distributions")

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
          signature(dist='Canonical_Mixture_Distribution'),
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

    #one row for each component, one column for each value of x
    component.densities = matrix(do.call(dist@density.function.name,
                           args=c(list(x=rep(x[in.bounds], each=dist@n.components), log=log), dist@parameters)),
                   nrow=dist@n.components)

    densities = apply(component.densities*dist@weights/dist@p.in.bounds, 2, sum)
    if (log)
    {
        log.derivative = dist@transformation@log.abs.transformation.derivative(orig.x[in.bounds])
        log.derivative[log.derivative==Inf] = 0

        log.densities = apply(component.densities + dist@log.weights - dist@log.p.in.bounds, 2, log.sum.exp)

        rv[in.bounds] = log.densities + log.derivative
    }
    else
    {
        derivative = abs(dist@transformation@transformation.derivative(orig.x[in.bounds]))
        derivative[derivative==Inf] = 0

        densities = colSums(component.densities * dist@weights / dist@p.in.bounds)

        rv[in.bounds] = densities * derivative
    }

    rv
})



setMethod('do.calculate.cdf',
          signature(dist='Canonical_Mixture_Distribution'),
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

    #one row for each component, one column for each value of q
    component.cdfs = matrix(do.call(dist@cdf.function.name,
                                         args=c(list(q=rep(q[in.bounds], each=dist@n.components),
                                                     lower.tail=lower.tail, log.p=log.p), dist@parameters)),
                                 nrow=dist@n.components)

    densities = apply(component.cdfs*dist@weights/dist@p.in.bounds, 2, sum)
    if (log.p)
        rv[in.bounds] = apply(component.cdfs + dist@log.weights - dist@log.p.in.bounds, 2, log.sum.exp)
    else
        rv[in.bounds] = colSums(component.cdfs * dist@weights / dist@p.in.bounds)

    rv
})


setMethod('do.get.quantiles',
          signature(dist='Canonical_Mixture_Distribution'),
def=function(dist, p, lower.tail=T, log.p=F, n.sim=1000)
{
    p = as.numeric(p)
    if (dist@is.improper)
        return(rep(NA, length(p)))

#    if (log.p)
#        p = exp(p)
#    if (!lower.tail)
#        p = 1-p

    component.qs = matrix(do.call(dist@quantile.function.name,
                           args=c(list(p=dist@p.lower.bound + rep(p, each=dist@n.components) * dist@p.in.bounds,
                                       lower.tail=lower.tail, log.p=log.p), dist@parameters)),
                          nrow=dist@n.components)

    rv = sapply(1:length(p), function(i){

        objective.fn = function(q){
            p.for.q = do.calculate.cdf(dist, q, n.sim=n.sim, lower.tail=lower.tail, log.p=log.p)
            (p[i] - p.for.q)^2
        }

        if (min(component.qs[,i]) == max(component.qs[,i]))
            component.qs[1,i]
        else
        {
            opt = optimize(f=objective.fn, interval=c(min(component.qs[,i]), max(component.qs[,i])))
            opt$minimum
        }
    })

    rv = dist@transformation@reverse.transform(rv)

    rv
})

setMethod('do.generate.random.samples',
          signature(dist='Canonical_Mixture_Distribution'),
def=function(dist, n)
{
    if (dist@is.improper)
        return(rep(NA, n))

    p = runif(n)

    indices = sample(1:dist@n.components, prob=dist@weights, replace=T)
    params = c(list(p=dist@p.lower.bound[indices] + p * dist@p.in.bounds[indices]),
               lapply(dist@parameters, function(one.param){
                   if (length(one.param)==1)
                       one.param
                   else
                       one.param[indices]
               }))

    rv = do.call(dist@quantile.function.name, params)

    rv = dist@transformation@reverse.transform(rv)
    rv
})



setMethod('do.get.covariance.matrix',
          signature(dist='Canonical_Mixture_Distribution'),
def=function(dist, n.sim=1000)
{
    if (any(is.na(dist@variances)) || any(is.na(dist@means)))
        callNextMethod(dist, n.sim=n.sim)
    else
        sum(dist@weights * (dist@variances + dist@means^2)) - do.get.means(dist, n.sim=n.sim)^2
})

setMethod('do.get.means',
          signature(dist='Canonical_Mixture_Distribution'),
def=function(dist, n.sim=200)
{
    if (any(is.na(dist@means)))
        callNextMethod(dist, n.sim=n.sim)
    else
        sum(dist@weights * dist@means)
})


