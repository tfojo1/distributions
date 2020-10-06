
##--------------------------------------##
##-- CLASS DEFINITION and CONSTRUCTOR --##
##--------------------------------------##

setClassUnion('character_or_null', members = c('NULL', 'character'))

#'@title The abstract superclass for all Distribution objects
#'
#'@description Distribution objects represent a distribution of one or more variables
#'
#'@seealso Specific subclasses which can be instantiated:
#' \itemize{
#'     \item{\link{Empiric_Distribution}}
#'     \item{\link{Univariate_Canonical_Distribution}}
#'     \item{\link{Joint_Independent_Distributions}}
#' }
#'
#'@slot var.names A character vector of the names of the variables in the distribution (optional)
#'@slot n.var The number of variables in the distribution
#'@slot support A \link{Support} object defining the support for each variable in the distribution
#'@slot is.improper A logical vectors indicating whether each variable in the distribution is improper
#'
#'
#'@name Distribution
#'@rdname Distribution
#'@aliases Distribution-class
#'@exportClass Distribution
#'@export
setClass("Distribution",
         representation=list(var.names='character_or_null',
                             n.var='integer',
                             support='Support',
                             is.discrete='logical',
                             is.improper='logical'))

#'@title A constructor
#'
#'@description This function is provided for use by developers, and should not be called directly. This function should be called by the subclass constructor. It performs consistency checks on the basic slots for a Distribution object
#'
#'@param subclass The name of the subclass to be instantiated with a call to 'new'
#'@param var.names,n.var The name and number of variables in the distribution. var.names may be NULL
#'@param lower.bounds,upper.bounds Vectors of lower and upper bounds for each of the variables in the distribution. May be NA
#'@param is.discrete,is.improper Logical vectors indicating whether each variable in the distribution is discrete or improper. If a scalar logical is passed, it is assumed to be the value for all variables in the distribution
#'@param ... Other arguments to be passed to the call to 'new'
#'
#'@keywords internal
#'@export
setMethod('initialize',
          signature(.Object='Distribution'),
          def=function(.Object,
                       var.names=NULL,
                       n.var=length(var.names),
                       support=Multivariate.Support(lapply(1:n.var, function(i){Continuous.Support()})),
                       is.improper=F)
{
    #Check names
    if (!is.null(var.names) && !is(var.names, 'character'))
        stop("'var.names' must be either a character vector or NULL")
    if (!is.null(var.names) && any(is.na(var.names) && !all(is.na(var.names))))
        stop("'var.names' cannot contain NA values")

    named.vars = !is.null(var.names) && length(var.names)>0 && all(!is.na(var.names))
    if (named.vars)
    {
        if (any(table(var.names)>2))
            stop("var.names must be unique")

        .Object@var.names = var.names
        .Object@n.var = length(var.names)
    }
    else
    {
        if (is.null(n.var) || (!is(n.var, 'integer') && !is(n.var, 'numeric')) ||
            length(n.var)!=1 || is.na(n.var) ||
            round(n.var != n.var) || n.var <= 0 || is.infinite(n.var) || round(n.var)!=n.var)
            stop("If var.names are not specified, n.var must be a positive integer")

        .Object@var.names = NULL
        .Object@n.var = as.integer(n.var)
    }

    #Check support

    if (!is(support, 'Support'))
        stop("'support' must be an object of class support")
    if (support@n.var == 1 && support@n.var < .Object@n.var)
        support = Multivariate.Support(lapply(1:.Object@n.var, function(i){support}))
    if (support@n.var != .Object@n.var)
        stop(paste0("'support' must describe n.var (", .Object@n.var, ") variables"))

    .Object@support = support
    .Object@is.discrete = support@is.discrete

    #Check lengths of is.improper

    if (length(is.improper)!=n.var)
    {
        if (length(is.improper)==1)
            is.improper = rep(is.improper, n.var)
        else
            stop(paste0("is.improper must have a value for each of the ", n.var, " variables in the distribution"))
    }
    .Object@is.improper = is.improper

    # Return
    .Object
})

##-------------------------------------------------##
##--        DENSITY and RANDOM GENERATION        --##
##-- (Must be implemented at the subclass level) --##
##-------------------------------------------------##

#'@title Calculate the density for a distribution
#'
#'@param dist An object of class Distribution or a subclass
#'@param x The values of the parameters at which to calculate the density. May be either a vector, if the density at a single point is desired, or a matrix where each row represents one point and each column represents a variable in the distribution
#'@param log A logical indicating whether to return the density on the log scale
#'@param n.sim If the values are to be calculated by random sample generation, how many random samples to generate
#'
#'@return If x is a vector, a scalar density. If x is a matrix, a numeric vector of densities, one for each row in x
#'
#'@details Note to developers: This function should not be overridden directly. Insteady override \code{\link{do.calculate.density}}
#'
#'@keywords internal
#'@export
setGeneric('calculate.density',
           def=function(dist, x, log=F, n.sim=1000){standardGeneric("calculate.density")})
setMethod('calculate.density',
          signature(dist='Distribution'),
def=function(dist, x, log=F, n.sim=1000)
{
    x = match.variables(dist, x)

    rv = do.calculate.density(dist=dist, x=x, log=log, n.sim=n.sim)

    as.numeric(rv)
})


#'@title Internal function for subclasses to calculate the density for a distribution
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{calculate.density}}
#'
#'@param x A matrix where each row represents one point and each column represents a variable in the distribution (the wrapper matches variable names before passing to this function)
#'
#'@inherit calculate.density
#'@details The wrapper \code{\link{calculate.density}} takes care of matching variable names and of formatting the return value
#'
#'@export
setGeneric('do.calculate.density',
           def=function(dist, x, log=F, n.sim=1000){standardGeneric("do.calculate.density")})
setMethod('do.calculate.density',
def=function(dist, x, log=F, n.sim=1000)
{
    stop(paste0("How to calculate density has not been defined for this ", class(dist)[1]))
})

#'@title Generate random samples from a distribution
#'
#'@inheritParams calculate.density
#'@param n The number of random samples to be generated
#'
#'@return If n==1, a numeric vector with one value for each variable in the distribution. If there is only one variable in the distribution, a numeric vector with n values. If n>1 and there is more than one variable in the distribution, a matrix with n rows, where each row represents one random sample, and each column represents a variable in the distribution
#'@details Note to developers: This function should not be overridden directly. Insteady override \code{\link{do.generate.random.samples}}
#'
#'@export
setGeneric('generate.random.samples',
           def=function(dist, n){standardGeneric("generate.random.samples")})
setMethod('generate.random.samples',
          signature(dist='Distribution'),
def=function(dist, n)
{
    print.improper.warnings(dist, description="generate random samples from")

    rv = do.generate.random.samples(dist=dist, n=n)

    dim(rv) = c(sample=n, variable=dist@n.var)
    dimnames(rv) = list(sample=NULL, variable=dist@var.names)

    if (n==1)
        rv[1,]
    else if (dist@n.var==1)
        rv[,1]
    else
        rv
})

#'@title Internal function for subclasses to generate random samples from a distribution
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{generate.random.samples}}
#'
#'@inherit generate.random.samples
#'@details The wrapper \code{\link{generate.random.samples}} takes care of formatting the return value and generating warnings for any improper variables
#'
#'@keywords internal
#'@export
setGeneric('do.generate.random.samples',
           def=function(dist, n){standardGeneric("do.generate.random.samples")})
setMethod('do.generate.random.samples',
          signature(dist='Distribution'),
def=function(dist, n)
{
    stop(paste0("How to generate random samples has not been defined for this ", class(dist)[1]))
})

##------------------------------------------------##
##--      MEDIANS, QUANTILES, and INTERVALS     --##
##-- (May be overwritten at the subclass level) --##
##------------------------------------------------##

##-- FUNCTION DEFINITIONS --##

#'@title Get quantiles for the variables in a distribution
#'
#'@inheritParams calculate.density
#'@param p A vector of probabilities for which to get the quantiles
#'@param lower.tail If true, the returned values are the greatest values x for each variable such that \eqn{P(X \leq x) \geq p}. If fale, the returned values are the smallest values x for each variable such that \eqn{P(X > x) \geq p}. This can be either a scalar (in which case the value is applied to all quantiles), or a vector with one logical value for each value of p
#'
#'@return If there is only one variable in the distribution, a vector of quantiles for each value of p. If p is length one, a vector of quantiles for each variable in the distribution. If there is more than one variable in the distribution and p>1, a matrix with one row for each value of p, where each column represents a variable in the distribution
#'
#'@details Note to developers: This function should not be overridden directly. Insteady override \code{\link{do.get.quantiles}}
#'
#'@export
setGeneric('get.quantiles',
           def=function(dist, p, lower.tail=T, log.p=F, n.sim=1000){standardGeneric("get.quantiles")})
setMethod('get.quantiles',
          signature(dist='Distribution'),
def=function(dist, p, lower.tail=T, log.p=F, n.sim=1000)
{

    print.improper.warnings(dist, description="calculate quantiles for")

    if (!is(p, 'numeric') && !is(p, 'integer'))
        stop("p must be a numeric vector")

    rv = do.get.quantiles(dist=dist, p=p, lower.tail=lower.tail, log.p=log.p, n.sim=n.sim)

    if (dist@n.var==1)
    {
        rv = as.numeric(rv)
        names(rv) = as.character(p)
    }
    else
    {
        if (length(p)==1)
        {
            rv = as.numeric(rv)
            names(rv) = dist@var.names
        }
        else
        {
            dim(rv) = c(quantile=length(p), variable=dist@n.var)
            dimnames(rv) = list(quantile=as.character(p), variable=dist@var.names)
        }
    }

    rv
})

#'@title Get medians of the variables in a distribution
#'
#'@inheritParams get.quantiles
#'
#'@return A numeric vector with a median for each variable in the distribution
#'@details Note to developers: This function should not be overridden directly. Insteady override \code{\link{do.get.medians}}
#'
#'@export
setGeneric('get.medians',
           def=function(dist, n.sim=1000){standardGeneric("get.medians")})
setMethod('get.medians',
          signature(dist='Distribution'),
def=function(dist, n.sim=1000)
{
    print.improper.warnings(dist, description="calculate medians for")

    rv = do.get.medians(dist=dist, n.sim=n.sim)

    rv = as.numeric(rv)
    names(rv) = dist@var.names

    rv
})

#'@title Get equal-tailed intervals (confidence/credible intervals) for the variables in a distribution
#'
#'@inherit get.intervals
#'
#'@details Note to developers: This function should not be overridden directly. Insteady override \code{\link{do.get.equal.tailed.intervals}}
#'
#'@export
setGeneric('get.equal.tailed.intervals',
           def=function(dist, coverage=0.95, n.sim=1000){standardGeneric("get.equal.tailed.intervals")})
setMethod('get.equal.tailed.intervals',
          signature(dist='Distribution'),
def=function(dist, coverage=0.95, n.sim=1000)
{
    get.intervals(dist, type='equal-tailed', coverage=coverage, n.sim=n.sim)
})

#'@title Get highest-density intervals (confidence/credible intervals) for the variables in a distribution
#'
#'@inherit get.intervals
#'
#'@details Note to developers: This function should not be overridden directly. Insteady override \code{\link{do.get.highest.density.intervals}}
#'@export
setGeneric('get.highest.density.intervals',
           def=function(dist, coverage=0.95, n.sim=1000){standardGeneric("get.highest.density.intervals")})
setMethod('get.highest.density.intervals',
          signature(dist='Distribution'),
def=function(dist, coverage=0.95, n.sim=1000)
{
    get.intervals(dist, type='highest-density', coverage=coverage, n.sim=n.sim)
})

#'@title Get intervals (confidence/credible intervals) for the variables in a distribution
#'
#'@inheritParams get.quantiles
#'@param type Either 'equal-tailed' or 'highest-density' depending on the type of interval desired (equal-tailed is generally easier to compute)
#'@param coverage A vector of probabilities (from 0 to 1) which the intervals should cover
#'
#'@return If length(coverage)==1, and the distribution has one variable, returns a numeric vector with two values: the lower and upper bound of the interval.
#'If length(coverage)==1 and the distribution has more than one variable, returns a matrix with one row for the lower bound and one row for the upper, and one column for each variable.
#'If length(coverage)>1 and the distribution has one variable, returns a matrix with one row for the lower bound and one row for the upper bound, and one column for each value in coverage.
#'If length(coverage)>1, and the distribution has more than one variable returns a 3-dimensional array; the first dimension has values 'lower' and 'upper'. The second dimension has one value for each variable. The third dimension has one value for each value of coverage.
#'
#'@seealso \code{\link{get.equal.tailed.intervals}}, \code{\link{get.highest.density.intervals}}
#'@details Note to developers: this function should not be overridden. Override \code{\link{do.get.equal.tailed.intervals}} and \code{\link{do.get.highest.density.intervals}} instead
#'
#'@export
setGeneric('get.intervals',
           def=function(dist, type=c('equal-tailed','highest-density')[1], coverage=0.95, n.sim=1000){standardGeneric('get.intervals')})
setMethod('get.intervals',
          signature(dist="Distribution"),
function(dist, type=c('equal-tailed','highest-density')[1], coverage=0.95, n.sim=1000)
{
    print.improper.warnings(dist, description="calculate intervals for")

    if (type=='equal-tailed' || type=='equal.tailed')
        rv = do.get.equal.tailed.intervals(dist=dist, coverage=coverage, n.sim=n.sim)
    else if (type=='highest-density' || type=='highest.density')
        rv = do.get.highest.density.intervals(dist=dist, coverage=coverage, n.sim=n.sim)
    else
        stop("type must be either 'equal-tailed' or 'highest-density")

    dim.names = list(bound=c('lower','upper'),
                     variable=dist@var.names,
                     coverage=as.character(coverage))

    dim(rv) = c(bound=2,  variable=dist@n.var, coverage=length(coverage))
    dimnames(rv) = dim.names

    if (length(coverage)==1)
        rv[,,1]
    else if (dist@n.var==1)
        rv[,1,]
    else
        rv

})

##-- DEFAULT IMPLEMENTATIONS --##

#'@title Internal function for subclasses to get quantiles for the variables in a distribution
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{get.quantiles}}
#'
#'@inherit get.quantiles
#'
#'@details The default implementation takes random samples and pulls the sample quantiles from them.
#'The wrapper \code{\link{get.quantiles}} takes care of formatting the return value and generating warnings for any improper variables
#'
#'@keywords internal
#'@export
setGeneric('do.get.quantiles',
           def=function(dist, p, lower.tail=T, log.p=F, n.sim=1000){standardGeneric('do.get.quantiles')})
setMethod('do.get.quantiles',
          signature(dist='Distribution'),
def=function(dist, p, lower.tail=T, log.p=F, n.sim=1000)
{
    rands = matrix(do.generate.random.marginal.samples(dist, n=n.sim), nrow=n.sim)

    if (log.p)
        p = exp(p)
    if (!lower.tail)
        p = 1-p

    rv = apply(rands, 2, quantile, probs=p)

    rv
})


#'@title Internal function for subclasses to get medians of the variables in a distribution
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{get.medians}}
#'@inherit get.medians
#'
#'@details The wrapper \code{\link{get.medians}} takes care of formatting the return value and generating warnings for any improper variables
#'
#'@keywords internal
#'@export
setGeneric('do.get.medians',
           def=function(dist, n.sim){standardGeneric('do.get.medians')})
setMethod('do.get.medians',
          signature(dist='Distribution'),
def=function(dist, n.sim=1000)
{
    do.get.quantiles(dist=dist, p=0.5, n.sim=n.sim)
})


#'@title Internal function for subclasses to get equal-tailed intervals (confidence/credible intervals) for the variables in a distribution
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{get.intervals}} or \code{\link{get.equal.tailed.intervals}}
#'
#'@inherit get.equal.tailed.intervals
#'
#'@details The default implementation calls \code{\link{get.quantiles}} on the alpha/2 and 1-alpha/2 corresponding to each covarage.
#'The wrapper \code{\link{get.equal.tailed.intervals}} takes care of formatting the return value and generating warnings for any improper variables
#'
#'@keywords internal
#'@export
setGeneric('do.get.equal.tailed.intervals',
           def=function(dist, coverage=0.95, n.sim=1000){standardGeneric("do.get.equal.tailed.intervals")})
setMethod('do.get.equal.tailed.intervals',
          signature(dist='Distribution'),
def=function(dist, coverage=0.95, n.sim=1000)
{
    alphas = (1-coverage)
    sapply(alphas, function(alpha){
        get.quantiles(dist=dist, p=c(alpha/2, 1-alpha/2), n.sim=n.sim)
    })
})

#'@title Internal function for subclasses to get equal-tailed intervals (confidence/credible intervals) for the variables in a distribution
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{get.intervals}} or \code{\link{get.highest.density.intervals}}
#'
#'@inherit get.highest.density.intervals
#'
#'@details The default implementation takes random samples and calculates an empiric highest-density interval.
#'The wrapper \code{\link{get.highest.density.intervals}} takes care of formatting the return value and generating warnings for any improper variables
#'
#'@keywords internal
#'@export
setGeneric('do.get.highest.density.intervals',
           def=function(dist, coverage=0.95, n.sim=1000){standardGeneric("do.get.highest.density.intervals")})
setMethod('do.get.highest.density.intervals',
          signature(dist='Distribution'),
def=function(dist, coverage=0.95, n.sim=1000)
{
    rands = matrix(do.generate.random.marginal.samples(dist, n=n.sim), nrow=n.sim)
    do.get.highest.density.interval.from.samples(samples=rands, coverage=coverage)
})

##------------------------------------------------##
##--          MEANS, VARIANCE, and CDF          --##
##-- (May be overwritten at the subclass level) --##
##------------------------------------------------##

#'@title Get means of the variables in a distribution
#'
#'@inheritParams get.quantiles
#'
#'@return A numeric vector with a mean for each variable in the distribution
#'@details Note to developers: this function should not be overridden. Override \code{\link{do.get.means}} instead
#'
#'@export
setGeneric('get.means',
           def=function(dist, n.sim=1000){standardGeneric("get.means")})
setMethod('get.means',
def=function(dist, n.sim=1000)
{
    print.improper.warnings(dist, description="calculate means for")

    rv = do.get.means(dist=dist, n.sim=n.sim)
    rv = as.numeric(rv)
    names(rv) = dist@var.names

    rv
})

#'@title Get the standard deviations of the variables in a distribution
#'
#'@inheritParams get.quantiles
#'
#'@return A numeric vector with a standard deviation for each variable in the distribution
#'@details Note to developers: this function should not be overridden. Override \code{\link{do.get.sds}} instead
#'
#'@export
setGeneric('get.sds',
           def=function(dist, n.sim=1000){standardGeneric("get.sds")})
setMethod('get.sds',
def=function(dist, n.sim=1000)
{
    print.improper.warnings(dist, description="calculate.sds.for")

    rv = do.get.sds(dist=dist, n.sim=n.sim)
    rv = as.numeric(rv)
    names(rv) = dist@var.names

    rv
})

#'@title Get the variance-covariance matrix for variables in a distribution
#'
#'@inheritParams get.quantiles
#'
#'@return A square matrix of dimension number variables in distribution representing the variance-covariance matrix for the distribution
#'@details Note to developers: this function should not be overridden. Override \code{\link{do.get.covariance.matrix}} instead
#'
#'@export
setGeneric('get.covariance.matrix',
           def=function(dist, n.sim=1000){standardGeneric("get.covariance.matrix")})
setMethod('get.covariance.matrix',
          signature(dist='Distribution'),
def=function(dist, n.sim=1000)
{
    print.improper.warnings(dist, description="calculate covariances for")

    rv = do.get.covariance.matrix(dist=dist, n.sim=n.sim)
    dim(rv) = c(dist@n.var, dist@n.var)

    dimnames(rv) = list(dist@var.names, dist@var.names)

    rv
})

#'@title Calculate the cumulative distribution function for a distribution
#'
#'@inheritParams calculate.density
#'@param x The values of the parameters at which to calculate the CDF. May be either a vector, if the CDF at a single point is desired, or a matrix where each row represents one point and each column represents a variable in the distribution
#'@param log A logical indicating whether to return the CDF on the log scale
#'
#'@details Note to developers: this function should not be overridden. Override \code{\link{do.calculate.cdf}} instead
#'@return If x is a vector, a scalar density. If x is a matrix, a numeric vector of densities, one for each row in x
#'
#'@export
setGeneric('calculate.cdf',
           def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000*dist@n.var){standardGeneric("calculate.cdf")})
setMethod('calculate.cdf',
          signature(dist='Distribution'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000*dist@n.var)
{
    print.improper.warnings(dist, description="calculate cumulative distributions for")

    q = match.variables(dist, q)

    rv = do.calculate.cdf(dist=dist, q=q, lower.tail=lower.tail, log.p=log.p, n.sim=n.sim)

    names(rv) = NULL
    rv
})

##-- DEFAULT IMPLEMENTATIONS --##

#'@title Internal function for subclasses to get means of the variables in a distribution
#'
#'@inherit get.means
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{get.means}}
#'
#'@details The default implementation generates random samples and takes the sample mean for each variable.
#'The wrapper \code{\link{get.means}} takes care of formatting the return value and generating warnings for any improper variables
#'
#'@keywords internal
#'@export
setGeneric('do.get.means',
           def=function(dist, n.sim=1000){standardGeneric("do.get.means")})
setMethod('do.get.means',
          signature(dist='Distribution'),
def=function(dist, n.sim=1000)
{
    rands = matrix(do.generate.random.marginal.samples(dist, n=n.sim), nrow=n.sim)
    colMeans(rands)
})

#'@title An internal function for subclasses to get the standard deviations of the variables in a distribution
#'
#'@inheritParams get.sds
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{get.sds}}
#'
#'@details The default implementation calls \code{\link{get.covariance.matrix}} and takes the square-root of the diagonal
#'The wrapper \code{\link{get.sds}} takes care of formatting the return value and generating warnings for any improper variables
#'
#'@keywords internal
#'@export
setGeneric('do.get.sds',
           def=function(dist, n.sim=1000){standardGeneric("do.get.sds")})
setMethod('do.get.sds',
signature(dist='Distribution'),
def=function(dist, n.sim=1000)
{
    cov.mat = get.covariance.matrix(dist=dist, n.sim=n.sim)
    sqrt(diag(cov.mat))
})

#'@title Internal function for subclasses to get the variance-covariance matrix for variables in a distribution
#'
#'@inheritParams get.covariance.matrix
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{get.covariance.matrix}}
#'
#'@details The default implementation generates random samples and takes the sample covariance.
#'The wrapper \code{\link{get.covariance.matrix}} takes care of formatting the return value and generating warnings for any improper variables
#'
#'@keywords internal
#'@export
setGeneric('do.get.covariance.matrix',
           def=function(dist, n.sim){standardGeneric('do.get.covariance.matrix')})
setMethod('do.get.covariance.matrix',
          signature(dist='Distribution'),
def=function(dist, n.sim=1000)
{
    rands = matrix(do.generate.random.samples(dist, n=n.sim), nrow=n.sim)
    if (dist@n.var==1)
        var(rands)
    else
        cov(rands)
})

#'@title Internal function for subclasses to calculate the cumulative density function for a distribution
#'
#'@inheritParams calculate.cdf
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{calculate.cdf}}
#'
#'@details The default implementation generates random samples and takes the mean number of samples where all variables are <= the given values of q
#'The wrapper \code{\link{calculate.cdf}} takes care of matching variables (so the 'q' passed to this function is clean) and formatting the return value and generating warnings for any improper variables
#'
#'@keywords internal
#'@export
setGeneric('do.calculate.cdf',
           def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000*dist@n.var){standardGeneric("do.calculate.cdf")})
setMethod('do.calculate.cdf',
          signature(dist='Distribution'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000*sqrt(dist@n.var))
{
    rands = generate.random.samples(dist, n=n.sim)
    if (is.null(dim(rands)))
        dim(rands) = c(length(rands),1)
    rv = apply(x, 1, function(one.q){
        mean(apply(rands, 1, function(one.rand){
            all(one.rand <= one.q)
        }))
    })

    if (!lower.tail)
        rv = 1-rv
    if (log.p)
        rv = log(rv)

    rv
})

##------------------------------------------------##
##--        MARGINAL DENSITIIES AND CDFs        --##
##-- (May be overwritten at the subclass level) --##
##------------------------------------------------##

#'@title Calculate a marginal density for each variable in a distribution
#'
#'@inheritParams calculate.density
#'
#'@return If x is a vector or a matrix with only one row, a vector of densities with one value for each variable in the distribution. If x is a matrix, a matrix of densities, with one for each row in x and one column for each variable in the distribution
#'
#'@details Note to developers: this function should not be overridden. Override \code{\link{do.calculate.marginal.densities}} instead
#'@export
setGeneric('calculate.marginal.densities',
           def=function(dist, x, log=F, n.sim=1000){standardGeneric('calculate.marginal.densities')})
setMethod('calculate.marginal.densities',
          signature(dist='Distribution'),
def=function(dist, x, log=F, n.sim=1000)
{
    x = match.variables(dist, x)

    rv = do.calculate.marginal.densities(dist=dist, x=x, log=log, n.sim=n.sim)
    dim(rv) = c(dim(x)[1], variable=dist@n.var)
    dimnames(rv) = list(NULL, variable=dist@var.names)

    if (dim(x)[1]==1)
        rv[1,]
    else
        rv
})

#'@title Calculate a marginal cumulative distribution function for each variable in a distribution
#'
#'@inheritParams calculate.cdf
#'
#'@return If q is a vector, a vector of probabilities with one value for each variable in the distribution. If q is a matrix, a matrix of probabilites, with one for each row in x and one column for each variable in the distribution
#'@details Note to developers: this function should not be overridden. Override \code{\link{do.calculate.marginal.cdfs}} instead
#'
#'@export
setGeneric('calculate.marginal.cdfs',
           def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000){standardGeneric('calculate.marginal.cdfs')})
setMethod('calculate.marginal.cdfs',
          signature(dist='Distribution'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000)
{
    print.improper.warnings(dist, description="calculate cumulative distributions for")

    q = match.variables(dist, q)

    rv = do.calculate.marginal.cdfs(dist=dist, q=q, lower.tail=lower.tail, log.p=log.p, n.sim=n.sim)
    dim(rv) = c(dim(q)[1], variable=dist@n.var)
    dimnames(rv) = list(NULL, variable=dist@var.names)

    if (dim(q)[1]==1)
        rv[1,]
    else
        rv
})

#'@title Generate marginal random samples for each variable in a distribution
#'
#'@description Generates random samples for each variable in the distribution, but does NOT preserve the correlations across variables within samples. For some distribution implementations, this can be substantially faster than generating true, multivariate random samples, so this method is preferred if random samples are only needed for one variable at a time
#'
#'@inherit generate.random.samples
#'
#'@details Note to developers: this function should not be overridden. Override \code{\link{do.generate.random.marginal.samples}} instead
#'
#'@export
setGeneric('generate.random.marginal.samples',
           def=function(dist, n){standardGeneric('generate.random.marginal.samples')})
setMethod('generate.random.marginal.samples',
          signature(dist='Distribution'),
def=function(dist, n)
{
  print.improper.warnings(dist, description="generate random samples from")

  rv = do.generate.random.marginal.samples(dist=dist, n=n)

  dim(rv) = c(sample=n, variable=dist@n.var)
  dimnames(rv) = list(sample=NULL, variable=dist@var.names)

  if (n==1)
      rv[1,]
  else if (dist@n.var==1)
      rv[,1]
  else
      rv
})

##-- DEFAULT IMPLEMENTATIONS --##

#'@title An internal function for subclasses to calculate a marginal density for each variable in a distribution
#'
#'@inheritParams calculate.density
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{calculate.marginal.densities}}
#'
#'@details By default, implemented only for distributions with only one variable
#'The wrapper \code{\link{calculate.marginal.densities}} takes care of matching variables (so the 'x' passed to this function is clean) and formatting the return value
#'
#'@keywords internal
#'@export
setGeneric('do.calculate.marginal.densities',
           def=function(dist, x, log=F, n.sim=1000){standardGeneric('do.calculate.marginal.densities')})
setMethod('do.calculate.marginal.densities',
          signature(dist='Distribution'),
def=function(dist, x, log=F, n.sim=1000)
{
    if (dist@n.var==1)
        calculate.density(dist, x=x, log=log, n.sim=n.sim)
    else
        stop(paste0("How to calculate marginal densities has not been defined for this ", class(dist)[1]))
})


#'@title Internal function for subclasses tocalculate a marginal cumulative distribution function for each variable in a distribution
#'
#'@inheritParams calculate.marginal.cdfs
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{calculate.marginal.cdfs}}
#'
#'@details By default, implemented only for distributions with only one variable
#'The wrapper \code{\link{calculate.marginal.cdfs}} takes care of matching variables (so the 'q' passed to this function is clean) and formatting the return value and generating warnings for any improper variables
#'
#'@keywords internal
#'@export
setGeneric('do.calculate.marginal.cdfs',
           def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000){standardGeneric('do.calculate.marginal.cdfs')})
setMethod('do.calculate.marginal.cdfs',
          signature(dist='Distribution'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000)
{
    if (dist@n.var==1)
        calculate.cdf(dist, q=q, lower.tail=lower.tail, log.p=log.p, n.sim=n.sim)
    else
        stop(paste0("How to calculate marginal cdfs has not been defined for this ", class(dist)[1]))
})

#'@title Internal function for subclasses to generate marginal random samples for each variable in a distribution
#'
#'@inheritParams generate.random.marginal.samples
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{generate.random.marginal.samples}}
#'
#'@details The default implementation just calls the standard \link{do.generate.random.samples}
#'The wrapper \code{\link{generate.random.marginal.samples}} takes care of formatting the return value and generating warnings for any improper variables
#'
#'@keywords internal
#'@export
setGeneric('do.generate.random.marginal.samples',
           def=function(dist, n){standardGeneric('do.generate.random.marginal.samples')})
setMethod('do.generate.random.marginal.samples',
signature(dist='Distribution'),
def=function(dist, n)
{
    do.generate.random.samples(dist, n)
})

##------------------------------------------------##
##--          SUBSETTING DISTRIBUTIONS          --##
##-- (May be overwritten at the subclass level) --##
##------------------------------------------------##


#'@title Get a distribution that describes a subset of the variables in a given distribution
#'
#'@inheritParams calculate.density
#'@param vars.to.keep Either a character vector with the names of the variables to keep, an integer/numeric vector with the indices of the variables to keep, or a logical vector of length n.var indicating which variables to keep
#'
#'@return An object of class \code{\link{Distribution-class}}
#'@details Note to developers: this function should not be overridden. Override \code{\link{do.subset.distribution}} instead
#'
#'@export
setGeneric('subset.distribution',
           def=function(dist, vars.to.keep){standardGeneric("subset.distribution")})
setMethod('subset.distribution',
          signature(dist='Distribution'),
def=function(dist, vars.to.keep)
{
    keep.indices = map.keep.indices(dist, vars.to.keep)
    if (length(keep.indices)==0)
        stop("In subsetting a distribution, you must keep at least one variable")
    else if (setequal(keep.indices, 1:dist@n.var))
        dist
    else
        do.subset.distribution(dist=dist, keep.indices=keep.indices)
})


#'@title Internal function for subclasses to get a distribution that describes a subset of the variables in a given distribution
#'
#'@inheritParams subset.distribution
#'@param keep.indices The numeric indices of which variables of the distribution (in the order they appear in var.names) to keep (this will be a subset of 1:dist@n.var)
#'
#'@description This function is provided to be overriden by developers. End users should instead call \code{\link{subset.distribution}}
#'
#'@details Not implemented by default
#'The wrapper \code{\link{subset.distribution}} takes care of matching variables to their numeric indices. This function is only called by the wrapper if a proper subset (ie, not all the variables) is desired
#'
#'@keywords internal
#'@export
setGeneric('do.subset.distribution',
           def=function(dist, keep.indices){standardGeneric('do.subset.distribution')})
setMethod('do.subset.distribution',
          signature(dist='Distribution'),
def=function(dist, keep.indices)
{
    stop(paste0("How to subset the distribution has not been defined for this ", class(dist)[1]))
})


##-----------------------------------------------------------------------##
##-- FUNCTIONS THAT DON'T NEED TO BE OVERWRITTEN AT THE SUBCLASS LEVEL --##
##-----------------------------------------------------------------------##

#'@title Get a function that calculates the density for the distribution
#'
#'@description A convenience function that packages up a call to \code{\link{calculate.density}} on this distribution into its own stand-alone function (which could be passed along as an argument to other functions)
#'
#'@inheritParams calculate.density
#'@param default.log Whether the returned function should calculate the density on the log scale (by default, the returned function can have a log parameter set explicitly)
#'
#'@return A function that takes two arguments. The first, mandatory, is 'x' - the values of the parameters at which to calculate the density. May be either a vector, if the density at a single point is desired, or a matrix where each row represents one point and each column represents a variable in the distribution
#'The second, optional, is 'log' - whether to calculate the density on the log scale. By default, log is set to default.log
#'
#'@export
setGeneric('get.density.function',
           def=function(dist, default.log=T, n.sim=1000){standardGeneric("get.density.function")})
setMethod('get.density.function',
          signature(dist='Distribution'),
def=function(dist, default.log=T, n.sim=1000)
{
    function(x, log=default.log)
    {
        calculate.density(dist, x=x, log=log, n.sim=n.sim)
    }
})

#'@title Get a function that generates random samples form the distribution
#'
#'@description A convenience function that packages up a call to \code{\link{generate.random.samples}} on this distribution into its own stand-alone function (which could be passed along as an argument to other functions)
#'
#'@inheritParams calculate.density
#'
#'@return A function that takes one argument, 'n' - the number of random samples to generate
#'
#'@export
setGeneric('get.random.generator.function',
           def=function(dist){standardGeneric("get.random.generator.function")})
setMethod('get.random.generator.function',
          signature(dist='Distribution'),
def=function(dist)
{
    function(n)
    {
        generate.random.samples(dist, n=n)
    }
})


##-------------------##
##-- OTHER METHODS --##
##-------------------##

setMethod('get.description',
          signature(object='Distribution'),
          function(object){
              if (object@n.var==1)
              {
                  if (object@is.discrete)
                      prefix = "A discrete univariate distribution with support on "
                  else
                      prefix = "A continuous univariate distribution with support on "
                  paste0(prefix, get.description(object@support))
              }
              else
              {
                  if (all(object@is.discrete))
                      paste0("A multivariate, discrete distribution over ", object@n.var, " variables")
                  else if (!any(object@is.discrete))
                      paste0("A multivariate, continuous distribution over ", object@n.var, " variables")
                  else
                      paste0("A multivariate, mixed discrete and continuous distribution over ", object@n.var, " variables")

              }
          })

setMethod('show',
          signature(object='Distribution'),
          function(object){
              print(get.description(object))
          })
