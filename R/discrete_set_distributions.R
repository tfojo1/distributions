
##--------------------------------------##
##-- CLASS DEFINITION and CONSTRUCTOR --##
##--------------------------------------##

#'@title The Discrete_Set_Distribution class
#'
#'@description A class that represents a discrete probability distribution over an arbitrary set of numbers
#'
#'@slot values The numeric set of values
#'@slot weights The probability of each value in the set. Sum to 1
#'
#'@seealso The \link{Constant.Distribution}
#'
#'@name Discrete_Set_Distribution
#'@rdname Discrete_Set_Distribution
#'@aliases Discrete_Set_Distribution-class
#'@exportClass Discrete_Set_Distribution
#'@export
setClass('Discrete_Set_Distribution',
         contains='Distribution',
         representation=list(values='numeric',
                             weights='numeric'))

setMethod('initialize',
          signature(.Object='Discrete_Set_Distribution'),
          def=function(.Object,
                       values,
                       weights=rep(1,length(values)),
                       var.name=NULL)
{
    if (!is.numeric(values) || length(values) == 0  || any(is.na(values)))
        stop("values must be a non-NA, numeric vector with at least one value")

    if (length(unique(values)) != length(values))
        stop("numbers can only appear in 'values' once")

    if (!is.numeric(weights) || length(weights)!=length(values)  || any(is.na(weights)))
        stop("weights must be a non-NA, numeric vector with the same length as values")

    o = order(vlaues)

    .Object = callNextMethod(.Object,
                           var.names=var.name,
                           support=Discrete.Set.Support(supported.values=values[o]),
                           is.improper=F)

    .Object@values = values[o]
    .Object@weights = weights[o] / sum(weights)

    #-- Return --#
    .Object
})

#'@title Create a distribution representing a Discrete Set
#'
#'@description Creates an object representing a discrete probability distribution over an arbitrary set of numbers
#'
#'@param values The numeric set of values
#'@param weights The weight given to the corresponding value. Need not sum to 1
#'@param var.name The name of the single variable in the distribution, or NULL if no name is specified
#'
#'@family Discrete Distribution Constructors
#'@family Distribution Constructors
#'
#'@seealso \link{Constant.Distribution}
#'
#'@export
Discrete.Set.Distribution <- function(values,
                                      weights,
                                      var.name=NULL)
{
    new('Discrete_Set_Distribution', values=values, weights=weights, var.name=var.name)
}


#'@title The Constant_Distribution class
#'
#'@description A class that represents a probility distribution with support for exactly one value
#'
#'@slot value The single value supported by this distribution
#'
#'@seealso The \link{Discrete.Set.Distribution}
#'
#'@name Constant_Distribution
#'@rdname Constant_Distribution
#'@aliases Constant_Distribution-class
#'@exportClass Constant_Distribution
#'@export
setClass('Constant_Distribution',
         contains='Distribution',
         representation=list(value='numeric'))

setMethod('initialize',
          signature(.Object='Constant_Distribution'),
          def=function(.Object,
                       value,
                       var.name=NULL)
{

    if (!is.numeric(value) || length(value) != 1 || is.na(value))
        stop("value must be a single, non-NA, numeric value")

    .Object = callNextMethod(.Object,
                           var.names=var.name,
                           n.var=1,
                           support=Discrete.Set.Support(supported.values=value),
                           is.improper=F)
    .Object@value = value


    #-- Return --#
    .Object
})

#'@title Create a distribution representing a Constant Distribution
#'
#'@description Creates an object representing a degenerate distribution that only supports a single, discrete value
#'
#'@param value The single value supported by the distribution
#'@param var.name The name of the single variable in the distribution, or NULL if no name is specified
#'
#'@family Discrete Distribution Constructors
#'@family Distribution Constructors
#'
#'@seealso \link{Discrete.Set.Distribution}
#'
#'@export
Constant.Distribution <- function(value, var.name=NULL)
{
    new('Constant_Distribution', value=value, var.name=var.name)
}


##----------------------------##
##-- METHODS IMPLEMENTATION --##
##----------------------------##


setMethod('do.calculate.density',
          signature(dist='Discrete_Set_Distribution'),
def=function(dist, x, log=F, n.sim=1000)
{
    sapply(x, function(one.x){
        mask = one.x == dist@values
        if (any(mask))
            dist@weights[mask]
        else
            0
    })
})

setMethod('do.calculate.density',
          signature(dist='Constant_Distribution'),
def=function(dist, x, log=F, n.sim=1000)
{
    if (log)
    {
        rv = rep(0, length(x))
        rv[x!=dist@value] = -Inf
    }
    else
    {
        rv = rep(1, length(x))
        rv[x!=dist@value] = 0
    }

    rv
})


setMethod('do.calculate.cdf',
          signature(dist='Discrete_Set_Distribution'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000)
{
    if (lower.tail)
        rv = sapply(q, function(one.q){
            sum(weights[dist@values <= one.q])
        })
    else
        rv = sapply(q, function(one.q){
            sum(weights[dist@values >= one.q])
        })

    if (log.p)
        log(rv)
    else
        rv
})

setMethod('do.calculate.cdf',
          signature(dist='Constant_Distribution'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000)
{
    if (log.p)
        sapply(q, function(one.q){
            if (one.q==dist@value)
                0
            else
                -Inf
        })
    else
        sapply(q, function(one.q){
            if (one.q==dist@value)
                1
            else
                0
        })
})


setMethod('do.get.quantiles',
          signature(dist='Discrete_Set_Distribution'),
def=function(dist, p, lower.tail=T, log.p=F, n.sim=1000)
{
    if (log.p)
        p = exp(p)

    if (lower.tail)
    {
        cum.prob = cumsum(dist@weights)
        sapply(p, function(one.p){
            mask = cum.prob >= one.p
            dist@values[mask]
        })
    }
    else
    {
        cum.prob = cumsum(rev(adt@weights))
        rev.values = rev(dist@values)
        sapply(p, function(one.p){
            mask = cum.prob >= one.p
            rev.values[mask]
        })
    }
})

setMethod('do.get.quantiles',
          signature(dist='Constant_Distribution'),
def=function(dist, p, lower.tail=T, log.p=F, n.sim=1000)
{
    rep(dist@value, length(p))

})

setMethod('do.generate.random.samples',
          signature(dist='Discrete_Set_Distribution'),
def=function(dist, n)
{
    sample(dist@values, size=n, replace=T, prob = dist@weights)
})


setMethod('do.generate.random.samples',
          signature(dist='Constant_Distribution'),
          def=function(dist, n)
{
    rep(dist@value, n)
})

setMethod('do.get.covariance.matrix',
          signature(dist='Discrete_Set_Distribution'),
def=function(dist, n.sim=1000)
{
    sum(dist@values^2 * dist@weights) - sum(dist@values * dist@weights)^2
})

setMethod('do.get.covariance.matrix',
          signature(dist='Constant_Distribution'),
def=function(dist, n.sim=1000)
{
    0
})

setMethod('do.get.means',
          signature(dist='Discrete_Set_Distribution'),
def=function(dist, n.sim=200)
{
    sum(dist@values * dist@weights)
})

setMethod('do.get.means',
          signature(dist='Constant_Distribution'),
def=function(dist, n.sim=200)
{
    dist@value
})

setMethod('get.description',
          signature(object='Discrete_Set_Distribution'),
def=function(object){
    rv = paste0("A discrete distribuution over ",
                length(dist@values), " value",
                ifelse(length(dist@values)==1, '', 's'))
    if (length(object@values)<=5)
        rv = paste0(rv, ": <", paste0(object@values, collapse=', '), ">")

    rv
})

setMethod('get.description',
          signature(object='Constant_Distribution'),
def=function(object){
    paste0("A constant distribution with value ", object@value)
})

