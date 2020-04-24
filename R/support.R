
##------------------------------------##
##-- THE MAIN CLASS AND SUB-CLASSES --##
##------------------------------------##

#'@title An object representing support over real numbers
#'
#'@name Support
#'@rdname Support
#'@aliases Support-class
#'@exportClass Support
#'@export
setClass('Support',
         representation=list(
             n.var='integer',
             is.discrete='logical',
             integers.only='logical',
             is.contiguous='logical'))

setMethod('initialize',
          signature(.Object='Support'),
def=function(.Object,
             n.var=1,
             is.discrete,
             integers.only,
             is.contiguous)
{
    if ((!is(n.var, 'integer') && (!is(n.var, 'numeric') || round(n.var)!=n.var)) ||
        length(n.var) != 1 || is.na(n.var) || n.var < 1)
        stop("'n.var' must be a scalar, non-NA, positive integer")

    if (length(is.discrete) != n.var)
        stop(paste0("'is.discrete' must the same length as n.var (", n.var, ")"))
    if (length(integers.only) != n.var)
        stop(paste0("'integers.only' must the same length as n.var (", n.var, ")"))
    if (length(is.contiguous) != n.var)
        stop(paste0("'is.contiguous' must the same length as n.var (", n.var, ")"))

    if (any(integers.only & !is.discrete))
        stop(paste0("is.discrete must be TRUE if integers.only is TRUE"))

    .Object@n.var = as.integer(n.var)
    .Object@is.discrete = is.discrete
    .Object@integers.only = integers.only
    .Object@is.contiguous = is.contiguous

    .Object
})

#'@title An object representing support over an arbitrary set of real numbers
#'
#'@name Discrete_Set_Support
#'@rdname Discrete_Set_Support
#'@aliases Discrete_Set_Support-class
#'@exportClass Discrete_Set_Support
#'@export
setClass('Discrete_Set_Support',
         contains='Support',
         representation=list(supported.values='numeric'))

setMethod('initialize',
          signature(.Object='Discrete_Set_Support'),
def=function(.Object, supported.values)
{
    .Object = callNextMethod(.Object, n.var=1, is.discrete=T, integers.only=F, is.contiguous=F)

    if (!is(supported.values, 'numeric') && !is(supported.values, 'integer'))
        stop("'supported.values' must be a numeric vector with length >= 1; NA values are not permitted")

    if (length(supported.values)==0 || any(is.na(supported.values)))
        stop("'supported.values' must be a numeric vector with length >= 1; NA values are not permitted")

    .Object@supported.values = supported.values

    .Object
})

#'@title An object representing support over a range of real numbers
#'
#'@name Range_Support
#'@rdname Range_Support
#'@aliases Range_Support-class
#'@exportClass Range_Support
#'@export
setClass('Range_Support',
         contains='Support',
         representation=list(
             lower = 'numeric',
             upper = 'numeric',
             lower.inclusive = 'logical',
             upper.inclusive = 'logical'
         ))
setMethod('initialize',
signature(.Object='Range_Support'),
def=function(.Object, lower, upper, lower.inclusive=integers.only, upper.inclusive=integers.only, integers.only=F)
{

    #Check data types
    if (!is(lower, 'numeric') && !is(lower, 'integer'))
        stop("'lower' and 'upper' must be non-empty, non-NA, numeric vectors of the same length")
    if (!is(upper, 'numeric') && !is(upper, 'integer'))
        stop("'lower' and 'upper' must be non-empty, non-NA, numeric vectors of the same length")

    if (length(lower)==0 || length(upper)==0)
        stop("'lower' and 'upper' must be non-empty, non-NA, numeric vectors of the same length")

    #Make sure lengths of lower and upper are the same

    if (length(lower)==0 || length(upper)==0 || length(upper) != length(lower))
        stop("'lower' and 'upper' must be non-empty, non-NA, numeric vectors of the same length")

    #Make sure inclusive markers are the appropriate length

    lower.inclusive = as.logical(lower.inclusive)
    upper.inclusive = as.logical(upper.inclusive)

    if (any(is.na(lower.inclusive)) || any(is.na(upper.inclusive)))
        stop("'lower.inclusive' and 'upper.inclusive' must be non-NA logical vectors")

    if (length(lower.inclusive)==1)
        lower.inclusive = rep(lower.inclusive, length(lower))
    if (length(upper.inclusive)==1)
        upper.inclusive = rep(upper.inclusive, length(upper))
    if (length(lower.inclusive) != length(lower) || length(upper.inclusive) != length(upper))
        stop("'lower.inclusive' and 'upper.inclusive' must have the same length as 'lower' and 'upper")

    # Make sure uppers are not below lowers
    if (any((!lower.inclusive | !upper.inclusive) & upper <= lower))
        stop("upper must be > lower for exclusive ranges")
    if (any(lower.inclusive & upper.inclusive & upper < lower))
        stop("upper must be >= lower for inclusive ranges")

    # Set it and return

    .Object = callNextMethod(.Object, n.var=1,
                             is.discrete=integers.only, integers.only=integers.only,
                             is.contiguous = length(lower)==1)
    .Object@lower = lower
    .Object@upper = upper
    .Object@lower.inclusive = lower.inclusive
    .Object@upper.inclusive = upper.inclusive

    .Object
})

#'@title An object representing the union of two or more supports over a range of real numbers
#'
#'@name Union_Support
#'@rdname Union_Support
#'@aliases Union_Support-class
#'@exportClass Union_Support
#'@export
setClass('Union_Support',
         contains='Support',
         representation=list(component.supports='list',
                             n.components='integer'))
setMethod('initialize',
          signature(.Object='Union_Support'),
def=function(.Object, ...)
{
    component.supports = list()

    for (elem in list(...))
    {
        if (is(elem, 'Support'))
            component.supports = c(component.supports, list(elem))
        else if (is(elem, 'list'))
        {
            if (!all(sapply(elem, function(sub.elem){is(sub.elem, 'Support')})))
                stop("... must contain only object of class 'Support' or lists which contain only objects of class 'Support'")

            component.supports = c(component.supports, elem)
        }
        else
            stop("... must contain only object of class 'Support' or lists which contain only objects of class 'Support'")
    }

    if (any(sapply(component.supports, function(supp){supp@n.var != 1})))
        stop("Can only create a Union of Supports for univariate supports (n.var=1)")

    .Object = callNextMethod(.Object,
                             is.contiguous=F,
                             is.discrete=all(sapply(component.supports, function(support){support@is.discrete})))
    .Object@component.supports = component.supports
    .Object@n.components = length(component.supports)

    .Object
})

#'@title An object representing support (on the real numbers) over across multiple variables
#'
#'@name Multivariate_Support
#'@rdname Multivariate_Support
#'@aliases Multivariate_Support-class
#'@exportClass Multivariate_Support
#'@export
setClass('Multivariate_Support',
         contains='Support',
         representation=list(component.supports='list'))
setMethod('initialize',
          signature(.Object='Multivariate_Support'),
def=function(.Object, ...)
{
    component.supports = list()
    for (elem in list(...))
    {
        if (is(elem, 'Multivariate_Support'))
            component.supports = c(component.supports, elem@component.supports)
        else if (is(elem, 'Support'))
            component.supports = c(component.supports, list(elem))
        else if (is(elem, 'list'))
        {
            for (sub.elem in elem)
            {
                if (is(sub.elem, 'Multivariate_Support'))
                    component.supports = c(component.supports, sub.elem@component.supports)
                else if (is(sub.elem, 'Support'))
                    component.supports = c(component.supports, list(sub.elem))
                else
                    stop("... must contain only object of class 'Support' or lists which contain only objects of class 'Support'")
            }
        }
        else
            stop("... must contain only object of class 'Support' or lists which contain only objects of class 'Support'")
    }

    is.discrete = unlist(sapply(component.supports, function(one.support){one.support@is.discrete}))
    integers.only = unlist(sapply(component.supports, function(one.support){one.support@integers.only}))
    is.contiguous = unlist(sapply(component.supports, function(one.support){one.support@is.contiguous}))

    .Object = callNextMethod(.Object, n.var=length(component.supports),
                             is.discrete=is.discrete, integers.only=integers.only, is.contiguous=is.contiguous)
    .Object@component.supports = component.supports
    names(.Object@component.supports) = NULL

    .Object

})

#'@title Subset a Support object
#'
#'@param support An object of class \link{Support}
#'@param keep Which dimensions of the support to keep. Either a set of integer indices or a logical vector
#'
#'@return A \link{Support} object
#'
#'@export
setGeneric('subset.support',
           function(support, keep){standardGeneric('subset.support')})
setMethod('subset.support',
          signature(support='Support'),
def=function(support, keep)
{
    if (!((is(keep, 'integer') || is(keep, 'numeric')) && length(setdiff(keep, 1:support@n.var))==0) &&
        !(is(keep, 'logical') && length(keep) == support@n.var))
        stop(paste("In subsetting this Support, keep must be either integer indices from 1 to ",
                    support@n.var, " or a logical vector of length ", support@n.var))

    if (any(is.na(keep)))
        stop("'keep' cannot have NA values")

    if (is(keep, 'logical'))
        keep = (1:support@n.var)[keep]

    if (setequal(keep, 1:support@n.var))
        support
    else
        stop(paste0("Subsetting is not implemented for this ", class(support)))
})

setMethod('subset.support',
          signature(support='Multivariate_Support'),
def=function(support, keep)
{
    if (!((is(keep, 'integer') || is(keep, 'numeric')) && length(setdiff(keep, 1:support@n.var))==0) &&
        !(is(keep, 'logical') && length(keep) == support@n.var))
        stop(paste("In subsetting this Support, keep must be either integer indices from 1 to ",
                   support@n.var, " or a logical vector of length ", support@n.var))

    if (any(is.na(keep)))
        stop("'keep' cannot have NA values")

    if (is(keep, 'logical'))
        keep = (1:support@n.var)[keep]
    else
        keep = unique(keep)

    if (setequal(keep, 1:support@n.var))
        support
    else if (length(keep)==1)
        support@component.supports[[keep]]
    else
        Multivariate.Support(support@component.supports[keep])
})

##------------------------------------------------------##
##-- IMPLEMENTATIONS OF MERGE (AT THE SUBCLASS LEVEL) --##
##------------------------------------------------------##

setGeneric('merge.support',
           def=function(s1, s2){standardGeneric('merge.support')})


##-----------------------------------------------------##
##-- IMPLEMENTATIONS OF SHOW (AT THE SUBCLASS LEVEL) --##
##-----------------------------------------------------##

setMethod('show',
          signature(object='Support'),
          def=function(object){
              print(paste0("Support on ", get.description(object)))
          })

setMethod('get.description',
          signature(object='Discrete_Set_Support'),
def=function(object)
{
    if (length(object@supported.values <= 5))
        paste0("the discrete set <",
               paste0(sort(object@supported.values), collapse=', '),
               ">")
    else
        paste0("a set of ",
               length(object@supported.values),
               " discrete values")
})

setMethod('get.description',
          signature(object='Range_Support'),
def=function(object)
{
    str.ranges = paste0(c('(','[')[1+as.numeric(object@lower.inclusive)],
                        object@lower,
                        ' to ',
                        object@upper,
                        c(')',']')[1+as.numeric(object@upper.inclusive)])

    all.str.ranges = paste0(str.ranges, collapse = ' U ')

    if (object@integers.only)
        paste0("integers in ", all.str.ranges)
    else
        paste0("real numbers in ", all.str.ranges)
})

setMethod('get.description',
          signature(object='Union_Support'),
def=function(object)
{
    paste0("any of ", object@n.components, " Support objects")
})

setMethod('get.description',
          signature(object='Multivariate_Support'),
def=function(object)
{
    if (!any(object@is.discrete))
        descriptor = 'continuous'
    else if (all(object@integers.only))
        descriptor = 'integer'
    else if (all(object@discrete))
        descriptor = 'discrete'
    else
        descriptor = "mixed continuous and discrete"
  paste0(descriptor, " values for ", object@n.var, " variables")
})

##-------------------------------------------------------------##
##-- IMPLEMENTATIONS OF IS.SUPPORTED (AT THE SUBCLASS LEVEL) --##
##-------------------------------------------------------------##

#'@title Test if numeric values fall within a defined support
#'
#'@param support An object of class \link{Support}
#'@param x If the support is univariate, a numeric vector. If the support is multivariate, a matrix with one column for each variable and one row for each point to test or a numeric vector with one value for each variable in the support
#'@param by.variable If this is a Multivariate_Support, if by.variable is TRUE, calculates support for each variable independently
#'
#'@return If the support is univariate or if by.variable==FALSE or if only a vector is passed to x, returns a logical vector, with one value for each value of x, corresponding to whether that value of x falls within the defined support. If the support is multivariate and by.variable==TRUE and a matrix with more than one row was passed to x, returns a matrix with one column for each variable and one row for each point
#'
#'@export
setGeneric('is.supported',
           def=function(support, x, by.variable=F){standardGeneric('is.supported')})

setMethod('is.supported',
          signature(support='Discrete_Set_Support'),
          def=function(support, x, by.variable=F)
          {
              sapply(x, function(one.x){
                  any(one.x == support@supported.values)
              })
          })

setMethod('is.supported',
          signature(support='Range_Support'),
          def=function(support, x, by.variable=F)
          {
              rv = sapply(x, function(one.x){
                  above.lower = one.x > support@lower | (support@lower.inclusive & one.x == support@lower)
                  below.upper = one.x < support@upper | (support@upper.inclusive & one.x == support@upper)
                  any(above.lower & below.upper)
              })

              if (support@integers.only)
                  rv = rv & (round(x)==x)

              rv
          })

setMethod('is.supported',
          signature(support='Union_Support'),
          def=function(support, x, by.variable=F)
          {
              by.component = sapply(support@component.supports, is.supported, x=x)
              apply(by.component, 1, any)
          })

setMethod('is.supported',
signature(support='Multivariate_Support'),
def=function(support, x, by.variable=F)
{
    if (is.null(dim(x)))
    {
        if (length(x) != support@n.var)
            stop(paste0("If x is a vector, then it must have one value for each variable in the support (", support@n.var, ")"))

        x = matrix(x, nrow=1)
    }
    else
    {
        if (ncol(x) != support@n.var)
            stop(paste0("If x is a matrix, then it must have one column for each variable in the support (", support@n.var, ")"))
    }

    by.component = sapply(1:support@n.var, function(i){
        is.supported(support@component.supports[[i]], x=x[,i])
    })

    if (by.variable)
        by.component
    else if (nrow(x)==1)
        all(by.component)
    else
        apply(by.component, 1, all)
})

##-------------------------------------------------------------------##
##-- IMPLEMENTATIONS OF GET.SUPPORT.BOUNDS (AT THE SUBCLASS LEVEL) --##
##-------------------------------------------------------------------##

#'@title Get the maximum and minimm values within a support over real numbers
#'
#'@param support An object of class \link{Support}
#'
#'@return If called on an object of class \link{Support}, returns a numeric vector of length 2, with lower and upper bounds. If called on a list containing only \link{Support} objects, returns a matrix with two rows (lower and upper bounds) and a column for each element in the list
#'
#'@export
setGeneric('get.support.bounds',
           def=function(support){standardGeneric('get.support.bounds')})

setMethod('get.support.bounds',
          signature(support='Discrete_Set_Support'),
def=function(support)
{
    c(lower=min(support@supported.values), upper=max(support@supported.values))
})

setMethod('get.support.bounds',
          signature(support='Range_Support'),
def=function(support)
{
    c(lower=min(support@lower), upper=max(support@upper))
})

setMethod('get.support.bounds',
          signature(support='Union_Support'),
def=function(support)
{
    by.component = sapply(support@component.supports, get.support.bounds)
    c(lower=min(by.component[1,]), upper=max(by.component[2,]))
})

setMethod('get.support.bounds',
          signature(support='list'),
def=function(support)
{
    if (!all(sapply(support, function(one.support){is(one.support, 'Support')})))
        stop("If 'support' is a list, it can only contain objects of class 'Support'")

    rv = sapply(support, get.support.bounds)
    dimnames(rv)[[2]] = names(support)
    rv
})


setMethod('get.support.bounds',
          signature(support='Multivariate_Support'),
def=function(support)
{
    rv = sapply(support@component.supports, get.support.bounds)
    rv
})

##---------------------------##
##-- SPECIFIC CONSTRUCTORS --##
##---------------------------##

#'@title Create an object which represents support over a range of integers
#'
#'@param lower,upper The lower and upper bounds of support. For a single range of support, these should both be of length 1. If support exists over multiple disjoint ranges, then lower and upper should be numeric vectors of the same length, with the lower and upper bounds of each of the disjoint ranges
#'@param lower.inclusive,upper.inclusive Logical vectors indicating whether the lower and upper bounds are inclusive (vs. exclusive). If length(lower)/length(upper) are >1, these can be either logical vectors of the same length, or logical scalars (in which case the value applies to each element of lower and upper)
#'
#'@return An object of class \link{Range_Support}, such that is.supported is true if x falls within any of the ranges lower[i] to upper[i]
#'
#'@export
Continuous.Support <- function(lower=-Inf, upper=Inf, lower.inclusive=F, upper.inclusive=F)
{
    new('Range_Support',
        lower=lower,
        upper=upper,
        lower.inclusive=lower.inclusive,
        upper.inclusive=upper.inclusive,
        integers.only=F)
}

#'@title Create an object which represents support over a range of integers
#'
#'@inherit Continuous.Support
#'
#'@export
Integer.Support <- function(lower, upper, lower.inclusive=T, upper.inclusive=T)
{
    new('Range_Support',
        lower=lower,
        upper=upper,
        lower.inclusive=lower.inclusive,
        upper.inclusive=upper.inclusive,
        integers.only=T)
}

#'@title Create a support over an arbitrary set of values
#'
#'@param supported.values A numeric vector of the allowed values in this support
#'
#'@return An object of class \link{Discrete_Set_Support}
#'
#'@export
Discrete.Set.Support <- function(supported.values)
{
    new('Discrete_Set_Support', supported.values=supported.values)
}

#'@title Join multiple univariate Support objects into one univariate Support
#'
#'@param ... Objects of class \link{Support}, or lists that contain only objects of class \link{Support}
#'
#'@return A \link{Union_Support} object, such that is.supported is true if is.supported is true for any of its component supports
#'
#'@export
Support.Union <- function(...)
{
    new('Union_Support', ...)
}

#'@title Join multiple Support objects into a multivariate Support
#'
#'@param ... Objects of class \link{Support}, or lists that contain only objects of class \link{Support}
#'
#'@return A \link{Multivariate_Support} object
#'
#'@export
Multivariate.Support <- function(...)
{
    args = list(...)
    if (length(args)==1)
    {
        if (is(args[[1]], 'Support'))
            return (args[[1]])
        else if (is(args[[1]], 'list') && length(args[[1]])==1 && is(args[[1]][[1]], 'Support'))
            return (args[[1]][[1]])
    }

    new('Multivariate_Support', ...)
}
