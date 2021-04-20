#'@title The transformation class
#'
#'@description A class that represents a functional transformation, along with how to reverse the transformation and take its derivative
#'
#'@slot transform.function,reverse.transform.function,transformation.derivative,log.abs.transformation.derivative Functions that compute the transformation, inverse transformation, derivative of the transformation function, and log of the absolute value of the derivative of the transformation function, directly. All should take one argument, 'x', which can be either a numeric vector or numeric scalar. log.abs.transformation.derivative may be supplied as NULL; the other three functions must be defined
#'@slot name A character value with a descriptive name of the transformation
#'
#'@name transformation
#'@rdname transformation
#'@aliases transformation-class
#'@exportClass transformation
#'@export
setClass('transformation',
         representation=list(transform='function',
                             reverse.transform='function',
                             transformation.derivative='function',
                             log.abs.transformation.derivative='function',
                             name='character'))

#'@title Create an object representing a transformation
#'
#'@param transform.function,reverse.transform.function,transformation.derivative,log.transformation.derivative Functions that compute the transformation, inverse transformation, derivative of the transformation function, and log of the derivative of the transformation function, directly. All should take one argument, 'x', which can be either a numeric vector or numeric scalar. log.transformation.derivative may be supplied as NULL; the other three functions must be defined
#'@param name A character value with a descriptive name of the transformation
#'
#'@export
create.transformation <- function(transform.function,
                                  reverse.transform.function,
                                  transformation.derivative,
                                  log.abs.transformation.derivative=NULL,
                                  name)
{
    if (is.null(log.abs.transformation.derivative))
        log.abs.transformation.derivative <- function(x){log(abs(transformation.derivative(x)))}

    new('transformation',
        transform=transform.function,
        reverse.transform=reverse.transform.function,
        transformation.derivative=transformation.derivative,
        log.abs.transformation.derivative=log.abs.transformation.derivative,
        name=name[1])
}

#'@title Get Pre-Defined Transformation Objects by Name
#'
#'@description Gets a pre-defined transformation object. The available transformations are 'log', 'logit', and 'identity'
#'
#'@return If transformation.names is a single character value, returns an object of class transformation. Otherwise, returns a list of transformation objects
#'@export
get.defined.transformation <- function(transformation.names,
                                       throw.error.if.no.match=F)
{
    if (!is(transformation.names, 'character'))
        stop("transformation.names must be a character vector")

    rv = DEFINED.TRANSFORMATIONS[transformation.names]
    if (throw.error.if.no.match && any(sapply(rv, is.null)))
    {
        unmatched.names = transformation.names[sapply(rv, is.null)]
        stop(paste0("The predefined transformations are: ",
                    paste0("'", names(DEFINED.TRANSFORMATIONS), "'", collapse=', '),
                    ". The following are not predefined transformations: ",
                    paste0("'", unmatched.names, "'", collapse=', ')))
    }

    if (length(rv)==1)
        rv[[1]]
    else
        rv
}

setMethod('show',
          signature(object='transformation'),
          def=function(object){
              cat(paste0(object@name, " transformation object"))
          })

DEFINED.TRANSFORMATIONS = list(log=create.transformation(transform.function = log,
                                                         reverse.transform.function = exp,
                                                         transformation.derivative = function(x){1/x},
                                                         log.abs.transformation.derivative = function(x){-log(x)},
                                                         name='log'),
                               logit=create.transformation(transform.function = function(x){log(x) - log(1-x)},
                                                           reverse.transform.function = function(x){1/(1+exp(-x))},
                                                           transformation.derivative = function(x){1/x/(1-x)},
                                                           log.abs.transformation.derivative = function(x){-log(x)-log(1-x)},
                                                           name='logit'),
                               identity=create.transformation(transform.function = function(x){x},
                                                              reverse.transform.function = function(x){x},
                                                              transformation.derivative = function(x){rep(1, length(x))},
                                                              log.abs.transformation.derivative = function(x){rep(0, length(x))},
                                                              name='identity'),
                               pnorm=create.transformation(transform.function = pnorm,
                                                           reverse.transform.function = qnorm,
                                                           transformation.derivative = dnorm,
                                                           log.abs.transformation.derivative = function(x){dnorm(x, log=T)},
                                                           name='pnorm'),
                               qnorm=create.transformation(transform.function = qnorm,
                                                           reverse.transform.function = pnorm,
                                                           transformation.derivative = function(x){
                                                               inv.rv = dnorm(qnorm(x))
                                                               rv = 1/inv.rv
                                                               rv[inv.rv==0] = 0
                                                               rv
                                                           },
                                                           log.abs.transformation.derivative = function(x){
                                                               inv.rv = dnorm(qnorm(x), log=T)
                                                               rv = -inv.rv
                                                               rv[inv.rv==-Inf] = -Inf
                                                               rv
                                                           },
                                                           name='qnorm')
                               )



match.transformations.to.variables <- function(transformations, var.names, n.var=length(var.names))
{
    parsed.transformations = lapply(1:n.var, function(name){get.defined.transformation('identity')})
    names(parsed.transformations) = var.names
    transformations = as.list(transformations)

    if (is.null(transformations) || length(transformations)==0)
    {}
    else if (is.null(names(transformations)) || is.null(var.names))
    {
        if (length(transformations)==1)
            parsed.transformations = lapply(1:n.var, function(name){transformations[[1]]})
        else if (length(transformations) == n.var)
        {
            parsed.transformations = transformations
            names(parsed.transformations) = var.names
        }
        else
            stop(paste0("If names are not specified for the 'transformations' parameter or variables are not named, ",
                        "'transformations' must have the same length as the number of variables (", n.var, ")"))
    }
    else
    {
        unused.transformation.names = setdiff(names(transformations), var.names)
        if (length(unused.transformation.names) > 0)
            warning(paste0("The following transformation",
                           ifelse(length(unused.transformation.names>1, "s were", 'was')),
                           " given but ",
                           ifelse(length(unused.transformation.names>1, "are not named variables: ", 'is not a named variable: ')),
                           paste0("'", unused.transformation.names, "'", collapse=", ")
            ))

        parsed.transformations[names(transformations)] = transformations[names(transformations)]
    }

    transformations = lapply(parsed.transformations, function(one.transformation){
        if (is.null(one.transformation))
            distributions::get.defined.transformation('identity')
        else if (is(one.transformation, 'transformation'))
            one.transformation
        else if (is(one.transformation, 'character'))
        {
            if (length(one.transformation)!=1)
                stop("Only single names (ie not vectors) are permitted as elements of 'transformations")

            get.defined.transformation(one.transformation, throw.error.if.no.match = T)
        }
        else
            stop("transformations must be either NULL, an instance of the 'transformation' class, or a character name of a defined transformation")
    })

    transformations
}

#x is a matrix with one column per variable
do.transform.variables <- function(transformations, x)
{
    transformed.x = sapply(1:length(transformations), function(i){
        transformations[[i]]@transform(x[,i])
    })

    dim(transformed.x) = dim(x)
    dimnames(transformed.x) = dimnames(x)

    transformed.x
}

do.reverse.transform.variables <- function(transformations, x)
{
    transformed.x = sapply(1:length(transformations), function(i){
        transformations[[i]]@reverse.transform(x[,i])
    })

    dim(transformed.x) = dim(x)
    dimnames(transformed.x) = dimnames(x)

    transformed.x
}

do.transformation.derivative.for.variables <- function(transformations, x)
{
    derivative = sapply(1:length(transformations), function(i){
        transformations[[i]]@transformation.derivative(x[,i])
    })

    dim(derivative) = dim(x)
    dimnames(derivative) = dimnames(x)

    derivative
}

do.log.abs.transformation.derivative.for.variables <- function(transformations, x)
{
    log.derivative = sapply(1:length(transformations), function(i){
        transformations[[i]]@log.abs.transformation.derivative(x[,i])
    })

    dim(log.derivative) = dim(x)
    dimnames(log.derivative) = dimnames(x)

    log.derivative
}
