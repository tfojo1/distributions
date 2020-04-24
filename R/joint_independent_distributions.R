
##--------------------------------------##
##-- CLASS DEFINITION and CONSTRUCTOR --##
##--------------------------------------##

#'@title The Joint_Independent_Distributions class
#'
#'@description A class that joins multiple independent sub-distributions into a joint distribution

#'@slot subdistributions A list of the Distribution objects representing the component distributions of this joint distribution
#'@slot indices.for.subdistributions A list whose elements are integer vectors, such that indices.for.subdistributions[[i]] is the indices into the overall distribution of the variables represented by subdistributions[[i]]
#'@slot n.subdistributions The number of component distributions
#'
#'@seealso \link{join.distributions}
#'
#'@name Joint_Independent_Distributions
#'@rdname Joint_Independent_Distributions
#'@aliases Joint_Independent_Distributions-class
#'@exportClass Joint_Independent_Distributions
#'@export
setClass('Joint_Independent_Distributions',
         contains='Distribution',
         representation=list(subdistributions='list',
                             indices.for.subdistributions='list',
                             n.subdistributions='integer'))

#'@title Join Multiple Distributions
#'
#'@description Joins multiple Distribution objects into a single Distribution object, on the premise that the variables in each sub-distribution are independent of the variables in any other sub-distribution
#'
#'@param ... One or more Distribution objects or lists of Distribution objects to be joined. If var.names are set on the component distributions, those names will be preserved. If not, and if the elements of ... are named (or elements in lists contained in ... are named), those names will be used as the var.names for the corresponding sub-distributions
#'
#'@return An object of class Joint_Independent_Distributions
#'
#'@family Distribution Constructors
#'
#'@export
join.distributions <- function(..., var.names=NULL)
{
    new('Joint_Independent_Distributions', ..., var.names=var.names)
}

setMethod('initialize',
          signature(.Object='Joint_Independent_Distributions'),
def=function(.Object, ..., var.names=NULL)
{
    components = list()

    arguments = list(...)
    for (i in 1:length(arguments))
    {
        if (is(arguments[[i]], 'Distribution'))
            components = c(components, arguments[i]) #passing list subset keeps the name
        else if (is(arguments[[i]], 'list'))
        {
            if (!all(sapply(arguments[[i]], function(sub.elem){is(sub.elem, 'Distribution')})))
                stop("The arguments in ... must all be either 'Distribution' objects or lists that contain only 'Distribution' objects")

            components = c(components, arguments[[i]])
        }
        else
            stop("The arguments in ... must all be either 'Distribution' objects or lists that contain only 'Distribution' objects")
    }

    # Arrange components

    listed.var.names = character()
    supports = list()
    indices.for.subdistributions = list()
    max.index.so.far=0

    for (i in 1:length(components))
    {
        dist = components[[i]]

        #get the name
        component.name = NULL
        if (!is.null(names(components)))
            component.name = names(components)[i]

        if (!is.null(component.name) && (is.na(component.name) || component.name==''))
            component.name = NULL

        # Pull var names
        if (!is.null(dist@var.names) &&
            (is.null(component.name) || dist@n.var>1))
            var.names.for.dist = dist@var.names
        else
        {
            if (is.null(component.name))
                var.names.for.dist = NULL
            else if (dist@n.var==1)
                var.names.for.dist = component.name
            else
                var.names.for.dist = paste0(component.name, '.', 1:dist@n.var)
        }
        listed.var.names = c(listed.var.names, var.names.for.dist)

        supports = c(supports, list(dist@support))

        # Set up indices
        dist.indices = max.index.so.far + 1:dist@n.var
        max.index.so.far = max.index.so.far + dist@n.var
        indices.for.subdistributions = c(indices.for.subdistributions, list(dist.indices))
    }

    # Check var.names if given

    n.var = sum(sapply(components, function(dist){dist@n.var}))
    if (!is.null(var.names) || length(var.names)>0)
    {
        if (!is(var.names, 'character') || length(var.names) != n.var)
            stop(paste0("var.names must be a character vector with one element for each variable"))
    }
    else
        var.names = listed.var.names


    .Object = callNextMethod(.Object,
                             var.names=var.names,
                             n.var=n.var,
                             support=Multivariate.Support(supports),
                             is.improper=unlist(sapply(components, function(dist){dist@is.improper})))

    .Object@subdistributions=components
    .Object@indices.for.subdistributions=indices.for.subdistributions
    .Object@n.subdistributions=length(components)

    .Object
})


##----------------------------##
##-- METHODS IMPLEMENTATION --##
##----------------------------##


setMethod('do.calculate.density',
          signature(dist='Joint_Independent_Distributions'),
def=function(dist, x, log=F, n.sim=1000)
{
    rv.components = sapply(1:dist@n.subdistributions, function(i){
        do.calculate.density(dist@subdistributions[[i]],
                          x=matrix(x[,dist@indices.for.subdistributions[[i]]], ncol=length(dist@indices.for.subdistributions[[i]])),
                          log=log)
    })
    dim(rv.components) = c(dim(x)[1], dist@n.subdistributions)

    if (log)
        rowSums(rv.components)
    else
        apply(rv.components, 1, prod)
})


setMethod('do.calculate.cdf',
          signature(dist='Joint_Independent_Distributions'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000)
{
    rv.components = sapply(1:dist@n.subdistributions, function(i){
        do.calculate.cdf(dist@subdistributions[[i]],
                         q=matrix(q[,dist@indices.for.subdistributions[[i]]], ncol=length(dist@indices.for.subdistributions[[i]])),
                         lower.tail=lower.tail, log.p=log.p, n.sim=n.sim)
    })
    dim(rv.components) = c(dim(q)[1], dist@n.subdistributions)

    if (log.p)
        rowSums(rv.components)
    else
        apply(rv.components, 1, prod)
})


setMethod('do.get.quantiles',
          signature(dist='Joint_Independent_Distributions'),
def=function(dist, p, lower.tail=T, log.p=F, n.sim=1000)
{
    rv = matrix(0, nrow=length(p), ncol=dist@n.var)
    for (i in 1:dist@n.subdistributions)
        rv[,dist@indices.for.subdistributions[[i]]] = do.get.quantiles(dist@subdistributions[[i]],
                                                                    p=p, lower.tail=lower.tail, log.p=log.p, n.sim=n.sim)

    rv
})

setMethod('do.generate.random.samples',
signature(dist='Joint_Independent_Distributions'),
def=function(dist, n)
{
    rv = matrix(0, nrow=n, ncol=dist@n.var)
    for (i in 1:dist@n.subdistributions)
        rv[,dist@indices.for.subdistributions[[i]]] = do.generate.random.samples(dist@subdistributions[[i]], n=n)

    rv
})


setMethod('do.get.covariance.matrix',
signature(dist='Joint_Independent_Distributions'),
def=function(dist, n.sim=1000)
{
    rv = matrix(0, dist@n.var, dist@n.var)
    for (i in 1:dist@n.subdistributions)
        rv[dist@indices.for.subdistributions[[i]],dist@indices.for.subdistributions[[i]]] =
            do.get.covariance.matrix(dist@subdistributions[[i]], n.sim=n.sim)

    rv
})

setMethod('do.get.means',
          signature(dist='Joint_Independent_Distributions'),
def=function(dist, n.sim=1000)
{
    unlist(sapply(dist@subdistributions, do.get.means, n.sim=n.sim))
})

setMethod('do.get.medians',
          signature(dist='Joint_Independent_Distributions'),
def=function(dist, n.sim=1000)
{
    unlist(sapply(dist@subdistributions, do.get.medians, n.sim=n.sim))
})

setMethod('do.get.sds',
          signature(dist='Joint_Independent_Distributions'),
def=function(dist, n.sim=1000)
{
    unlist(sapply(dist@subdistributions, do.get.sds, n.sim=n.sim))
})

setMethod('do.get.equal.tailed.intervals',
          signature(dist='Joint_Independent_Distributions'),
def=function(dist, coverage=0.95, n.sim=1000)
{
    rv = array(0, dim=c(2, dist@n.var, length(coverage)))

    for (i in 1:dist@n.subdistributions)
        rv[,dist@indices.for.subdistributions[[i]],] = do.get.equal.tailed.intervals(dist@subdistributions[[i]],
                                                                                  coverage=coverage, n.sim=n.sim)

    rv
})

setMethod('do.get.highest.density.intervals',
          signature(dist='Joint_Independent_Distributions'),
def=function(dist, coverage=0.95, n.sim=1000)
{
    rv = array(0, dim=c(2, dist@n.var, length(coverage)))

    for (i in 1:dist@n.subdistributions)
        rv[,dist@indices.for.subdistributions[[i]],] = do.get.highest.density.intervals(dist@subdistributions[[i]],
                                                                                     coverage=coverage, n.sim=n.sim)

    rv
})

setMethod('do.calculate.marginal.densities',
          signature(dist='Joint_Independent_Distributions'),
def=function(dist, x, log=F, n.sim=1000)
{
    rv = matrix(0, nrow=nrow(x), ncol=dist@n.var)
    for (i in 1:dist@n.subdistributions)
        rv[,dist@indices.for.subdistributions[[i]]] = do.calculate.marginal.densities(dist@subdistributions[[i]],
                                                                    x=matrix(x[,dist@indices.for.subdistributions[[i]]], ncol=length(dist@indices.for.subdistributions[[i]])),
                                                                    log=log, n.sim=n.sim)

    rv
})

setMethod('do.calculate.marginal.cdfs',
          signature(dist='Joint_Independent_Distributions'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000)
{
    rv = matrix(0, nrow=nrow(q), ncol=dist@n.var)
    for (i in 1:dist@n.subdistributions)
        rv[,dist@indices.for.subdistributions[[i]]] = do.calculate.marginal.cdfs(dist@subdistributions[[i]],
                                                                              q=q, lower.tail=lower.tail, log.p=log.p, n.sim=n.sim)

    rv
})

setMethod('do.generate.random.marginal.samples',
signature(dist='Joint_Independent_Distributions'),
def=function(dist, n)
{
    rv = matrix(0, nrow=n, ncol=dist@n.var)
    for (i in 1:dist@n.subdistributions)
        rv[,dist@indices.for.subdistributions[[i]]] = do.generate.random.marginal.samples(dist@subdistributions[[i]], n=n)

    rv
})

setMethod('do.subset.distribution',
          signature(dist='Joint_Independent_Distributions'),
def=function(dist, keep.indices)
{
    components = lapply(1:dist@n.subdistributions, function(i){
        keep.indices.for.subdir = intersect(keep.indices, dist@indices.for.subdistributions[[i]])
        if (length(keep.indices.for.subdir)==0)
            NULL
        else
            subset.distribution(dist@subdistributions[[i]],
                                1 + keep.indices.for.subdir - min(dist@indices.for.subdistributions[[i]]))
    })
    components = components[!sapply(components, is.null)]

    if (length(components)==1)
        components[[1]]
    else
        join.distributions(components, var.names = dist@var.names[keep.indices])
})
