
##--------------------------------------##
##-- CLASS DEFINITION and CONSTRUCTOR --##
##--------------------------------------##

#'@title The Smoothed_Empiric_Distribution Class
#'
#'@description A class that defines a distribution based on a set of samples from that distribution, but applies Gaussian kernel smoothers. For an unsmoothed representation, see \link{Empiric_Distribution}
#'
#'@slot points A matrix containing the samples
#'@slot weights A numeric vector with the weight given to each row in points. Sum to 1
#'@slot n.points The number of rows in points
#'
#'@seealso  \link{Empiric_Distribution}, \link{Smoothed.Empiric.Distribution}
#'
#'@name Smoothed_Empiric_Distribution
#'@rdname Smoothed_Empiric_Distribution
#'@aliases Smoothed_Empiric_Distribution-class
#'@exportClass Smoothed_Empiric_Distribution
#'@export
setClass("Smoothed_Empiric_Distribution",
         contains="Empiric_Distribution",
         representation = list(marginal.smoothed.distributions='Joint_Independent_Distributions'))


#'@title Create an smoothed empiric distribution from a set of samples
#'
#'@param points A matrix of values, where each row represents one sample and each column represents one variable. If a numeric vector is passed, assumed to be samples from a univariate distribution
#'@param weights A numeric vector of weights, one for each row in points
#'@param var.names The names of the variables in the distribution. If NULL, uses the column names from points
#'
#'@family Distribution Constructors
#'
#'@export
Smoothed.Empiric.Distribution <- function(points, weights=1, var.names=NULL,
                                          bandwidth.multipliers=1)
{
    if (is.null(dim(points)))
        points = matrix(points, ncol=1)

    n.points = dim(points)[1]

    if (length(weights)==1)
        weights = rep(1, n.points)

    if (length(weights)!=n.points)
        stop("weights must have one value for each row in points")

    if (any(is.na(points)))
        stop("points cannot contain NA values")

    if (any(is.na(weights)))
        stop("weights cannot contain NA values")
    weights = weights / sum(weights)

    if (is.null(var.names) && !is.null(dim(points)))
        var.names = dimnames(points)[[2]]

    create.distribution(subclass='Smoothed_Empiric_Distribution',
                        var.names=var.names,
                        n.var=dim(points)[2],
                        lower.bounds = apply(points, 2, min),
                        upper.bounds = apply(points, 2, max),
                        is.discrete=T,
                        is.improper=F,
                        points=points,
                        weights=weights,
                        n.points=n.points)
}

##-----------------------------##
##-- METHODS IMPLEMENTATIONS --##
##-----------------------------##

setMethod('do.calculate.density',
          signature(dist='Smoothed_Empiric_Distribution'),
def=function(dist, x, log=F, n.sim=1000)
{
    rv = apply(x, 1, function(one.x){
        row.equals = apply(dist@points, function(one.points){
            all(one.points==one.x)
        })
        sum(row.equals * dist@weights)
    })

    if (log)
        log(rv)
    else
        rv
})


setMethod('do.calculate.cdf',
          signature(dist='Smoothed_Empiric_Distribution'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000)
{
    rv = apply(q, 1, function(one.q){
        include.row = apply(dist@points, function(one.points){
            if (lower.tail)
                all(one.points<=one.q)
            else
                all(one.points>=one.q)
        })
        sum(include.row * dist@weights)
    })

    if (log.p)
        log(rv)
    else
        rv
})


setMethod('do.get.quantiles',
signature(dist='Smoothed_Empiric_Distribution'),
def=function(dist, p, lower.tail=T, log.p=F, n.sim=1000)
{
    if (log.p)
        p = exp(p)

    apply(dist@points, 2, function(x){
        o = order(x, decreasing = !lower.tail)
        sorted.points = x[o]
        cum.weight = cumsum(dist@weights[o])

        sapply(p, function(one.p){
            if (one.p < 0 || one.p > 1)
                NaN
            else
                sorted.points[one.p <= cum.weight][1]
        })
    })
})

setMethod('do.generate.random.samples',
          signature(dist='Smoothed_Empiric_Distribution'),
def=function(dist, n)
{
    indices = sample(1:dist@n.points, n, replace = T)
    dist@points[indices,]
})

setMethod('do.get.covariance.matrix',
          signature(dist='Smoothed_Empiric_Distribution'),
def=function(dist, n.sim=1000)
{
    sample.means = get.means(dist)
    normalizing.constant = 1 / (1 - sum(dist@weights^2))
    sapply(1:dist@n.points, function(i){
        sapply(1:dist@n.points, function(j){
            sum(dist@weights * (dist@points[,i]-sample.means[i]) * (dist@points[,j]-sample.means[j]))
        })
    }) * normalizing.constant
})

setMethod('do.get.sds',
          signature(dist='Smoothed_Empiric_Distribution'),
def=function(dist, n.sim=1000)
{
    sample.means = get.means(dist)
    normalizing.constant = 1 / (1 - sum(dist@weights^2))
    variances = sapply(1:dist@n.points, function(i){
        sum(dist@weights * (dist@points[,i]-sample.means[i])^2)
    }) * normalizing.constant
    sqrt(variances)
})

# do.get.medians and do.get.equal.tailed.intervals use the default implementations

setMethod('do.get.means',
signature(dist='Smoothed_Empiric_Distribution'),
def=function(dist, n.sim=200)
{
    colSums(dist@points * dist@weights)
})

setMethod('do.get.highest.density.intervals',
          signature(dist='Smoothed_Empiric_Distribution'),
def=function(dist, coverage=0.95, n.sim=1000)
{
    do.get.highest.density.interval.from.samples(samples=dist@points, coverage=coverage, weights=dist@weights)
})

setMethod('do.calculate.marginal.densities',
          signature(dist='Smoothed_Empiric_Distribution'),
def=function(dist, x, log=F, n.sim=1000)
{
    rv = apply(1:dist@n.var, function(i){
        sapply(x, function(one.x){
            sum(dist@weights * (dist@points[,i]==one.x))
        })
    })

    if (log)
        log(rv)
    else
        rv
})

setMethod('do.calculate.marginal.cdfs',
          signature(dist='Smoothed_Empiric_Distribution'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000)
{
    rv = apply(1:dist@n.var, function(i){
        sapply(q, function(one.q){
            if (lower.tail)
                sum(dist@weights * (dist@points[,i]<=one.q))
            else
                sum(dist@weights * (dist@points[,i]>=one.q))
        })
    })

    if (log.p)
        log(rv)
    else
        rv
})

setMethod('do.subset.distribution',
          signature(dist='Smoothed_Empiric_Distribution'),
def=function(dist, keep.indices)
{
    create.empiric.distribution(points=dist@points[,keep.indices],
                                weights=dist@weights,
                                var.names=dist@var.names[keep.indices])
})

