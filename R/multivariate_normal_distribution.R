


##--------------------------------------##
##-- CLASS DEFINITION and CONSTRUCTOR --##
##--------------------------------------##

#'@title The Multivariate_Normal_Distribution Class
#'
#'@description Defines a multivariate normal distribution
#'
#'@slot mu,sigma The mean vector and covariance matrix of the distribution
#'
#'@name Multivariate_Normal_Distribution
#'@rdname Multivariate_Normal_Distribution
#'@aliases Multivariate_Normal_Distribution-class
#'@exportClass Multivariate_Normal_Distribution
#'@export
setClass("Multivariate_Normal_Distribution",
         contains="Joint_Independent_Distributions",
         representation = list(mu='numeric',
                               sigma='matrix',
                               bounded='logical',
                               transformations='list',
                               transformed.lower.bounds='numeric',
                               transformed.upper.bounds='numeric'))


setMethod('initialize',
          signature(.Object='Multivariate_Normal_Distribution'),
def=function(.Object,
             mu, sigma=diag(rep(1, length(mu))),
             support=NULL, var.names=NULL,
             lower=-Inf, upper=Inf,
             transformations=NULL)
{
    # Check mu
    if (!is(mu, 'numeric') && !is(mu, 'integer') && !is(mu, 'matrix') && !is(mu, 'array'))
        stop('mu must be a non-NA, numeric vector or matrix')
    if (length(mu)==0 || any(is.na(mu)))
        stop('mu must be a non-NA, numeric vector')

    n.var = length(mu)

    # Check sigma
    if (is.null(dim(sigma)) || length(dim(sigma))!=2 ||
        dim(sigma)[1] != n.var || dim(sigma)[2] != n.var)
        stop(paste0("For mu with length ", n.var, ", sigma must be a ", n.var, "x", n.var, " matrix"))
    if (any(is.na(sigma)))
        stop("sigma cannot contain NA values")
    if (!is(sigma[1,1], 'numeric') && !is(sigma[1,1], 'integer'))
        stop("sigma can only contain numeric values")
    if (!matrixcalc::is.positive.semi.definite(sigma))
        stop("sigma must be positive semi-definite")

    # Check support
    if (is.null(support))
    {
        if (length(lower)==1)
            lower = rep(lower, n.var)
        if (length(upper)==1)
            upper = rep(upper, n.var)

        if (length(lower) != n.var || length(upper) != n.var)
            stop("lower and upper must either be length 1 or have one value for each variable in the distribution")

        support = Multivariate.Support(lapply(1:n.var, function(i){Continuous.Support(lower[i], upper[i])}))
    }

    if (support@n.var==1 && n.var!=1)
        support = Multivariate.Support(lapply(1:n.var, function(i){support}))

    if (any(!support@is.contiguous) || any(support@is.discrete))
        stop("Multivariate_Normal_Distributions must have continuous, contiguous support")

    # Check names
    if (is.null(var.names))
    {
        if (!is.null(names(mu)))
            var.names = names(mu)
        else if (!is.null(dimnames(sigma)[[1]]))
            var.names = dimnames(sigma)[[1]]
        else if (!is.null(dimnames(sigma)[[2]]))
            var.names = dimnames(sigma)[[2]]
    }
    else
    {
        if (length(var.names) != n.var)
            stop(paste0("mu has length ", n.var, " but var.names has length ", length(var.names)))
    }

    names(mu) = var.names
    dimnames(sigma) = list(var.names, var.names)

    # Match up transformations
    transformations = match.transformations.to.variables(transformations,
                                                         var.names=var.names, n.var=n.var)

    # Get the bounds
    bounds = get.support.bounds(support)
    transformed.bounds = do.transform.variables(transformations, bounds)

    # Get the marginal distributions
    marginal.distributions = lapply(1:n.var, function(i){
        Transformed.Normal.Distribution(mean=mu[i], sd=sqrt(sigma[i,i]),
                                        lower=bounds['lower',i], upper=bounds['upper',i],
                                        var.name=var.names[i],
                                        transformation = transformations[[i]])
    })

    .Object = callNextMethod(.Object,
                             marginal.distributions,
                             var.names=var.names)

    .Object@mu = as.numeric(mu)
    .Object@sigma = sigma
    .Object@bounded = apply(transformed.bounds, 2, function(bounds){
        bounds[1] > -Inf || bounds[2] < Inf
    })
    .Object@transformations = transformations
    .Object@transformed.lower.bounds = transformed.bounds['lower',]
    .Object@transformed.upper.bounds = transformed.bounds['upper',]

    .Object
})

##-----------------------------##
##-- METHODS IMPLEMENTATIONS --##
##-----------------------------##

setMethod('do.calculate.density',
signature(dist='Multivariate_Normal_Distribution'),
def=function(dist, x, log=F, n.sim=1000)
{
    if (log)
        rv = rep(-Inf, dim(x)[1])
    else
        rv = rep(0, dim(x)[1])

    in.bounds = is.supported(dist@support, x)
    x.in.bounds = matrix(x[in.bounds,], nrow=sum(in.bounds))

    if (any(in.bounds))
    {
        transformed.x = do.transform.variables(dist@transformations, x.in.bounds)

        if (any(dist@bounded))
        {
            rv[in.bounds] = tmvtnorm::dtmvnorm(x=transformed.x, mean=dist@mu, sigma=dist@sigma,
                               lower=dist@transformed.lower.bounds, upper=dist@transformed.upper.bounds,
                               log=log)
        }
        else
        {
            rv[in.bounds] = mvtnorm::dmvnorm(x=transformed.x, mean=dist@mu, sigma=dist@sigma, log=log)
        }

        if (log)
        {
            log.derivative = rowSums(do.log.abs.transformation.derivative.for.variables(dist@transformations, x.in.bounds))
            rv[in.bounds] = rv[in.bounds] + log.derivative
        }
        else
        {
            derivative = apply(do.transformation.derivative.for.variables(dist@transformations, x.in.bounds), 1, prod)
            rv[in.bounds] = rv[in.bounds] * derivative
        }
    }

    rv
})


setMethod('do.calculate.cdf',
          signature(dist='Multivariate_Normal_Distribution'),
def=function(dist, q, lower.tail=T, log.p=F, n.sim=1000)
{
    transformed.q = do.transform.variables(dist@transformations, x=q)
    if (any(dist@bounded))
    {
        rv = apply(transformed.q, 1, function(one.q){
            if (lower.tail)
                tmvtnorm::ptmvnorm(lowerx=dist@transformed.lower.bounds, upperx=one.q,
                                   mean=dist@mu, sigma=dist@sigma,
                                   lower=dist@transformed.lower.bounds, upper=dist@transformed.upper.bounds)
            else
                tmvtnorm::ptmvnorm(lowerx=one.q, upperx=dist@transformed.upper.bounds,
                                   mean=dist@mu, sigma=dist@sigma,
                                   lower=dist@transformed.lower.bounds, upper=dist@transformed.upper.bounds)
        })
    }
    else
    {
        rv = apply(transformed.q, 1, function(one.q){
            if (lower.tail)
                mvtnorm::pmvnorm(upper=one.q, mean=dist@mu, sigma=dist@sigma)
            else
                mvtnorm::pmvnorm(lower=one.q, mean=dist@mu, sigma=dist@sigma)
        })
    }

    if (log.p)
        log(rv)
    else
        rv
})

#Handles (as marginals) by the superclass
# do.get.quantiles
# do.get.sds
# do.get.medians
# do.get.equal.tailed.intervals
# do.get.means
# do.get.highest.density.intervals
# do.calculate.marginal.densities
# do.calculate.marginal.cdfs

setMethod('do.generate.random.samples',
          signature(dist='Multivariate_Normal_Distribution'),
def=function(dist, n)
{
    if (any(dist@bounded))
        transformed.rv = tmvtnorm::rtmvnorm(n, mean=dist@mu, sigma=dist@sigma,
                                            lower=dist@transformed.lower.bounds,
                                            upper=dist@transformed.upper.bounds)
    else
        transformed.rv = mvtnorm::rmvnorm(n, mean=dist@mu, sigma=dist@sigma)

    dim(transformed.rv) = c(n, dist@n.var)
    do.reverse.transform.variables(dist@transformations, transformed.rv)
})

setMethod('do.get.covariance.matrix',
          signature(dist='Multivariate_Normal_Distribution'),
def=function(dist, n.sim=1000)
{
    pull.covariance.from.sigma = !dist@bounded & sapply(dist@transformations, function(tf){tf@name=='identity'})
    browser()
    if (all(pull.covariance.from.sigma))
        dist@sigma
    else
    {
        rands = generate.random.samples(dist, n=n.sim)
        cov.mat = cov(rands)
        cov.mat[pull.covariance.from.sigma, pull.covariance.from.sigma] = dist@sigma[pull.covariance.from.sigma, pull.covariance.from.sigma]
        sds = get.sds(dist, n.sim=n.sim)
        corr.mat = cov2cor(cov.mat)
        sds %*% t(sds) * corr.mat
    }
})



setMethod('do.subset.distribution',
signature(dist='Multivariate_Normal_Distribution'),
def=function(dist, keep.indices)
{
    if (length(keep.indices)==1)
    {
        bounds = get.support.bounds(dist@support)[,keep.indices]
        Transformed.Normal.Distribution(mean=dist@mu[keep.indices],
                                        sd=sqrt(dist@sigma[keep.indices, keep.indices]),
                                        tranformation=transformations[[keep.indices]],
                                        lower = bounds[1], upper = bounds[2],
                                        var.names = var.names[keep.indices])
    }
    else
    {
        new('Multivariate_Normal_Distribution',
            mu=dist@mu[keep.indices],
            sigma=dist@sigma[keep.indices, keep.indices],
            support=subset.support(dist@support, keep.indices),
            var.names=dist@var.names[keep.indices],
            transformations=dist@transformations[keep.indices])
    }
})


