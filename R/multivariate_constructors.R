
#'@title Create a Multivariate Normal Distribution Object
#'
#'@param mu,sigma The mean vector and covariance matrix of the distribution
#'@param var.names The names of the variables in the distribution. If NULL, uses the names of mu, the row names of sigma, or the column names of sigma, in that order
#'@param lower,upper The upper and lower bounds for each variable, if bounded. If only one value is given, assumed to apply to all variables in the distribution
#'
#'@family Distribution Constructors
#'
#'@export
Multivariate.Normal.Distribution <- function(mu=rep(0,nrow(sigma)),
                                             sigma=diag(rep(1, length(mu))),
                                             lower=-Inf,
                                             upper=Inf,
                                             var.names=NULL)
{
    new('Multivariate_Normal_Distribution',
        mu=mu,
        sigma=sigma,
        lower=lower,
        upper=upper,
        var.names=var.names,
        transformations='identity')
}

#'@title Create a Multivariate Log-Normal Distribution Object
#'
#'@param mu,sigma The mean vector and covariance matrix of the distribution (after log transformation)
#'@param lower,upper The upper and lower bounds for each variable (on the original scale, NOT the log scale). If only one value is given, assumed to apply to all variables in the distribution
#'@inheritParams Multivariate.Normal.Distribution
#'
#'@family Distribution Constructors
#'
#'@export
Multivariate.Lognormal.Distribution <- function(mu=rep(0,nrow(sigma)),
                                                sigma=diag(rep(1, length(mu))),
                                                lower=0,
                                                upper=Inf,
                                                var.names=NULL)
{
    if (any(lower<0))
        stop("A log-normal distribution cannot have a lower bound less than zero")

    if (length(mu) < dim(sigma)[1])
        mu = rep(mu, dim(sigma)[1])

    new('Multivariate_Normal_Distribution',
        mu=mu,
        sigma=sigma,
        lower=lower,
        upper=upper,
        var.names=var.names,
        transformations='log')
}

#'@title Create a Multivariate Logit-Normal Distribution Object
#'
#'@param mu,sigma The mean vector and covariance matrix of the distribution (after logit transformation)
#'@param lower,upper The upper and lower bounds for each variable (on the original scale, NOT the logit scale). If only one value is given, assumed to apply to all variables in the distribution
#'@inheritParams Multivariate.Normal.Distribution
#'
#'@family Distribution Constructors
#'
#'@export
Multivariate.Logitnormal.Distribution <- function(mu=rep(0,nrow(sigma)),
                                                sigma=diag(rep(1, length(mu))),
                                                lower=0,
                                                upper=1,
                                                var.names=NULL)
{
    if (any(lower<0))
        stop("A logit-normal distribution cannot have a lower bound less than zero")
    if (any(upper>1))
        stop("A logit-normal distribution cannot have an upper bound greater than one")

    new('Multivariate_Normal_Distribution',
        mu=mu,
        sigma=sigma,
        lower=lower,
        upper=upper,
        var.names=var.names,
        transformations='logit')
}

#'@title Create a transformed multivariate normal distribution object
#'
#'@description Creates a distribution where all the variables, after arbitrary transformations, follow a normal distribution
#'
#'@param mu,sigma The mean vector and covariance matrix of the distribution (after transformations)
#'@param lower,upper The upper and lower bounds for each variable (on the original scale, NOT the transformed scale). If only one value is given, assumed to apply to all variables in the distribution
#'@param transformations Objects of class \link{transformation} to be applied to each variable. Can be either an object of class \link{transformation}, the character name of a defined transformation, or a list of \link{transformation} objects.   If only one transformation is supplied, assumed to apply to all variables. If this is a named list, any variables missing transformation are assmed to not be transformed
#'@inheritParams Multivariate.Normal.Distribution
#'
#'@family Distribution Constructors
#'
#'@export
Transformed.Multivariate.Normal.Distribution <- function(mu=rep(0,nrow(sigma)),
                                                              sigma=diag(rep(1, length(mu))),
                                                              lower=-Inf,
                                                              upper=Inf,
                                                              var.names=NULL,
                                                              transformations=NULL)
{
    if (is.null(transformations))
        transformation.names = 'identity'
    else if (is(transformations, 'transformation'))
        transformation.names = transformations@name
    else if (is(transformations, 'character'))
        transformation.names = transformations
    else if (is(transformations, 'list'))
        transformation.names = sapply(transformations, function(transformation){
            if (is.null(transformation))
                'identity'
            else if (is(transformation, 'transformation'))
                transformation@name
            else if (is(transformation, 'character'))
                transformation
            else
                stop("'transformations' must be an object of class 'transformation' or the name of a pre-defined transformation, or a list whose elements are either an object of class 'transformation' or the name of a pre-defined transformation")
        })
    else
        stop("'transformations' must be an object of class 'transformation' or the name of a pre-defined transformation, or a list whose elements are either an object of class 'transformation' or the name of a pre-defined transformation")

    if (all(transformation.names=='identity'))
        Multivariate.Normal.Distribution(mu=mu, sigma=sigma, lower=lower, upper=upper, var.names = var.names)
    else if (all(transformation.names=='log'))
        Multivariate.Lognormal.Distribution(mu=mu, sigma=sigma, lower=lower, upper=upper, var.names=var.names)
    else if (all(transformation.names=='logit'))
        Multivariate.Logitnormal.Distribution(mu=mu, sigma=sigma, lower=lower, upper=upper, var.names=var.names)
    else
        new('Multivariate_Normal_Distribution',
            mu=mu,
            sigma=sigma,
            lower=lower,
            upper=upper,
            var.names=var.names,
            transformations=transformations)
}


#'@title Create a multivariate normal distribution with a Continuous Autoregressive 1 (CAR-1) correlation structure
#'
#'@description Creates a distribution where all the variables, after arbitrary transformations, follow a normal distribution with a Continuous Autoregressive 1 (CAR-1) correlation structure
#'
#'@param rho The autoregressive correlation coefficient
#'@param times The times for each observation in the distribution. The correlation between any two observations at times t1 and t2 is rho^(abs(t1-t2))
#'@param mu,sds The mean and standard deviation vectors. If given a scalar value, assumes that the mean/sd is the same for every variable in the distribution
#'@inheritParams Transformed.Multivariate.Normal.Distribution
#'
#'@family Distribution Constructors
#'
#'@export
Autoregressive.Multivariate.Normal.Distribution <- function(rho,
                                                            times,
                                                            mu=0,
                                                            sds=1,
                                                            lower=-Inf,
                                                            upper=Inf,
                                                            var.names=NULL,
                                                            transformations=NULL)
{
    if (length(mu) == 1)
        mu = rep(mu, length(times))

    if (length(sds) == 1)
        sds = rep(sds, length(times))

    if (rho < 0 || rho > 1)
        stop("rho must be between 0 and 1")
    sigma = rho ^ matrix(abs(rep(times, each=length(times)) - rep(times, length(times))),
                         nrow=length(times), ncol=length(times)) *
        sds %*% t(sds)

    Transformed.Multivariate.Normal.Distribution(mu=mu,
                                                 sigma=sigma,
                                                 lower=lower,
                                                 upper=upper,
                                                 var.names=var.names,
                                                 transformations=transformations)
}

#'@title Create a multivariate normal distribution with a Compound Symmetry correlation structure
#'
#'@description Creates a distribution where all the variables, after arbitrary transformations, follow a normal distribution with a Compound Symmetry (AKA Exchangeable) correlation structure
#'
#'@param rho The correlation coefficient
#'@param n The number of variables in the distribution
#'@param mu,sds The mean and standard deviation vectors. If given a scalar value, assumes that the mean/sd is the same for every variable in the distribution
#'@inheritParams Transformed.Multivariate.Normal.Distribution
#'
#'@family Distribution Constructors
#'
#'@export
Compound.Symmetry.Multivariate.Normal.Distribution <- function(rho,
                                                               n=max(length(mu), length(sds)),
                                                               mu=0,
                                                               sds=1,
                                                               lower=-Inf,
                                                               upper=Inf,
                                                               var.names=NULL,
                                                               transformations=NULL)
{
    if (length(mu) == 1)
        mu = rep(mu, n)

    if (length(sds) == 1)
        sds = rep(sds, n)

    if (rho < 0 || rho > 1)
        stop("rho must be between 0 and 1")

    sigma = matrix(rho, nrow=n, ncol=n)
    diag(sigma) = 1
    sigma = sigma * sds %*% t(sds)

    Transformed.Multivariate.Normal.Distribution(mu=mu,
                                                 sigma=sigma,
                                                 lower=lower,
                                                 upper=upper,
                                                 var.names=var.names,
                                                 transformations=transformations)
}

#'@title Create a Multivariate Uniform Distribution Object Where the Variables are Correlated
#'
#'@description Uses a rough approximation to generate a multivariate normal distribution such that the correlations between variables are approximately equal to the given correlation.matrix. Specifically, generates random samples from a multivariate normal distribution with mean=0 and covariance matrix = cov2cor(correlation.matrix + 0.04*I)
#'
#'@param min,max The vectors of minima and maxima for each variable in the distribution. If given a scalar value, assumes that value applies for all variables
#'@param correlation.matrix The matrix specifying how variables should be correlated
#'@param var.names The names of the variables in the distribution. If NULL, uses the names of mu, the row names of sigma, or the column names of sigma, in that order
#'
#'@family Distribution Constructors
#'
#'@export
Multivariate.Correlated.Uniform.Distribution <- function(correlation.matrix,
                                                         min=rep(0,dim(correlation.matrix)[1]),
                                                         max=rep(1,dim(correlation.matrix)[1]),
                                                         var.names=NULL)
{
    # Check covariance matrix
    if (is.null(correlation.matrix) ||
        (!is(correlation.matrix, 'matrix') && !is(correlation.matrix, 'array')) ||
        length(dim(correlation.matrix))!=2 || dim(correlation.matrix)[1] != dim(correlation.matrix)[2])
        stop("correlation matrix must be a square matrix")

    correlation.matrix = cov2cor(correlation.matrix)
    if (!matrixcalc::is.positive.definite(correlation.matrix))
        stop("correlation matrix must be positive definite")

    n = dim(correlation.matrix)[1]
    correlation.matrix = cov2cor(correlation.matrix + matrix(.04, nrow=n, ncol=n))

    # Check min and max
    if (length(min)==1)
        min = rep(min, n)
    if (length(max)==1)
        max = rep(max, n)

    if (length(min) != n || length(max) != n)
        stop("min and max must be either length one, or have the same length as the dimension of correlation.matrix")
    if (any(min>=max))
        stop("min must be < max")

    # Set up the transformations
    transformations = lapply(1:n, function(i){
        create.transformation(transform.function = function(x){qnorm( (x-min[i]) / (max[i]-min[i]) )},
                              reverse.transform.function = function(x){pnorm(x) * (max[i]-min[i]) + min[i]},
                              transformation.derivative = function(x){
                                  x = (x-min[i]) / (max[i]-min[i])
                                  inv.rv = dnorm(qnorm(x))
                                  rv = 1/inv.rv
                                  rv[inv.rv==0] = 0
                                  rv
                              },
                              log.abs.transformation.derivative = function(x){
                                  x = (x-min[i]) / (max[i]-min[i])
                                  inv.rv = dnorm(qnorm(x), log=T)
                                  rv = -inv.rv
                                  rv[inv.rv==-Inf] = -Inf
                                  rv
                              },
                              name=paste0('qnorm_shift_to_', min[i], "_", max[i]))
    })
    Transformed.Multivariate.Normal.Distribution(mu=rep(0,n),
                                                 sigma=correlation.matrix,
                                                 lower=min,
                                                 upper=max,
                                                 var.names=var.names,
                                                 transformations=transformations)
}
