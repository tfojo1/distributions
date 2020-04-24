library(distributions)

dists = join.distributions(norm=normal.distribution(0,1),
                           lognorm=lognormal.distribution(0,1))

calculate.density(dists, c(1,2))
dnorm(1) * dlnorm(2)

calculate.cdf(dists, c(1,2))
pnorm(1) * plnorm(2)


get.means(dists)


mu=0
sigma=1
dist = normal.distribution(mu,sigma)
calculate.density(dist, 0)
calculate.cdf(dist, 1)
get.quantiles(dist, c(0,.2,.6,.9))

check = data.frame(density=c(calculate.density(dist, 0), dnorm(0, mu, sigma)),
                   cdf=c(calculate.cdf(dist, 1), pnorm(1, mu, sigma)),
                   quantile=c(get.quantiles(dist, 0.4), qnorm(0.4, mu, sigma)))

mu=0
sigma=1
dist = lognormal.distribution(mu,sigma)
calculate.density(dist, 5)
calculate.cdf(dist, 1)
get.quantiles(dist, c(0,.2,.6,.9))

check = data.frame(density=c(calculate.density(dist, 5), dlnorm(5, mu, sigma)),
                   cdf=c(calculate.cdf(dist, 1), plnorm(1, mu, sigma)),
                   quantile=c(get.quantiles(dist, 0.4), qlnorm(0.4, mu, sigma)))
check

lower=0
upper=2
dist = uniform.distribution(lower,upper)
calculate.density(dist, .2)
calculate.cdf(dist, .8)
get.quantiles(dist, c(0,.2,.6,.9))

check = data.frame(density=c(calculate.density(dist, .2), dunif(.2, lower, upper)),
                   cdf=c(calculate.cdf(dist, .8), punif(.8, lower, upper)),
                   quantile=c(get.quantiles(dist, 0.4), qunif(0.4, lower, upper)))
check

dist = uniform.distribution(-Inf, Inf)
dist = loguniform.distribution(0, Inf)
