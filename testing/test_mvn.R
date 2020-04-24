
library(mvtnorm)
library(tmvtnorm)
library(reshape2)
library(ggplot2)

N.DIM = 5
mu = rnorm(N.DIM, 0, 20)
sigma = rWishart(1, df=N.DIM+1, Sigma=diag(rep(1, length(mu))))[,,1]

dist = Multivariate.Normal.Distribution(mu=mu, sigma=sigma, var.names=paste0('v',1:N.DIM))

test.rands = generate.random.samples(dist, 4)

#Density
x=cbind(
    dist=calculate.density(dist, test.rands),
    check=dmvnorm(test.rands, mean=mu, sigma=sigma));x

all(x[,1]==x[,2])

#CDF
x=cbind(
    dist=calculate.cdf(dist, test.rands[,1]),
    check=pmvnorm(upper=test.rands[,1], mean=mu, sigma=sigma));x

all(x[,1]==x[,2])


#Random

test.rands = generate.random.samples(dist, 1000)
df = melt(test.rands)
ggplot(df, aes(x=value)) + geom_histogram() + facet_wrap(~variable, scales='free')

#CIs

cbind(dist=get.intervals(dist)[1,],
      check=mu + qnorm(.025)*sqrt(diag(sigma)) )


## Log
dist = Multivariate.Lognormal.Distribution(mu=mu, sigma=sigma, var.names=paste0('v',1:N.DIM))


test.rands = generate.random.samples(dist, 4)
x=cbind(
    dist=calculate.density(dist, test.rands, log=T),
    check=dmvnorm(log(test.rands), mean=mu, sigma=sigma, log=T) + log(apply(1/test.rands, 1, prod)));x
all(x[,1]==x[,2])

cbind(calculate.marginal.densities(dist, test.rands)[,1])

#Random

test.rands = generate.random.samples(dist, 1000)
df = melt(test.rands)
ggplot(df, aes(x=value)) + geom_histogram() + facet_wrap(~variable, scales='free')
