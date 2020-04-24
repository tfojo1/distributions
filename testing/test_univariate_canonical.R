
library(distributions)

dist = Lognormal.Distribution()

test.rands = rlnorm(6)

cbind(test=calculate.density(dist, test.rands),
      check=dlnorm(test.rands))
