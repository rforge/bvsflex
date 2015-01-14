rtruncBeta <- function(n,shape1,shape2,  a = -Inf, b = Inf){
  u <- runif(n, min = 0, max = 1)
  x <- qbeta(pbeta(a,shape1,shape2) + u*(pbeta(b,shape1,shape2) - pbeta(a,shape1,shape2)),shape1,shape2)
  return(x)
}

# n=1e6
# a=rbeta(n=n,shape1 = 1,shape2 = 3)
# plot(density(a))
# b=rtruncBeta(n=n,a=0,b=.5,shape1 = 1,shape2 = 3)
# lines(density(b))
# 
# a=rbeta(n=n,shape1 = 1,shape2 = 1)
# plot(density(a))
# p=3
# (g0uni=pi^2/(3*p)) #uniform
# abPi=0.5   #shape==scale parameter for beta distribution on P(y=1|g); 1 for uniform, .5 for jeffreys
# (g0jef=2*trigamma(abPi)/p)
# buni=rtruncBeta(n=n,shape1 = 1,shape2 = 1,a=0,b=g0uni/(g0uni+1))
# bjef=rtruncBeta(n=n,shape1 = 1,shape2 = 1,a=0,b=g0jef/(g0jef+1))
# lines(density(buni),col="red")
# lines(density(bjef),col="blue")
