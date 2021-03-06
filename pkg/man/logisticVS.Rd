\name{logisticVS}
\alias{logisticVS}
\title{Logistic Bayesian variable selection}
\usage{
logisticVS(X, Y, b, h0, g = 1, block = NULL, 
    mu = NULL, phi = NULL, m0 = NULL, g0 = NULL, 
    aPhi = 1, bPhi = 0.1, aG = NULL, bG = NULL,
    MCMC, thinn = 1, seed = 1234, outdir = NULL, 
    Piupdate = FALSE, MuPhiUpdate = FALSE, gupdate = "none", wlsgprior=FALSE, 
    sigmaG = 1, sigmaMu = 1, sigmaPhi = 1)
}

\arguments{
 \item{X}{The input matrix with n rows (no. samples) and p columns (no. variables).}
 \item{Y}{The binary response vector of length n with elements in \{0, 1\}.}
 \item{b}{The prior mean vector for beta, of length p.}
 \item{h0}{The part of the prior precision matrix for beta, which is to be multiplied by 1/g. This can be a vector of length p to be interpreted as the diagonal of a diagonal matrix, or it can be a matrix of dimension p x p.}
 \item{g}{The hyper-parameter g for the prior covariance matrix for beta (a scalar). The default is 1. Ignored if gupdate <> "none".}
 \item{block}{A (sparse) p x p matrix with entries <> 0, if the corresponding two variables should be updated together (optional). The default is diag(1,p).}
 \item{mu}{Vector of prior mean mu[i] of Beta prior distribution of Pi[i] (i = 1,...,p). The default is rep(0.5,p). Note, that together with the default for phi this implies a Uniform[0,1] prior distribution for all Pi[i], which might not be appropriate in very high-dimensional situations with p>>n. If Piupdate=FALSE this implies that a fixed Pi[i]=0.5 is assumed for i = 1,...,p.}
 \item{phi}{Vector of prior precision parameters phi[i] of Beta prior distribution of Pi[i] (i = 1,...,p). The default is rep(2,p).}
 \item{aPhi, bPhi}{Hyper-prior shape and RATE parameters for Gamma prior of phi. Default values are 1 and 0.1. (DEFAULTS PROBABLY NOT GOOD)}
 \item{m0, g0}{Hyper-prior mean and precision parameters for Beta prior of mu. This can be made informative by including external information. Default values are rep(0.5,p) and rep(2,p) - see help for mu and phi parameters for more details.}
 \item{aG}{Hyper-prior parameter for g (for gupdate in \{"IG", "hyperg"\}). Shape parameter in IG(aG,bG) (aG>0, default=0.5), or in hyperg(aG,bG) (aG>0, default=1).}
 \item{bG}{Hyper-prior parameter for g (for gupdate in \{"IG", "hyperg"\}). Scale parameter in IG(aG,bG) or in hyperg(aG,bG) (bG>0, default=0.5 in both cases).}
 \item{MCMC}{Number of MCMC samples.}
 \item{thinn}{Thinning parameter. Every 'thinn'-ed MCMC iteration will be kept.}
 \item{seed}{Random seed for reproducibility.}
 \item{outdir}{A character string specifying the directory for saving MCMC output files.}
 \item{Piupdate}{Should Pi be sampled according to a prior Beta distribution? The default is FALSE, which implies that Pi is assumed to be fixed as aBeta/(aBeta+bBeta).}
 \item{MuPhiUpdate}{Should mu and phi be sampled according to a prior Beta and Gamma distribution? The default is FALSE.}
 \item{gupdate}{Class of prior distributions for the hyper-parameter g. Possible values are "IG" (inverse-gamma distribution), "hyperg" (hyper-g prior, see Liang et al. 2008 and Sheng and Li 2014), and "none". The default is "none", which implies that g is assumed to be fixed with the value specified.}
 \item{wlsgprior}{Should the 'weighted least squares g-prior' be used? If TRUE, then h[gamma] is computed as 1/g * 1/n t(x[gamma]) \%*\% solve(Lambda) \%*\% x[gamma] in each iteration, instead of using h = 1/g * h0 (i.e. h0 is not used in this case).}
 \item{sigmaG}{Standard deviation of proposal density for log(g) in gupdate='hyperg' (sigmaG>0, default=1).}
 \item{sigmaMu}{Standard deviation of proposal density for log(mu) in gupdate='hyperg' (sigmaMu>0, default=1).}
 \item{sigmaPhi}{Standard deviation of proposal density for log(phi) in gupdate='hyperg' (sigmaPhi>0, default=1).}
}
\description{
The BVSflex package implements efficient Bayesian variable selection models for high-dimensional input data. A flexible selection prior allows the incorporation of additional information, e.g. a second data source including all sources of variation. 

In the logisticVS() function this is implemented for a logistic regression model.
}
\details{
The model does not contain an intercept. Therefore, the input variables (columns of X) as well as the outcome variable Y are centered automatically to zero mean, before the model is run. 

Note that the data are not scaled automatically to unit variance.
}
\value{
outdir (The directory containing all the results from the MCMC run.)
}
\references{
Holmes C, Held L (2006). ``Bayesian auxiliary variable models for binary and multinomial regression'', Bayesian Analysis 1:145-168.

Liang F, Paulo R, Molina G, Clyde MA, Berger JO (2008). ``Mixtures of g priors for Bayesian variable selection. Journal of the American Statistical Association'' 103:410-423.

Shang Z, Li P (2014). ``Bayesian high-dimensional screening via MCMC''. Journal of Statistical Planning and Inference 155:54-78.

Zucknick M, Richardson S (2014). ``MCMC algorithms for Bayesian variable selection in the logistic regression model for large-scale genomic applications''.
}
\author{Manuela Zucknick}
\examples{
# Low-dimensional example (this will take several seconds):
set.seed(1234)
B <- 1000
n <- 100
p <- 18
d <- 3
X <- matrix(rnorm(n*p), ncol=p)
beta <- c(rep(1,d),rep(0,p-d))
expeta <- exp(X\%*\%beta)
Y <- rbinom(n, 1, prob=expeta/(expeta + 1))
res <- logisticVS(X, Y, b=rep(0,p), h0=1/n * t(X)\%*\%X,
       mu=rep(d/p,p), phi=rep(p,p), aG=0.5, bG=0.5*2.5^2,
       MCMC=B, thinn=1, seed=1234, 
       outdir=paste(getwd(), "example", sep="/"), 
       Piupdate=TRUE, MuPhiUpdate=FALSE, gupdate="IG", wlsgprior=FALSE)

# Low-dimensional example with hyper-priors mu and phi (this will take several more seconds):
set.seed(1234)
B <- 5000
res <- logisticVS(X, Y, b=rep(0,p), h0=1/n * t(X)\%*\%X,
       mu=rep(d/p,p), phi=rep(p,p), aG=0.5, bG=0.5*2.5^2,
       m0=rep(d/p,p),  g0=rep(2*p,p), MCMC=B, thinn=1, seed=1234, 
       outdir=paste(getwd(), "example", sep="/"), 
       Piupdate=TRUE, MuPhiUpdate=TRUE, gupdate="IG", wlsgprior=FALSE)

\dontrun{
require(Matrix)
tmp <- read.table(paste(res, "GAM.txt", sep="/"))
GAM <- as.matrix(sparseMatrix(i = tmp[,1], j = tmp[,2], x = tmp[,3]))
image(y=1:p, x=1:B, z=t(GAM), col=c("white","black"), ylab="Variables", 
      xlab="MCMC iteration", main="Trace plot for GAM")}

\dontrun{
# High-dimensional example (this will take several minutes):
set.seed(1234)
B <- 100000
n <- 100
p <- 500
d <- 3
X <- matrix(rnorm(n*p), ncol=p)
beta <- c(rep(2,d),rep(0,p-d))
expeta <- exp(X\%*\%beta)
Y <- rbinom(n, 1, prob=expeta/(expeta + 1))
res <- logisticVS(X, Y, b=rep(0,p), h0=rep(1,p), g=1, block=diag(1,p), 
       mu=rep(d/p,p), phi=rep(p/2,p),
       MCMC=B, thinn=1, seed=1234, 
       outdir=paste(getwd(), "example", sep="/"), 
       Piupdate=FALSE, MuPhiUpdate=FALSE, gupdate="none")
require(Matrix)
tmp <- read.table(paste(res, "GAM.txt", sep="/"))
GAM <- as.matrix(sparseMatrix(i = tmp[,1], j = tmp[,2], x = tmp[,3]))
plot(apply(GAM,2,sum), ylab="Model size")
plot(apply(GAM,1,mean), ylab="Posterior inclusion estimates")}
}
