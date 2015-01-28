\name{logisticVS}
\alias{logisticVS}
\title{Logistic Bayesian variable selection}
\usage{
logisticVS(X, Y, b, v, g = 1, block = NULL, 
    aBeta = NULL, bBeta = NULL, aG = NULL, bG = NULL,
    MCMC, thinn = 1, seed = 1234, outdir = NULL, Piupdate = FALSE, gupdate = "none")
}
\arguments{
 \item{X}{The input matrix with n rows (no. samples) and p columns (no. variables).}
 \item{Y}{The binary response vector of length n with elements in \{0, 1\}.}
 \item{b}{The prior mean vector for beta, of length p.}
 \item{v}{The prior covariance matrix for beta. This can be a vector of length p to be interpreted as the diagonal of a diagonal matrix, or it can be a matrix of dimension p x p.}
 \item{g}{The hyper-parameter g for the prior covariance matrix for beta (a scalar). The default is 1. Ignored if gupdate <> "none".}
 \item{block}{A (sparse) p x p matrix with entries <> 0, if the corresponding two variables should be updated together (optional). The default is diag(1,p).}
 \item{aBeta}{Vector of prior parameters a[i] in Beta(a[i],b[i]) which is the prior distribution of Pi[i] (i = 1,...,p). The default is rep(1,p). Note, that together with the default for bBeta this implies a Uniform[0,1] prior distribution for all Pi[i], which might not be appropriate in very high-dimensional situations with p>>n. If Piupdate=FALSE this implies that a fixed Pi[i]=0.5 is assumed for i = 1,...,p.}
 \item{bBeta}{Vector of prior parameters b[i] in Beta(a[i],b[i]) which is the prior distribution of Pi[i] (i = 1,...,p). The default is rep(1,p).}
 \item{aG}{Hyper-prior parameter for g (for gupdate in \{"IG", "hyperg"\}). Shape parameter in IG(aG,bG).}
 \item{bG}{Hyper-prior parameter for g (for gupdate = "IG" only). Scale parameter in IG(aG,bG).}
 \item{MCMC}{Number of MCMC samples.}
 \item{thinn}{Thinning parameter. Every 'thinn'-ed MCMC iteration will be kept.}
 \item{seed}{Random seed for reproducibility.}
 \item{outdir}{A character string specifying the directory for saving MCMC output files.}
 \item{Piupdate}{Should Pi be sampled according to a prior Beta distribution? The default is FALSE, which implies that Pi is assumed to be fixed as aBeta/(aBeta+bBeta).}
 \item{gupdate}{Class of prior distributions for the hyper-parameter g. Possible values are "IG" (inverse-gamma distribution), "hyperg" (hyper-g prior, see Liang et al. 2008), and "none". The default is "none", which implies that g is assumed to be fixed with the value specified.}
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
Holmes, C. and Held, L. (2006). ``Bayesian auxiliary variable models for binary and multinomial regression'', Bayesian Analysis 1:145-168.

Zucknick M, Richardson S (2014). ``MCMC algorithms for Bayesian variable selection in the logistic regression model for large-scale genomic applications''.
}
\author{Manuela Zucknick}
\examples{
set.seed(1234)
B <- 500
n <- 100
p <- 20
X <- matrix(rnorm(n*p), ncol=p)
beta <- c(1,1,1,rep(0,p-3))
expeta <- exp(X\%*\%beta)
Y <- rbinom(n, 1, prob=expeta/(expeta + 1))
res <- logisticVS(X, Y, b=rep(0,p), v=rep(1,p), g=1, block=diag(1,p), 
       aBeta=rep(1,p), bBeta=rep(1,p), aG=0.5, bG=0.5,
       MCMC=B, thinn=1, seed=1234, 
       outdir=paste(getwd(), "example", sep="/"), 
       Piupdate=TRUE, gupdate="IG")

\dontrun{
require(Matrix)
tmp <- read.table(paste(res, "GAM.txt", sep="/"))
GAM <- as.matrix(sparseMatrix(i = tmp[,1], j = tmp[,2], x = tmp[,3]))
image(y=1:p, x=1:B, z=t(GAM), col=c("white","black"), ylab="Variables", 
      xlab="MCMC iteration", main="Trace plot for GAM")
}
}
