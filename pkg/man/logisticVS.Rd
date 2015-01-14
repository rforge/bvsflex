\name{logisticVS}
\alias{logisticVS}
\title{Logistic Bayesian variable selection}
\usage{
logisticVS(X, Y, b, v, block = NULL, aBeta = NULL, bBeta = NULL, 
    MCMC, thinn = 1, seed = 1234, outdir = NULL, Piupdate = FALSE, 
    burnin = 0)
}
\arguments{
 \item{X}{The input matrix of dimension n (no. samples) x p (no. variables).}
 \item{Y}{The binary response vector of length n with elements in \{0, 1\}.}
 \item{b}{The prior mean vector of length p for beta.}
 \item{v}{A vector of length p with the diagonal elements of the (diagonal) prior covariance matrix for beta.}
 \item{block}{A (sparse) p x p matrix with entries <>0, if the corresponding two variables should be updated together (optional). The default is diag(1,p).}
 \item{aBeta}{Vector of prior parameters a[i] in Beta(a[i],b[i]) which is the prior distribution of Pi[i] (i = 1,...,p). The default is rep(1,p). Note, that together with the default for bBeta this implies a Uniform[0,1] prior distribution for all Pi[i], which might not be appropriate in very high-dimensional situations with p>>n. If Piupdate=FALSE this implies that a fixed Pi[i]=0.5 is assumed for i = 1,...,p.}
 \item{bBeta}{Vector of prior parameters b[i] in Beta(a[i],b[i]) which is the prior distribution of Pi[i] (i = 1,...,p). The default is rep(1,p).}
 \item{MCMC}{Number of MCMC samples.}
 \item{thinn}{Thinning parameter. Every 'thinn'-ed MCMC iteration will be kept.}
 \item{seed}{Random seed for reproducibility.}
 \item{outdir}{A character string specifying the directory for saving MCMC output files.}
 \item{Piupdate}{Should Pi be sampled according to a prior Beta distribution? The default is FALSE, which implies that Pi is assumed to be fixed as aBeta/(aBeta+bBeta).}
 \item{burnin}{Number of burnin samples. Default is 0.}
}
\description{
The BVSflex package implements efficient Bayesian variable selection models for high-dimensional input data. A flexible selection prior allows the incorporation of additional information, e.g. a second data source including all sources of variation. 

In the logisticVS() function this is implemented for a logistic regression model.
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
res <- logisticVS(X, Y, b=rep(0,p), v=rep(1,p),
       block=diag(1,p), MCMC=B, thinn=1, seed=1234, 
       outdir=paste(getwd(), "example", sep="/"), Piupdate=TRUE)

\dontrun{
require(Matrix)
tmp <- read.table(paste(res, "GAM.txt", sep="/"))
GAM <- as.matrix(sparseMatrix(i = tmp[,1], j = tmp[,2], x = tmp[,3]))
image(y=1:p, x=1:B, z=t(GAM), col=c("white","black"), ylab="Variables", 
      xlab="MCMC iteration", main="Trace plot for GAM")
}
}
