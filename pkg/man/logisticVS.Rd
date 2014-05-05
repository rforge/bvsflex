\name{logisticVS}
\alias{logisticVS}
\title{Logistic Bayesian variable selection}
\usage{
logisticVS(X, Y, b, v, k, d, block=NULL, MCMC, thinn=1, seed=1234, outdir=NULL)
}
\arguments{
 \item{X}{The input matrix of dimension n (no. samples) x p (no. variables).}
 \item{Y}{The binary response vector of length n with elements in \{0, 1\}.}
 \item{b}{The prior mean vector of length p for beta.}
 \item{v}{A vector of length p with the diagonal elements of the (diagonal) prior covariance matrix for beta.}
 \item{k}{A scalar specifying the overall expected number of variables in a model (E(Pi[i] = k/p).}
 \item{d}{A vector of length p containing the distances for variables i=1,...,p.}
 \item{block}{A (sparse) p x p matrix with entries <>0, if the corresponding two variables should be updated together (optional). The default is diag(1,p).}
 \item{MCMC}{Number of MCMC samples.}
 \item{thinn}{Thinning parameter. Every 'thinn'-ed MCMC iteration will be kept.}
 \item{seed}{Random seed for reproducibility.}
 \item{outdir}{A character string specifying the directory for saving MCMC output files.}
}
\description{
The BVSflex package implements efficient Bayesian variable selection models for high-dimensional input data. A flexible selection prior allows the incorporation of additional information, e.g. a second data source including all sources of variation. 

In the logisticVS() function this is implemented for a logistic regression model.
}
\value{
outdir 
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
res <- logisticVS(X, Y, b=rep(0,p), v=rep(1,p), k=3, d=rep(1/p,p), block=diag(1,p), 
       MCMC=B, thinn=1, seed=1234, outdir=paste(getwd(), "example", sep="/"))

\dontrun{
require(Matrix)
tmp <- read.table(paste(res, "GAM.txt", sep="/"))
GAM <- as.matrix(sparseMatrix(i = tmp[,1], j = tmp[,2], x = tmp[,3]))
image(y=1:p, x=1:B, z=t(GAM), col=c("white","black"), ylab="variables", 
      xlab="MCMC iteration", main="Trace plot for GAM")
}
}
