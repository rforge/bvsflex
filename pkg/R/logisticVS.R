#Logistic Bayesian variable selection (with 'singleGibbs' block Gibbs sampler for updating GAM and beta), with flexible variable selection prior Pi, allow input of (empirical) distribution for Pi in form of a matrix 
#Manuela Zucknick & Ana Corberan, 21-05-2011, last modified 05-05-2014

logisticVS <- function(X, Y, b, v, k, d, block=NULL, MCMC, thinn=1, seed=1234, outdir=NULL){
# k = a scalar specifying the overall expected number of variables in a model (E(Pi[i] = k/p)
# d = a vector of length p containing the distances for variables i=1,...,p
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  set.seed(seed);
  
  startdir <- getwd();
  if(is.null(outdir)){
    outdir <- startdir;
  }else{
    dir.create(outdir);
  }
  setwd(outdir);  
  
  n = nrow(X);
  p = ncol(X);
  
  if(is.null(block)) block <- diag(1,p);
  
  if(length(Y) != n) stop("Length of Y has to be the same as nrow(X) = n.\n");
  if(all(sort(unique(Y)) != c(0,1))) stop("The class vector Y should always contain both 0's and 1's (and only these values).\n");
  if(length(b) != p) stop("Length of b has to be the same as ncol(X) = p.\n");
  if(length(v) != p) stop("Length of v has to be the same as ncol(X) = p.\n");
  if(ncol(block) != p | nrow(block) != p) stop("Dimension of 'block' has to be p x p, where p = ncol(X).\n");
  
  if(!is.wholenumber(MCMC) | MCMC<=0) stop("'MCMC' is the number of MCMC iterations and hence has to be a positive integer.\n")	
  if(!is.wholenumber(thinn) | MCMC<=thinn) stop("Every 'thinn'-ed MCMC iteration will be kept (i.e. if 'thinn'=1 then all iterations are kept, if 'thinn'=10 then only every tenth iteration is returned for posterior inference. Consequently, 'thinn' has to be a positive integer.\n")
  
  if(any(d < 0)) stop("'d' is a distance and hence no value can be negative.\n")
  if(!is.wholenumber(k) | k<=0 | k>p) stop("'k' has to be a positive integer smaller than 'p'.\n")
  
  cat("Starting MCMC...\n");
  .C("logisticVS", 
      X = as.double(X), Y = as.double(Y), n = as.integer(n), p = as.integer(p),
      b = as.double(b), v = as.double(v),
      k = as.integer(k), d = as.double(d),
      block = as.integer(block), MCMC = as.integer(MCMC), 
      thinn = as.integer(thinn),
      PACKAGE = "bvsflex");   
   
  setwd(startdir);  
  return(outdir);
}

