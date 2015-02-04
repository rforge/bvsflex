#Logistic Bayesian variable selection (with 'singleGibbs' block Gibbs sampler for updating GAM and beta), with flexible variable selection prior Pi
#v0.1: Manuela Zucknick, 2011-05-21, last modified 2013-12-12
#v0.2: Manuela Zucknick & Ana Corberan, 2014-05-05, last modified 2014-06-01
#v0.24/v0.25: Manuel Wiesenfarth, 2014-11-27, added to R-forge 2015-01-14
#v0.3: Manuela Zucknick, 2015-01-26 (add prior distributions for g)
#v0.31/v0.32: Manuela Zucknick, 2015-01-28 (add non-diagonal v), debugging 2015-02-04

logisticVS <- function(X, Y, b, v, g=1, block=NULL, 
                       aBeta=NULL, bBeta=NULL,
                       aG=NULL, bG=NULL,
                       MCMC, thinn=1, seed=1234, outdir=NULL, 
                       Piupdate=FALSE, gupdate="none"){
 
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
  
  if(length(Y) != n) stop("Length of Y has to be n = nrow(X).\n");
  if(all(sort(unique(Y)) != c(0,1))) stop("The class vector Y should always contain both 0's and 1's (and only these values).\n");
  if(length(b) != p) stop("Length of b has to be p = ncol(X).\n");
  if(is.vector(v)){
    if(length(v) != p) stop("If v is a vector, length of v has to be p = ncol(X).\n");
  }else{
    if(!is.matrix(v) | ncol(v) != p | nrow(v) != p){
      stop("If v is not a vector, it has to be a matrix of dimension p x p, where p = ncol(X).\n");
    }
  }
  if(length(g) != 1) stop("g should be a scalar.\n");
  if(ncol(block) != p | nrow(block) != p) stop("Dimension of 'block' has to be p x p, where p = ncol(X).\n");
  
  if(!is.wholenumber(MCMC) | MCMC<=0) stop("'MCMC' is the number of MCMC iterations and hence has to be a positive integer.\n")	
  if(!is.wholenumber(thinn) | MCMC<=thinn) stop("Every 'thinn'-ed MCMC iteration will be kept (i.e. if 'thinn'=1 then all iterations are kept, if 'thinn'=10 then only every tenth iteration is returned for posterior inference. Consequently, 'thinn' has to be a positive integer.\n")
  
  if(is.null(aBeta)) aBeta <- rep(1, p);
  if(is.null(bBeta)) bBeta <- rep(1, p);
  
  if(is.na(Piupdate)) stop("'Pipudate' should be either FALSE or TRUE (no missing values allowed).\n") 
  Piupdate <- ifelse(Piupdate==TRUE, 1, 0);
  
  if(gupdate!="none" & gupdate!="IG" & gupdate!="hyperg") stop("'gupdate' should be either 'none', 'IG', or 'hyperg'.\n")
  gupdate.int <- 0;
  
  if(gupdate == "IG"){
    gupdate.int <- 1;
    
    #default: Zellner-Siow prior (Zellner & Siow 1980)
    if(is.null(aG)) aG <- 0.5;
    if(is.null(bG)) bG <- 0.5;
    if(aG <=0 | bG <= 0) stop("The parameters 'aG' and 'bG' for gupdate='IG' have to be positive.\n")
  }
  if(gupdate == "hyperg"){
    gupdate.int <- 2;
    
    if(is.null(aG)) aG <- 3; #default: prior suggested by Liang et al. (2008)
    if(aG <= 2) stop("The parameter 'aG' for gupdate='hyperg' has to be larger than 2.\n")
  }
  
  #if v is a vector, make a diagonal matrix out of it:
  if(is.vector(v)) v <- diag(v);
  
  #Center X[,j] (j=1,...,p) and Y (because we don't have an intercept in the model):
  Xs <- scale(X, scale=FALSE);
  Ys <- scale(Y, scale=FALSE);

  cat("Starting MCMC...\n");
  .C("logisticVS", 
      X = as.double(Xs), Y = as.double(Ys), n = as.integer(n), p = as.integer(p),
      b = as.double(b), v = as.double(v), g = as.double(g),
      aBeta=as.double(aBeta), bBeta=as.double(bBeta),
      aG=as.double(aG), bG=as.double(bG),
      block = as.integer(block), MCMC = as.integer(MCMC), 
      thinn = as.integer(thinn),
      Piupdate = as.integer(Piupdate),
      gupdate = as.integer(gupdate.int),
      PACKAGE = "bvsflex");   
  
  setwd(startdir);  
  return(outdir);
}

