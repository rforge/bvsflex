#Logistic Bayesian variable selection (with 'singleGibbs' block Gibbs sampler for updating GAM and beta), with flexible variable selection prior Pi
#v0.1: Manuela Zucknick, 2011-05-21, last modified 2013-12-12
#v0.2: Manuela Zucknick & Ana Corberan, 2014-05-05, last modified 2014-06-01
#v0.24/v0.25: Manuel Wiesenfarth, 2014-11-27, added to R-forge 2015-01-14
#v0.3: MZ, 2015-01-26 (add prior distributions for g)
#v0.31-v0.33: MZ, 2015-01-28 (add non-diagonal v), debugging 2015-02-04
#v0.35: MZ, 2015-02-20, replace v0 by h0 (prior precision matrix)
#v0.4: MZ, 2015-02-27, bug fix (was introduced in version 0.35), add wlsgprior = "weighted least squares g-prior" as option 
#      (then h is computed as 1/g * t(Xgam)%*%invLAMXgam in each iteration instead of using a prespecified h0. 
#       h0 is ignored in this case)
#v0.41-v0.45: MZ, MW, 2015-03-02 to 2015-03-10, bug fixes (in g_update), h in wlsgprior is now computed as 1/g * 1/n t(Xgam)%*%invLAMXgam

logisticVS <- function(X, Y, b, h0, g=1, block=NULL, 
                       aBeta=NULL, bBeta=NULL,
                       aG=NULL, bG=NULL,
                       MCMC, thinn=1, seed=1234, outdir=NULL, 
                       Piupdate=FALSE, gupdate="none", 
                       wlsgprior=FALSE, sigmaG=NULL){
 
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
  if(!(wlsgprior %in% c(TRUE,FALSE))){
    stop("'wlsgprior' is a binary, only values TRUE and FALSE are allowed.\n");
  }
  if(wlsgprior == FALSE){
    if(is.null(h0)) stop("Please provide h0.\n");
    if(is.vector(h0)){
      if(length(h0) != p) stop("If h0 is a vector, length of h0 has to be p = ncol(X).\n");
    }else{
      if(!is.matrix(h0) | ncol(h0) != p | nrow(h0) != p){
        stop("If h0 is not a vector, it has to be a matrix of dimension p x p, where p = ncol(X).\n");
      }
    }
    #if h0 is a vector, make a diagonal matrix out of it:
    if(is.vector(h0)) h0 <- diag(h0);
  }else{
    h0 <- 0.0; #if h0 is a scalar and equal to 0.0, in C code this is interpreted as wlsgprior==TRUE
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
    if(aG <0 | bG < 0) stop("The parameters 'aG' and 'bG' for gupdate='IG' have to be non-negative\n")
  }
  if(gupdate == "hyperg"){
    gupdate.int <- 2;
    
    #default: prior suggested by Liang et al. (2008), also see Shang and Li (2014)
    if(is.null(aG)) aG <- 1;
    if(is.null(bG)) bG <- 0.5;
    if(aG <=0 | bG <= 0) stop("The parameters 'aG' and 'bG' for gupdate='hyperg' have to be positive.\n")
  }
  if(is.null(sigmaG)) sigmaG <- 1;
  
  #Center X[,j] (j=1,...,p) and Y (because we don't have an intercept in the model):
  Xs <- scale(X, scale=FALSE);
  Ys <- scale(Y, scale=FALSE);

  cat("Starting MCMC...\n");
  .C("logisticVS", 
      X = as.double(Xs), Y = as.double(Ys), n = as.integer(n), p = as.integer(p),
      b = as.double(b), h0 = as.double(h0), g = as.double(g),
      aBeta=as.double(aBeta), bBeta=as.double(bBeta),
      aG=as.double(aG), bG=as.double(bG), sigmaG=as.double(sigmaG),
      block = as.integer(block), MCMC = as.integer(MCMC), 
      thinn = as.integer(thinn),
      Piupdate = as.integer(Piupdate),
      gupdate = as.integer(gupdate.int),
      PACKAGE = "bvsflex");   
  
  setwd(startdir);  
  return(outdir);
}

