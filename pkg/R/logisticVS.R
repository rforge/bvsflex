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
#v0.46: MZ, 2015-06-17 (change to mu and phi parametrisation of prior for pi)

logisticVS <- function(X, Y, b, h0, g=1, block=NULL, 
                       mu=NULL, phi=NULL,
                       m0=NULL, g0=NULL,
                       aPhi=1, bPhi=0.1,
                       aG=NULL, bG=NULL,
                       MCMC, thinn=1, seed=1234, outdir=NULL, 
                       Piupdate=FALSE, MuPhiUpdate=FALSE, gupdate="none", 
                       wlsgprior=FALSE, sigmaG=1, sigmaMu=1, sigmaPhi=1){
 
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
  
  if(is.na(MuPhiUpdate)) stop("'MuPhiUpdate' should be either FALSE or TRUE (no missing values allowed).\n") 
  MuPhiUpdate <- ifelse(MuPhiUpdate==TRUE, 1, 0);
  
  
  if (!is.null(phi) &length(phi)==1) {phi=rep(phi,p); warning("phi should be a vector of lenght p. Continue with phi=rep(phi,p).")}
  
  if(MuPhiUpdate != 1){
    if(is.null(mu)) mu <- rep(0.5, p); #this is the same as assuming aBeta[i]=1...
    if(is.null(phi)) phi <- rep(2, p); #...and bBeta[i]=1 for all i
  }else{
    if(is.null(m0)) m0 <- rep(0.5, p); #this is the same as assuming aBeta[i]=1...
    if(is.null(g0)) g0 <- rep(2, p);   #...and bBeta[i]=1 for all i
    if(aPhi <=0 | bPhi <= 0) stop("The parameters 'aPhi' and 'bPhi' have to be positive.\n")
    if(any(m0<=0) | any(m0>=1) | any(g0<=0)) stop("The parameters 'm0[i]' must be in (0,1) and 'g0[i]' must be positive.\n")
  }
  
  if(is.na(Piupdate)) stop("'Piupdate' should be either FALSE or TRUE (no missing values allowed).\n") 
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
  
  #Center X[,j] (j=1,...,p) and Y (because we don't have an intercept in the model):
  Xs <- scale(X, scale=FALSE);
  Ys <- scale(Y, scale=FALSE);

  cat("Starting MCMC...\n");
  .C("logisticVS", 
      X = as.double(Xs), Y = as.double(Ys), n = as.integer(n), p = as.integer(p),
      b = as.double(b), h0 = as.double(h0), g = as.double(g),
      mu = as.double(mu), phi = as.double(phi),
      m0 = as.double(m0), g0 = as.double(g0), sigmaMu = as.double(sigmaMu),
      aPhi = as.double(aPhi), bPhi = as.double(bPhi), sigmaPhi = as.double(sigmaPhi),
      aG = as.double(aG), bG = as.double(bG), sigmaG = as.double(sigmaG),
      block = as.integer(block), MCMC = as.integer(MCMC), thinn = as.integer(thinn),
      Piupdate = as.integer(Piupdate),
      MuPhiUpdate = as.integer(MuPhiUpdate),
      gupdate = as.integer(gupdate.int),
      PACKAGE = "bvsflex");   
  
  setwd(startdir);  
  return(outdir);
}

