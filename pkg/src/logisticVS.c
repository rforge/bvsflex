#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Print.h>

void diagonalize(double *y, int n, double *dy){
	//y is a vector of length n. dy is corresponding diagonal matrix of dim nxn.
	
	for (int i=0; i<(n*n); i++) {
		dy[i] = 0;
	}
	for (int i=0; i<n; i++) {
		dy[(n+1) * i] = y[i];
	}
}

void makeSymmetric(double *V, int p){
	//V is a matrix of dim pxp. This function mirrors the upper triangular part 
	//of the matrix at the diagonal resulting in a symmetric matrix.
	
	for (int i=1; i<p; i++) {
		for (int j=0; j<i; j++) {
			V[p*j + i] = V[p*i + j];
		}
	}
}

void matrixInverse(double *V, int p){
	//V is a symmetric positive definite matrix of dim pxp. This function computes the inverse
	//of V and returns it in place of V.
	const char Upper = 'U';
	int info;
	
	F77_NAME(dpotrf)(&Upper, &p, V, &p, &info);
	if(info != 0) error("Error in Lapack function dpotrf, while in function matrixInverse.\n");
	F77_NAME(dpotri)(&Upper, &p, V, &p, &info);	
	if(info != 0) error("Error in Lapack function dpotri, while in function matrixInverse.\n");
	
	makeSymmetric(V, p);
	//dpotrf/dpotri only deal with values in the upper triangular part of the matrix. 
	//makeSymmetric fills the rest.
}

void LowerTri_fillZeros(double *V, int p){
	//V is a matrix of dim pxp. This function keeps the lower triangular part 
	//of the matrix but fills the upper triangular (above diagonal) with zeros.
	
	for (int i=1; i<p; i++) {
		for (int j=0; j<i; j++) {
			V[p*i + j] = 0;
		}
	}
}

void UpperTri_fillZeros(double *V, int p){
	//V is a matrix of dim pxp. This function keeps the upper triangular part 
	//of the matrix but fills the lower triangular (below diagonal) with zeros.
	
	for (int i=1; i<p; i++) {
		for (int j=0; j<i; j++) {
			V[p*j + i] = 0;
		}
	}
}

void matrixRow(double *X, int n, int p, int *j, double *Xj){
	//j is assumed to be counting from 1, while within this function array indices are counted from 0.
	//X is matrix of dim nxp
	
	for (int i=0; i<p; i++) {
		Xj[i] = X[j[0] + n*i];
    }
}

void matrixSubset(double *X, int *GAM, int n, int p, double *Xgam){
	//GAM is binary indicator vector of length p, indicating which columns of X to keep in Xgam.
	//X is matrix of dim nxp.
	
	int igam = 0;	//for counting through elements of GAM that are 1
	for (int i=0; i<p; i++) {
		if (GAM[i] == 1) {
			for (int j=0; j<n; j++){
				Xgam[j + n*igam] = X[j + n*i];
			}
			igam = igam + 1;
		}
	}
}

void matrixSubsetRows(double *X, int *GAM, int p, int n, double *Xgam){
  //GAM is binary indicator vector of length p, indicating which rows of X to keep in Xgam.
	//X is matrix of dim pxn.
		  
  int pgam = 0;
  for (int k = 0; k < p; k++) {
		pgam += GAM[k];
	}
    
  for (int i=0; i<n; i++){
    int jgam = 0;  //for counting through elements of GAM that are 1 (start again with each new column)
	  for (int j=0; j<p; j++) {
		  if (GAM[j] == 1) {
				Xgam[jgam + pgam*i] = X[j + p*i];
        jgam = jgam + 1;
			}
		}
	}
}

void vectorSubset(double *y, int *GAM, int p, double *ygam){
	//GAM is binary indicator vector of length p, indicating which elements of y to keep in ygam.
	//y is vector of length p.
	
	int igam = 0;	//for counting through those elements of GAM that are 1
	for (int i=0; i<p; i++) {
		if (GAM[i] == 1) {
			ygam[igam] = y[i];
			igam = igam + 1;
		}
	}
}

double dmvnorm(double *Z, int n, double *Zmean, double *Zprec, char Log, int pgam){
	//Returns value of multivariate normal density belonging to Z (vector of length n)
	//Zmean is the mean vector, Zprec the precision matrix (=inverse covariance matrix)
	//If Log='L' then return log density, else return density itself.
	//(For internal notation, see dmvnorm() in R package mvtnorm.)
	
	const int one = 1;
	const char Upper = 'U', Lower = 'L', Trans = 'T', NoDiag = 'N';
	
	double Zcentered[n];
	for (int j = 0; j < n; j++) {
		Zcentered[j] = Z[j] - Zmean[j];
	}
  
  //input Zprec, output Zchol (lower Cholesky factorisation matrix)
  double Zchol[(n * n)];
	for (int i = 0; i < (n * n); i++) {
		Zchol[i] = Zprec[i];
	}
  
	int info;
	F77_NAME(dpotrf)(&Upper, &n, Zchol, &n, &info);
	if(info != 0){
    error("Error in Lapack function dpotrf, while computing Zchol.\n");
	}
	
	//input Zcentered, output Zchol %*% Zcentered
	F77_NAME(dtrmv)(&Lower, &Trans, &NoDiag, &n, Zchol, &n, Zcentered, &one);
	
	double distval;
	distval = F77_NAME(ddot)(&n, Zcentered, &one, Zcentered, &one);
	
	//log determinant of Zchol = 0.5 * log determinant of Zprec
	double logdet2 = 0;
	for (int i = 0; i < n; i++) {
		logdet2 += log(Zchol[i * (n+1)]);
	}
	
	//M_LN_SQRT_2PI = 0.5 * log(2 * PI)
	double retval = -n * M_LN_SQRT_2PI + logdet2 - 0.5*distval;
    if (Log != 'L') retval = exp(retval);	
	
	return(retval);
}

void Vgam_Bgam_update(double *Vgam, double *Bgam, int pgam, int n,
					  double *Xgam, double *v0gam, double g, double *bgam, double *LAM, double *Z){
	
	const int one = 1;
	const double oned = 1.0, zerod = 0.0, moned = -1.0;
	const char Upper = 'U', NoTrans = 'N', Trans = 'T';
	
	double sqrtinvLAM[n];
	for (int j = 0; j < n; j++) {
		sqrtinvLAM[j] = 1.0/sqrt(LAM[j]);
	}
	double dsqrtinvLAM[(n * n)];
	diagonalize(sqrtinvLAM, n, dsqrtinvLAM);
	
	double dLAM[(n * n)];
	diagonalize(LAM, n, dLAM);
  
  double invvgam[(pgam * pgam)];
  for (int i=0; i<(pgam * pgam); i++) {
    invvgam[i] = g * v0gam[i]; 
	}
  matrixInverse(invvgam, pgam);

	if(pgam > n){
		//inversion of nxn matrix: 
		//Vtmp = vgam - vgam%*%t(Xgam)%*%solve(LAM + Xgam%*%vgam%*%t(Xgam))%*%Xgam%*%vgam;
    
		double Xvgam[(n * pgam)]; //Xgam %*% (g * v0gam)
		F77_NAME(dgemm)(&NoTrans, &NoTrans, &n, &pgam, &pgam, &g, Xgam, &n, v0gam, &pgam, &zerod, Xvgam, &n);
 
		double LAM_XvXgam[(n * n)]; //input = dLAM, output = LAMXvXgam
		diagonalize(LAM, n, LAM_XvXgam);
		F77_NAME(dgemm)(&NoTrans, &Trans, &n, &n, &pgam, &oned, Xgam, &n, Xvgam, &n, &oned, LAM_XvXgam, &n);
		
		int info;
		double invLAMXvX_Xvgam[(n * pgam)]; //input = Xvgam, output = invLAMXvX_Xvgam
		for (int i=0; i<(n * pgam); i++) {
			invLAMXvX_Xvgam[i] = Xvgam[i];
		}
		F77_NAME(dposv)(&Upper, &n, &pgam, LAM_XvXgam, &n, invLAMXvX_Xvgam, &n, &info);
		if(info != 0) error("Error in Lapack function dposv, while computing invLAMXvX. Error %d\n",info);
    
  	for (int i=0; i<(pgam * pgam); i++) {
			Vgam[i] = g * v0gam[i]; //input = vgam, output = Vtmp
		}
		F77_NAME(dgemm)(&Trans, &NoTrans, &pgam, &pgam, &n, &moned, Xvgam, &n, invLAMXvX_Xvgam, &n, &oned, Vgam, &pgam);
	
	}else{
		//inversion of dxd matrix (d=dim(GAM)):
		//Vtmp = solve(t(Xgam)%*%invLAM%*%Xgam + solve(vgam));    
		
		double sqrtinvLAMXgam[(n * pgam)];
		F77_NAME(dgemm)(&NoTrans, &NoTrans, &n, &pgam, &n, &oned, dsqrtinvLAM, &n, Xgam, &n, &zerod, sqrtinvLAMXgam, &n);

		// 1) input: invvgam, output: t(X)%*%invLAM%*%X + invvgam, 
    for (int i=0; i<(pgam * pgam); i++) {
  	  Vgam[i] = invvgam[i];
	  }		    
		F77_NAME(dgemm)(&Trans, &NoTrans, &pgam, &pgam, &n, &oned, sqrtinvLAMXgam, &n, sqrtinvLAMXgam, &n, &oned, Vgam, &pgam);		
    
    // 2) compute inverse (first cholesky decomposition, then inverse from that).
		matrixInverse(Vgam, pgam);
	}
  
	//Bgam = Vgam %*% (invvgam%*%bgam + tXgam%*%invLAM%*%Z);
	double invvbgam_XinvLAMZ[pgam];
  F77_NAME(dgemv)(&NoTrans, &pgam, &pgam, &oned, invvgam, &pgam, bgam, &one, &zerod, invvbgam_XinvLAMZ, &one);  
  // up to here "invvbgam_XinvLAMZ" is actually only invvgam%*%bgam
  
	double invLAMZ[n];
	for (int j=0; j<n; j++) {
		invLAMZ[j] = (1.0/LAM[j]) * Z[j];
	}
	F77_NAME(dgemv)(&Trans, &n, &pgam, &oned, Xgam, &n, invLAMZ, &one, &oned, invvbgam_XinvLAMZ, &one); 
	// here invvbgam_XinvLAMZ is updated to include everything that is multiplied with Vgam resulting in Bgam
	
	F77_NAME(dgemv)(&Trans, &pgam, &pgam, &oned, Vgam, &pgam, invvbgam_XinvLAMZ, &one, &zerod, Bgam, &one);	
}

double Pgam_update(int pgam, int n, double *Xgam, double *bgam, 
				   double *Z, double *invLAM, double *dinvLAM, double *Vgam, char NoLog){
	
	const int one = 1;
	const double oned = 1.0, zerod = 0.0, moned = -1.0;
	const char Lower = 'L', NoTrans = 'N', Trans = 'T';
	
	//Zprecgam = invLAM - invLAM%*%Xgam%*%Vgam%*%t(Xgam)%*%invLAM;
	double Lgam[(pgam * pgam)];
	for (int i=0; i<(pgam * pgam); i++) {
		Lgam[i] = Vgam[i]; //input Vgam, output Lgam (lower Cholesky factorisation matrix)
	}
  
	int info;
	F77_NAME(dpotrf)(&Lower, &pgam, Lgam, &pgam, &info);
	if(info != 0) error("Error in Lapack function dpotrf, while computing Lgam. info = %d\n", info);
	LowerTri_fillZeros(Lgam, pgam);
	
	double XgamL[(n * pgam)];
	F77_NAME(dgemm)(&NoTrans, &NoTrans, &n, &pgam, &pgam, &oned, Xgam, &n, Lgam, &pgam, &zerod, XgamL, &n);
	
	double invLAMXgamL[(n * pgam)];
	F77_NAME(dgemm)(&NoTrans, &NoTrans, &n, &pgam, &n, &oned, dinvLAM, &n, XgamL, &n, &zerod, invLAMXgamL, &n);
	
	//input dinvLAM, output Zprecgam
	double Zprecgam[(n * n)];
	diagonalize(invLAM, n, Zprecgam);
	F77_NAME(dgemm)(&NoTrans, &Trans, &n, &n, &pgam, &moned, invLAMXgamL, &n, invLAMXgamL, &n, &oned, Zprecgam, &n);
  
	//Pgam = dmvnorm(Z, mean=as.vector(Xgam%*%bgam), sigma=chol2inv(chol(Zprecgam)), log=FALSE);
	double Zmeangam[n];
	F77_NAME(dgemv)(&NoTrans, &n, &pgam, &oned, Xgam, &n, bgam, &one, &zerod, Zmeangam, &one);
  
	double Pgam = dmvnorm(Z, n, Zmeangam, Zprecgam, NoLog, pgam);
	return(Pgam);
}

void betaGAM_update(int n, int p, double *X, 
					double *b, double *v0, double g, double *Pi,
					double *beta, int *GAM, double *LAM, double *Z,  
					int *block, double T, int *numneigh, int *select,
					double *logprior, double *logpost)
{	
	//const double K1 = 0, K2 = 0; other options are not yet implemented
	const int one = 1;
	const double oned = 1.0;
	const char Upper = 'U', Trans = 'T', NoLog = 'N', Log = 'L';
  
	int pgam = 0;
	for (int k = 0; k < p; k++) {
		pgam += GAM[k];
	}
	
	//multiply T to LAM
	if(T != 1.0){
		double T1 = T - 1.0;
		F77_NAME(daxpy)(&n, &T1, LAM, &one, LAM, &one); 
		//trick&Trans: use fact that (y = alpha y) is the same as (y = (alpha-1)y + y).
	}
		
	double invLAM[n];
	for (int j = 0; j < n; j++) {
		invLAM[j] = 1.0/LAM[j];
	}
	double dinvLAM[(n * n)];
	diagonalize(invLAM, n, dinvLAM);
	
	//(0) Compute values corresponding to current GAM 
	//(the final Xgam after passing through betaGAM_update will be returned by reference for use in LAMZ_update)
	
	double Xgam[(pgam * n)];
	matrixSubset(X, GAM, n, p, Xgam);
  
  double v0tmpgam[(p * pgam)], v0gam[(pgam * pgam)];
  matrixSubset(v0, GAM, p, p, v0tmpgam); //matrixSubset subsets wrt. columns
  matrixSubsetRows(v0tmpgam, GAM, p, pgam, v0gam); //matrixSubsetRows subsets wrt. rows
  
	double bgam[pgam];
	vectorSubset(b, GAM, p, bgam);
	
	double Vgam[pgam * pgam];
	double Bgam[pgam];
	Vgam_Bgam_update(Vgam, Bgam, pgam, n, Xgam, v0gam, g, bgam, LAM, Z);
  
	double Pgam = Pgam_update(pgam, n, Xgam, bgam, Z, invLAM, dinvLAM, Vgam, NoLog);
	
	//(1) Proposal for GAM
	//(select gamma_i at random, find 'mates' and for all of these,
	//draw from p(GAM_i|Z) [Lee et al. (2003)]: draw by computing density
	//for both possible values (1 and 0). As one value is available from 
	//GAMstar in previous iteration, only compute the other one = GAMtry)
	
	//find variables to update together with GAM(select) (= neighbours)
  select[0] = p * unif_rand(); //random sample from 0:(p-1)
	int neighbour[p];
	numneigh[0] = 0;
	for (int i = (select[0] * p); i < ((select[0]+1) * p); i++) {
		if (block[i] != 0) {
			neighbour[numneigh[0]] = i % p; //modulo p
			numneigh[0] = numneigh[0] + 1;
		}
	}
	//for 'neighbour' to include 'select', the diagonal of 'block' should be !=0.
	
	int GAMstar[p], GAMtry[p];
	for (int k = 0; k < p; k++) {
		GAMstar[k] = GAM[k];
	}
	
	double Pstar = Pgam;
	int i, ptry, pstar;
	for (int iter = 0; iter < numneigh[0]; iter++) {
	
		i = neighbour[iter];
		
		//in each iteration compare GAMtry to currently accepted gamma:		
		for (int k = 0; k < p; k++) {
			GAMtry[k] = GAMstar[k];
		}
		GAMtry[i] = (GAMstar[i] == 0 ? 1 : 0);
		ptry = 0;
		for (int k = 0; k < p; k++) ptry += GAMtry[k];
		
		if(ptry > 0){	
			//if GAMtry is all zeros then keep GAMstar as it is
			
			double Xtry[(ptry * n)];
			matrixSubset(X, GAMtry, n, p, Xtry);
			
      double v0tmptry[(p * ptry)], v0try[(ptry * ptry)];
      matrixSubset(v0, GAMtry, p, p, v0tmptry); //matrixSubset subsets wrt. columns
      matrixSubsetRows(v0tmptry, GAMtry, p, ptry, v0try); //matrixSubsetRows subsets wrt. rows
      
			double btry[ptry];
			vectorSubset(b, GAMtry, p, btry);
			
			double Vtry[(ptry * ptry)];
			double Btry[ptry];
			Vgam_Bgam_update(Vtry, Btry, ptry, n, Xtry, v0try, g, btry, LAM, Z);
			
			double Ptry = Pgam_update(ptry, n, Xtry, btry, Z, invLAM, dinvLAM, Vtry, NoLog);
      
			//only multiply by p(GAM(i)) because all the rest cancels out anyway
			double probTRY = Ptry * pow(Pi[i], GAMtry[i]) * pow((1-Pi[i]), (1-GAMtry[i]));
			double probSTAR = Pstar * pow(Pi[i], GAMstar[i]) * pow((1-Pi[i]), (1-GAMstar[i]));
      
			if (probTRY==0 && probSTAR==0){  //to reduce computational problems with precision
				if (GAMtry[i]==0){           //-> favour more sparse GAM
					GAMstar[i] = GAMtry[i];  //else: GAMstar(i) = GAM(i)
					Pstar = Ptry;           
				}
			}else{
				double U = runif(0.0, 1.0);
				if (U <= probTRY/(probSTAR + probTRY)){
					GAMstar[i] = GAMtry[i];  //else: GAMstar(i) = GAM(i)
					Pstar = Ptry;			 //else: Pstar stays the same
				}
			}
		}
	}
	
	//(2) compute values corresponding to proposed GAMstar = new GAM
	//(3) GAM = GAMstar: acceptance/rejection step is not needed since alpha = 1 always
	
	pstar = 0;
	for (int k = 0; k < p; k++) pstar += GAMstar[k];
	
	double Xstar[(pstar * n)];
	matrixSubset(X, GAMstar, n, p, Xstar);

	double v0tmpstar[(p * pstar)], v0star[(pstar * pstar)];
	matrixSubset(v0, GAMstar, p, p, v0tmpstar); //matrixSubset subsets wrt. columns
  matrixSubsetRows(v0tmpstar, GAMstar, p, pstar, v0star); //matrixSubsetRows subsets wrt. rows
	
	double bstar[pstar];
	vectorSubset(b, GAMstar, p, bstar);
	
	double Vstar[(pstar * pstar)];
	double Bstar[pstar];
	Vgam_Bgam_update(Vstar, Bstar, pstar, n, Xstar, v0star, g, bstar, LAM, Z);
	
	double logPgam = Pgam_update(pstar, n, Xstar, bstar, Z, invLAM, dinvLAM, Vstar, Log);
	
	double Lstar[(pstar * pstar)];
	for (int i=0; i<(pstar * pstar); i++) {
		Lstar[i] = Vstar[i]; //input Vgam, output Lgam (lower Cholesky factorisation matrix)
	}
	int info;
	F77_NAME(dpotrf)(&Upper, &pstar, Lstar, &pstar, &info);
	if(info != 0) error("Error in Lapack function dpotrf, while computing Lgam. info = %d\n", info);
	UpperTri_fillZeros(Lstar, pstar);
	
	double Tvec[pstar];
	for (int i = 0; i < pstar; i++) {
		Tvec[i] = rnorm(0.0, 1.0);
	}
	
	double betastar[pstar]; 
	//Create betastar (as local copy of beta) as array with exactly pstar elements (for computations within this function).
	//Contrary to that beta has *p elements (for flexibility when transfering values outside of this function, because I can't change the dimensionality of an array.) 
	for (int i = 0; i < pstar; i++) {
		betastar[i] = Bstar[i]; //input to dgemv is Bgam=Bstar, output will be betastar
	}
	F77_NAME(dgemv)(&Trans, &pstar, &pstar, &oned, Lstar, &pstar, Tvec, &one, &oned, betastar, &one);	
	
	int cnt = 0;
	for (int i = 0; i < p; i++) {
		GAM[i] = GAMstar[i];
		if(GAM[i] == 1) {
			beta[i] = betastar[cnt]; //fill in the elements of beta with new betastar values
			cnt++;
		}
	}
	
	//log probabilities log(P(beta,GAM|Z,LAM)) and log(P(beta,GAM))
	double logpriorGAM = 0.0;
	for (int i = 0; i < p; i++) {
		logpriorGAM += log(Pi[i]) * GAM[i] + log(1 - Pi[i]) * (1 - GAM[i]);
	}	
	
	double logpriorbeta = 0.0, logpostbeta = 0.0; //P(beta=empty|GAM,Z,LAM) = 1 (if pgam = 0)
	if (pstar > 0){
		double Vprec[(pstar * pstar)], vprec[(pstar * pstar)];
		for (int i = 0; i < (pstar * pstar); i++) {
			Vprec[i] = Vstar[i];
		}
		matrixInverse(Vprec, pstar);  //input Vstar, output inverse(Vstar)
    
    for (int i = 0; i < (pstar * pstar); i++) {
  		vprec[i] = g * v0star[i];
		}
		matrixInverse(vprec, pstar);  //input vstar = g * v0star, output inverse(vstar)
		
		logpostbeta = dmvnorm(betastar, pstar, Bstar, Vprec, Log, pstar);
		logpriorbeta = dmvnorm(betastar, pstar, bstar, vprec, Log, pstar);
	}
  
	logprior[0] = logpriorGAM + logpriorbeta;
	logpost[0] = logPgam + logpostbeta;  //logpriorGAM cancels out anyway
}

//////////////////////////////////////////////////////////////////////////////////////////////

double lefttrunclogis(double mu, double std, double lower)
{
	double lowerProb = plogis((lower - mu)/std, 0.0, 1.0, 1, 0);
	//x, location, scale, lower_tail, give_log
	
	double ru = runif(0.0, 1.0);
	double u = lowerProb + (1.0 - lowerProb) * ru;
	
	double ql = qlogis(u - DBL_EPSILON, 0.0, 1.0, 1, 0); 
	//p, location, scale, lower_tail, log_p
	double f = mu + ql * std;
	
	return(f);
}

double righttrunclogis(double mu, double std, double upper)
{
	double upperProb = plogis((upper - mu)/ std, 0.0, 1.0, 1, 0); 
	//x, location, scale, lower_tail, give_log
	
	double ru = runif(0.0, 1.0);
	double u = upperProb * ru;
	
	double ql = qlogis(u, 0.0, 1.0, 1, 0);
	//p, location, scale, lower_tail, log_p
	
	double f = mu + ql * std;
	
	return(f);
}

double KS_pdf(double x)
{	
	double pdf = 0;
	
	double x2 = x * x;
	double pdf_tmp = x * exp(-2 * x2);
	
	int n = 1;
	int n2 = 1;
	while (fabs(pdf_tmp) > DBL_EPSILON){
		pdf = pdf + pdf_tmp;
		n = n + 1;
		n2 = n * n;
		pdf_tmp = pow(-1, n+1) * n2 * x * exp(-2 * n2 * x2);
	}
	pdf = 8 * pdf;
	
	return(pdf);
}

int local_right(double U, double LAM, double T)
{
	// LOCAL_RIGHT Help function for function LAMBDADIST 
	double lH = 0.5 * LAM * (T-1) - log(T);
	double lU = log(U);
	double Z = 1.0;
	double X = exp(-0.5 * LAM);
	int j = 0;
	int j2 = 0;
	int STOP = 0;
	int OK = 0;
	while(!STOP){
		j = j + 1;
		j2 = (j + 1) * (j + 1);
		Z = Z - j2 * pow(X, (j2 - 1));
		if(lH + log(Z) > lU){
			STOP = 1;
			OK = 1;
			break;
		}
		j = j + 1;
		j2 = (j + 1) * (j + 1);
		Z = Z + j2 * pow(X, (j2 - 1));
		if(lH + log(Z) < lU){
			STOP = 1;
			OK = 0;
			break;
		}
	}
    return(OK);
}

int local_left(double U, double LAM, double T)
{
	// LOCAL_LEFT Help function for function LAMBDADIST
	double PI2 = M_PI * M_PI;
	double lH = 0.5 * log(2.0) + 2.5 * log(M_PI) - 2.5 * log(LAM) - PI2/(2.0 * LAM) + 0.5 * LAM * T - log(T);
	double lU = log(U);
	double X = exp(-PI2/(2 * LAM));
	double K = LAM/PI2;
	double Z = 1.0;
	int j = 0;
	int j2 = 0;
	int STOP = 0;
	int OK = 0;
	while(!STOP){
		j = j + 1;
		Z = Z - K * pow(X, (j * j - 1));
		if (lH + log(Z) > lU){
			STOP = 1;
			OK = 1;
			break;
		}
		j = j + 1;
		j2 = (j + 1) * (j + 1);
		Z = Z + j2 * pow(X, (j2 - 1));
		if (lH + log(Z) < lU){
			STOP = 1;
			OK = 0;
			break;
		}
	}
	return(OK);
}

double lambdadist_temp(double r, double T)
{
	/* LAMBDADIST_TEMP Sampling from tempered full conditional distribution of LAM
	 # in Bayesian auxiliary variable model for logistic regression
	 # using joint update to {z,beta}
	 # Tempering: p(LAM|Z,beta,Xgam) = C * N(Xgam*beta, T*LAM) * p(LAM) 
	 # Based on pseudocode A4 by Chris Holmes and Leonhard Held (2006)
	 #
	 # Reference: 
	 # Holmes, C. and Held, L. (2006). Bayesian Auxiliary Variable Models for 
	 # Binary and Multinomial Regression, Bayesian Analysis 1:145-168
	 */
	
	// Global constants that are needed often,
	// compute only once here to save computing time
	double Tr = (1.0/T) * r;
	
	double LAMlocal = 0;
	double Y;
	double U;
	int OK = 0;
	int iter = 0;
	
	while(!OK || iter<100){
		Y = runif(0.0, 1.0);
		if(Y != 0){
			Y = Y * Y;
			Y = 1 + (Y - sqrt(Y*(4*r + Y)))/(2*r);
			U = runif(0.0, 1.0);
			if(U <= 1.0/(1.0 + Y)){
				if(Y == 0){
					Y = DBL_EPSILON;
				}
				LAMlocal = Tr/Y;
			}else{
				LAMlocal = Tr*Y;
			}
			// Now, LAM = LAMu/T (LAMu = original LAM) where
			// f_LAMu(LAMu) = GIG(0.5, 1, r^2) = rejection sampling distribution
			
			U = runif(0.0, 1.0);
			if(LAMlocal > 4/3){
				OK = local_right(U, LAMlocal, T);
			}else{
				OK = local_left(U, LAMlocal, T);
			}
		}
		iter = iter + 1;
	}
	return(LAMlocal);
}

void LAMz_update(int n, int p, double *X, double *Y,
				 double *beta, int *GAM, double *LAM, double *Z,
				 double T, double *logprob)
{
	/* LAMZ_UPDATE Gibbs sampling update for LAM and Z in Bayesian auxiliary 
	 # variable model for logistic regression with variable selection
	 # Help function for LOGISTICVS
	 */
	
	const int one = 1;
	const double oned = 1.0, zerod = 0.0;
	const char NoTrans = 'N';
	
	int pgam = 0;
	for (int k = 0; k < p; k++) pgam += GAM[k];
	
	double Xgam[(n * pgam)];
	matrixSubset(X, GAM, n, p, Xgam);
	
	double betagam[pgam];  //only the first pgam elements of beta are relevant.
	for (int k = 0; k < pgam; k++) betagam[k] = beta[k];
	
	double m[n]; //m = Xgam %*% beta
	F77_CALL(dgemv)(&NoTrans, &n, &pgam, &oned, Xgam, &n, betagam, &one, &zerod, m, &one);
	
	double logZ = 0.0;
	double logpriorLAM = 0.0;
	for(int j = 0; j < n; j++)
	{
		if(Y[j] > 0){
			Z[j] = lefttrunclogis(m[j], sqrt(T), 0);
		}else{
			Z[j] = righttrunclogis(m[j], sqrt(T), 0);
		}
		
		double R = fabs(Z[j] - m[j]);
    
		if(isinf(R) || isnan(R) || R>1e10){
			error("Object R is (close to) infinite or NaN.\n");
		}else{
			LAM[j] = lambdadist_temp(R, T);
			logpriorLAM -= 0.5*log(16* LAM[j]) + log(KS_pdf(0.5*sqrt(LAM[j])));
			logZ += dnorm(Z[j], m[j], sqrt(T * LAM[j]), 1);
		}
	}
  
	logprob[0] = logZ + logpriorLAM;
}


//////////////////////////////////////////////////////////////////////////////////////////////

void Pi_update(double *Pi, int *GAM, double *aBeta, double *bBeta, int p){
   //Direct sampling:
  for (int i = 0; i < p; i++) {
    Pi[i] = rbeta(aBeta[i]+GAM[i], bBeta[i]+1.0-GAM[i]);
  }
}

// class of inverse gamma priors for g:
void g_update_IG(double *g, double aG, double bG){
   *g = 1.0 / rgamma(aG, 1.0 / bG); 
   //Note: Rmath.h rgamma uses shape and rate parameters, while 
   //rgamma() in R uses shape and scale parameters per default and 
   //aG and bG are meant as shape and scale parameters.
}
// hyper-g prior for g according to Liang et al. (2008) 
// (use g/(1+g) ~ Beta formulation):
void g_update_hyperg(double *g, double aG){
   double gratio = rbeta(1.0, 0.5 * aG - 1.0);
   *g = 1.0 / (1.0 / gratio - 1.0);
}


//////////////////////////////////////////////////////////////////////////////////////////////
//New version (with aBeta and bBeta and allowing for hyper-prior for g)
void logisticVS(double *X, double *Y, int *n, int *p,
  			double *b, double *v0, double *g, 
        double *aBeta, double *bBeta,
        double *aG, double *bG,
				int *block, int *MCMC, int *thins, 
        int *Piupdate, int *gupdate){
	
	GetRNGstate();	
	
	int thin[1];
	thin[0] = thins[0];
	double T = 1.0;
	
	int K2;
	int numneigh[1];
	int select[1];
	double logpriorGAMbeta[1];
	double logpostGAMbeta[1];
	double logLAMz[1];
  
  double Pi[*p];
  for (int i = 0; i < *p; i++){
    Pi[i] = (aBeta[i] / (aBeta[i] + bBeta[i]));
  }
  
 	int GAM[*p];
	double beta[*p];
	for (int i = 0; i < *p; i++) {
		double U = runif(0.0, 1.0);
		GAM[i] = (U < Pi[i] ? 1 : 0);
		beta[i] = 0.0;
	}
  
	double LAM[*n];
	double Z[*n];
	for (int j = 0; j < *n; j++) {
		double U2 = fabs(rnorm(0.0, 1.0));
		Z[j] = (Y[j] > 0.0 ? U2 : -U2);
		LAM[j] = 1;
	}
	
	//save in sparse matrix form (col1=rowID, col2=colID, col3=beta):
	FILE *fidbeta;
	fidbeta = fopen("beta.txt", "w");
	fprintf(fidbeta, "%d %d %.10e\n", *p, (*MCMC/ *thin), 0.0);    //store dimensions
	
	FILE *fidGAM;
	fidGAM = fopen("GAM.txt", "w");
	fprintf(fidGAM, "%d %d %d\n", *p, (*MCMC/ *thin), 0);    //store dimensions
	
  FILE *fidlogprob, *fidnumneigh, *fidselect, *fidg; //*fidpi;
    fidlogprob = fopen("logprobs.txt", "w");
    fidnumneigh = fopen("numneighs.txt", "w");
    fidselect = fopen("selects.txt", "w");
    fidg = fopen("g.txt", "w");
    //fidpi = fopen("pi.txt", "w"); 
  
	for (int K = 0; K < *MCMC; K++) {
    
		// Allow R interrupts; check every 10 iterations
		if (!(K % 10)) R_CheckUserInterrupt();
  	if (!((K+1) % 1000)) Rprintf("iteration %d\n", K+1);    
    
    if(*gupdate == 1){
      //class of inverse-gamma priors
      g_update_IG(g, *aG, *bG);
    }
    if(*gupdate == 2){
      //class of inverse-gamma priors
      g_update_hyperg(g, *aG);
    }
    
    if(*Piupdate == 1){
      Pi_update(Pi, GAM, aBeta, bBeta, *p);
    }
  
		betaGAM_update(*n, *p, X,
					   b, v0, *g, Pi,
					   beta, GAM, LAM, Z,
					   block, T, numneigh, select, logpriorGAMbeta, logpostGAMbeta);
    
		LAMz_update(*n, *p, X, Y,
					beta, GAM, LAM, Z,
					T, logLAMz);
    
		if (!(K % *thin)) {
			K2 = K / *thin;
		      
          //fprintf(fidpi, "%.4f ", Pi[0]);
          fprintf(fidg, "%.4f ", g[0]);
        	fprintf(fidlogprob, "%.4f ", logpriorGAMbeta[0] + logLAMz[0]);
        	fprintf(fidnumneigh, "%d ", numneigh[0]);
        	fprintf(fidselect, "%d ", select[0]);
		
			for (int i = 0; i < *p; i++) {
        
				if (GAM[i] == 1){
					fprintf(fidbeta, "%d %d %.10e\n", i+1, K2+1, beta[i]);
					fprintf(fidGAM, "%d %d %d\n", i+1, K2+1, 1);
				}
			}
		}
	}
	
  //fclose(fidpi);
  fclose(fidg);
	fclose(fidbeta);
	fclose(fidGAM);
	fclose(fidlogprob);
	fclose(fidnumneigh);
	fclose(fidselect);
	
	PutRNGstate();	
}



