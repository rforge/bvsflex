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
	//of the matrix but fills the upper triangular (below diagonal) with zeros.
	
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
	//GAM is binary indicator vector of length p, indicating which columns of X go keep in Xgam.
	//X is matrix of dim nxp.
	
	int j;
	int igam = 0;	//for counting through those elements of GAM that are 1
	for (int i=0; i<p; i++) {
		if (GAM[i] == 1) {
			for(j=0; j<n; j++){
				Xgam[j + n*igam] = X[j + n*i];
			}
			igam = igam + 1;
		}
	}
}

void vectorSubset(double *y, int *GAM, int p, double *ygam){
	//GAM is binary indicator vector of length p, indicating which elements of y go keep in ygam.
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
	const char Lower = 'L', Trans = 'T', NoDiag = 'N';
	
	double Zcentered[n];
	for (int j = 0; j < n; j++) {
		Zcentered[j] = Z[j] - Zmean[j];
	}
	
	//input Zprec, output Zchol (lower Cholesky factorisation matrix)
	double Zchol[(n * n)];
	for (int i = 0; i < (n * n); i++) {
		Zchol[i] = Zprec[i];
	}
	
	/*
	if(pgam > n){
	Rprintf("Zprecgam in dmvnorm:\n");
	for (int i=0; i<(n*n); i++) {
		Rprintf("%6.4f ", Zchol[i]);
		if ((i+1)%n == 0) {
			Rprintf("\n");
		}
	}
	error("dmvnorm: Okay up to here.\n");
	}
	*/
	
	int info;
	F77_NAME(dpotrf)(&Lower, &n, Zchol, &n, &info);
	if(info != 0) error("Error in Lapack function dpotrf, while computing Zchol.\n");
	
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
					  double *Xgam, double *vgam, double *bgam, double *LAM, double *Z){
	
	const int one = 1;
	const double oned = 1.0, zerod = 0.0, moned = -1.0;
	const char Upper = 'U', NoTrans = 'N', Trans = 'T';
	
	double sqrtinvLAM[n];
	for (int j = 0; j < n; j++) {
		sqrtinvLAM[j] = 1.0/sqrt(LAM[j]);
	}
	double dsqrtinvLAM[(n * n)];
	diagonalize(sqrtinvLAM, n, dsqrtinvLAM);
	
	// Up to here, LAM and vgam are still vectors. Convert to diagonal matrices.
	double dvgam[(pgam * pgam)];
	diagonalize(vgam, pgam, dvgam);
	
	double dLAM[(n * n)];
	diagonalize(LAM, n, dLAM);
	
	double Vtmp[(pgam * pgam)];

	if(pgam > n){
		//inversion of nxn matrix: 
		//Vtmp = vgam - vgam%*%t(Xgam)%*%solve(LAM + Xgam%*%vgam%*%t(Xgam))%*%Xgam%*%vgam;
		
		double Xvgam[(n * pgam)];
		F77_NAME(dgemm)(&NoTrans, &NoTrans, &n, &pgam, &pgam, &oned, Xgam, &n, dvgam, &pgam, &zerod, Xvgam, &n);
 
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
		
		diagonalize(vgam, pgam, Vtmp); //input = vgam, output = Vtmp
		
		F77_NAME(dgemm)(&Trans, &NoTrans, &pgam, &pgam, &n, &moned, Xvgam, &n, invLAMXvX_Xvgam, &n, &oned, Vtmp, &pgam);
		
	}else{

		//inversion of dxd matrix (d=dim(GAM)):
		//Vtmp = solve(t(Xgam)%*%invLAM%*%Xgam + solve(vgam));    
		
		double sqrtinvLAMXgam[(n * pgam)];
		F77_NAME(dgemm)(&NoTrans, &NoTrans, &n, &pgam, &n, &oned, dsqrtinvLAM, &n, Xgam, &n, &zerod, sqrtinvLAMXgam, &n);
		
		// 1) compute t(X)%*%LAM%*%X, 2) add 1/vgam[i] to diagonal elements, 
		// 3) compute inverse (first cholesky decomposition, then inverse from that).
		F77_NAME(dgemm)(&Trans, &NoTrans, &pgam, &pgam, &n, &oned, sqrtinvLAMXgam, &n, sqrtinvLAMXgam, &n, &zerod, Vtmp, &pgam);		
		
		for (int i=0; i<pgam; i++) {
			Vtmp[(pgam+1) * i] += 1.0/vgam[i];
		}
		matrixInverse(Vtmp, pgam);
	}
	for (int i=0; i<(pgam * pgam); i++) {
		Vgam[i] = Vtmp[i];
	}
	
	//Bgam = Vgam %*% (invvgam%*%bgam + tXgam%*%invLAM%*%Z);
	double Btmp[pgam];
	double invvbgam_XinvLAMZ[pgam];	//for now this is only invvgam%*%bgam
	for (int i=0; i<pgam; i++) {
		invvbgam_XinvLAMZ[i] = (1.0/vgam[i]) * bgam[i];
	}
	double invLAMZ[n];
	for (int j=0; j<n; j++) {
		invLAMZ[j] = (1.0/LAM[j]) * Z[j];
	}
	F77_NAME(dgemv)(&Trans, &n, &pgam, &oned, Xgam, &n, invLAMZ, &one, &oned, invvbgam_XinvLAMZ, &one); 
	// here invvbgam_XinvLAMZ is updated to include everything that is multiplied with Vgam resulting in Bgam
	
	F77_NAME(dgemv)(&Trans, &pgam, &pgam, &oned, Vgam, &pgam, invvbgam_XinvLAMZ, &one, &zerod, Btmp, &one);	
	for (int i=0; i<pgam; i++) {
		Bgam[i] = Btmp[i];
	}
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
					double *b, double *v, double *Pi,
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
	
	double vgam[pgam];
	vectorSubset(v, GAM, p, vgam);
	
	double bgam[pgam];
	vectorSubset(b, GAM, p, bgam);
	
	double Vgam[pgam * pgam];
	double Bgam[pgam];
	Vgam_Bgam_update(Vgam, Bgam, pgam, n, Xgam, vgam, bgam, LAM, Z);
	
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
			
			double vtry[ptry];
			vectorSubset(v, GAMtry, p, vtry);
			
			double btry[ptry];
			vectorSubset(b, GAMtry, p, btry);
			
			double Vtry[(ptry * ptry)];
			double Btry[ptry];
			Vgam_Bgam_update(Vtry, Btry, ptry, n, Xtry, vtry, btry, LAM, Z);
			
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
	
	double vstar[pstar];
	vectorSubset(v, GAMstar, p, vstar);
	
	double bstar[pstar];
	vectorSubset(b, GAMstar, p, bstar);
	
	double Vstar[(pstar * pstar)];
	double Bstar[pstar];
	Vgam_Bgam_update(Vstar, Bstar, pstar, n, Xstar, vstar, bstar, LAM, Z);
	
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
		double Vprec[(pstar * pstar)], dvprec[(pstar * pstar)];
		for (int i = 0; i < (pstar * pstar); i++) {
			Vprec[i] = Vstar[i];
		}
		matrixInverse(Vprec, pstar);  //input Vstar, output inverse(Vstar)
		
		double vprec[pstar];
		for (int i = 0; i < pstar; i++) {
			vprec[i] = 1.0/vstar[i];
		}
		diagonalize(vprec, pstar, dvprec);
		
		logpostbeta = dmvnorm(betastar, pstar, Bstar, Vprec, Log, pstar);
		logpriorbeta = dmvnorm(betastar, pstar, bgam, dvprec, Log, pstar);
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
	double Tr = (1/T) * r;
	
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
			if(U <= 1/(1 + Y)){
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
void C_update(double *C, double *Pi, int *GAM, double *d,
              int p, int k, double CpropRange, int *accept_count){

  double sum_d = 0;
  for (int i = 0; i < p; i++){
    sum_d += d[i];
  }

  //Propose new C' (if Cprop.range==0 sample from uniform(0,1), otherwise sample from (max{0, C' - 0.2}, min{C' + 0.2, 1})):
  double Cprop;
  Cprop = runif(fmax(0.0, *C-CpropRange), fmin(*C+CpropRange, 1.0));

  double aprop[p];
  double bprop[p];
  double a[p];
  double b[p];
  for (int i = 0; i < p; i++){
    aprop[i] = Cprop * k * p * d[i] + (1.0-Cprop) * k * sum_d;
    bprop[i] = p * sum_d - aprop[i];
    
    a[i] = *C * k * p * d[i] + (1.0 - *C) * k * sum_d;
    b[i] = p * sum_d - a[i];
  }
  
  //Compute density for proposed C':
  double logdens_Cprop = 0.0;
  for (int i = 0; i < p; i++){
    logdens_Cprop += dbeta(Pi[i], aprop[i]+GAM[i], bprop[i]+1.0-GAM[i], 1);
  }
  
  //Compute density for actual C:
  double logdens_C = 0.0;
  for (int i = 0; i < p; i++){
    logdens_C += dbeta(Pi[i], a[i]+GAM[i], b[i]+1.0-GAM[i], 1);
  }
  
  double dens_ratio = exp(logdens_Cprop - logdens_C);
  double accept_prob = fmin(1.0, dens_ratio);
  
  double U = runif(0.0, 1.0);
  if(U < accept_prob){
    *C = Cprop;
    
    *accept_count = *accept_count + 1;
  }
  
  //Rprintf("%.5e %.5e %.5e\n", *C, Cprop, accept_prob);
  //Rprintf("%.5e %.5e\n", logdens_C, logdens_Cprop);
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
// alte version von Pi_update
/*
void Pi_update(double *Pi, int *GAM, double C, double *d,  
               int p, int k){
  
  double sum_d = 0.0;
  for (int i = 0; i < p; i++){
    sum_d += d[i];
  }
  
  //Direct sampling:
  double a[p];
  double b[p];
  for (int i = 0; i < p; i++) {
    a[i] = C * k * p * d[i] + (1.0-C) * k * sum_d;
    b[i] = p * sum_d - a[i];
    
    Pi[i] = rbeta(a[i]+GAM[i], b[i]+1.0-GAM[i]);
  }
}
*/


// neue version von Pi_update
void Pi_update(double *Pi, int *GAM, double *aBeta, double *bBeta, int p){
   //Direct sampling:
  for (int i = 0; i < p; i++) {
    Pi[i] = rbeta(aBeta[i]+GAM[i], bBeta[i]+1.0-GAM[i]);
  }
}


/* Old version (with d, k and C):
void logisticVS(double *X, double *Y, int *n, int *p,
      	double *b, double *v, int *k, double *d,
				int *block, int *MCMC, int *thins, double *CpropRange, 
        int *Piupdate, int *Cupdate, int *burnin){
	
	GetRNGstate();	
	
	int thin[1];
	thin[0] = thins[0];
	double T = 1.0;
  int accept_count[1];
  accept_count[0] = 0;
	
	int K2;
	int numneigh[1];
	int select[1];
	double logpriorGAMbeta[1];
	double logpostGAMbeta[1];
	double logLAMz[1];
  
  double sum_d = 0;
  for (int i = 0; i < *p; i++){
    sum_d += d[i];
  }
	
  double C[1];
  C[0] = 1.0;
  if(*Cupdate==1){
    C[0] = runif(0.0, 1.0);
  }
  
  double Pi[*p];
  for (int i = 0; i < *p; i++){
    Pi[i] = (float)(*k) * d[i]/sum_d;
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
	
  FILE *fidPi, *fidC, *fidlogprob, *fidnumneigh, *fidselect, *fidaccept;
    fidC = fopen("C.txt", "w");
    fidPi = fopen("Pi.txt", "w");
    fidlogprob = fopen("logprobs.txt", "w");
    fidnumneigh = fopen("numneighs.txt", "w");
    fidselect = fopen("selects.txt", "w");
    fidaccept = fopen("acceptrates.txt", "w");
	
	for (int K = 0; K < *MCMC; K++) {
		
		// Allow R interrupts; check every 10 iterations
		if (!(K % 10)) R_CheckUserInterrupt();
  	if (!((K+1) % 1000)) Rprintf("iteration %d\n", K+1);
			
    //Sample C and update the px1 vector Pi
    if(*Piupdate == 1){
      Pi_update(Pi, GAM, *C, d, *p, *k);
    
      if(*Cupdate == 1){
        C_update(C, Pi, GAM, d, *p, *k, *CpropRange, accept_count);
	      //Rprintf("%.5e\n", *C);
    
        if(!(K % 20) & (K < *burnin) & (K > 0)) {
          if(accept_count[0] / 20.0 < 0.3){
            *CpropRange = *CpropRange * 0.9;
          }
          if(accept_count[0] / 20.0 > 0.5){
            *CpropRange = *CpropRange * 1.1;
         }
          //Rprintf("accept_rate %.3f\n", accept_count[0] / 20.0);
          //Rprintf("CpropRange %.3f\n", *CpropRange);
          fprintf(fidaccept, "%.4f ", accept_count[0] / 20.0);
      
          accept_count[0] = 0;
        }
      }
    }
  
		betaGAM_update(*n, *p, X,
					   b, v, Pi,
					   beta, GAM, LAM, Z,
					   block, T, numneigh, select, logpriorGAMbeta, logpostGAMbeta);
		
		LAMz_update(*n, *p, X, Y,
					beta, GAM, LAM, Z,
					T, logLAMz);
		
		if (!(K % *thin)) {
			K2 = K / *thin;
		
        	fprintf(fidlogprob, "%.4f ", logpriorGAMbeta[0] + logLAMz[0]);
        	fprintf(fidnumneigh, "%d ", numneigh[0]);
        	fprintf(fidselect, "%d ", select[0]);
		
			for (int i = 0; i < *p; i++) {
        fprintf(fidPi, "%.10e\t", Pi[i]);
        
				if (GAM[i] == 1){
					fprintf(fidbeta, "%d %d %.10e\n", i+1, K2+1, beta[i]);
					fprintf(fidGAM, "%d %d %d\n", i+1, K2+1, 1);
				}
			}
      fprintf(fidPi, "\n");
      fprintf(fidC, "%.6e\n", *C);
		}
	}
	
  fclose(fidC);
  fclose(fidPi);
	fclose(fidbeta);
	fclose(fidGAM);
	fclose(fidlogprob);
	fclose(fidnumneigh);
	fclose(fidselect);
  fclose(fidaccept);
	
	PutRNGstate();	
}
*/

//New version (with aBeta and bBeta)
void logisticVS(double *X, double *Y, int *n, int *p,
  			double *b, double *v, double *aBeta, double *bBeta,
				int *block, int *MCMC, int *thins, 
        int *Piupdate){
	
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
  
  Rprintf("Starting Pi with prior mean, e.g. Pi[1] = %.3e\n", Pi[0]);
  
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
	
  FILE *fidlogprob, *fidnumneigh, *fidselect;
    fidlogprob = fopen("logprobs.txt", "w");
    fidnumneigh = fopen("numneighs.txt", "w");
    fidselect = fopen("selects.txt", "w");
	
	for (int K = 0; K < *MCMC; K++) {
		
		// Allow R interrupts; check every 10 iterations
		if (!(K % 10)) R_CheckUserInterrupt();
  	if (!((K+1) % 1000)) Rprintf("iteration %d\n", K+1);
			
    if(*Piupdate == 1){
      Pi_update(Pi, GAM, aBeta, bBeta, *p);
    }
  
		betaGAM_update(*n, *p, X,
					   b, v, Pi,
					   beta, GAM, LAM, Z,
					   block, T, numneigh, select, logpriorGAMbeta, logpostGAMbeta);
		
		LAMz_update(*n, *p, X, Y,
					beta, GAM, LAM, Z,
					T, logLAMz);
		
		if (!(K % *thin)) {
			K2 = K / *thin;
		
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
	
	fclose(fidbeta);
	fclose(fidGAM);
	fclose(fidlogprob);
	fclose(fidnumneigh);
	fclose(fidselect);
	
	PutRNGstate();	
}



