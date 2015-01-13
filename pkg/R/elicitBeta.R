#,xlim=c(rangeModelsize[1],1.5*rangeModelsize[2])/p

plot.elicit <- function(object, which=1, ylim = NULL,xlim=NULL,dpbeta=T) {
  a=object$beta[[which]]$a
  b=object$beta[[which]]$b

  x <- seq(0, 1, length = 10025)
  
  
  fx <- cbind(dbeta(x, a,b), pbeta(x, a,b))
  f <- fx; f[fx == Inf] <- 1e100
 
  if (is.null(xlim)) xlim=c(object$beta[[which]]$rangePi[1],1.5*object$beta[[which]]$rangePi[2])
  if (is.null(ylim)) ylim=c(0,max(f[,1])+.2)
  
  if(dpbeta) par(mfrow=c(2,1))
  plot(x, f[,1], ylab="", type="l", ylim=ylim,xlim=xlim, main = sprintf("dbeta(x, a=%g, b=%g)", a,b))
  abline(h = 0, col="gray", lty=3)
  abline(v=object$beta[[which]]$priorPi,col="red")
  abline(v=a/(a+b),col="black",lty=3)
  abline(v=(a-1)/(a+b-2),col="blue",lty=3)
  abline(v=qbeta(0.5,a,b),col="lightblue",lty=3)
  abline(v=object$beta[[which]]$rangePi,col="green",lty="dashed")
  legend("topright", c("mean","median","mode","crit","rangeModelsize"),
         col=c("black","lightblue","blue","red","green"), lty=c(3,3,3,1,2), bty = "n")
  
  if(dpbeta){
    plot(x, f[,2], ylab="", type="l", ylim=c(0,1.1),xlim=xlim, main = sprintf("pbeta(x, a=%g, b=%g)", a,b))
    abline(h = 0:1, col="gray", lty=3)
    abline(h = 0, col="gray", lty=3)
    abline(v=object$beta[[which]]$priorPi,col="red")
    abline(v=a/(a+b),col="black",lty=3)
    abline(v=(a-1)/(a+b-2),col="blue",lty=3)
    abline(v=qbeta(0.5,a,b),col="lightblue",lty=3)
    abline(v=object$beta[[which]]$rangePi,col="green",lty="dashed")
  }
  invisible(cbind(x, fx))
}


#p: inverse weight of one gene (=ncols(X) in case of noinformation)
#modelsize: expected number of selected genes (modelsize)  (previously ==k)
#rangeModelsize: range of modelsizes that should be covered with error 1-coverage,
#upperBound.aBeta): Upper bound for possible values of aBeta
#coverage; coverage probability
#crit: either "mode" or "mean"


elicitBeta=function(p=NULL,d=NULL,modelsize=10,rangeModelsize=c(0,3*modelsize),coverage=0.9, crit="mode",lowerBound.aBeta=1, upperBound.aBeta=2*modelsize^2){
  
  if (lowerBound.aBeta<1 & (crit=="mode" | crit=="median")) {warning("if crit='mode' or 'median', lowerBound.aBeta has to be >=1. Setting lowerBound.aBeta=1");lowerBound.aBeta=1}
  if (modelsize <rangeModelsize[1] | modelsize >rangeModelsize[2]) stop("rangeModelsize must include modelsize")
  if (is.null(p) & is.null(d)) stop("either p or d must be given.")
  
  if (crit=="mean") bwork=function(a,priorPi) a/priorPi-a
  else if (crit=="mode") bwork=function(a,priorPi) (a-1)/priorPi+2-a
  else if (crit=="median") bwork=function(a,priorPi) a/priorPi-a+(2-1/priorPi)/3   # (a-1/3)/priorPi-a+2/3
  else stop("crit must be eather 'mode', 'mean' or 'mean'")

  
  search=function(a,coverage,priorPi,range,lowerBound.aBeta){
    #in order to force a to be >1 (condition such that density is uni modal (if a>1,b>1)), use a=a+1
    #for mode, possible values of a have to include 0 in uniroot, but for mode to be defined, aBeta has to be >1
    a=a+lowerBound.aBeta 
    #suppressWarnings(integrate(function(r) (1/beta(a,bwork(a,priorPi)))*r^(a-1) *(1-r)^(bwork(a,priorPi)-1),lower=range[1],upper=range[2])$value-coverage)
    suppressWarnings(pbeta(range[2],a,bwork(a,priorPi))-pbeta(range[1],a,bwork(a,priorPi))-coverage)
  }
  prior.upperBound.aBeta=upperBound.aBeta
  chooseBeta=function(priorPi){ 
    aa=list()
    class(aa)="try-error"
    while(class(aa)=="try-error"){
      aa <- try(uniroot(search, 
                        c(0, upperBound.aBeta),          #possible values of aBeta
                        tol = 0.001, 
                        coverage= coverage,              #coverage probability
                        priorPi=priorPi,             #prior mean or mode of pi   
                        range=rangePi,                 #range of pi's that should be covered with error 1-coverage, 
                        #e.g. range=c(0,5*priorPi) and priorPi=k/p would induce that we allow prior model sizes in the interval of (0,5k)
                        lowerBound.aBeta=lowerBound.aBeta), silent = TRUE)
      upperBound.aBeta=upperBound.aBeta-1
      if(upperBound.aBeta<0) stop("not converged")
    }
    #str(aa)
    a=aa$root+lowerBound.aBeta
    if(upperBound.aBeta+1<prior.upperBound.aBeta) cat("upperBound.aBeta ",upperBound.aBeta+1,"\n")
    b=bwork(a, priorPi)
    
    summaryPi=data.frame(
      mean=a/(a+b),
      mode=(a-1)/(a+b-2),
      median=qbeta(0.5,a,b),    #median=(a-1/3)/(a+b+2/3),
      variance=(a*b)/((a+b)^2*(a+b+1))
    )
    summaryModelsize=summaryPi*p
    return(list(a=a,b=b, priorPi=priorPi,rangePi=rangePi,summaryPi=summaryPi,summaryModelsize=summaryModelsize))
  }
  
  if (!is.null(d)) {
    if (!is.null(p)) cat("Spezification of p disregarded when d vector is given\n")
    p=length(d)
    sum.d=sum(d)
    pp=sum.d/d 
  } else pp=p
  
  
  priorPi=modelsize/pp
  rangePi=rangeModelsize/p
  if (length(unique(priorPi))==1) {
    res=chooseBeta(priorPi=priorPi[1])
    res=replicate(p,res,simplify=F)
  }
  else res=lapply(priorPi,chooseBeta)
  
 result=list(beta=res,spec=list(p=p,d=d,modelsize=modelsize,rangeModelsize=rangeModelsize,coverage=coverage, crit=crit,lowerBound.aBeta=lowerBound.aBeta, upperBound.aBeta=upperBound.aBeta))
 class(result)="elicit"
 return(result)
 
}
