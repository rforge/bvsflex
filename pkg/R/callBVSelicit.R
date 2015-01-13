callBVSelicit=function(ind,shrink,crit="mode",ylim=c(0,55)){
  library(bvsflex)
  #  dir <- paste("BVS_",dvec[i],"_",format(Sys.time(), "%Y-%m-%d--%Hh%M"), sep="")
  dir <- paste("BVS_",dvec[ind], "_p",p,"_k",k,"_shrink",shrink,sep="")
  #  dirvec[i,2]=dir
  print(dir)
  if (!file.exists(dir)) dir.create(dir)
  setwd(dir)
  gc()
  sink(file("iteration.txt", open = "wt"))
  tryerr <- try({
    elicit<-elicitBeta(d=eval(parse(text =dvec[ind]))+shrink,modelsize=k,rangeModelsize=rangeModelsize,coverage=coverage, crit=crit,lowerBound.aBeta=1, upperBound.aBeta=10000)
    pdf("betadists.pdf")
    for(iii in 1:p) {
      plot(object=elicit,which=iii,dpbeta=F,ylim=ylim,xlim=c(0,0.2))
      abline(v=k/p,col="orange",lty="dashed",lwd=3)
    }
    dev.off()
    time <- system.time(logisticVSbeta(X=X, Y=Y, b=b, v=v, block=block, 
                                       Piupdate=T,aBeta=sapply(elicit$beta,function(x) x$a),bBeta=sapply(elicit$beta,function(x) x$b),
                                       MCMC=BVS.TOT, thinn=BVS.thin, burnin=burnin, seed=1234, outdir=NULL ))
  },     silent=FALSE)
  sink()
  if(class(tryerr) == "try-error"){ 
    message="Call to logisticVS failed."
    try({
      save(message,tryerr,file="errorlog.RData")
      write(tryerr,file="errorlog.txt")
      save(annoX,file="annoX.Rdata")
      save(X,file="X.Rdata")
      d=lapply(dvec,function(x) eval(parse(text =x)))
      names(d)=dvec
      save(b,v,k,burnin,BVS.TOT,BVS.thin,d,file="parameters.Rdata")
    })
    #   stop("Call to logisticVS failed.\n") 
  } else {
    try({
      save(time,tryerr,file="succeeded.RData")
      save(annoX,file="annoX.Rdata")
      save(X,file="X.Rdata")
      d=lapply(dvec,function(x) eval(parse(text =x)))
      names(d)=dvec
      save(b,v,k,burnin,BVS.TOT,BVS.thin,d,file="parameters.Rdata")
    })
  }
  setwd("../")
  #write.table(dirvec,file="currentdir.txt")
}
