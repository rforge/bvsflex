callBVS=function(i,shrink){
  library(bvsflex)
  #  dir <- paste("BVS_",dvec[i],"_",format(Sys.time(), "%Y-%m-%d--%Hh%M"), sep="")
  dir <- paste("BVS_",dvec[i], "_p",p,"_k",k,"_shrink",shrink,sep="")
  #  dirvec[i,2]=dir
  print(dir)
  if (!file.exists(dir)) dir.create(dir)
  setwd(dir)
  gc()
  sink(file("iteration.txt", open = "wt"))
  tryerr <- try(time <- 
                  system.time(logisticVS(X=X, Y=Y, b=b, v=v, k=k, d=eval(parse(text =dvec[i]))+shrink, block=block, 
                                         Cupdate=FALSE, Piupdate=TRUE, 
                                         MCMC=BVS.TOT, thinn=BVS.thin, burnin=burnin, seed=1234, outdir=NULL )), 
                silent=FALSE)
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
