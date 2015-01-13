getBVSflex=function(path,burnin,thin,extraThin,m,geneSymbols){ 
  require(Matrix)
  
  wd=getwd()
  setwd(path)
  
  #Get the top m models:
    logprobs <- scan("logprobs.txt",na.strings="1.#INF")
    lp.after <- logprobs[-(1:(burnin/thin))]
    logprobs=logprobs[seq(1,length(logprobs),by=extraThin)]
    lp.after=lp.after[seq(1,length(lp.after),by=extraThin)]
    
    lp.max <- sort(lp.after, decreasing=TRUE)[m]
    lp.top <- which(lp.after >= lp.max) + (burnin/(thin*extraThin))

  gam <- read.table("GAM.txt",fill=T)
  gam= gam[complete.cases(gam),]
  gam <- sparseMatrix(i=gam[,1], j=gam[,2], x=gam[,3])
  gam <- as.matrix(gam)
  gam=gam[,seq(1,ncol(gam),by=extraThin)]
  modelsize <- apply(gam, 2, sum)


  bet <- read.table("beta.txt",fill=T)
  bet= bet[complete.cases(bet),]
  bet <- sparseMatrix(i=bet[,1], j=bet[,2], x=bet[,3])
  bet <- as.matrix(bet)
  bet=bet[,seq(1,ncol(bet),by=extraThin)]
  
  setwd(wd)

#   #(1) "0.05-quantile probability model":
#   gam.mean <- apply(gam[,-(1:(burnin/(thin*extraThin))),drop=FALSE], 1, mean)
#   #plot(gam.mean)
#   gam.qpm <- which(gam.mean > 0.05)
#   bet.qpm <- apply(bet[gam.qpm,-(1:(burnin/(thin*extraThin))),drop=FALSE], 1, mean)
  
  #(2) Top models based on logprobs:
  top.beta <- bet[,lp.top,drop=FALSE]
  top.list <- table(unlist(apply(top.beta, 2, function(x){which(x!=0)})))


  gamma.p <- t(gam)
  beta.p <- t(bet)
  bet.mean <- apply(beta.p[-(1:(burnin/(thin*extraThin))),], 2, function(x){mean(x[x!=0])})
  bet.mean[is.na(bet.mean)] <- 0
  bet.sd <- apply(beta.p[-(1:(burnin/(thin*extraThin))),], 2, function(x){sd(x[x!=0])})
  bet.sd[is.na(bet.sd)] <- 0
  #bet.mean=bet.mean[bet.sd!=0]
  #bet.sd=bet.sd[bet.sd!=0]
  
  gam.mean <- apply(gamma.p[-(1:(burnin/(thin*extraThin))),], 2, mean)
  gam.sd <- apply(gamma.p[-(1:(burnin/(thin*extraThin))),]
                  , 2, sd)
  #gam.mean=gam.mean[gam.sd!=0]
  #gam.sd=gam.sd[gam.sd!=0]

  print(summary(logprobs))
  print(summary(modelsize))

  results=list(
    logprobs=logprobs,lp.top=lp.top, modelsize=modelsize,top.beta=top.beta,top.list=top.list,
    gamma=gamma.p, beta=beta.p,beta.mean=bet.mean,beta.sd=bet.sd,gamma.mean=gam.mean,gamma.sd=gam.sd,
    burnin=burnin,thin=thin,extraThin=extraThin ,geneSymbols=geneSymbols
    )
class(results)="bvsflex"
results
}


plotPaths.bvsflex=function(object,burnin,thin){
  require(ggplot2)
  require(grid)
  grid.newpage();pushViewport(viewport(layout = grid.layout(2,1)))
  with(object, print(qplot((1:length(logprobs))[!is.na(logprobs)&logprobs>-1e6],logprobs[!is.na(logprobs)&logprobs>-1e6], geom="line",alpha = I(1/1))+geom_vline(xintercept=(burnin/(thin)),colour=I("red"))#+scale_x_discrete(breaks=NULL)
        ,vp=viewport(layout.pos.row =1,layout.pos.col =1)))
  with(object,print(qplot(1:length(modelsize),modelsize, geom="line", alpha = I(1/1))+geom_vline(xintercept=(burnin/(thin)),colour=I("red"))#+scale_x_discrete(breaks=NULL)
        ,vp=viewport(layout.pos.row =2,layout.pos.col =1)))
}


plotACF.bvsflex=function(object,burnin,thin,...) with(object, plot(acf(modelsize[-(1:(burnin/(thin)))],lag=50),main="modelsize",...))


predict.bvsflex=function(object,X2,sdX) with(object, X2 %*% t(t(top.beta)*sdX))


ROC.bvsflex=function(object,X2,Y2,sdX,...){
  require(pROC)
  require(ROCR)
  
  linpred <- predict.bvsflex(object,X2,sdX)
  
  err <- rep(NA, ncol(linpred))
  for(i in 1:ncol(linpred)){
    tab <- table(Y2, linpred[,i]>0)
    # tab <- table(Y2, factor(linpred[,i]>0,levels=c(FALSE,TRUE)))
    err[i] <- (tab[1,2]+tab[2,1])/sum(tab)
  }
  
  auc <- rep(NA, ncol(linpred))
  roc <- roc(Y2, linpred[,1])
  for(i in 1:ncol(linpred)){
    roc <- roc(Y2, linpred[,i]) 
    auc[i] <- roc$auc
  }
  pred <- prediction(linpred, matrix(rep(Y2, m),ncol=m))
  perf <- performance(pred,"tpr","fpr")
  
  op<-  par(mar=c(4,4,2,2), mgp=c(2.5,1,0))
  plot(perf, col="gray", ...)
  plot(perf, avg="vertical", spread.estimate="stderror",
       show.spread.at=seq(min(linpred),max(linpred),0.1), add=TRUE)
  abline(0,1)
  par(font=2)
  legend("bottomright", legend=c(paste("mean ERR =",round(mean(err),3)),
                                 paste("mean AUC =",round(mean(auc),3))), bty="n",cex=1.3)
  par(op)
}


plotBeta.bvsflex=function(object,chromosomes=NULL,...){
  require(gplots)
  op<-par(mgp=c(2.5,1,0), mar=c(4,4,1,1))
  #  locs <- match(unique(annoX$Chr), annoX$Chr)
  with(object,plotCI(beta.mean, uiw=beta.sd, ylab=as.list(expression(paste("Coefficient ",beta[gamma][i]))), 
                     xlab="Chromosome", xaxt="n", pch=20, gap=0, barcol="gray", cex=0.5, 
                     ylim=c(min(beta.mean-beta.sd), max(beta.mean+beta.sd)),...))
  #     axis(1, at=match(1:22, suppressWarnings(as.numeric(annoX$Chr))), labels=1:22, las=2)
  if (!is.null(chromosomes)) axis(1, at=match(1:22, suppressWarnings(as.numeric(chromosomes))), labels=1:22, las=2)
  abline(h=0)
  par(op)
}

summary.bvsflex=function(object,cutoff=rev(c(0.01,0.05,0.1,0.2,0.5))){ #cutoff may be a vector
  cat("Number of selected:","\n")
  gmean=with(object,data.frame(gamma.mean,row.names=geneSymbols))
  gmean=gmean[order(gmean$gamma.mean,decreasing=T),,drop=F]
  numincl=as.matrix(sapply(cutoff,function(x) sum(gmean>x)))
  rownames(numincl)=paste0(">",cutoff)
  incl=(lapply(cutoff,function(x) (gmean)[gmean>x,,drop=F]))
  names(incl)=rownames(numincl)
  print(numincl)
  cat("\n\n")
  print(incl)
  cat("\n\n")
}



plotGamma.bvsflex=function(object,cutoff=c(0.05,0.1),chromosomes=NULL,names=NULL,ylim=c(0,0.24),... ) {#cutoff vector of cutoffs (maximum length 3); text vector of same length specifying whether associated gene symbols should be printed
  require(gplots)
  if (is.null(names)) names=rep(FALSE,length(cutoff))
  op<-par(mgp=c(2.5,1,0), mar=c(4,4,1,1))
  means=with(object,gamma.mean[gamma.sd!=0])
  sds=with(object,gamma.sd[gamma.sd!=0])
  #plotCI(means, li=ifelse(((gam.mean[gam.sd!=0]>2*Pi0)& (means-sds >0)), means-sds, NA), ui=ifelse(((gam.mean[gam.sd!=0]>2*Pi0)& (means-sds >0)), means+sds, NA), 
  plotCI(means, li=ifelse(means-sds >0, means-sds, 0), ui=means+sds, 
         ylab=expression(paste("Posterior inclusion probability p(",gamma[i]," |data)")), 
         xlab="Location (chromosome)", xaxt="n", main="", pch=20, col="black", gap=0, barcol="gray", ylim=ylim,...)
  #axis(1, at=match(1:22, suppressWarnings(as.numeric(annoX$Chr))), labels=1:22, las=2)
  if (!is.null(chromosomes)) axis(1, at=match(1:22, suppressWarnings(as.numeric(chromosomes))), labels=1:22, las=2)
  abline(h=0, col="gray")
  #Pi0 <- mean(Pi)
  
  col=c("green","orange","red")
  pch=c(20,17,15)
  
 # leg=NULL
  for (i in 1:length(cutoff)){
    whichfreq <- which(means>cutoff[i])
    points(whichfreq, means[whichfreq], pch=pch[i], col=col[i]) 
    if(names[i] & length(whichfreq)>0){
      text(whichfreq, means[whichfreq], labels=object$geneSymbols[whichfreq], pos=3, col=col[i], cex=1.2)
    }
  #  leg=c(leg,substitute(expression('p(',gamma[i],' |data) > ',cutoff), list(cutoff=cutoff[i])))
  }
 # legend(x="topright", legend=(unlist(leg)),  
 #        pch=pch[1:length(cutoff)], col=col[1:length(cutoff)], bty="n")
  legend(x="topright", legend=paste("p(gamma|data)>",cutoff),  
         pch=pch[1:length(cutoff)], col=col[1:length(cutoff)], bty="n")
 par(op)

}


plotBetaGamma.bvsflex=function(object,cutoff,ylim=c(0,0.24),... ){
  require(gplots)
  op<- par(mgp=c(2.5,1,0), mar=c(4,4,1,1))
  with(object,plotCI(x=beta.mean, y=gamma.mean, li=ifelse(gamma.mean-gamma.sd > 0, gamma.mean-gamma.sd, 0), ui=gamma.mean+gamma.sd, 
         xlab=expression(paste("Coefficient ",beta[gamma][i])), 
         ylab=expression(paste("Posterior inclusion probability p(",gamma[i]," |data)")),
         pch=20, ylim=ylim,  gap=0, barcol="gray", err="y",...))
  with(object,plotCI(x=beta.mean, y=gamma.mean, uiw=beta.sd, pch=20, ylim=ylim, gap=0, barcol="gray", err="x", add=TRUE))
  abline(v=0, col="gray")
  abline(h=0, col="gray")
  abline(h=cutoff, col="blue")
  whichfreq0 <- which(object$gamma.mean>cutoff)
  if(length(whichfreq0)>0){
    with(object, points(beta.mean[whichfreq0], gamma.mean[whichfreq0], pch=20, col="blue"))
    with(object,text(x=beta.mean[whichfreq0], y=gamma.mean[whichfreq0], labels=geneSymbols[whichfreq0], pos=3, col="blue", cex=1))
    #  text(x=beta.mean[whichfreq0], y=gamma.mean[whichfreq0], labels=annoX[whichfreq0,"Symbol"], pos=3, col="blue", cex=1)
    #  text(x=bet.mean[whichfreq0], y=gam.mean[whichfreq0], labels=colnames(X)[whichfreq0], pos=3, col="blue", cex=1)
    legend(x="topleft", legend=c(expression(),
                                 substitute(paste("p(",gamma[i]," |data) > ",cutoff), list(cutoff=cutoff))), 
           pch=20, col=c("blue"), bty="n", cex=1.0)
  }
  par(op)
}


