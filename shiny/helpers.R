library(sharpeRratio)
library(data.table)

source('mcrollapplr.R')
options(mc.cores=detectCores())


estimateSharpeRatio = function(myts,N=252,numPerm=100){
  #print(paste(start(myts),end(myts)))
  mySNR=mcrollapplyr(myts,N,mylapply=mclapply,FUN=function(x){
    myres=as.data.table(estimateSNR(x))
    mu=mean(x)
    mu3=mean((x-mu)^3)
    mu4=mean((x-mu)^4)
    sigma3=mean((x-mu)^2)^(3/2)
    sigma4=mean((x-mu)^2)^2
    gamma3=mu3/sigma3
    gamma4=mu4/sigma4
    
    SNR0=mean(x)/sd(x)
    SNR0sd=(1+SNR0^2/4*(mu4/sigma4-1)-SNR0*mu3/sigma3)/sqrt(length(x))
    DT=cbind(myres,data.table(SNR0=SNR0,SNR0sd=SNR0sd,N=length(x)))
    return(DT)
  })
  mySNR=rbindlist(mySNR)
  mySNR[,SNR0:=SNR0*sqrt(252)]
  mySNR[,SNR:=SNR*sqrt(252)]
  mySNR=cbind(date=tail(index(myts),nrow(mySNR)),mySNR)

  return(mySNR)
}




