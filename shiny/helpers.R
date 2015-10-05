library(ghyp)
library(Rcpp)
sourceCpp("R0bar.cpp")

source('mcrollapplr.R')
options(mc.cores=detectCores())

allsplines=readRDS("allsplines_N252_Navg1000.rds")



estimateSharpeRatio = function(myts,N=252,numPerm=100){
  print(paste(start(myts),end(myts)))
  mySNR=mcrollapplyr(myts,N,mylapply=mclapply,FUN=function(x){
    print(rnorm(1))
    myfit=tryCatch(fit.tuv(x,silent=TRUE,nu=4.5),error=function(e){print(e);return(NA)})
    nuOrig=coef(myfit)$nu
    nu=round(10*coef(myfit)$nu)/10
    nu=min(10,nu)
    nu=max(2.5,nu)
    R0bar=R0bar(coredata(x),numPerm = numPerm)
    if(abs(R0bar)<1){
      R0bar=0
    }
    
    SNR_raw=allsplines[[as.character(nu)]](abs(R0bar))    #splines have been calibrated for positive SNRs
    SNR_raw=SNR_raw*  (sign(SNR_raw)>0)                   #if(spline(abs(x)))<0, then set SNR to 0 (this occurs when SNR is very small)     
    SNR=SNR_raw*sign(R0bar)                               #then apply sign of R0bar 
    
    mu=mean(x)
    mu3=mean((x-mu)^3)
    mu4=mean((x-mu)^4)
    sigma3=mean((x-mu)^2)^(3/2)
    sigma4=mean((x-mu)^2)^2
    gamma3=mu3/sigma3
    gamma4=mu4/sigma4
    
    SNR0=mean(x)/sd(x)
    SNR0sd=(1+SNR0^2/4*(mu4/sigma4-1)-SNR0*mu3/sigma3)/sqrt(length(x))
    return(list(nuOrig=nuOrig,nu=nu,SNR=SNR,SNR0=SNR0,SNR0sd=SNR0sd,R0bar=R0bar,N=length(x)))
  })
  return(mySNR)
}




