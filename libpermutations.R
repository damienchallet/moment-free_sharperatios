library(parallel,quietly = TRUE)
library(boot,quietly = TRUE)
library(BEST,quietly = TRUE)
library(Deducer,quietly=TRUE)


source('libdc.R')
source('SIGN.test.R')

if(Sys.info()[["nodename"]]=="stanislao-Precision-T5610"){
options(mc.cores=detectCores()/3)
}else{
options(mc.cores=detectCores())
}

estimateFromDeviations=function(N,Nperms,snr){
  vec=rnorm(N,mean = snr,sd = 1)
  
  ups_perm=unlist(mclapply(1:Nperms,function(x) num_records_up(cumsum(sample(vec)))))
  snrUp=estimateSNRFromNumRecordsGauss(N,ups_perm)
  downs_perm=unlist(mclapply(1:Nperms,function(x) num_records_down(cumsum(sample(vec)))))
  snrDown=-estimateSNRFromNumRecordsGauss(N,downs_perm)
  return(list(up=ups_perm,down=downs_perm,snrUp=snrUp,snrDown=snrDown,snrDirect=mean(vec)/sd(vec),vec=vec))
}

snr.boot.func=function(x,sel){
  x=x[sel]
  return(mean(x)/sd(x))
}

median.test <- function(x, y){
  z <- c(x, y)
  g <- rep(1:2, c(length(x), length(y)))
  m <- median(z)
  fisher.test(z < m, g)$p.value
}


sampleEstimates=function(N,Navg=100,numPerm=10*N,SNR=0.001,dtype="Gauss",withBoot=FALSE,nu=3,mylapply=mclapply,nreps=1,withBEST=FALSE){
  allres=mylapply(1:Navg,function(i){
    #    print(i)
    if(dtype=="Gauss"){
      vec0=rnorm(N,mean=SNR,sd=1)
    }else if(dtype=="Student"){
      if(nu>1){
        sdStudent=T_moments(2,nu)
      }else{
        stStudent=1
      }
      vec0=rt(N,nu)/sdStudent+SNR
    }
    else if(dtype=="Exponential"){
      vec0=sign(runif(N)-0.5)*rexp(N)/sqrt(2)+SNR
    }else if(dtype=="Uniform"){
      vec0=(runif(N)-0.5)*2*sqrt(3)+SNR
    }else if(dtype=="Cauchy"){
      vec0=rcauchy(N)+SNR
    }
    vec0=rep(vec0,nreps)
    myres=computeUpsDownsRandomPerms(vec0,Navg)
    #     
    #     myres=lapply(1:Nsamples,function(i){
    #       vec=cumsum(sample(vec0,replace = FALSE))
    #       up=num_records_up(vec)
    #       down=num_records_down(vec)
    #       return(list(up=up,down=down))
    #     })
    ups=myres$up
    downs=myres$down
    pH0_ups=sapply(ups,function(R) prob_num_records(R,N))
    pH0_downs=sapply(downs,function(R) prob_num_records(R,N))
    
    R0=rstat(vec0)
    RT=computeRTbar(vec0,numPerm)
    
    snr=mean(vec0)/sd(vec0)*sqrt(length(vec0))
    rbis=rstatBIS(vec0,numPermutations = numPerm)
    wilcox=wilcox.test(vec0)$p.value
    SIGN=SIGN.test(vec0)$statistic
    if(withBEST){
      BEST.res=BESTmcmc(vec0,verbose = FALSE)
      BEST.summary=summary(BEST.res)
      BEST=BEST.summary["mu","HDIup"]-BEST.summary["mu","HDIlo"]
    }else{
      BEST=NULL
    }
    if(withBoot){
      snr.boot=boot(vec0,snr.boot.func,N)
      snr.ci=boot.ci(snr.boot,type="bca")
      snr_boot=sum(snr.boot$t>0)/length(snr.boot$t)
      #      snr_sd=tail(diff(as.vector(snr.ci$bca)),1)/2
    }else{
      snr.boot=NULL
      snr.ci=NULL
      snr_sd=NA
      snr_boot=NA
    }
    return(list(ups=ups,downs=downs,
                up_avg=mean(ups),up_sd=sd(ups),
                down_avg=mean(downs),down_sd=sd(downs),
                snr=snr,snr_boot=snr_boot,ups=ups,downs=downs,
                pH0_ups=mean(pH0_ups),pH0_downs=mean(pH0_downs),
                pH0_ups_sd=sd(pH0_ups),pH0_downs_sd=sd(pH0_downs),
                wilcox=wilcox,
                BEST=BEST,
                rbis=rbis,
                R0=R0,
                RT=RT,
                SIGN=SIGN
    ))
  })
  
  
  R0s=sapply(allres,function(x) x$R0)
  pH0_ups=sapply(allres,function(x) x$pH0_ups)
  pH0_ups_sd=sapply(allres,function(x) x$pH0_ups_sd)
  pH0_downs=sapply(allres,function(x) x$pH0_downs)
  pH0_downs_sd=sapply(allres,function(x) x$pH0_downs_sd)
  
  wilcox=sapply(allres,function(x) x$wilcox)
  SIGN=sapply(allres,function(x) x$SIGN)
  
  BEST=sapply(allres,function(x) x$BEST)
  
  
  snrs=sapply(allres,function(x) x$snr)
  snr_avg_avg=mean(snrs)
  snr_avg_sd=sd(snrs)
  snr_boot=(sapply(allres,function(x) x$snr_boot))
  snr_boot_sd=sd(sapply(allres,function(x) x$snr_boot))
  
  rbis=sapply(allres,function(x) x$rbis)
  
  RT=sapply(allres,function(x) x$RT)
  
  return(list(R0s=R0s,
              snr_avg_avg=snr_avg_avg,snr_avg_sd=snr_avg_sd,
              snr_boot=snr_boot,snr_boot_sd=snr_boot_sd,
              snrs=snrs,
              #              ups_all=ups_all,downs_all=downs_all,
              N=N,Nsamples=numPerm,Navg=Navg,
              pH0_ups=pH0_ups,pH0_downs=pH0_downs,
              wilcox=wilcox,
              BEST=BEST,
              rbis=rbis,
              RT=RT,
              SIGN=SIGN
  ))
}

twosampleEstimates=function(N,Navg=100,numPerm=10*N,numSamplesH0=10000,SNR=0.001,dtype="Gauss",withBoot=FALSE,nu=3,mylapply=mclapply,volRatio=1,lengthRatio=1,diffType="records"){
  allres=mylapply(1:Navg,function(i){
    #    print(i)
    N1=N
    N2=floor(N*lengthRatio)
    if(dtype=="Gauss"){
      vec2=rnorm(N2,mean=SNR,sd=1)
      vec1=rnorm(N1)*volRatio
    }else if(dtype=="Student"){
      if(nu>1){
        sdStudent=T_moments(2,nu)
      }else{
        stStudent=1
      }
      vec2=rt(N2,nu)/sdStudent+SNR
      vec1=rt(N1,nu)/sdStudent*volRatio
    }
    else if(dtype=="Exponential"){
      vec2=sign(runif(N2)-0.5)*rexp(N2)/sqrt(2)+SNR
      vec1=sign(runif(N1)-0.5)*rexp(N1)/sqrt(2)*volRatio
    }else if(dtype=="Uniform"){
      vec2=(runif(N2)-0.5)*2*sqrt(3)+SNR
      vec1=(runif(N1)-0.5)*2*sqrt(3)*volRatio
    }
    myres=computeUpsDownsTwoSamplesRandomPerms(vec2,vec1,numPerm)
    ups=myres$up1-myres$up2
    downs=myres$down1-myres$down2
    
    myres=computeUpsDownsDiffTwoSamplesRandomPerms(vec2,vec1,numPerm)
    ups_v=myres$up
    downs_v=myres$down
    #     
    #     myres=lapply(1:Nsamples,function(i){
    #       vec=cumsum(sample(vec0,replace = FALSE))
    #       up=num_records_up(vec)
    #       down=num_records_down(vec)
    #       return(list(up=up,down=down))
    #     })
    
    #     pH0_ups=sapply(ups,function(R) prob_num_records(R,T))
    #     pH0_downs=sapply(downs,function(R) prob_num_records(R,T))
    #     
    snr=t.test(vec2,vec1)$statistic
    wilcox=wilcox.test(vec1,vec2)$statistic
    rtest2_list=rtest2_pval(vec1,vec2)
    rtest2=min(unlist(rtest2_list))
    
    permt=perm.t.test(vec1,vec2)$statistic
    
    if(withBoot){
      snr.boot=boot(vec1,snr.boot,numPerm)
      snr.ci=boot.ci(snr.boot,type="bca")
      snr_sd=tail(diff(as.vector(snr.ci$bca)),1)/2
    }else{
      snr.boot=NULL
      snr.ci=NULL
      snr_sd=NA
    }
    return(list(up_avg=mean(ups),up_sd=sd(ups),
                down_avg=mean(downs),down_sd=sd(downs),
                snr=snr,snr_sd=snr_sd,
                ups=ups,downs=downs,
                up_avg_v=mean(ups_v),down_avg_v=mean(downs_v),
                #                 pH0_ups=mean(pH0_ups),pH0_downs=mean(pH0_downs),
                #                 pH0_ups_sd=sd(pH0_ups),pH0_downs_sd=sd(pH0_downs),
                wilcox=wilcox,rtest2=rtest2,permt=permt
    ))
  })
  
  ups=sapply(allres,function(x) x$up_avg)
  downs=sapply(allres,function(x) x$down_avg)
  
  ups_v=sapply(allres,function(x) x$up_avg_v)
  downs_v=sapply(allres,function(x) x$down_avg_v)
  
  
  ups_all=sapply(allres,function(x) x$ups)
  downs_all=sapply(allres,function(x) x$downs)
  
  #   pH0_ups=sapply(allres,function(x) x$pH0_ups)
  #   pH0_ups_sd=sapply(allres,function(x) x$pH0_ups_sd)
  #   pH0_downs=sapply(allres,function(x) x$pH0_downs)
  #   pH0_downs_sd=sapply(allres,function(x) x$pH0_downs_sd)
  
  wilcox=sapply(allres,function(x) x$wilcox)
  
  rtest2=sapply(allres,function(x) x$rtest2)

  permt=sapply(allres,function(x) x$permt)
  
  up_avg_avg=mean(ups)
  up_avg_sd=sd(ups)
  up_sd_avg=mean(sapply(allres,function(x) x$up_sd))
  up_sd_sd=sd(sapply(allres,function(x) x$up_sd))
  
  down_avg_avg=mean(downs)
  down_avg_sd=sd(downs)
  down_sd_avg=mean(sapply(allres,function(x) x$down_sd))
  down_sd_sd=sd(sapply(allres,function(x) x$down_sd))
  
  snrs=sapply(allres,function(x) x$snr)
  snr_avg_avg=mean(snrs)
  snr_avg_sd=sd(snrs)
  snr_sd_avg=mean(sapply(allres,function(x) x$snr_sd))
  snr_sd_sd=sd(sapply(allres,function(x) x$snr))
  return(list(up_avg_avg=up_avg_avg,up_avg_sd=up_avg_sd,
              up_sd_avg=up_sd_avg,up_sd_sd=up_sd_sd,
              down_avg_avg=down_avg_avg,down_avg_sd=down_avg_sd,
              down_sd_avg=down_sd_avg,down_sd_sd=down_sd_sd,
              snr_avg_avg=snr_avg_avg,snr_avg_sd=snr_avg_sd,
              snr_sd_avg=snr_sd_avg,snr_sd_sd=snr_sd_sd,
              ups=ups,downs=downs,
              ups_v=ups_v,downs_v=downs_v,
              snrs=snrs,
              #              ups_all=ups_all,downs_all=downs_all,
              N=N,Nsamples=numPerm,Navg=Navg,
              #               pH0_ups=pH0_ups,pH0_downs=pH0_downs,
              wilcox=wilcox,rtest2=rtest2,permt=permt
  ))
}

