source('libH0.R')

Ns=252
Navg=10000
numPerm=1000
SNRmin=0.001
Npoints=30
SNRs=SNRmin*(1000^(1/Npoints))^c(0:Npoints)  #from SNRmin to 100*SNRmin
SNRmax=max(SNRs)

dtypes=c("Student")#,"Student")
nus=seq(2.5,10,0.1) #3

forceRUN=FALSE
dirSave="precomputed/efficiency"
dir.create(dirSave,recursive = TRUE,showWarnings = FALSE)

mylapply=mclapply

allR0s=list()
for(N in Ns){
  for(dtype in dtypes){
      if(dtype!="Student"){
        oldnus=nus
        nus=""
      }else{
        if(length(nus)==1 && nus==""){
          nus=oldnus
        }
      }
    for(nu in nus){
      
      allres=lapply(toList(SNRs),function(SNR){
        print(paste(dtype,SNR))
        filename=paste(dirSave,"/myres_","N",N,"_Navg",Navg,"_numPerm",numPerm,"_dtype:",dtype,nu,"_SNR",SNR,".rds",sep="")
        print(filename)
        if(file.exists(filename) && file.info(filename)$size>4000 && !forceRUN){
          print("  already computed, skipping")
          R0s=readRDS(filename)
          return(mean(R0s))
        }else{
          R0s=unlist(mylapply(seq_len(Navg),function(x){
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
            myres=computeR0bar(vec0,numPerm = numPerm)
            return(myres)
          }))
          saveRDS(R0s,file=filename)
          return(mean(R0s))
        }
      })
      filename=paste(dirSave,"/allres_","N",N,"_Navg",Navg,"_numPerm",numPerm,"_dtype:",dtype,nu,"_SNR",SNRmin,"-",SNRmax,".rds",sep="")
      saveRDS(allres,file=filename)
      
      
      allR0s=unlist(allres)

      myspline=splinefun(c(0,allR0s),c(0,SNRs))
      
      filename=paste(dirSave,"/spline_","N",N,"_Navg",Navg,"_numPerm",numPerm,"_dtype:",dtype,nu,"_SNR",SNRmin,"-",SNRmax,".rds",sep="")
      saveRDS(myspline,file=filename)
#      plot(allR0s,SNRs,col=1,log='xy')
#      lines(allR0s,myspline(allR0s),col=2)
    }
  }
  #mypattern=paste0(dirSave,"/allres_","N",N,"_Navg",Navg,"_numPerm",numPerm,"_Gauss",Gauss,"_",i,"_SNR.*")
  
}

  