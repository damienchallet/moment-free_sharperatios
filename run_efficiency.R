source('libdc.R')
source('libpermutations.R')

Ns=126
Navg=10000
numPerm=1000
SNRmin=0.001
SNRs=SNRmin*(100^(1/20))^c(0:15)  #from 0.01 to 0.1
GaussS=FALSE#rev(c(TRUE,FALSE))

forceRUN=FALSE
dirSave="precomputed/efficiency"
dir.create(dirSave,recursive = TRUE,showWarnings = FALSE)

for(N in Ns){
  for(Gauss in GaussS){
    if(Gauss){
      what="Gauss"
    }else{
      what="Student"
    }
    for(i in 1:1){
      allres=lapply(toList(SNRs),function(SNR){
        print(paste(what,SNR))
        filename=paste(dirSave,"/allres_","N",N,"_Navg",Navg,"_numPerm",numPerm,"_Gauss",Gauss,"_",i,"_SNR",SNR,".rds",sep="")
        print(filename)
        if(file.exists(filename) && file.info(filename)$size>4000 && !forceRUN){
          print("  already computed, skipping")
          return(NULL)
        }
        myres=sampleEstimates(N=N, Navg=Navg, numPerm=numPerm, SNR=SNR, dtype=what,mylapply=mclapply)
        saveRDS(myres,file=filename)
        return(NULL)
      })
      save(allres,file=filename)
      rm(allres)
      gc()
    }
  }
}