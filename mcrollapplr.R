library(parallel)
library(digest)

mcrollapplyr <- function (x, width, fill=NA, FUN=NULL, FUN_postTreatment=NULL,saveDir=NA,mylapply=mclapply,skipMissingResults=FALSE,...) {
  if(is.null(dim(x))){ #x is a vector
    x=xts(as.matrix(x),index(x))
  }
  if(nrow(x)<width){
    return(NULL)
  }
  starting.indexes <- seq(width,nrow(x),1)
  numdata = ncol(x)
  if(!is.na(saveDir)){
    dir.create(saveDir,recursive = TRUE,showWarnings = FALSE)
  }
  
  doCompute=TRUE
  
  result = mylapply(starting.indexes, function(it) {
#    print(index(x)[it])
    myx=x[seq(it-width+1, it)]
    if(!is.na(saveDir)){
      myhash=digest(myx)
      filename=paste(saveDir,"/",myhash,".rda",sep="")
      if(file.exists(filename) && file.info(filename)$size>5000){
        load(filename)
        doCompute=FALSE
      }else{
        doCompute=TRUE
      }
    }
    if(doCompute){
      if(skipMissingResults){
        return(NULL)
      }
      myres=FUN(coredata(myx), ...)
      if(!is.na(saveDir)){
        save(myres,file=filename)
      }
      
      if(!is.null(FUN_postTreatment)){
        myres=FUN_postTreatment(myres)
      }
    }
    return(myres)
    
  })
  
  areNull=unlist(lapply(result,is.null))
  if(all(areNull)){
    return(NULL)
  }
  if(all(!unlist(lapply(result,is.numeric)))){
    names(result)=index(x)[starting.indexes]
    return(result)
  }else{
    result = do.call(rbind, result)
    
    result = xts(
      rbind(
        matrix(rep(fill,width-1),ncol=numdata,nrow=width-1),
        result )
      ,order.by=time(x))
  }
  
  return(result)
}



