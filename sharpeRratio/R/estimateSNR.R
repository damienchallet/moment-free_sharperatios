#' @importFrom stats var quantile
#' @importFrom ghyp fit.tuv coef
#' @import Rcpp
#' @importFrom methods is
#' @useDynLib sharpeRratio


a_full=function(r0){
  val=(r0<0.7)*a(r0)+(r0>=0.7)*(r0<=0.997)*a_medium(r0)+(r0>0.997)*exp(2.877)/sqrt(252)/(1-r0)^0.163
  return(val)
}


f_full=function(x){    # splines are by definition unstable
  return(ifelse(x<3,f(x),f(3)))
}

b=function(r0,N){
  r0Nalpha=r0*N^.42
  return(exp(f_full(r0Nalpha))*r0^1.6)
}
  

theta=function(r0,N,nu,nu_fixed){

  myEstimate=a_full(abs(r0))-b(abs(r0),N)*nu^(-1.5)
  myEstimate=myEstimate * (sign(myEstimate)>0)                   #if(spline(abs(r0)))<0, then set SNR to 0 (this may occur when SNR is very small)
  myEstimate=myEstimate * sign(r0)                               #then apply sign of r0
  
  if(nu_fixed){
    myEstimate=myEstimate*(1-(exp(3.7726-6.2661*log(nu))-exp(-0.36538*nu-1.58686)))*(1-0.009248)
  }else{
    myEstimate=myEstimate*((1-(exp(4.4635-6.9240*log(nu))-exp(-0.18441*nu-3.18909)))*(1+exp(-log(nu)*3.2+0.9)-0.009248))
  }
  
  return(myEstimate)
}

Test.N = function(x) {    #Jelito, D., & Pitera, M. (2018). arXiv preprint arXiv:1811.05464.
  x=x[is.finite(x)]
  n = length(x)
  q1 = quantile(x,0.2 )
  q2 = quantile(x,0.8 )
  low = x [ x <= q1 ]
  med = x [ x > q1 & x < q2 ]
  high = x [ x >= q2 ]
  N = var(low)+var(high)-2*var(med)
  N = N*sqrt(n)/(var(x)*1.8)
  return(N) 
}


#'computes the signal-to-noise ratio
#' @export
#' @param x A (non-empty) numeric vector of data values.
#' @param numPerm The number of permutations (or shuffling) of the order of the sample values. By default set to \code{min(100,3 log(length(x)))}.
#' @param nu the Student t-distribution tail exponent of the sample data (if know). By default: NA. If set to NA, the tail exponent of the data is obtained from fit to a Student t-distribution. If NA, nu is estimated.
#' @param quantiles a vector of the lower and upper quantile needed to compute the confidence interval (use only if nu is known).
#' @return a list element \itemize{
#' \item SNR The signal-to-noise ratio. To have something comparable with a t-statistics, multiply by \code{sqrt(length(x))}. To have a Sharpe ratio, multiply by the correct factor (\code{sqrt(252)}) for daily returns)
#' \item SNR.ci The 95% confidence interval if $\nu$ is known. If not, do not trust this.
#' \item nu The fitted Student t-distribution tail exponent.
#' \item R0bar The number of upper records minus the number of lower records of the cumulated sum of \code{x}.
#' \item N The length of the vector \code{x}. It may be smaller than the input length if x contains NAs.}
#' @examples
#'   x <- rt(100,3)/sqrt(3)+0.05  #some Student-t distributed synthetic price log-returns
#'   estimateSNR(x)    

estimateSNR=function(x,numPerm=NA,nu=NA,quantiles=c(0.05,0.95)){
  if(is.na(nu))
    nu_fixed=FALSE
  else
    nu_fixed=TRUE
  
  if(length(x)<1){
    return(NULL)
  }
  x=x[is.finite(x)]
  
  N=length(x)   # note: this is the length of x[is.finite(x)], which removes NAs and Inf
  if(N<1){
    return(NULL)
  }
  
  if(!is.finite(numPerm)){
    numPerm=min(100,ceiling(3*log(N)))
  }

  x=as.numeric(x)
  if(is.na(nu) || nu<0 ){    # nu not given or the given nu has an absurd value
    if(abs(Test.N(x))<3 || is.infinite(nu)){    # Gaussian
      nu=1e13
      
    }else{
      myfit=tryCatch(ghyp::fit.tuv(x,silent=TRUE,nu=6),error=function(e){return(NA)})
      if(is(myfit,"ghyp")){
        nu=1e13  #fit.tuv fails on Gaussian data. 
      }else{
        nu=coef(myfit)$nu
      }
    }
  }
  R0bar_res=computeR0bar(x,numPerm = numPerm)   #computes R_+ - R_- over numPerm permutations

  R0bar=R0bar_res$mean
  if(abs(R0bar)<1){
    R0bar=0
    return(list(nu=nu,SNR=0,R0bar=0,N=length(x)))
  }

  SNR=theta(R0bar/N,N,nu,nu_fixed)                #splines have been calibrated for positive SNRs, hence abs(R0bar)

  SNR_q1=theta(R0bar_res$q1/N,N,nu,nu_fixed)
  SNR_q2=theta(R0bar_res$q2/N,N,nu,nu_fixed)

  
  return(list(SNR=SNR,SNR.ci=c(SNR_q1,SNR_q2),nu=nu,R0bar=R0bar,N=length(x)))
}

.onUnload <- function (libpath) {
  library.dynam.unload("sharpeRratio", libpath)
}
