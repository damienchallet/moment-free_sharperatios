#include <algorithm>    // std::shuffle
#include <Rcpp.h>
#include <random>

using namespace Rcpp;

//' Computes the number of upper records of the cumulative sum of \code{x}
//'
//' @param x a vector of sample values
//' @return the number of upper records 
// [[Rcpp::export]]
int num_records_up(NumericVector x) {
  int num_records=1;
  double xmax=x[0];
  int i;
  for(i=1;i<x.size();i++){
    if(x[i]>xmax){
      num_records++;
      xmax=x[i];
    }
  }
  return num_records;
}

//' Computes the number of lower records of the cumulative sum of \code{x}
//'
//' @param x a vector of sample values
//' @return the number of lower records 
// [[Rcpp::export]]
int num_records_down(NumericVector x){
  int num_records=1;
  double xmin=x[0];
  int i;
  for(i=1;i<x.size();i++){
    if(x[i]<xmin){
      num_records++;
      xmin=x[i];
    }
  }
  return num_records;
}


//' Computes the average difference between the number of upper and lower records of the cumulative sum of the sample values.
//'
//' @param x a vector of sample values
//' @param numPerm the number of random permutations (or shuffles) of the sample value order
//' @param q1 a real number for computing the lower confidence interval  
//' @param q2 a real number for computing the upper confidence interval  
//' @return a list \itemize{
//' \item mean the average difference of upper and lower records of the cumulative sum of \code{x}
//' \item q1 the q1 quantile of the difference
//' \item q1 the q2 quantile of the difference}
//' 
// [[Rcpp::export]]
List computeR0bar(NumericVector x, int numPerm=100, double q1=0.025, double q2=0.975){
  std::random_device rng;
  std::mt19937 urng(rng());
  
  int N=x.size();
  NumericVector cumvec(N);
  double R0bar=0.;
  double R0;
  NumericVector R0s(numPerm);
  
  for(int perm=0;perm<numPerm;perm++){
    cumvec[0]=x[0];
    for(int j=1;j<N;j++){
      cumvec[j]=cumvec[j-1]+x[j];
    }
    R0=num_records_up(cumvec)-num_records_down(cumvec);
    R0s[perm]=R0;
    R0bar+=R0;
    std::shuffle (x.begin(), x.end(),urng);
  }
  
  std::sort(R0s.begin(), R0s.end());
  int idx;
  
  idx = std::floor(q1 * numPerm*1.);
  double R0_q1 = R0s[idx];
  idx = std::floor(q2 * numPerm*1.);
  double R0_q2 = R0s[idx];
  
  R0bar/=numPerm;
  List L=List::create(_["mean"]=R0bar, _["q1"]=std::min(R0_q1,R0_q2), _["q2"]=std::max(R0_q1,R0_q2));
  return L;
}

