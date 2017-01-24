#include <algorithm>    // std::random_shuffle
#include <Rcpp.h>
// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

using namespace Rcpp;


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


// [[Rcpp::export]]

double computeR0bar(NumericVector vec,int numPerm=100){
  int N=vec.size();
  NumericVector cumvec(N);
  double R0bar=0.;
  for(int perm=0;perm<numPerm;perm++){
    cumvec[0]=vec[0];
    for(int j=1;j<N;j++){
      cumvec[j]=cumvec[j-1]+vec[j];
    }
    R0bar+=num_records_up(cumvec)-num_records_down(cumvec);
    std::random_shuffle (vec.begin(), vec.end(),randWrapper);
  }
  R0bar/=numPerm;
  return R0bar;
}

