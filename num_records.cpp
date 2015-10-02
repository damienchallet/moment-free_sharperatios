// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]
// 
#include <algorithm>    // std::random_shuffle
//#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

Function printR("print");
Function pasteR("paste");
Function meanR("mean");
Function wilcox("wilcox.test");

// [[Rcpp::depends(RcppArmadillo)]]


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

// [[Rcpp::export]]

IntegerVector records_up(NumericVector x){
  int num_records=1;
  IntegerVector cumNumRecs(x.size());
  double xmax=x[0];
  cumNumRecs[0]=1;
  int i;
  for(i=1;i<x.size();i++){
    if(x[i]>xmax){
      num_records++;
      xmax=x[i];
    }
    cumNumRecs[i]=num_records;
  }
  return cumNumRecs;
}

// [[Rcpp::export]]

int age_record_max(NumericVector x){
  int num_records=1;
  int age_record=0;
  int age_record_max=0;
  double xmax=x[0];
  for(int i=1;i<x.size();i++){
    if(x[i]>xmax){
      xmax=x[i];
      age_record=0;
    }else{
      age_record++;
      age_record_max=std::max(age_record,age_record_max);
    }
  }
  return age_record_max;
}



// [[Rcpp::export]]

IntegerVector records_down(NumericVector x){
  int num_records=1;
  IntegerVector cumNumRecs(x.size());
  double xmin=x[0];
  cumNumRecs[0]=1;
  int i;
  for(i=1;i<x.size();i++){
    if(x[i]<xmin){
      num_records++;
      xmin=x[i];
    }
    cumNumRecs[i]=num_records;
  }
  return cumNumRecs;
}


// [[Rcpp::export]]

List computeUpsDownsRandomPerms(NumericVector vec,int numPerm){
  IntegerVector ups(numPerm);
  IntegerVector downs(numPerm);
  int N=vec.size();
  NumericVector cumvec(N);
  double R0bar=0.,Rupbar=0.,Rdownbar=0.;
  for(int perm=0;perm<numPerm;perm++){
    std::random_shuffle (vec.begin(), vec.end() );
    cumvec[0]=vec[0];
    for(int j=1;j<N;j++){
      cumvec[j]=cumvec[j-1]+vec[j];
    }
    ups[perm]=records_up(cumvec)[N-1];
    downs[perm]=records_down(cumvec)[N-1];
    Rupbar+=ups[perm];
    Rdownbar+=downs[perm];
  }
  Rupbar/=numPerm;
  Rdownbar/=numPerm;
  R0bar=Rupbar-Rdownbar;
  return List::create(Named("up")=ups,Named("down")=downs,
		      Named("Rupbar")=Rupbar,Named("Rdownbar")=Rdownbar,Named("R0bar")=R0bar);
}



// [[Rcpp::export]]

List computeUpsDownsDiffTwoSamplesRandomPerms(NumericVector vec1, NumericVector vec2, int numPerm){
  IntegerVector ups(numPerm);
  IntegerVector downs(numPerm);
  int N=min(NumericVector::create(vec1.size(),vec2.size()));
  NumericVector cumvec(N);
  for(int i=0;i<numPerm;i++){
    std::random_shuffle (vec1.begin(), vec1.end() );
    std::random_shuffle (vec2.begin(), vec2.end() );
    cumvec[0]=vec1[0]-vec2[0];
    for(int j=1;j<N;j++){
      cumvec[j]=cumvec[j-1]+vec1[j]-vec2[j];
    }
    ups[i]=records_up(cumvec)[N-1];
    downs[i]=records_down(cumvec)[N-1];
  }
//   IntegerVector R0bar=ups-downs;
  return List::create(Named("up")=ups,Named("down")=downs);//,Named("R0bar")=R0bar);
}

// List computeUpsDownsDiffTwoSamplesSameRandomPerms(NumericVector vec1, NumericVector vec2, int numPerm){
//   IntegerVector ups(numPerm);
//   IntegerVector downs(numPerm);
//   int N=min(NumericVector::create(vec1.size(),vec2.size()));
//   
//   NumericVector cumvec(N);
//   for(int i=0;i<numPerm;i++){
//     std::random_shuffle (vec1.begin(), vec1.end() );
//     std::random_shuffle (vec2.begin(), vec2.end() );
//     cumvec[0]=vec1[0]-vec2[0];
//     for(int j=1;j<N;j++){
//       cumvec[j]=cumvec[j-1]+vec1[j]-vec2[j];
//     }
//     ups[i]=records_up(cumvec)[N-1];
//     downs[i]=records_down(cumvec)[N-1];
//   }
//   double R0bar=mean(ups-downs);
//   return List::create(Named("up")=ups,Named("down")=downs,Named("R0bar")=R0bar);
// }
// 

// [[Rcpp::export]]

List computeUpsDownsTwoSamplesRandomPerms(NumericVector vec1, NumericVector vec2, int numPerm){
  IntegerVector ups1(numPerm);
  IntegerVector downs1(numPerm);
  IntegerVector ups2(numPerm);
  IntegerVector downs2(numPerm);
  
  int N=std::min(vec1.size(),vec2.size());
  
  NumericVector cumvec1(N),cumvec2(N);
  for(int perm=0;perm<numPerm;perm++){
    std::random_shuffle (vec1.begin(), vec1.end() );
    std::random_shuffle (vec2.begin(), vec2.end() );

    cumvec1[0]=vec1[0];
    for(int j=1;j<N;j++){
      cumvec1[j]=cumvec1[j-1]+vec1[j];
    }
    cumvec2[0]=vec2[0];
    for(int j=1;j<N;j++){
      cumvec2[j]=cumvec2[j-1]+vec2[j];
    }

    ups1[perm]=records_up(cumvec1)[N-1];
    downs1[perm]=records_down(cumvec1)[N-1];
    ups2[perm]=records_up(cumvec2)[N-1];
    downs2[perm]=records_down(cumvec2)[N-1];
  }
  double R0bar1=0;
  double R0bar2=0;
  for(int perm=0;perm<numPerm;perm++){
    R0bar1+=ups1[perm]-downs1[perm];
    R0bar2+=ups2[perm]-downs2[perm];
  }
  R0bar1/=numPerm;
  R0bar2/=numPerm;
  double R0bar12=R0bar1-R0bar2;
  return List::create(Named("up1")=ups1,Named("down1")=downs1,Named("R0bar")=R0bar1,
                            Named("up2")=ups2,Named("down2")=downs2,Named("R0bar")=R0bar2,
                                  Named("R0bar12")=R0bar12);
}




// [[Rcpp::export]]

double computeR0bar(NumericVector vec,int numPerm=100){
  int N=vec.size();
  NumericVector cumvec(N);
  double R0bar=0.;
  for(int perm=0;perm<numPerm;perm++){
    std::random_shuffle (vec.begin(), vec.end() );
    cumvec[0]=vec[0];
    for(int j=1;j<N;j++){
      cumvec[j]=cumvec[j-1]+vec[j];
    }
    R0bar+=num_records_up(cumvec)-num_records_down(cumvec);
  }
  R0bar/=numPerm;
  return R0bar;
}

// [[Rcpp::export]]


double computeRTbarOld(NumericVector vec,int numPerm){
  int N=vec.size();
  NumericVector cumvec(N);
  double RTbar=0.;
  for(int perm=0;perm<numPerm;perm++){
    std::random_shuffle (vec.begin(), vec.end() );
    cumvec[0]=vec[0];
    for(int j=1;j<N;j++){
      cumvec[j]=cumvec[j-1]+vec[j];
    }
    RTbar+=age_record_max(cumvec);
  }
  RTbar/=numPerm;
  return RTbar;
}

// [[Rcpp::export]]


int computeRTbar(NumericVector vec,int numPerm){
  int N=vec.size();
  NumericVector cumvec(N);
  int RTbar=0.;
  for(int perm=0;perm<numPerm;perm++){
    std::random_shuffle (vec.begin(), vec.end() );
    cumvec[0]=vec[0];
    for(int j=1;j<N;j++){
      cumvec[j]=cumvec[j-1]+vec[j];
    }
    RTbar=std::max(RTbar,age_record_max(cumvec));
  }
//  RTbar/=numPerm;
  return RTbar;
}

// 
// double RplusTest(NumericVector x, int numPermutations){
//   List mylist=computeUpsDownsRandomPerms(x,numPermutations);
//   int N=x.size();
//   double mystat=((double)(mylist["Rupbar"])-sqrt(4/M_PI*N))/sqrt((2-4/M_PI)*N)/(1-1/1.25/sqrt(N));
//   NumericVector mystatV(1);
//   mystatV[0]=mystat;  
//   NumericVector pvalue=pnorm(mystatV,0.,1.);
//   return pvalue[0];
//  
// }

// [[Rcpp::export]]

double RminusTest(NumericVector x, int numPermutations){
  List mylist=computeUpsDownsRandomPerms(x,numPermutations);

  int N=x.size();
  double mystat=((double)(mylist["Rdownbar"])-sqrt(4/M_PI*N))/sqrt((2-4/M_PI)*N)/(1-1/1.25/sqrt(N));

  NumericVector mystatV(1);
  mystatV[0]=mystat;  
  NumericVector pvalue=pnorm(mystatV,0.,1.);
  return pvalue[0]; 
}

// [[Rcpp::export]]
double rstat(NumericVector x,int numPermutations=100){

  int N=x.size();
  double mystat=computeR0bar(x,numPermutations)/sqrt((2-4/M_PI)*N)/(1.658*(1-0.88/sqrt(N)));

//   NumericVector mystatV(1);
//   mystatV[0]=mystat;  
//   NumericVector pvalue=pnorm(mystatV,0.,1.);
//  return pvalue[0];
  return mystat;
}


// [[Rcpp::export]]
double rstatBIS(NumericVector x,int numPermutations=1000){
  List Rs=computeUpsDownsRandomPerms(x,numPermutations);
  IntegerVector ups=Rs["up"];
  IntegerVector downs=Rs["down"];
  int numPos=0;
  for(int i=0;i<numPermutations;i++){
    numPos+=(ups[i]>downs[i]?1:0);
  }
  return numPos*1./numPermutations;
}


// [[Rcpp::export]]

List computeDistrR0RandomPerms(NumericVector vec,int numPerm){
  IntegerVector R0s(numPerm);
  int N=vec.size();
  NumericVector cumvec(N);
  double R0bar=0.,Rupbar=0.,Rdownbar=0.;
  for(int perm=0;perm<numPerm;perm++){
    std::random_shuffle (vec.begin(), vec.end() );
    cumvec[0]=vec[0];
    for(int j=1;j<N;j++){
      cumvec[j]=cumvec[j-1]+vec[j];
    }
    R0s[perm]=records_up(cumvec)[N-1]-records_down(cumvec)[N-1];
  }
  
  return List::create(Named("R0s")=R0s);
}

// [[Rcpp::export]]
double rstat_nonstat(NumericVector x,double pvalmin=0.05,double numPerm=100) {
  int Ndiv2=x.size()/2;
  List R0sList1=computeDistrR0RandomPerms(x[seq(0,Ndiv2-1)],numPerm);          // split vector into two halves
  List R0sList2=computeDistrR0RandomPerms(x[seq(Ndiv2,x.size()-1)],numPerm);
  
  List wilcoxRes=wilcox(R0sList1["R0s"],R0sList2["R0s"]);              // same distrib of R0 ?
  double rstatWhole=rstat(x);
  double wilcoxDiff=wilcoxRes["p.value"];
  if(wilcoxDiff<pvalmin){
    rstatWhole=0.;
  }
  return rstatWhole;
}



///////////////////////////// two-sample part


arma::colvec conc(const arma::colvec & x, const arma::colvec & y) {
  return arma::join_cols(x, y); 
}

IntegerVector cppcp(NumericVector samp, NumericVector ref, IntegerVector ord) {
  int nobs = samp.size();
  IntegerVector ans(nobs);
  for (int i = 0, j = 0; i < nobs; ++i) {
    int ind(ord[i] - 1); // C++ uses 0-based indices
    double ssampi(samp[ind]);
    while (ref[j] < ssampi && j < ref.size()) ++j;
    ans[ind] = j;     // j is the 1-based index of the lower bound
  }
  return ans;
}

// [[Rcpp::export]]

//permutation test
List rtest2_pval(NumericVector x,NumericVector y,int numSamplesH0=10000,int numPermutations=100){
  //build H0
  NumericVector xy=as<NumericVector>(wrap(conc(x,y)));
  NumericVector RplusH0_X(numSamplesH0);
  NumericVector RplusH0_Y(numSamplesH0);
  NumericVector RminusH0_X(numSamplesH0);
  NumericVector RminusH0_Y(numSamplesH0);
  
  int xSize=x.size();
  int ySize=y.size();

  List records_X,records_Y;
  int i;
  
//   omp_set_num_threads(1);
// #pragma omp parallel for default(none) firstprivate(xy) private(i,records_X,records_Y) shared(numPermutations,numSamplesH0,RplusH0_X,RplusH0_Y,RminusH0_X,RminusH0_Y,xSize,ySize)
  for(i=0;i<numSamplesH0;i++){
    std::random_shuffle(xy.begin(),xy.end());
    records_X=computeUpsDownsRandomPerms(xy[seq(0,xSize-1)],numPermutations);             // the first xSize elements go to X
    records_Y=computeUpsDownsRandomPerms(xy[seq(ySize,xSize+ySize-1)],numPermutations);   // the rest goes to Y
    RplusH0_X[i]=records_X["Rupbar"];
    RplusH0_Y[i]=records_Y["Rupbar"];
    RminusH0_X[i]=records_X["Rdownbar"];
    RminusH0_Y[i]=records_Y["Rdownbar"];
  }
  
  std::sort(RplusH0_X.begin(), RplusH0_X.end());
  std::sort(RplusH0_Y.begin(), RplusH0_Y.end());
  std::sort(RminusH0_X.begin(), RminusH0_X.end());
  std::sort(RminusH0_Y.begin(), RminusH0_Y.end());
  
  records_X=computeUpsDownsRandomPerms(x,numPermutations);
  records_Y=computeUpsDownsRandomPerms(y,numPermutations);
  

  double pval_plus_X= (std::lower_bound(RplusH0_X.begin(),RplusH0_X.end(), (double)records_X["Rupbar"]) - RplusH0_X.begin())*1./numSamplesH0;
  double pval_plus_Y= (std::lower_bound(RplusH0_Y.begin(),RplusH0_Y.end(), (double)records_Y["Rupbar"]) - RplusH0_Y.begin())*1./numSamplesH0;
  
  double pval_minus_X= (std::lower_bound(RminusH0_X.begin(),RminusH0_X.end(), (double)records_X["Rdownbar"]) - RminusH0_X.begin())*1./numSamplesH0;
  double pval_minus_Y= (std::lower_bound(RminusH0_Y.begin(),RminusH0_Y.end(), (double)records_Y["Rdownbar"]) - RminusH0_Y.begin())*1./numSamplesH0;
  
    return List::create(Named("pval_plus_X")=pval_plus_X,
                        Named("pval_plus_Y")=pval_plus_Y,
                        Named("pval_minus_X")=pval_minus_X,
                        Named("pval_minus_Y")=pval_minus_Y);
}



