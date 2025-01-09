\name{sharpeRratio-package}
\alias{sharpeRratio-package}
\alias{sharpeRratio}
\docType{package}
\title{A moment-free estimator of the Sharpe (signal-to-noise) ratio
}
\description{
This package implements a new estimator of Sharpe ratios that does not rely on the computation of any moment, despite the fact that its usual definition involves at least the first two moments (average and standard deviation).  
}
\details{

An implementation of an alternative method to measure Sharpe ratios, i.e. signal-to-noise ratios in time series with heavy-tailed increments. The method itself does not require the computation of any moment as it is based on counting the number of records of the cumulative sum of the increments. When increments are known to be Gaussian, the usual estimator has to be used. However, when the increments are heavy-tailed, the new estimator is more precise (efficient). Note that the increments are assumed to be i.i.d. Note also that the new estimator is almost as efficient as the usual Sharpe ratio for Gaussian variables.
}
\author{
Damien Challet
Maintainer: Damien Challet
}
\references{Challet, D. (2017). Sharper asset ranking from total drawdown durations. Applied Mathematical Finance, 24(1), 1-22.}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
%\seealso{
%}
%\examples{
%% ~~ simple examples of the most important functions ~~
%}
