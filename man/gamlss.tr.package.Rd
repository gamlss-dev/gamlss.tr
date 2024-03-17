\name{gamlss.tr-package}
\alias{gamlss.tr-package}
\alias{gamlss.tr}
\docType{package}

\title{
\packageTitle{gamlss.tr}
}
\description{
\packageDescription{gamlss.tr}
}

\details{
The DESCRIPTION file:
\packageDESCRIPTION{gamlss.tr}
\packageIndices{gamlss.tr}
}

\author{
\packageAuthor{gamlss.tr}

Maintainer: \packageMaintainer{gamlss.tr}
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}
\keyword{ package }
\keyword{ distribution }
\keyword{ regression }

\examples{
# generating a t-distribution from 0 to 100  	
gen.trun(par=c(0,100),family="TF", name="0to100", type="both")
op<-par(mfrow=c(2,2))
plot(function(x) dTF0to100(x, mu=80 ,sigma=20, nu=5), 0, 100, ylab="pdf")
plot(function(x) pTF0to100(x, mu=80 ,sigma=20, nu=5), 0, 100, ylab="cdf")
plot(function(x) qTF0to100(x, mu=80 ,sigma=20, nu=5), 0.01, .999, ylab="invcdf")
hist(s1<-rTF0to100(1000, mu=80 ,sigma=20, nu=5), ylab="hist", xlab="x", main="generated data")
par(op)
}
