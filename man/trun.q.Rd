\name{trun.q}
\alias{trun.q}

\title{Truncated Inverse Cumulative Density Function of a gamlss.family Distribution}
\description{
 Creates a function to produce the inverse of a truncated cumulative density function generated from a current GAMLSS family distribution.
  
  For continuous distributions left truncation at  3 means that the random variable can take the value 3.
For discrete distributions left truncation at 3 means that the random variable can take values from 4 onwards.
This is the same for right truncation. Truncation at 15 for a discrete variable means that 
15 and greater values are not allowed but for continuous variable it mean values greater that 15 are not allowed 
(so 15 is a possible value). 
}
\usage{
trun.q(par, family = "NO", type = c("left", "right", "both"),
        varying = FALSE, ...)
}

\arguments{
   \item{par}{ a vector with one (for \code{"left"} or \code{"right"} truncation) or two elements for \code{"both"}. When the argument \code{varying = TRUE} then \code{par} can be a vector or a matrix with two columns respectively.}
  \item{family}{ a \code{gamlss.family} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlss()} can be found in \code{gamlss.family}.
   Functions such as BI() (binomial) produce a family object. }
  \item{type}{ whether \code{left}, \code{right} or in \code{both} sides truncation is required, (left is the default)  }
    \item{varying}{whether the truncation  varies for diferent observations. This can be usefull in regression analysis. If \code{varying = TRUE} then \code{par} should be an n-length vector for \code{type} equal \code{"left"} and \code{"right"} and an n by 2 matrix for \code{type="both"}}
  \item{\dots}{for extra arguments }
}

\value{
Returns a q family function 
}
\references{Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. \doi{10.1201/9780429298547}
	An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \doi{10.18637/jss.v023.i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC. \doi{10.1201/b21973}

Stasinopoulos, M. D., Rigby, R. A., and De Bastiani F., (2018) GAMLSS: a distributional regression approach, \emph{Statistical Modelling}, Vol. \bold{18}, pp, 248-273, SAGE Publications Sage India: New Delhi, India.


Stasinopoulos, M.D., Kneib, T., Klein, N., Mayr, A. and Heller, G.Z., (2024). \emph{Generalized Additive Models for Location, Scale and Shape: A Distributional Regression Approach, with Applications} (Vol. \bold{56}). Cambridge University Press.

(see also  \url{https://www.gamlss.com/}).}
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@gre.ac.uk} and Bob Rigby}

\seealso{  \code{\link{trun.d}}, \code{\link{trun.q}}, \code{\link{trun.r}}, \code{\link{gen.trun}} }
\examples{
# trucated q  continuous function
# continuous
#----------------------------------------------------------------------------------------
# left
test1<-trun.q(par=c(0), family="TF", type="left")
test1(.6)
qTF(pTF(0)+0.6*(1-pTF(0)))
#----------------------------------------------------------------------------------------
# right
test2 <- trun.q(par=c(10), family="BCT", type="right")
test2(.6)
qBCT(0.6*pBCT(10))
#----------------------------------------------------------------------------------------
# both
test3<-trun.q(par=c(-3,3), family="TF", type="both")
test3(.6)
qTF(0.6*(pTF(3)-pTF(-3))+pTF(-3))
#----------------------------------------------------------------------------------------
# varying  par
#----------------------------------------------------------------------------------------
# left
test7<-trun.q(par=c(0,1,2), family="TF", type="left", varying=TRUE)
test7(c(.5,.5,.6))
qTF(pTF(c(0,1,2))+c(.5,.5,.6)*(1-pTF(c(0,1,2))))
#---------------------------------------------------------------------------------------
# right
test9 <- trun.q(par=c(10,11,12), family="BCT", type="right", varying=TRUE)
test9(c(.5,.5,.6))
qBCT(c(.5,.5,.6)*pBCT(c(10,11,12)))
#----------------------------------------------------------------------------------------
# both
test10<-trun.q(par=cbind(c(0,1,2), c(10,11,12)), family="TF", type="both", varying=TRUE)
test10(c(.5, .5, .7))
qTF(c(.5, .5, .7)*(pTF(c(10,11,12))-pTF(c(0,1,2)))+pTF(c(0,1,2)))
#----------------------------------------------------------------------------------------
# FOR DISCRETE DISTRIBUTIONS
# trucated q function
# left
test4<-trun.q(par=c(0), family="PO", type="left")
test4(.6)
qPO(pPO(0)+0.6*(1-pPO(0)))
# varying
test41<-trun.q(par=c(0,1,2), family="PO", type="left", varying=TRUE)
test41(c(.6,.4,.5))
qPO(pPO(c(0,1,2))+c(.6,.4,.5)*(1-pPO(c(0,1,2))))
#----------------------------------------------------------------------------------------
# right
test5 <- trun.q(par=c(10), family="NBI", type="right")
test5(.6)
qNBI(0.6*pNBI(10))
test5(.6, mu=10, sigma=2)
qNBI(0.6*pNBI(10, mu=10, sigma=2), mu=10, sigma=2)
# varying
test51 <- trun.q(par=c(10, 11, 12), family="NBI", type="right", varying=TRUE)
test51(c(.6,.4,.5))
qNBI(c(.6,.4,.5)*pNBI(c(10, 11, 12)))
test51(c(.6,.4,.5), mu=10, sigma=2)
qNBI(c(.6,.4,.5)*pNBI(c(10, 11, 12), mu=10, sigma=2), mu=10, sigma=2)
#----------------------------------------------------------------------------------------
# both
test6<-trun.q(par=c(0,10), family="NBI", type="both")
test6(.6)
qNBI(0.6*(pNBI(10)-pNBI(0))+pNBI(0))
# varying 
test61<-trun.q(par=cbind(c(0,1,2), c(10,11,12)), family="NBI", type="both", varying=TRUE)
test61(c(.6,.4,.5))
qNBI(c(.6,.4,.5)*(pNBI(c(10,11,12))-pNBI(c(0,1,2)))+pNBI(c(0,1,2)))
#----------------------------------------------------------------------------------------
}
\keyword{distribution}
\keyword{regression}
