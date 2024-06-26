\name{trun.d}
\alias{trun.d}

\title{Truncated Probability Density Function of a gamlss.family Distribution}
\description{
  Creates a truncated probability density function version from a current GAMLSS family distribution
  
  For continuous distributions left truncation at  3 means that the random variable can take the value 3.
For discrete distributions left truncation at 3 means that the random variable can take values from 4 onwards.
This is the same for right truncation. Truncation at 15 for a discrete variable means that 
15 and greater values are not allowed but for continuous variable it mean values greater that 15 are not allowed 
(so 15 is a possible value). 
}
\usage{
trun.d(par, family = "NO", type = c("left", "right", "both"),
       varying = FALSE, ...)
}

\arguments{
  \item{par}{ a vector with one (for \code{"left"} or \code{"right"} truncation) or two elements for \code{"both"}. When the argument \code{varying = TRUE} then \code{par} can be a vector or a matrix with two columns respectively.}
  \item{family}{a \code{gamlss.family} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlss()} can be found in \code{gamlss.family}.
   Functions such as BI() (binomial) produce a family object. }
  \item{type}{whether \code{left}, \code{right} or in \code{both} sides truncation is required, (left is the default).}
   \item{varying}{whether the truncation  varies for diferent observations. This can be usefull in regression analysis. If \code{varying = TRUE} then \code{par} should be an n-length vector for \code{type} equal \code{"left"} and \code{"right"} and an n by 2 matrix for \code{type="both"}}
  \item{\dots}{ for extra arguments }
}

\value{
Returns a d family function 
}
\references{ Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. \doi{10.1201/9780429298547}
	An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \doi{10.18637/jss.v023.i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC. \doi{10.1201/b21973}

Stasinopoulos, M. D., Rigby, R. A., and De Bastiani F., (2018) GAMLSS: a distributional regression approach, \emph{Statistical Modelling}, Vol. \bold{18}, pp, 248-273, SAGE Publications Sage India: New Delhi, India.


Stasinopoulos, M.D., Kneib, T., Klein, N., Mayr, A. and Heller, G.Z., (2024). \emph{Generalized Additive Models for Location, Scale and Shape: A Distributional Regression Approach, with Applications} (Vol. \bold{56}). Cambridge University Press.

(see also  \url{https://www.gamlss.com/}). }
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@gre.ac.uk} and Bob Rigby}

\seealso{  \code{\link{trun.p}}, \code{\link{trun.q}}, \code{\link{trun.r}}, \code{\link{gen.trun}} }
\examples{
#------------------------------------------------------------------------------------------
# continuous distribution 
# left truncation 
test1<-trun.d(par=c(0), family="TF", type="left")
test1(1)
dTF(1)/(1-pTF(0))
if(abs(test1(1)-(dTF(1)/pTF(0)))>0.00001) stop("error in left trucation")
test1(1, log=TRUE)
log(dTF(1)/(1-pTF(0)))
if(abs(test1(1, log=TRUE)-log(dTF(1)/pTF(0)))>0.00001) 
                   stop("error in left trucation")
integrate(function(x) test1(x, mu=-2, sigma=1, nu=1),0,Inf) 
# the pdf is defined even with negative mu
integrate(function(x) test1(x, mu=0, sigma=10, nu=1),0,Inf) 
integrate(function(x) test1(x, mu=5, sigma=5, nu=10),0,Inf)
plot(function(x) test1(x, mu=-3, sigma=1, nu=1),0,10)
plot(function(x) test1(x, mu=3, sigma=5, nu=10),0,10)
#----------------------------------------------------------------------------------------
# right truncation
test2<-trun.d(par=c(10), family="BCT", type="right")
test2(1)
dBCT(1)/(pBCT(10))
#if(abs(test2(1)-(dBCT(1)/pBCT(10)))>0.00001) stop("error in right trucation")
test2(1, log=TRUE)
log(dBCT(1)/(pBCT(10)))
if(abs(test2(1, log=TRUE)-log(dBCT(1)/(pBCT(10))))>0.00001) 
                   stop("error in right trucation")
integrate(function(x) test2(x, mu=2, sigma=1, nu=1),0,10) 
integrate(function(x) test2(x, mu=2, sigma=.1, nu=1),0,10) 
integrate(function(x) test2(x, mu=2, sigma=.1, nu=10),0,10) 
plot(function(x) test2(x, mu=2, sigma=.1, nu=1),0,10)
plot(function(x) test2(x, mu=2, sigma=1, nu=1),0,10)
#-----------------------------------------------------------------------------------------
# both left and right truncation
test3<-trun.d(par=c(-3,3), family="TF", type="both")
test3(0)
dTF(0)/(pTF(3)-pTF(-3))
if(abs(test3(0)-dTF(0)/(pTF(3)-pTF(-3)))>0.00001) 
              stop("error in right trucation")
test3(0, log=TRUE)
log(dTF(0)/(pTF(3)-pTF(-3)))
if(abs(test3(0, log=TRUE)-log(dTF(0)/(pTF(3)-pTF(-3))))>0.00001) 
            stop("error in both trucation")
plot(function(x) test3(x, mu=0, sigma=1, nu=1),-3,3)
integrate(function(x) test3(x, mu=2, sigma=1, nu=1),-3,3)
#-----------------------------------------------------------------------------------------
# discrete distribution
# left 
# Poisson truncated at zero means zero is excluded
test4<-trun.d(par=c(0), family="PO", type="left")
test4(1)
dPO(1)/(1-pPO(0))
if(abs(test4(1)-dPO(1)/(1-pPO(0)))>0.00001) stop("error in left trucation")
test4(1, log=TRUE)
log(dPO(1)/(1-pPO(0)))
if(abs(test4(1, log=TRUE)-log(dPO(1)/(1-pPO(0))))>0.00001) 
               stop("error in left trucation")
sum(test4(x=1:20, mu=2)) # 
sum(test4(x=1:200, mu=80)) #
plot(function(x) test4(x, mu=20), from=1, to=51, n=50+1, type="h") # pdf 
# right
# right truncated at 10 means 10 is excluded
test5<-trun.d(par=c(10), family="NBI", type="right")
test5(2)
dNBI(2)/(pNBI(9))
if(abs(test5(1)-dNBI(1)/(pNBI(9)))>0.00001) stop("error in right trucation")
test5(1, log=TRUE)
log(dNBI(1)/(pNBI(9)))
if(abs(test5(1, log=TRUE)-log(dNBI(1)/(pNBI(9))))>0.00001) stop("error in right trucation")
sum(test5(x=0:9, mu=2,   sigma=2)) # 
sum(test5(x=0:9, mu=300, sigma=5)) # can have mu > parameter
plot(function(x) test5(x, mu=20, sigma=3), from=0, to=9, n=10, type="h") # pdf
plot(function(x) test5(x, mu=300, sigma=5), from=0, to=9, n=10, type="h") # pdf
#----------------------------------------------------------------------------------------
# both
test6<-trun.d(par=c(0,10), family="NBI", type="both")
test6(2)
dNBI(2)/(pNBI(9)-pNBI(0))
if(abs(test6(2)-dNBI(2)/(pNBI(9)-pNBI(0)))>0.00001) 
        stop("error in right trucation")
test6(1, log=TRUE)
log(dNBI(1)/(pNBI(9)-pNBI(0)))
if(abs(test6(1, log=TRUE)-log(dNBI(1)/(pNBI(9)-pNBI(0))))>0.00001) 
  stop("error in right trucation")
sum(test6(x=1:9, mu=2,   sigma=2)) # 
sum(test6(x=1:9, mu=100, sigma=5)) # can have mu > parameter
plot(function(x) test6(x, mu=20, sigma=3), from=1, to=9, n=9, type="h") # pdf
plot(function(x) test6(x, mu=300, sigma=.4), from=1, to=9, n=9, type="h") # pdf
#------------------------------------------------------------------------------------------
# now try when the trucated points varies for each observarion
# this will be appropriate for regression models only 
# continuous
#----------------------------------------------------------------------------------------
# left truncation
test7<-trun.d(par=c(0,1,2), family="TF", type="left", varying=TRUE)
test7(c(1,2,3))
dTF(c(1,2,3))/(1-pTF(c(0,1,2)))
test7(c(1,2,3), log=TRUE)
#----------------------------------------------------------------------------------------
# right truncation
test8<-trun.d(par=c(10,11,12), family="BCT", type="right", varying=TRUE)
test8(c(1,2,3))
dBCT(c(1,2,3))/(pBCT(c(10,11,12)))
test8(c(1,2,3), log=TRUE)
#----------------------------------------------------------------------------------------
# both left and right truncation
test9<-trun.d(par=cbind(c(0,1,2),c(10,11,12) ), family="TF", type="both", 
             varying=TRUE)
test9(c(1,2,3))
dTF(c(1,2,3))/ (pTF(c(10,11,12))-pTF(c(0,1,2)))
test3(c(1,2,3), log=TRUE)
#----------------------------------------------------------------------------------------
# discrete
# left
test10<-trun.d(par=c(0,1,2), family="PO", type="left", varying=TRUE)
test10(c(1,2,3))
dPO(c(1,2,3))/(1-pPO(c(0,1,2)))
# right
test11<-trun.d(par=c(10,11,12), family="NBI", type="right", varying=TRUE)
test11(c(1,2,3))
dNBI(c(1,2,3))/pNBI(c(9,10,11))
# both
test12<-trun.d(par=rbind(c(0,10), c(1,11), c(2,12)), family="NBI", type="both", varying=TRUE)
test12(c(2,3,4))
dNBI(c(2,3,4))/(pNBI(c(9,10,11))-pNBI(c(0,1,2)))
}
\keyword{distribution}
\keyword{regression}
