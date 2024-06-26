\name{trun.p}
\alias{trun.p}

\title{Truncated Cumulative Density Function of a gamlss.family Distribution}
\description{
 Creates a truncated cumulative density function version from a current GAMLSS family distribution.
 
 For continuous distributions left truncation at  3 means that the random variable can take the value 3.
For discrete distributions left truncation at 3 means that the random variable can take values from 4 onwards.
This is the same for right truncation. Truncation at 15 for a discrete variable means that 
15 and greater values are not allowed but for continuous variable it mean values greater that 15 are not allowed 
(so 15 is a possible value). 
}
\usage{
trun.p(par, family = "NO", type = c("left", "right", "both"),
        varying = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{par}{ a vector with one (for \code{"left"} or \code{"right"} truncation) or two elements for \code{"both"}. When the argument \code{varying = TRUE} then \code{par} can be a vector or a matrix with two columns respectively.}
  \item{family}{ a \code{gamlss.family} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlss()} can be found in \code{gamlss.family}.
   Functions such as BI() (binomial) produce a family object. }
  \item{type}{ whether \code{left}, \code{right} or in \code{both} sides truncation is required, (left is the default)}
    \item{varying}{whether the truncation  varies for diferent observations. This can be usefull in regression analysis. If \code{varying = TRUE} then \code{par} should be an n-length vector for \code{type} equal \code{"left"} and \code{"right"} and an n by 2 matrix for \code{type="both"}}
  \item{\dots}{for extra arguments }
}

\value{
   Return a p family function 
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

(see also  \url{https://www.gamlss.com/}). }
\author{  Mikis Stasinopoulos \email{d.stasinopoulos@gre.ac.uk} and Bob Rigby}


\seealso{  \code{\link{trun.d}}, \code{\link{trun.q}}, \code{\link{trun.r}}, \code{\link{gen.trun}}}
\examples{
# trucated p  continuous function
# continuous
#----------------------------------------------------------------------------------------
# left
test1<-trun.p(par=c(0), family="TF", type="left")
test1(1)
(pTF(1)-pTF(0))/(1-pTF(0))
if(abs(test1(1)-(pTF(1)-pTF(0))/(1-pTF(0)))>0.00001) 
                      stop("error in left trucation of p")
plot(function(x) test1(x, mu=2, sigma=1, nu=2),0,10)
#----------------------------------------------------------------------------------------
# right
test2 <- trun.p(par=c(10), family="BCT", type="right")
test2(1)
pBCT(1)/pBCT(10)
if(abs(test2(1)-pBCT(1)/pBCT(10))>0.00001) stop("error in right trucation")
test2(1, lower.tail=FALSE)
1-pBCT(1)/pBCT(10)
if(abs(test2(1, lower.tail=FALSE)-(1-pBCT(1)/pBCT(10)))>0.00001) 
                stop("error in right trucation")
test2(1, log.p=TRUE)
log(pBCT(1)/pBCT(10))
if(abs(test2(1, log.p=TRUE)-log(pBCT(1)/pBCT(10)))>0.00001) 
               stop("error in right trucation")
plot(function(x) test2(x, mu=2, sigma=1, nu=2, tau=2),0,10)
plot(function(x) test2(x, mu=2, sigma=1, nu=2, tau=2, 
                       lower.tail=FALSE),0,10)
#----------------------------------------------------------------------------------------
# both 
test3<-trun.p(par=c(-3,3), family="TF", type="both")
test3(1)
(pTF(1)-pTF(-3))/(pTF(3)-pTF(-3))
if(abs(test3(1)-(pTF(1)-pTF(-3))/(pTF(3)-pTF(-3)))>0.00001) 
             stop("error in right trucation")
test3(1, lower.tail=FALSE)
1-(pTF(1)-pTF(-3))/(pTF(3)-pTF(-3))
if(abs(test3(0,lower.tail=FALSE)-
            (1-(pTF(0)-pTF(-3))/(pTF(3)-pTF(-3))))>0.00001) 
      stop("error in right trucation")
plot(function(x) test3(x, mu=2, sigma=1, nu=2, ),-3,3)
plot(function(x) test3(x, mu=2, sigma=1, nu=2, lower.tail=FALSE),-3,3)
#----------------------------------------------------------------------------------------
# Discrete
#----------------------------------------------------------------------------------------
# trucated p function
# left
test4<-trun.p(par=c(0), family="PO", type="left")
test4(1)
(pPO(1)-pPO(0))/(1-pPO(0))
if(abs(test4(1)-(pPO(1)-pPO(0))/(1-pPO(0)))>0.00001) 
               stop("error in left trucation of p")
plot(function(x) test4(x, mu=2), from=1, to=10, n=10, type="h")
cdf <- stepfun(1:40, test4(1:41, mu=5), f = 0)
plot(cdf, main="cdf", ylab="cdf(x)", do.points=FALSE )
#----------------------------------------------------------------------------------------
# right
test5<-trun.p(par=c(10), family="NBI", type="right")
test5(2)
pNBI(2)/(pNBI(9))
if(abs(test5(2)-(pNBI(2)/(pNBI(9))))>0.00001) 
          stop("error in right trucation of p")
plot(function(x) test5(x, mu=2), from=0, to=9, n=10, type="h")
cdf <- stepfun(0:8, test5(0:9, mu=5), f = 0)
plot(cdf, main="cdf", ylab="cdf(x)", do.points=FALSE )
#----------------------------------------------------------------------------------------
# both 
test6<-trun.p(par=c(0,10), family="NBI", type="both")
test6(2)
(pNBI(2)-pNBI(0))/(pNBI(9)-pNBI(0))
if(abs(test6(2)-(pNBI(2)-pNBI(0))/(pNBI(9)-pNBI(0)))>0.00001) 
               stop("error in the both trucation")
test6(1, log=TRUE)
log((pNBI(1)-pNBI(0))/(pNBI(9)-pNBI(0)))
if(abs(test6(1, log=TRUE)-log((pNBI(1)-pNBI(0))/(pNBI(9)-pNBI(0))))>0.00001) 
             stop("error in both trucation")
plot(function(y) test6(y, mu=20, sigma=3), from=1, to=9, n=9, type="h") 
plot(function(y) test6(y, mu=300, sigma=.4), from=1, to=9, n=9, type="h")
cdf <- stepfun(1:8, test6(1:9, mu=5), f = 0)
plot(cdf, main="cdf", ylab="cdf(x)", do.points=FALSE )
#----------------------------------------------------------------------------------------
# varying  truncation
#----------------------------------------------------------------------------------------
# coninuous
# left
test6<-trun.p(par=c(0,1,2), family="TF", type="left", varying=TRUE)
test6(c(2,3,4))
(pTF(c(2,3,4))-pTF(c(0,1,2)))/(1-pTF(c(0,1,2)))
test6(c(2,3,4), log.p=TRUE)
#----------------------------------------------------------------------------------------
# right
test7 <- trun.p(par=c(10,11,12), family="BCT", type="right", varying=TRUE)
test7(c(1,2,3))
pBCT(c(1,2,3))/pBCT(c(10,11,12))
test7(c(1,2,3), lower.tail=FALSE)
1-pBCT(c(1,2,3))/pBCT(c(10,11,12))
test7(c(1,2,3), log.p=TRUE)
#--------------------------------------------------------------------------------------- 
# both 
test8<-trun.p(par=cbind(c(0,1,2), c(10,11,12)), family="TF", 
           type="both", varying=TRUE)
test8(c(1,2,3))
(pTF(c(1,2,3))-pTF(c(0,1,2)))/(pTF(c(10,11,12))-pTF(c(0,1,2)))
test8(c(1,2,3), lower.tail=FALSE)
1-(pTF(c(1,2,3))-pTF(c(0,1,2)))/(pTF(c(10,11,12))-pTF(c(0,1,2)))
#--------------------------------------------------------------------------------------
# discrete
#--------------------------------------------------------------------------------------
# left
test9<-trun.p(par=c(0,1,2), family="PO", type="left", varying=TRUE)
test9(c(1,2,3))
(pPO(c(1,2,3))-pPO(c(0,1,2)))/(1-pPO(c(0,1,2)))
#--------------------------------------------------------------------------------------
# right
test10<-trun.p(par=c(10,11,12), family="NBI", type="right", varying=TRUE)
test10(c(2,3,4))
pNBI(c(2,3,4))/(pNBI(c(9,10,11)))
#-------------------------------------------------------------------------------------
# both
test11<-trun.p(par=rbind(c(0,10), c(1,11), c(2, 12)), family="NBI", 
              type="both", varying=TRUE)
test11(c(2,3,4))
(pNBI(c(2,3,4))-pNBI(c(0,1,2)))/(pNBI(c(9,10,11))-pNBI(c(0,1,2)))
#-------------------------------------------------------------------------------------
}
\keyword{distribution}
\keyword{regression}
