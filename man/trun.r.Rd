\name{trun.r}
\alias{trun.r}

\title{Generates Random Values from a Truncated Density Function of a gamlss.family Distribution}
\description{
 Creates a function to generate randon values from a truncated probability density function created from a current GAMLSS family distribution

For continuous distributions left truncation at  3 means that the random variable can take the value 3.
For discrete distributions left truncation at 3 means that the random variable can take values from 4 onwards.
This is the same for right truncation. Truncation at 15 for a discrete variable means that 
15 and greater values are not allowed but for continuous variable it mean values greater that 15 are not allowed 
(so 15 is a possible value). 
}
\usage{
trun.r(par, family = "NO", type = c("left", "right", "both"), 
          varying = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
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
Returns a r family function  
}
\references{ Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}). }
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@gre.ac.uk} and Bob Rigby}

\seealso{ \code{\link{trun.p}}, \code{\link{trun.q}}, \code{\link{trun.d}}, \code{\link{gen.trun}} }
\examples{
# trucated r function
# continuous
#----------------------------------------------------------------------------------------
# left
test1<-trun.r(par=c(0), family="TF", type="left")
rr<-test1(1000)
hist(rr)
#----------------------------------------------------------------------------------------
# right
test2 <- trun.r(par=c(10), family="BCT", type="right")
rr<-test2(1000)
hist(rr)
#----------------------------------------------------------------------------------------
# both
test3<-trun.r(par=c(-3,3), family="TF", type="both")
rr<-test3(1000)
hist(rr)
#----------------------------------------------------------------------------------------
# discrete 
# trucated r function
# left
test4<-trun.r(par=c(0), family="PO", type="left")
tN <- table(Ni <- test4(1000))
r <- barplot(tN, col='lightblue')
#----------------------------------------------------------------------------------------
# right
test5 <- trun.r(par=c(10), family="NBI", type="right")
tN <- table(Ni <- test5(1000))
r <- barplot(tN, col='lightblue')
tN <- table(Ni <- test5(1000,mu=5))
r <- barplot(tN, col='lightblue')
tN <- table(Ni <- test5(1000,mu=10, sigma=.1))
r <- barplot(tN, col='lightblue')
#----------------------------------------------------------------------------------------
# both 
test6<-trun.r(par=c(0,10), family="NBI", type="both")
tN <- table(Ni <- test6(1000,mu=5))
r <- barplot(tN, col='lightblue')
#----------------------------------------------------------------------------------------
# varying = TRUE
#----------------------------------------------------------------------------------------
# continuous
#----------------------------------------------------------------------------------------
# left
test7<-trun.r(par=c(0,1,2), family="TF", type="left", varying=TRUE)
test7(3)

#----------------------------------------------------------------------------------------
# right
test8 <- trun.r(par=c(10,11,12), family="BCT", type="right", varying=TRUE)
test8(3)
#----------------------------------------------------------------------------------------
# both
test9<-trun.r(par=rbind(c(-3,3), c(-1,5), c(0,6)), , family="TF", type="both", varying=TRUE)
test9(3)
#----------------------------------------------------------------------------------------
# discrete 
# trucated r function
# left
test10<-trun.r(par=c(0,1,2), family="PO", type="left", varying=TRUE)
test10(3)
#----------------------------------------------------------------------------------------
# right
test11 <- trun.r(par=c(10,11,12), family="NBI", type="right", varying=TRUE)
test11(3)
test11(3, mu=10, sigma=.1)
#----------------------------------------------------------------------------------------
# both 
test12<-trun.r(par=rbind(c(0,10), c(1,11), c(2,12)), family="NBI", type="both", varying=TRUE)
test12(3,mu=5)
}
\keyword{distribution}
\keyword{regression}
