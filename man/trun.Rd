\name{trun}
\alias{trun}

\title{Fits a Truncate Distribution from a gamlss.family}
\description{
This function can be used to fit truncated distributions. It takes as an argument an existing GAMLSS family distribution and 
a parameter vector,  of the type c(left.value, right.value), and  generates a \code{gamlss.family} object which then can be used to fit 
 a truncated distribution.    
}
\usage{
trun(par = c(0), family = "NO",  type = c("left", "right", "both"), name = "tr", 
        local = TRUE, delta=NULL, varying = FALSE, ...)
}

\arguments{
  \item{par}{ a vector with one (for \code{"left"} or \code{"right"} truncation) or two elements for \code{"both"}. When the argument \code{varying = TRUE} then \code{par} can be a vector or a matrix with two columns respectively.}
  \item{family}{an existing \code{gamlss.family} distribution}
  \item{type}{what type of truncation is required, \code{left}, \code{right} or \code{both}. If \code{both} the \code{par} should be a vector of length two. (the default is left truncation)}
  \item{name}{a character string to be added to name of the created object i.e. with \code{family=TF} and \code{name=trZero} the \code{gamlss.family} object will be called \code{TFtrZero}}
  \item{local}{if TRUE the function will try to find the environment of \code{gamlss} to generate the \code{d} and \code{p} functions required for the fitting, if FALSE the functions will be generated in the global environment}
  \item{delta}{the delta increment used in the numerical derivatives} 
   \item{varying}{whether the truncation  varies for diferent observations. This can be usefull in regression analysis. If \code{varying = TRUE} then \code{par} should be an n-length vector for \code{type} equal \code{"left"} and \code{"right"} and an n by 2 matrix for \code{type="both"}}
  \item{\dots}{for extra arguments }
}
\details{
This function is created to help the user to fit a truncated form of existing \code{gamlss} distribution.
It does this by taking an existing \code{gamlss.family} and changing some of the components of the distribution to help the fitting process. 
It particular it i) creates a pdf (\code{d}) and a cdf (\code{p}) function within \code{gamlss}, 
ii) changes  the global deviance function \code{G.dev.incr}, the first derivative functions (see note below) and the quantile residual function.   
}
\value{
It returns a \code{gamlss.family} object which has all the components needed for fitting a distribution in \code{gamlss}. 
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

(see also  \url{https://www.gamlss.com/}).}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@gre.ac.uk} and Bob Rigby}

\note{This function is experimental and could be changed. The function \code{trun}  changes 
      the first derivatives of the original gamlss family \code{d} function to numerical derivatives 
      for the new truncated \code{d} function. The default increment \code{delta},
      for this numerical derivatives function, is \code{eps * pmax(abs(x), 1)} where 
      \code{eps<-sqrt(.Machine$double.eps)}.  The default \code{delta} could be inappropriate 
      for specific applications and can be overwritten by using the argument \code{delta}.}
\seealso{\code{\link{trun.d}}, \code{\link{trun.p}}, \code{\link{trun.q}}, \code{\link{trun.r}},  \code{\link{gen.trun}} }
\examples{
# generate a left truncated zero t family
gen.trun(0,family="TF")
# take a random sample of 1000 observations
sam<-rTFtr(1000,mu=10,sigma=5, nu=5 )
hist(sam)
# fit the distribution to the data
mod1<-gamlss(sam~1, family=trun(0,TF))
mod1
# now create a gamlss.family object before the fitting 
Ttruc.Zero<- trun(par=0,family=TF, local=FALSE)
mod2<-gamlss(sam~1, family=Ttruc.Zero)
# now check the sensitivity of delta 
Ttruc.Zero<- trun(par=0,family=TF, local=FALSE, delta=c(0.01,0.01, 0.01))
mod3<-gamlss(sam~1, family=Ttruc.Zero)
}
\keyword{distribution}
\keyword{regression}
