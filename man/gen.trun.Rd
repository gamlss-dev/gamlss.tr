\name{gen.trun}
\alias{gen.trun}

\title{Generates a truncated distribution from a gamlss.family}
\description{
  The \code{gen.trun()} function allows the user to generate \code{d}, \code{p},  \code{q}, and \code{r} distribution functions plus  an extra
  \code{gamlss.family} function for fitting a truncated  distribution with \code{gamlss}. 

For continuous distributions left truncation at 3 means that the random variable can take the value 3.
For discrete distributions left truncation at 3 means that the random variable can take values from 4 onwards.
This is the same for right truncation. Truncation at 15 for a discrete variable means that 
15 and greater values are not allowed but for continuous variable it mean values greater that 15 are not allowed 
(so 15 is a possible value).    

If the user want a different link (rather the default) for any of the parameters she/he has to declare at the generation of the functions, see example.   
}
\usage{
gen.trun(par = c(0), family = "NO", name = "tr", 
         type = c("left", "right", "both"), 
         varying = FALSE,  print=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{ a vector with one (for \code{"left"} or \code{"right"} truncation) or two elements for \code{"both"}. When the argument \code{varying = TRUE} then \code{par} can be a vector or a matrix with two columns respectively.}
  \item{family}{ a \code{gamlss.family} object, which is used to define the distribution and the link functions of the various parameters. 
   The distribution families supported by \code{gamlss()} can be found in \code{gamlss.family}.
   }
  \item{name}{the extra characters to be added to the name of new truncated distribution, by default it adds \code{tr}}
  \item{type}{whether \code{"left"}, \code{"right"} or in \code{"both"} sides truncation is required  }
 \item{varying}{whether the truncation  varies for different observations. This can be useful in regression analysis. If \code{varying = TRUE} then \code{par} should be an n-length vector for \code{type} equal \code{"left"} and \code{"right"} and an n by 2 matrix for \code{type="both"}}
 \item{print}{whether to print the names of the created distribution}
  \item{\dots}{for extra arguments}
}

\value{
 Returns  the \code{d}, the \code{p}, the \code{q}, the \code{r} and the fitting functions of a truncated \code{gamlss.family} distribution. 
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

\seealso{ \code{\link{trun.d}}, \code{\link{trun.p}}, \code{\link{trun.q}}, \code{\link{trun.r}} }
\examples{
# generating a t-distribution from 0 to 100  	
gen.trun(par=c(0,100),family="TF", name="0to100", type="both")
op<-par(mfrow=c(2,2))
plot(function(x) dTF0to100(x, mu=80 ,sigma=20, nu=5), 0, 100, ylab="pdf")
plot(function(x) pTF0to100(x, mu=80 ,sigma=20, nu=5), 0, 100, ylab="cdf")
plot(function(x) qTF0to100(x, mu=80 ,sigma=20, nu=5), 0.01, .999, ylab="invcdf")
hist(s1<-rTF0to100(1000, mu=80 ,sigma=20, nu=5), ylab="hist", xlab="x", 
            main="generated data")
par(op)
m1<-histDist(s1, family=TF0to100, xlim=c(0,100))# fitting the data
# using the argumnt varying 
# left part varies right part equal 100
leftPAR <- rPO(100)
gen.trun(par=cbind(leftPAR,rep(100, 100)),family="TF", name="0to100Varying", 
            type="both", varying=TRUE)
YY<- rTF0to100Varying(100, mu=80, sigma=20, nu=5)
m1<-gamlss(YY~1, family=TF0to100Varying)
m1
}
\keyword{distribution}
\keyword{regression}
