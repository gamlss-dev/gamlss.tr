\name{fitTail}
\alias{fitTail}
\alias{fitTailAll}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{For fitting truncated distribution to the tails of data 
}
\description{
There are two functions here. The function \code{fitTail()} which fits a truncated distribution to certain percentage of the tail of a response variable and the function  \code{fitTailAll()} which does a sequence of truncated fits. Plotting  the results from those fits is analogous to the Hill plot, Hill (1975).
}
\usage{
fitTail(y, family = "WEI3", percentage = 10, howmany = NULL, 
      type = c("right", "left"), ...)
   
fitTailAll(y, family = "WEI3", percentage = 10, howmany = NULL, 
      type = c("right", "left"), plot = TRUE, 
      print = TRUE, save = FALSE, start = 5, trace = 0,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The variable of interest
}
  \item{family}{a \code{gamlsss.family} distribution}
  \item{percentage}{what percentage of the tail need to be modelled, default is 10\%}
  \item{howmany}{how many observations in the tail needed. This is an alternative to \code{percentage}. If it specified it take over from the \code{percentage} argument otherwise \code{percentage} is used.}
  \item{type}{which tall needs checking the right (default) of the left}
  \item{plot}{whether to plot with default equal \code{TRUE} }
  \item{print}{whether to print the coefficients with default equal \code{TRUE}}
  \item{save}{whether to save the fitted linear model with default equal \code{FALSE}}
   \item{start}{where to start fitting from the tail of the data}
   \item{trace}{0: no output 1: minimal 2: print estimates}
  \item{\dots}{
for further argument to the fitting function
}
}
\details{
The idea here is to fit a truncated distribution to the tail of the data.
Truncated log-normal and Weibull distributions could be appropriate distributions. More details can be found in Chapter 6 of  "The Distribution Toolbox of GAMLSS" book  which can be found in   \url{https://www.gamlss.com/}).
}
\value{
A fitted gamlss model
}
\references{

Hill B. M. (1975) A Simple General Approach to Inference About the Tail of a Distribution
\emph{Ann. Statist.} Volume \bold{3}, Number 5, pp 1163-1174.

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
\author{
Bob Rigby \email{r.rigby@gre.ac.uk}, Mikis Stasinopoulos \email{d.stasinopoulos@gre.ac.uk} and Vlassios Voudouris 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[gamlss]{loglogSurv}}, \code{\link[gamlss]{logSurv}}
}
\examples{
data(film90)
F90 <- exp(film90$lborev1)# original scale
# trucated plots
# 10%
w403<- fitTail(F90, family=WEI3)
qqnorm(resid(w403))
abline(0,1, col="red")

\dontrun{
# hill -sequential plot 10% weibull
w1<-fitTailAll(F90)
# plot sigma
plot(w1[,2])
#-----------------
#LOGNO
l403<- fitTail(F90, family=LOGNO)
plot(l403)
qqnorm(resid(l403))
abline(0,1)
#  hill -sequential plot 10% log normal
l1<-fitTailAll(F90, family=LOGNO)
plot(l1[,2])
#-------------------------
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{regression}% __ONLY ONE__ keyword per line
