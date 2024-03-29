################################################################################
################################################################################
################################################################################
################################################################################
## last change to allow binomial distribution fitting 
## 28-6-18
## NOTE the for one parameter binomial type distribution the order 
#  of the agrument is
## x, bd , mu
## while for all the rest is 
## x, mu, sigma,...,bd
## this meens that the order in line 92 has to different from the rest
## If a new binomial type distribution is created it has to have the 
## binomial in order 
################################################################################
################################################################################
################################################################################
################################################################################
## This function creates a truncated gamlss.family object which can 
##  be used for
## fitting a GAMLSS model
## to correct the residuals for four parameter distributions
## change 30-10-2015 MS
# to make sure that the function creates a function in which the links 
# can be modified
################################################################################
################################################################################
################################################################################
################################################################################
## A problem arise when we need different link function 
## the solution at the moment is to as for the different link in the begining
## This hopefully is fixed now 22-2-
################################################################################
################################################################################
################################################################################
################################################################################
# source('/Users/stasinom/Dropbox/gamlss/library/gamlss.tr/R/trun.d-27-12-12.R')
# source('/Users/stasinom/Dropbox/gamlss/library/gamlss.tr/R/trun.p-27-12-12.R')
# source('/Users/stasinom/Dropbox/gamlss/library/gamlss.tr/R/trun.q-27-12-12.R')
# source('/Users/stasinom/Dropbox/gamlss/library/gamlss.tr/R/trun.r-28-12-12.R')
# source('/Users/stasinom/Dropbox/gamlss/library/gamlss.tr/R/gen.trun-21-06-13.R')
# require(gamlss)
################################################################################
################################################################################
################################################################################
################################################################################
trun <-function ( par = c(0), 
               family = "NO", 
                 type = c("left", "right", "both"),
                 name = "tr", 
                local = TRUE,
                delta = NULL, 
              varying = FALSE,
                ...)
{
#-------------------------------------------------------------------------------
     TEST <- "TEST" # dummy name
     type <- match.arg(type)
    fname <- if (is.name(family)) as.character(family)
             else if (is.character(family)) family
             else if (is.call(family)) as.character(family[[1]])
             else if (is.function(family)) deparse(substitute(family))
             else if (is(family, "gamlss.family"))  family$family[1]
             else stop("the family must be a character or a gamlss.family name")
  ifBinomial <- fname%in%gamlss::.gamlss.bi.list
        fam1 <- eval(parse(text=fname)) # the family to output
         fam <- as.gamlss.family(family) # this is created so I can get things
      family <- c("None", "None")  
      dorfun <- paste("d",fname,sep="")
      porfun <- paste("p",fname,sep="")
        dfun <- paste(paste("d",fname,sep=""), name, sep="")
        pfun <- paste(paste("p",fname,sep=""), name, sep="")
       nopar <- fam$nopar # or fam1$nopar
#------------------------------------------------------------------------------
if (local)
 {
#--trying to get gamlss sys.frame--  
     rexpr<-regexpr("gamlss",sys.calls())
for (i in 1:length(rexpr)){ 
    position <- i 
    if (rexpr[i]==1) break}
gamlss.environment <- sys.frame(position)      
 } else gamlss.environment <- sys.frame(0)

if (varying) assign("PAR_", par, envir=gamlss.environment)
#   generate d within gamlss
    eval(dummy <- trun.d(par, family = fname, type = type, varying = varying, ...))
    eval(call("<-",as.name(dfun),dummy), envir=gamlss.environment)# parent.frame(n = 1)
# generate p within gamlss
    eval(dummy <- trun.p(par, family = fname, type = type, varying = varying, ...))
    eval(call("<-",as.name(pfun),dummy), envir=gamlss.environment)# parent.frame(n = 1)
#-------------------------------------------------------------------------------
# rename the family 
   family[[1]] <- paste(paste(fname, name, sep=""))
   family[[2]] <- paste(type, "truncated",fam$family[[2]])
     fam$family <- family # in fam 
body(fam1)[[nopar+2]][[2]]$family <- family # and in fam1
# Global deviance increment  
                 sGD <- gsub(dorfun, dfun, deparse(body(fam$G.dev.incr)))
body(fam$G.dev.incr) <- parse(text=sGD)
body(fam1)[[nopar+2]][[2]]$G.dev.incr <- fam$G.dev.incr
# get the no of parameters  
         nopar <- fam$nopar
# check for the delta
 if (length(delta)==0) delta <- rep(NA,nopar) 
 if (length(delta)==1) delta <- rep(delta,nopar)
 if (length(delta)!=nopar)  stop("delta should be the same length the parameters in the family ")
#-------------------------------------------------------------------------------
# now change the first derivatives
switch(nopar,  
  { 
# 1 parameter-------------------------------------------------------------------
# dldm
  fam$dldm <- if (ifBinomial) 
  {
    function(y, bd, mu) as.vector(attr(gamlss::numeric.deriv(TEST(y, bd, mu, log=TRUE), "mu", delta=NULL), "gradient")) 
  } else  
  {
    function(y,mu)    as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, log=TRUE), "mu", delta=NULL), "gradient"))  
  }  
      sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1])) 
      sMU <- sub("NULL",  as.character(delta[1]), sMU) 
body(fam$dldm) <- parse(text=sMU[length(sMU)])
body(fam1)[[nopar+2]][[2]]$dldm  <- fam$dldm
# residuals
         sres <- gsub(porfun, pfun,  deparse(fam$rqres[[1]]))
if  (fam$type == "Discrete")
    {
      if (varying==FALSE)
        {
         sres <-  if (type=="left"|type=="both")  
           gsub("ymin = 0",  paste("ymin =",  par[1]+1),  sres) 
                  else sres  
        }
          else
        {
  if (type=="left")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
  if (type=="both")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
        }  
    }
         #  sres <- gsub("expression", "",  sres)
      fam$rqres <- parse(text=sres)
body(fam1)[[nopar+2]][[2]]$rqres <- fam$rqres
  },
  {
# 2 parameters -------------------------------------------------------  
      # dldm and dldd
      fam$dldm <- if(ifBinomial)
      {
    function(y,mu,sigma,bd) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, bd, log=TRUE), "mu", delta=NULL), "gradient"))  
      } else
      {
    function(y,mu,sigma)    as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, log=TRUE),     "mu", delta=NULL), "gradient"))
      }  
      fam$dldd <- if(ifBinomial) 
      {
        function(y,mu,sigma, bd) as.vector(attr(gamlss::numeric.deriv(TEST(y,mu, sigma, bd, log=TRUE), "sigma", delta=NULL), "gradient"))
      } else 
      {
        function(y,mu,sigma) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, log=TRUE), "sigma", delta=NULL), "gradient"))
      }  # mu
       sMU <- sub("TEST", dfun, body(fam$dldm))
      if (!is.na(delta[1])) 
        sMU <- sub("NULL",  as.character(delta[1]), sMU)  
       body(fam$dldm) <- parse(text=sMU[length(sMU)])
       body(fam1)[[nopar+2]][[2]]$dldm  <- fam$dldm
       # sigma   
        sSIGMA <- sub("TEST", dfun, body(fam$dldd))
      if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA)
body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)]) 
body(fam1)[[nopar+2]][[2]]$dldd  <- fam$dldd
# residuals
      sres <- gsub(porfun, pfun,  deparse(fam$rqres[[1]]))
    if  (fam$type == "Discrete")
      {
      if (varying==FALSE)
        {
                sres <-  if (type=="left"|type=="both")  gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
                else sres  
        }
          else
        {
    if (type=="left")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
    if (type=="both")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
        }  
      }
            #sres <- gsub("expression", "",  sres)
       fam$rqres <- parse(text=sres) 
body(fam1)[[nopar+2]][[2]]$rqres <- fam$rqres
          }, # 3 parameters----------------------------------------------------
          # dldm dldd dldv 
          {
   fam$dldm <- if(ifBinomial) function(y,mu,sigma,nu,bd) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, bd, log=TRUE), "mu", delta=NULL), "gradient")) 
                         else function(y,mu,sigma,nu)    as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE),     "mu", delta=NULL), "gradient"))
            

   fam$dldd <- if(ifBinomial) function(y,mu,sigma,nu,bd) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, bd, log=TRUE), "sigma", delta=NULL), "gradient")) 
                        else  function(y,mu,sigma,nu)    as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu,     log=TRUE), "sigma", delta=NULL), "gradient"))
     
     
   fam$dldv <- if(ifBinomial) function(y,mu,sigma,nu,bd) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, bd, log=TRUE), "nu", delta=NULL), "gradient")) 
                         else function(y,mu,sigma,nu)    as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu,     log=TRUE), "nu", delta=NULL), "gradient"))
   
     
      sMU <- sub("TEST", dfun, body(fam$dldm))
  if (!is.na(delta[1]))sMU <- sub("NULL",  as.character(delta[1]), sMU)          
  body(fam$dldm) <- parse(text=sMU[length(sMU)])  
  body(fam1)[[nopar+2]][[2]]$dldm  <- fam$dldm
    sSIGMA <- sub("TEST", dfun, body(fam$dldd))
  if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA) 
  body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)])  
  body(fam1)[[nopar+2]][[2]]$dldd  <- fam$dldd
      sNU <- sub("TEST", dfun, body(fam$dldv))
  if (!is.na(delta[3])) sNU <- sub("NULL",  as.character(delta[3]), sNU)
  body(fam$dldv) <- parse(text=sNU[length(sNU)])
  body(fam1)[[nopar+2]][[2]]$dldv  <- fam$dldv
  sres <- gsub(porfun, pfun,  deparse(fam$rqres[[1]]))
  if  (fam$type == "Discrete")
  {
    if (varying==FALSE)
    {
      sres <-  if (type=="left"|type=="both")  gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
      else sres  
    }
    else
    {
      if (type=="left")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
      if (type=="both")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
    }  
  }
  #sres <- gsub("expression", "",  sres)
  fam$rqres <- parse(text=sres)   
  body(fam1)[[nopar+2]][[2]]$rqres <- fam$rqres  
          }, # 4 paramers------------------------------------------------------
          # dldm dldd dldv dldt
          {
  fam$dldm <- if(ifBinomial) function(y,mu,sigma,nu,tau,bd) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau, bd, log=TRUE), "mu", delta=NULL), "gradient")) 
                        else function(y,mu,sigma,nu,tau)    as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau,     log=TRUE), "mu", delta=NULL), "gradient"))
            
            
  fam$dldd <- if(ifBinomial) function(y,mu,sigma,nu,tau,bd)  as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau, bd, log=TRUE), "sigma", delta=NULL), "gradient")) 
                        else function(y,mu,sigma,nu,tau)     as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau,     log=TRUE), "sigma", delta=NULL), "gradient"))
            
            
  fam$dldv <- if(ifBinomial) function(y,mu,sigma,nu,tau,bd) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau, bd, log=TRUE), "nu", delta=NULL), "gradient")) 
                       else  function(y,mu,sigma,nu,tau)    as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE),     "nu", delta=NULL), "gradient"))
            
  fam$dldt <- if(ifBinomial) function(y,mu,sigma,nu,tau,bd) as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau, bd, log=TRUE), "tau", delta=NULL), "gradient")) 
                       else  function(y,mu,sigma,nu,tau)    as.vector(attr(gamlss::numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE),     "tau", delta=NULL), "gradient"))

       sMU <- sub("TEST", dfun, body(fam$dldm))
  if (!is.na(delta[1])) sMU <- sub("NULL",  as.character(delta[1]), sMU)      
  body(fam$dldm) <- parse(text=sMU[length(sMU)])   
  body(fam1)[[nopar+2]][[2]]$dldm  <- fam$dldm
     sSIGMA <- sub("TEST", dfun, body(fam$dldd))
  if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA) 
  body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)])  
  body(fam1)[[nopar+2]][[2]]$dldd  <- fam$dldd
       sNU <- sub("TEST", dfun, body(fam$dldv))
  if (!is.na(delta[3])) sNU <- sub("NULL",  as.character(delta[3]), sNU)           
  body(fam$dldv) <- parse(text=sNU[length(sNU)]) 
  body(fam1)[[nopar+2]][[2]]$dldv  <- fam$dldv
  sTAU <- sub("TEST", dfun, body(fam$dldt))
  if (!is.na(delta[4])) sTAU <- sub("NULL",  as.character(delta[4]), sTAU)
  body(fam$dldt) <- parse(text=sTAU[length(sTAU)]) 
  body(fam1)[[nopar+2]][[2]]$dldt  <- fam$dldt
  sres <- gsub(porfun, pfun,  deparse(fam$rqres[[1]]))
  if  (fam$type == "Discrete")
  {
    if (varying==FALSE)
    {
      sres <-  if (type=="left"|type=="both")  gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
      else sres  
    }
    else
    {
      if (type=="left")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
      if (type=="both")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
    }  
  }
  #sres <- gsub("expression", "",  sres)
  fam$rqres <- parse(text=sres) 
  body(fam1)[[nopar+2]][[2]]$rqres <- fam$rqres  
})

# create a function and save it (rather than a list)
#                 fun <- function() {}
#body(fun)[[nopar+2]] <- fam  
#  {Body[2:3] 
#   structure(fam)}
# body(fam1)[[nopar+2]][[2]]$dldm
fam1
}
 
