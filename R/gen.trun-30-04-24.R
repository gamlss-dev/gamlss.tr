gen.trun <-function(par = c(0), 
                     family = "NO", 
                   extra.name = "tr", 
                     type = c("left", "right", "both"),
                  varying = FALSE,
                    print = TRUE,
                     ...)
{
   type <- match.arg(type)
   fam  <- as.gamlss.family(family) #ds Monday, March 10, 2008 at 10:07
  fname <- fam$family[[1]] 
   dfun <- paste(paste("d",fname,sep=""), extra.name, sep="")
   pfun <- paste(paste("p",fname,sep=""), extra.name, sep="")
   qfun <- paste(paste("q",fname,sep=""), extra.name, sep="")
   rfun <- paste(paste("r",fname,sep=""), extra.name, sep="")
    fun <- paste(fname, extra.name, sep="")
   alldislist <-c(dfun,pfun,qfun,rfun,fun)
# generate d 
eval(dummy <- trun.d(par, family = fname, type = type, varying = varying, ...))
eval(call("<-",as.name(dfun),dummy), envir=parent.frame(n = 1))
# generate p
eval(dummy <- trun.p(par, family = fname, type = type, varying = varying, ...))
eval(call("<-",as.name(pfun),dummy), envir=parent.frame(n = 1))
# generate q
eval(dummy <- trun.q(par, family = fname, type = type, varying = varying, ...))
eval(call("<-",as.name(qfun),dummy), envir=parent.frame(n = 1))
# generate r
eval(dummy <- trun.r(par, family = fname, type = type, varying = varying, ...))
eval(call("<-",as.name(rfun),dummy), envir=parent.frame(n = 1))
# generate the fitting distribution
eval(dummy <- trun(par, family = substitute(family), type = type, 
                   extra.name=extra.name, 
                   local=FALSE, varying = varying, ...))
eval(call("<-",as.name(fun),dummy), envir=parent.frame(n = 1))
if (print)
{
   cat("A truncated family of distributions from",  fname, "has been generated \n", 
       "and saved under the names: ", "\n",paste(alldislist,sep=","),"\n")#
   cat("The type of truncation is", type, "\n",
       "and the truncation parameter is", if(varying==FALSE) par else "varying", " \n") 
}  
 }
