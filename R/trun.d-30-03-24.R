################################################################################
################################################################################
################################################################################
################################################################################
# this is to create a truncating distribution with varying truncation 
# according to par 
# TO DO : outside truncation set to zero
################################################################################
################################################################################
################################################################################
################################################################################ 
trun.d <-function(par, family = "NO", 
                         type = c("left", "right", "both"), 
                      varying = FALSE, ...)
  {
     type <- match.arg(type) # the truncation type 
    fname <- family          # the family      
if (mode(family) != "character" && mode(family) != "name")
    fname <- as.character(substitute(family))
  distype <- eval(call(family))$type
     dfun <- paste("d",fname,sep="")
     pfun <- paste("p",fname,sep="")
      pdf <- eval(parse(text=dfun))
      cdf <- eval(parse(text=pfun))  
if (!varying) ### not varying ##################################################
{
if (type=="both" && length(par)!= 2)  
  stop(paste("the length of par should be 2 \n"))
if (type!="both" && length(par)!= 1)  
  stop(paste("the length of par should be 1 \n"))
#--
fun <- if (type=="left")  ############# LEFT ###################################
       function(x, log = FALSE, ...)
    {
      dfun <- pdf(x,log = TRUE,...)-log(1-cdf(par,...))
       if (log == TRUE)
       {
        dfun <- dfun
        dfun <- if (distype=="Discrete") ifelse(x <= par, NA, dfun)
              else  ifelse(x < par, NA, dfun)
      } else 
      {
        dfun <- exp(dfun)
        dfun <- if (distype=="Discrete") ifelse(x <= par, 0, dfun)
                else  ifelse(x < par, 0, dfun)
      }
      dfun
    }
    else if (type=="right")  ############# Right ###############################
       function(x, log = FALSE, ...)
    {     
     dfun <-  if (distype=="Discrete") 
                 pdf(x, log = TRUE,...)-log(cdf(par-1,...))
            else pdf(x, log = TRUE,...)-log(cdf(par,...))     
         if (log == TRUE)
         {
           dfun <- dfun 
           dfun <- if (distype=="Discrete") ifelse(x >= par, NA, dfun)
                    else ifelse(x > par, NA, dfun) 
         } else 
         {
           dfun <-  exp(dfun)  
           dfun <- if (distype=="Discrete") ifelse(x >= par, 0, dfun)
                   else ifelse(x > par, 0, dfun)
         }   
      dfun
    } 
  else if (type=="both")  ############# BOTH ###################################  
    function(x, log = FALSE, ...)
    {
      dfun <- if (distype=="Discrete") pdf(x, log = TRUE,...) - 
                 log(cdf(par[2]-1,...)-cdf(par[1],...)) 
            else pdf(x, log = TRUE,...) - log(cdf(par[2],...)-
                  cdf(par[1],...))   
      if (log == TRUE)
      {
        dfun  <- dfun
        dfun <- if (distype=="Discrete")                      
          ifelse( (x <= par[1] | x >= par[2]), NA, dfun)
        else ifelse( (x < par[1] | x > par[2]), NA, dfun)
      } else 
      {
        dfun  <-  exp(dfun)
        dfun <- if (distype=="Discrete")                      
               ifelse( (x <= par[1] | x >= par[2]), 0, dfun)
              else ifelse( (x < par[1] | x > par[2]), 0, dfun)
      }
      dfun
    }   
}######### end of not varying ##################################################
  else # this is for varying truncation only 
{#########  varying ############################################################
    if (type=="both" && dim(par)[2]!= 2)  
   stop(paste("the rows of par should be 2 \n"))
fun <- if (type=="left") ############# LEFT ####################################
      function(x, log = FALSE, ...)
  {
  if (length(x)!=length(par)) 
warning(paste("The length of x must be equal to the length of varying parameters which is ", length(par), "\n", "please procceed with caution because the two vectors are forced to have the same length by expanding the smalest", "\n"))
      if (length(x)>length(par)) par <- rep(par, length=length(x))
      if (length(x)<length(par))   x <- rep(x, length=length(par))
          dfun <- pdf(x,log = TRUE,...)-log(1-cdf(par,...))
  if (log == TRUE)
    {
          dfun <-  dfun
          dfun <- if (distype=="Discrete") ifelse(x <= par, NA, dfun)
                  else  ifelse(x < par, NA, dfun)
    }
  else 
    {
         dfun <-   exp(dfun)
         dfun <- if (distype=="Discrete") ifelse(x <= par, 0, dfun)
                 else  ifelse(x < par, 0, dfun)
    }
        dfun
  }
      else if (type=="right") ############# RIGHT  #############################
      function(x, log = FALSE, ...)
      {
        if (length(x)!=length(par)) 
  warning(paste("The length of x must be equal to the length of varying parameters which is ", length(par), "\n", "please procceed with caution because the two vectors are forced to have the same length by expanding the smalest", "\n"))
        if (length(x)>length(par)) par <- rep(par, length=length(x))
        if (length(x)<length(par))   x <- rep(x, length=length(par))
        dfun <- if (distype=="Discrete") 
                     pdf(x, log = TRUE,...)-log(cdf(par-1,...))
                else pdf(x, log = TRUE,...)-log(cdf(par,...))

        if (log == TRUE)
          {
          dfun <-  dfun 
          dfun <- if (distype=="Discrete") ifelse(x >= par, NA, dfun)
                 else ifelse(x > par, NA, dfun) 
        } else 
        {
          dfun <- exp(dfun)  
          dfun <- if (distype=="Discrete") ifelse(x >= par, 0, dfun)
                  else ifelse(x > par, 0, dfun)
        }  
        dfun
      } 
    else if (type=="both")############# BOTH ###################################
      function(x, log = FALSE, ...)
      {
        if (length(x)!=length(par[,1])) 
          warning(paste("The length of x must be equal to the length of varying parameters which is ", length(par[,1]), "\n", "please procceed with caution because the two objects are forced to have the same length by expanding the smalest", 
                        "\n"))
        if (length(x)>dim(par)[1]) par[,1] <- rep(par[,1], length=length(x))
        if (length(x)>dim(par)[2]) par[,2] <- rep(par[,2], length=length(x))
        if (length(x)<dim(par)[1])       x <- rep(x, length=length(par[,1]))  
        dfun <-   if (distype=="Discrete") pdf(x, log = TRUE,...) - 
                       log(cdf(par[,2]-1,...)-cdf(par[,1],...)) 
                  else pdf(x, log = TRUE,...) - log(cdf(par[,2],...)-
                       cdf(par[,1],...)) 
        if (log == TRUE)
        {
          dfun <- dfun
          dfun <- if (distype=="Discrete")                      
                      ifelse( (x <= par[,1] | x >= par[,2]), NA, dfun)
                 else ifelse( (x < par[,1] | x > par[,2]), NA, dfun)
        }
        else
        {
         dfun <- exp(dfun) 
         dfun <- if (distype=="Discrete")                      
                      ifelse( (x <= par[,1] | x >= par[,2]), 0, dfun)
                 else ifelse( (x < par[,1] | x > par[,2]), 0, dfun)
        }  
          
        dfun
      }  
}####### end varying ###########################################################
fun    
}
