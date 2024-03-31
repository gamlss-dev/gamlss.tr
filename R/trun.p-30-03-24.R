################################################################################
################################################################################
################################################################################
################################################################################
trun.p <- function(par, family = "NO", 
                          type = c("left", "right", "both"), 
                       varying = FALSE, ...)
  {
    type <- match.arg(type)
   fname <- family
if (mode(family) != "character" && mode(family) != "name")
   fname <- as.character(substitute(family))
 distype <- eval(call(family))$type
    pfun <- paste("p",fname,sep="")
     cdf <- eval(parse(text=pfun))
if (!varying)
 {############### NOT varying ##################################################  
if (type=="both" && length(par)!= 2)  
  stop(paste("the length of par should be 2 \n"))
if (type!="both" && length(par)!= 1)  
  stop(paste("the length of par should be 1 \n"))  
fun <- if (type=="left")  
      function(q, lower.tail = TRUE, log.p = FALSE, ...)
     {
      cof <- (cdf(q,...)-cdf(par,...))/(1-cdf(par,...))
      cof <- if(lower.tail == TRUE) cof  else 1-cof   
      cof <- if(log.p==FALSE)  cof else  log(cof) 
      cof <- if (distype=="Discrete") ifelse(q <= par, 0, cof)
      else  ifelse(q < par, 0, cof)
      cof
      }
     else if (type=="right")
       function(q, lower.tail = TRUE, log.p = FALSE, ...)
      {
      cof <- if (distype=="Discrete") cdf(q,...)/cdf(par-1,...)
             else cdf(q,...)/cdf(par,...) # added Friday, February 26, 2010 
      cof <- if(lower.tail == TRUE) cof  else 1-cof   
      cof <- if(log.p==FALSE)  cof else  log(cof) 
      cof <- if (distype=="Discrete") ifelse(q >= par, 0, cof)
               else ifelse(q > par, 0, cof)
      cof
      }
     else if (type=="both")    
      function(q, lower.tail = TRUE, log.p = FALSE, ...)
      {
     cof <- if (distype=="Discrete") 
       (cdf(q,...)-cdf(par[1],...))/(cdf(par[2]-1,...)-cdf(par[1],...)) 
          else (cdf(q,...)-cdf(par[1],...))/(cdf(par[2],...)-cdf(par[1],...))   
     # added Friday, February 26, 2010     
      cof <- if(lower.tail == TRUE) cof  else 1-cof   
      cof <- if(log.p==FALSE)  cof else  log(cof) 
      cof <- if (distype=="Discrete")                      
                ifelse( (q <= par[1] | q >= par[2]), 0, cof)
              else 
                ifelse( (q < par[1] | q > par[2]), 0, cof)
      cof
      }
}################ end for not varying ##########################################
else ############ this is for varying truncation only ##########################
{
  	if (type=="both" && dim(par)[2]!= 2)  
  	   stop(paste("the rows of par should be 2 \n"))
  	#--
  	fun <- if (type=="left")  
       function(q, lower.tail = TRUE, log.p = FALSE, ...)
      {
      if (length(q)!=length(par)) 
          stop(paste("length of q must be equal to ", length(par), "\n","" ))	
      if (distype=="Discrete" &&  any(q <= par))  
          stop(paste("q must be greater than ", par, "\n", ""))
      if (distype!="Discrete" && any(q < par))
          stop(paste("q must be greater or equal than ", par, "\n", ""))
      cof <- (cdf(q,...)-cdf(par,...))/(1-cdf(par,...))
      cof <- if(lower.tail == TRUE) cof  else 1-cof   
      cof <- if(log.p==FALSE)  cof else  log(cof) 
      cof
      }
     else if (type=="right")
       function(q, lower.tail = TRUE, log.p = FALSE, ...)
      {
      if (length(q)!=length(par)) 
          stop(paste("length of q must be equal to ", length(par), "\n" ))	
      if (distype=="Discrete" &&  any(q >= par))  
          stop(paste("q must be less than ", par, "\n", ""))
        if (distype!="Discrete" && any(q > par))
          stop(paste("q must be less or equal than ", par, "\n", ""))
      cof <- if (distype=="Discrete") cdf(q,...)/cdf(par-1,...)
             else cdf(q,...)/cdf(par,...) # added Friday, February 26, 2010 
      cof <- if(lower.tail == TRUE) cof  else 1-cof   
      cof <- if(log.p==FALSE)  cof else  log(cof) 
      cof
      }
     else if (type=="both")    
      function(q, lower.tail = TRUE, log.p = FALSE, ...)
      {
     if (dim(par)[1]!= length(q))  
          stop(paste("the length of q should be the equal to ", dim(par)[1],  
                     "\n")) 	
     if (distype=="Discrete" &&  (any(q <= par[,1]) || any(q >= par[,2])) )  
       stop(paste("q must be between the lower and uppper of the par argument", "\n", ""))
        if (distype!="Discrete" && (any(q < par[,1]) || any(q > par[,2])) )
          stop(paste("q must be between the lower and uppper of the par argument", "\n", ""))  
      cof <- if (distype=="Discrete") (cdf(q,...)-
                  cdf(par[,1],...))/(cdf(par[,2]-1,...)-cdf(par[,1],...)) 
             else (cdf(q,...)-cdf(par[,1],...))/(cdf(par[,2],...)-
                  cdf(par[,1],...))   # added Friday, February 26, 2010     
      cof <- if(lower.tail == TRUE) cof  else 1-cof   
      cof <- if(log.p==FALSE)  cof else  log(cof) 
      cof
      }########### end varying #################################################
}	  
  fun
 }
################################################################################
################################################################################
################################################################################
################################################################################