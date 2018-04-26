# prelimSearch <- function(mean, 
#                          var, 
#                          eps = 1e-16, 
#                          maxIter,
#                          B,
#                          refitFun,
#                          is_congruent,
#                          samplingFun,
#                          refPoint,
#                          dir,
#                          U = NULL
#                          )
# {
#   
#   qs <- qnorm(c(0.01, 0.99))*sqrt(var)
#   
#   range <- c(min(qs[1], qs[1] + mean), 
#              max(qs[2], qs[2] + mean))
#   
#   # use a range, which results in weights bigger eps
#   admisvals <- sqrt(-2 * var * log(eps * sqrt(2 * pi * var)))
#   if(mean > 0) range[1] <- max(range[1],-admisvals) else
#     range[2] <- min(range[2],admisvals)
#   
#   # roughly check region 
#   use_sd <- sqrt((range[2] - mean(range))/qnorm(0.99))
#   
#   # initialize values
#   lowB <- Inf
#   upB <- -Inf
#   nc <- TRUE
#   cn <- 0
#   logvals <- FALSE
#   
#   while(nc)
#   {
#     
#     ss <- generateSamples(refitFun = refitFun, 
#                           is_congruent = is_congruent, 
#                           samplingFun = function(B) 
#                             samplingFun(B, sd = use_sd), 
#                           B = B, 
#                           refPoint = refPoint, 
#                           dir = dir, U = U)
#     
#     upB <- max(ss$survive)
#     
#     if(length(upB)!=0 && (upB > mean & upB < Inf)) nc <- FALSE
#     
#     use_sd <- use_sd * 0.9
#     cn <- cn + 1
#     
#     # print(cn)
#     if(cn == maxIter){ 
#       
#       warning("Couldn't find an upper bound for sampling")
#       upB <- range[2]
#       break
#       
#     }
#     
#   }
#   
#   use_sd <- sqrt((range[2] - mean(range))/qnorm(0.99))
#   nc <- TRUE
#   cn <- 0
#   logvals <- FALSE
#   
#   while(nc)
#   {
#     
#     ss <- generateSamples(refitFun = refitFun, 
#                           is_congruent = is_congruent, 
#                           samplingFun = function(B) 
#                             samplingFun(B, sd = use_sd), 
#                           B = B, 
#                           refPoint = refPoint, 
#                           dir = dir)
#     
#     lowB <- min(ss$survive)
#     
#     if(length(lowB)!=0 && (lowB < mean & lowB > -Inf)) nc <- FALSE
#     
#     use_sd <- use_sd * 0.9
#     cn <- cn + 1
#     
#     if(cn == maxIter){ 
#       
#       warning("Couldn't find a lower bound for sampling")
#       lowB <- range[1]
#       break
#       
#     }
#     
#   }
#   
#   return(c(lowB, upB))
#   
# }


prelimSearch <- function(mean, 
                         var, 
                         eps = 1e-16, 
                         upper = TRUE,
                         maxIter,
                         B,
                         refitFun,
                         is_congruent,
                         samplingFun,
                         refPoint,
                         dir,
                         U = NULL
)
{
  qs <- qnorm(c(0.01, 0.99))*sqrt(var)
  
  range <- c(min(qs[1], qs[1] + mean), 
             max(qs[2], qs[2] + mean))
  
  # use a range, which results in weights bigger eps
  admisvals <- sqrt(-2 * var * log(eps * sqrt(2 * pi * var)))
  if(mean > 0) range[1] <- max(range[1], -admisvals) else
    range[2] <- min(range[2], admisvals)
  
  # roughly check region 
  use_sd <- sqrt((range[2] - mean(range))/qnorm(0.99))
  
  # initialize values
  if(upper) val <- Inf else val <- -Inf
  
  nc <- TRUE
  cn <- 0
  logvals <- FALSE
  
  while(nc)
  {
    
    ss <- generateSamples(refitFun = refitFun, 
                          is_congruent = is_congruent, 
                          samplingFun = function(B) 
                            samplingFun(B, sd = use_sd), 
                          B = B, 
                          refPoint = refPoint, 
                          dir = dir, U = U)
    
    val <- max(ss$survive)
    
    if(length(val)!=0 && (val > mean & val < Inf)) nc <- FALSE
    
    use_sd <- use_sd * 0.9
    cn <- cn + 1
    
    # print(cn)
    if(cn == maxIter){ 
      
      warning("Couldn't find an upper bound for sampling")
      val <- if(upper) range[2] else range[1]
      break
      
    }
    
  }
 
  return(val)
   
}
