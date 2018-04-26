impsamp <- function(refitFun, y, vT, is_congruent, B, var, ...)
{

  m <- nrow(vT)
  
  if(m==1) vT <- vT / as.numeric(sqrt(tcrossprod(vT)))
  R <- vT %*% y
  Z <- sqrt(sum(R^2))
  u <- t(vT) %*% R / Z
  yperp <- y - u * Z
  
  
  
  sigma1 <- sqrt(m - 2 * (gamma((m+1)/2) / gamma(m/2))^2 )
  
  # Do Monte Carlo (Importance Sampling)
  ss <- generateSamples(refitFun = refitFun, is_congruent = is_congruent, 
                        samplingFun = function(B){ 
                          
                          rBs <- Z + rnorm(B) * sigma1 * sqrt(var)
                          rBs[rBs>0]
                          
                        }, 
                        B = B, refPoint = yperp, dir = u)
  
  rBs <- ss$rBs
  survive <- rBs[ss$logvals]

  weights <- function(var = var){
    
    Z <- Z/sqrt(var)
    s <- survive/sqrt(var)
    (s^(m-1)) * (exp(-s^2/2)) / dnorm(s, mean = Z, sd = sigma1)
    
  }
  
  return(list(rBs = rBs, logvals = ss$logvals, obsval = Z, weights = weights))

}

normalsamp <- function(refitFun, y, vT, is_congruent, B, var, 
                       min_nr_res = 50,  
                       eps = 1e-7, maxIter = 1, ...)
{
  
  
  vTv = as.numeric(tcrossprod(vT))
  Z = as.numeric(vT%*%y)
  u = as.numeric(t(vT)) / vTv
  yperp = y - u*Z
  this_var <- var * vTv
  this_mean <- Z/2
  
  quant <- pnorm(this_mean, mean = 0, sd = sqrt(this_var))
  
  while(quant==1 | quant==0){
    
    this_mean <- this_mean*0.9
    quant <- pnorm(this_mean, mean = 0, sd = sqrt(this_var))
    
  }  
  
  weights <- rep(0, B)
  search <- TRUE
  counter <- 0
  
  while(search & counter <= maxIter){
    
    ss <- generateSamples(refitFun = refitFun, is_congruent = is_congruent, 
                          samplingFun = function(B) rnorm(n = B, mean = this_mean, sd = sqrt(this_var)), 
                          B = B, refPoint = yperp, dir = u)
    
    weights <- exp(-1/(2*this_var) * (2 * ss$survive * this_mean - this_mean^2))
    
    zeroW <- sum(weights < eps)
    lenW <- length(weights)
    
    if(lenW - zeroW > min_nr_res){
      
      search <- FALSE
      
    }else{
      
      # this_var = this_var*4
      this_mean <- this_mean * 0.9
      warning(paste0(zeroW, " out of ", 
                     lenW, 
                     " weights are zero. Trying sampling distribution ",
                     "with smaller mean value and repeating process."))
      
    }
    
    counter <- counter + 1
    
  }
  
  weights <- function(var = var) exp(-1/(2*var*vTv) * 
                                            (2 * ss$survive * this_mean - this_mean^2))
  
  return(list(rBs = ss$rBs, logvals = ss$logvals, obsval = Z, weights = weights))

}

unifsamp <- function(refitFun, y, vT, is_congruent, B, var, 
                     min_nr_res = 50,  
                     maxIter = 100, nInit = 20, ...)
{
  
  vTv = as.numeric(tcrossprod(vT))
  Z = as.numeric(vT%*%y)
  u = as.numeric(t(vT)) / vTv
  yperp = y - u*Z
  this_var <- var * vTv
  this_mean <- Z
  
  # get preliminary search region
  lo <- prelimSearch(mean = this_mean, 
                     var = this_var, 
                     upper = FALSE,
                     maxIter = maxIter,
                     refitFun = refitFun,
                     is_congruent = is_congruent,
                     samplingFun = function(B, sd) 
                       Z - abs(rnorm(B, mean = 0, sd = sd)),
                     B = nInit,
                     refPoint = yperp,
                     dir = u)
  
  up <- prelimSearch(mean = this_mean, 
                     var = this_var, 
                     maxIter = maxIter,
                     refitFun = refitFun,
                     is_congruent = is_congruent,
                     samplingFun = function(B, sd) 
                       Z + abs(rnorm(B, mean = 0, sd = sd)),
                     B = nInit,
                     refPoint = yperp,
                     dir = u)
  
  # now draw the B random numbers from a uniform distribution
  ss <- generateSamples(refitFun = refitFun, is_congruent = is_congruent, 
                        samplingFun = function(B) 
                          runif(n = B, min = lo, max = up), 
                        B = B, refPoint = yperp, dir = u)
  
  weights <- function(var = var) dnorm(ss$survive, mean = 0, 
                                            sd = sqrt(var * vTv))
  
  return(list(rBs = ss$rBs, logvals = ss$logvals, obsval = Z, 
              weights = weights))
  
}