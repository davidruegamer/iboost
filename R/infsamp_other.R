findtrunclimsGLM <- function(refitFun, n, r, rsd, B, betas, corB, ncore, confun, family, desMat, resids)
{
  
  # define checkfun
  
  if(family$family=="gaussian"){
    
    checkfun <- function(val){ 
      
      eta <- as.numeric(desMat %*% 
                          (betas + (val-r)*corB) + resids)
      confun(refitFun(eta))
      
    }
    
  }else if(family$family=="poisson"){
    
    checkfun <- function(val){ 
      
      eta <- as.numeric(desMat %*% 
                          (betas + (val-r)*corB) + resids)
      confun(refitFun(rpois(n = n, lambda = family$linkinv(eta))))
    }
    
  }else if(family$family=="binomial"){
    
    checkfun <- function(val){ 
      
      eta <- as.numeric(desMat %*% 
                          (betas + (val-r)*corB) + resids)
      confun(refitFun(factor(rbinom(n = n, size = 1, prob = family$linkinv(eta)))))
    }
    
  }else stop("Not implemented yet.")
  
  vlo <- binsearch(r - 3^c(3:-8) * rsd, nrIter = B, 
                   checkfun = checkfun, ncore = ncore, x = r)
  vup <- binsearch(r + 3^c(-8:3) * rsd, nrIter = B, 
                   checkfun = checkfun, ncore = ncore, x = r, lower = F)
  
  return(c(vlo,vup))
  
}



findtrunclimsAM <- function(refitFun, r, rsd, v, B, ncore, confun, resids)
{
  
  # define checkfun
  
  vTv <- as.numeric(crossprod(v))
  
  checkfun <- function(val){ 
    
    eta <- as.numeric(val * v / vTv + resids)
    confun(refitFun(eta))
    
  }
  
  
  vlo <- binsearch(r - 3^c(3:-8) * rsd, nrIter = B, 
                   checkfun = checkfun, ncore = ncore, x = r)
  vup <- binsearch(r + 3^c(-8:3) * rsd, nrIter = B, 
                   checkfun = checkfun, ncore = ncore, x = r, lower = F)
  
  return(c(vlo,vup))
  
}
