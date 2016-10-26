boot_inf <- function(mod, refit.mboost = NULL, evT, nrBL, confun, B, bootType, var, ncore, alpha)
{
  
  if(is.null(refit.mboost)) stopifnot(bootType == "nonparam")
  
  if(bootType == "nonparam"){
    
    ind <- 1:length(mod$`(weights)`)
    weightMatrix <- sapply(1:B, function(b){
      
      tt <- table(sample(ind, replace = T))
      sapply(ind, function(x) ifelse(x %in% names(tt), tt[names(tt) == x], 0))
      
    })
      
    ret <- mclapply(1:ncol(weightMatrix), function(i){ 
      
      mnew <- update(mod, weights = weightMatrix[,i])
      if(confun(mnew)) sapply(1:nrBL, function(w) evT(mnew, which = w) %*% mnew$response) else #/ 
        #as.numeric(crossprod(evT(mnew, which = w)))) else
        NA
       
    }, mc.cores = ncore)
    
  }else{
    
    mu <- getDesignmat(mod) %*% evT(mod, which = NULL) %*% mod$response
    newdat <- sapply(1:B, function(b) rnorm(length(mu), mu, sqrt(var)))
    
    ret <- mclapply(1:ncol(newdat), function(i){ 
      
      mnew <- refit.mboost(newdat[,i])
      if(confun(mnew)) 
        sapply(1:nrBL, function(w) evT(mod, which = w) %*% mnew$response) else # / 
        #as.numeric(crossprod(evT(mnew, which = w)))) else
          NA
      
    }, mc.cores = ncore)
    
  }
  
  ret <- ret[!sapply(ret, function(x) any(is.na(x)))]
  retDF <- lapply(1:nrBL, function(i) sapply(ret, "[[", i))
  names(retDF) <- 1:nrBL
  
  return(retDF)
  
}