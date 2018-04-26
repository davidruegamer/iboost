boot_inf <- function(mod, refit.mboost = NULL, which, confun, B = 1000, 
                     bootType = c("param", "nonparam"), var, ncore, alpha)
{
  
  nrBL <- length(unique(selected(mod)))
  
  # get testvector function
  evT <- function(mod, which){
    
    # TODO - CHECK: is order correct? 
    X <- getDesignmat(mod)[, order(sel)]
    W <- mod$`(weights)`
    if(is.null(which))
      solve(crossprod(X * W, X)) %*% t(X) %*% diag(W)  else
        ( solve(crossprod(X * W, X)) %*% t(X) %*% diag(W) )[which, ]
    
  } 
    
  
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
    
  }else if(bootType == "param"){
    
    mu <- getDesignmat(mod) %*% evT(mod, which = NULL) %*% mod$response
    newdat <- sapply(1:B, function(b) rnorm(length(mu), mu, sqrt(var)))
    
    ret <- mclapply(1:ncol(newdat), function(i){ 
      
      mnew <- refit.mboost(newdat[,i])
      if(confun(mnew)) 
        sapply(1:nrBL, function(w) evT(mod, which = w) %*% mnew$response) else # / 
        #as.numeric(crossprod(evT(mnew, which = w)))) else
          NA
      
    }, mc.cores = ncore)
    
  }else if(bootType == "paramCond"){
    
    vTmat <- evT(mod, which = NULL)
    p <- nrow(vTmat)
    n <- length(mod$response)
    cmod <- coef(mod)
    coefnames <- names(cmod)
    retDF <- list(length=p)
    
    for(j in 1:p){
    
      beta <- unlist(coef(mod))
      beta[j] <- 0
      mu <- getDesignmat(mod) %*% beta
      newdat <- sapply(1:B, function(b) rnorm(n, mu, sqrt(var)))
      cp <- coefnames[j]
      
      ret <- unlist(mclapply(1:ncol(newdat), function(i){ 
        
        mnew <- refit.mboost(newdat[,i])
        
        w <- which(names(coef(mnew))==cp)
        
        if(cp %in% names(coef(mnew))) 
          return(abs(cmod[[cp]])<abs(coef(mnew)[[w]])) else
            return(FALSE)
#         X <- getDesignmat(mnew)
#         coefMat <- solve(crossprod(X)) %*% t(X)
#         coefOLS <- coefMat %*% mnew$response

#         if(cp %in% names(coef(mnew))) return(coefOLS[w,]) else return(0)
        
      }, mc.cores = ncore))
     
      retDF[[j]] <- ret
      
    }
    
  }
  
  if(bootType != "paramCond"){
  
    ret <- ret[!sapply(ret, function(x) any(is.na(x)))]
    retDF <- lapply(1:nrBL, function(i) sapply(ret, "[[", i))
    names(retDF) <- 1:nrBL

  }
  
  return(retDF)
  
}