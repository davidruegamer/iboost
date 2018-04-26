getZtilde <- function(Ups, mod){
  
  nu <- mod$control$nu
  
  selpath <- selected(mod)
  
  sel <- sort(unique(selpath))
  
  vT <- mod$hatvalues()

  Zlist <- vector("list", length(selpath))
  
  Zstart <- lapply(vT, function(x) x*0)[sel]
  
  Zlist[1] <- (list(Zstart))
  
  names(Zlist[[1]]) <- sel
  
  for(i in 1:length(selpath))
  {
    
    thisZ <- nu * vT[[selpath[i]]] %*% Ups[[i]] #thisUps[[order(stssp)[stssp==j]]]
    if(i!=1) Zlist[[i]] <- Zlist[[i-1]]
    Zlist[[i]][[which(selpath[i]==sel)]] <- thisZ + Zlist[[i]][[which(selpath[i]==sel)]]
    names(Zlist[[i]]) <- sel
    
    
  }
    
  return(Zlist)
  
}



getZ <- function(Ups, mod){
  
  nu <- mod$control$nu
  
  
  selpath <- selected(mod)
  
  sel <- sort(unique(selpath))
  
  X <- getDesignmat(mod)
  # Xplus <- solve(crossprod(X)) %*% t(X)
  vT <- lapply(1:ncol(X), function(i) 
    t(X[,i])*as.numeric(solve(crossprod(X[,i])))) 
  
  Zlist <- vector("list", length(selpath))
  
  Zlist[[1]] <- (list(matrix(0, ncol = ncol(Ups[[1]]),
                             nrow = 1)))[rep(1,length(sel))]
  
  names(Zlist[[1]]) <- sel
  
  for(i in 1:length(selpath))
  {
    
    thisZ <- nu * vT[[selpath[i]]] %*% Ups[[i]] #thisUps[[order(stssp)[stssp==j]]]
    if(i!=1) Zlist[[i]] <- Zlist[[i-1]]
    Zlist[[i]][[selpath[i]]] <- thisZ + Zlist[[i]][[selpath[i]]]
    names(Zlist[[i]]) <- sel
    
    
  }
  
  return(Zlist)
  
}



# getUpsLast <- function(Ups, selpath)
# {
#   
#   cs <- sort(unique(selpath))
#   lapply(cs, function(i) ups[[max(which(selpath==i))]])
#   
# }
