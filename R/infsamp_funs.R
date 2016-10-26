infsamp <- function(refitFun, y, vT, confun, B, var, ncore)
{
 
  Pv <- lapply(vT, function(x) crossprod(x) / as.numeric(tcrossprod(x)))
  yperp <- lapply(Pv, function(p) y - p %*% y) 
  r <- lapply(vT, function(x) as.numeric(x %*% y) / as.numeric(tcrossprod(x)))
  v <- lapply(vT, t)
  
  lims <- list()
  
  for(j in 1:length(vT)){
    
    lims[[j]] <- findtrunclims(refitFun, r = r[[j]], yperp = yperp[[j]], B = B, 
                               var = var, v = v[[j]], ncore = ncore, confun = confun)
    
  
  }
   
  return(lims)
  
}


findtrunclims <- function(refitFun, r, yperp, B, var, v, ncore, confun)
{

  # define checkfun
  checkfun <- function(val) confun(refitFun(as.numeric(yperp + val * v)))

  vlo <- binsearch(r - 3^c(3:-8) * sqrt(var), nrIter = B, 
                   checkfun = checkfun, ncore = ncore, x = r)
  vup <- binsearch(r + 3^c(-8:3) * sqrt(var), nrIter = B, 
                   checkfun = checkfun, ncore = ncore, x = r, lower = F)
  
  return(c(vlo,vup))
  
}

binsearch <- function(grid, nrIter, checkfun, ncore, x, lower = T, tol = 1e-6)
{
  
  # search for the (lower / !lower = upper) limit, that satisfies the checkfun
  
  count <- 1
  res <- x
  
  while(count < nrIter){
    
    logvals <- unlist(mclapply(grid, checkfun, mc.cores = ncore))
    
    if(all(logvals)){ 
      
      if(count == 1){ 
        
        res <- ifelse(lower, -Inf, Inf)
        
      }else{
        
        res <- ifelse(lower, min(grid), max(grid))  
        
        } 
      
      break
      
    }
    
    if(all(!logvals) & count == nrIter) break
    if(count == nrIter) res <- ifelse(lower, grid[which.min(logvals)], 
                                      grid[which.max(logvals)])
    
    i0 <- ifelse(lower, max(which(!logvals)), min(which(!logvals)))
    
    # last step was too small, so that search jumped over the boundary
    if( (i0==12 & lower) | (i0==1 & !lower) ) break
    
    left <- ifelse(lower, grid[i0], grid[i0 - 1] + tol/10)
    right <- ifelse(lower, grid[i0 + 1] - tol/10, grid[i0])
    
    grid <- seq(left, right, length = 12)
        
    if(diff(grid)[1] < tol){
      
      res <- ifelse(lower, grid[i0 + 1], grid[i0 - 1])
      break
      
    }
    
    count <- count + 1  
      
    }
    
  return(res)
  
}
