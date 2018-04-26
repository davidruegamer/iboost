infsamp <- function(refitFun, y, vT, confun, B, var, ...)
{
  
  if(!is.list(vt) && nrow(vt)==1){
    
    Pv <- crossprod(vt) / as.numeric(tcrossprod(vt))
    yperp <- y - Pv %*% y
    r <- as.numeric(vt %*% y)
    vnorm <- t(vt) / as.numeric(tcrossprod(vt))
    lin <- TRUE
    
  }else{
    
    R <- vt %*% y
    r <- sqrt(sum(R^2))
    vnorm <- t(vt) %*% R / r
    yperp <- y - vnorm * r
    
    lin <- FALSE
    
  }
  
  lims <- findtrunclims(refitFun, 
                        r = r, 
                        yperp = yperp, 
                        B = B, 
                        var = var, 
                        vnorm = vnorm, 
                        confun = confun,
                        lin = lin, ...)
  attr(lims, "r") <- r
  
  return(lims)
  
}




findtrunclims <- function(refitFun, r, yperp, B, var, 
                          vnorm, ncore, confun, lin = TRUE,
                          gridvals = 3:-8)
{

  # define checkfun
  checkfun <- function(val) confun(refitFun(as.numeric(yperp + val * vnorm)))

  if(lin) lowEnd <- r - 3^gridvals * sqrt(var) else
    lowEnd <- seq(0, r - 3^min(gridvals), l = 12)
  
  vlo <- binsearch(lowEnd, nrIter = B, 
                   checkfun = checkfun, ncore = ncore, x = r)
  vup <- binsearch(r + 3^(-gridvals) * sqrt(var), nrIter = B, 
                   checkfun = checkfun, ncore = ncore, x = r, lower = F)
  
  return(c(vlo,vup))
  
}




binsearch <- function(grid, nrIter, checkfun, ncore, x, lower = T, tol = 1e-6)
{
  
  # search for the (lower / !lower = upper) limit, that satisfies the checkfun
  
  count <- 1
  res <- x
  
  while(count < nrIter){
    
    logvals <- unlist(mclapply(round(grid, -log(tol,10)), checkfun, mc.cores = ncore))
    
    if(all(logvals)){ 
      
      if(count == 1){ 
        
        res <- ifelse(lower, -Inf, Inf)
        if(lower && grid[1]==0) res <- 0
        
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
    if( (i0==length(grid) & lower) | (i0==1 & !lower) ) break
    
    left <- ifelse(lower, grid[i0], grid[i0 - 1] + tol/10)
    right <- ifelse(lower, grid[i0 + 1] - tol/10, grid[i0])
    
    grid <- seq(left, right, length = length(grid))
        
    if(diff(grid)[1] < tol){
      
      # better safe than sorry -> select with some margin
      # -> with this tolerance it can happen that logval == TRUE
      # but the result later will be FALSE
      res <- ifelse(lower, grid[pmin(length(grid), i0 + 1 + 1)], 
                    grid[pmax(1, i0 - 1 - 1)])
      break
      
    }
    
    count <- count + 1  
      
    }
    
  return(res)
  
}
