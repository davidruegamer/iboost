polyh_inf <- function(obj, vT, alpha, Ups = NULL, ..., returnComps = FALSE)
{
  
  ### get stuff from the model
  Y <- obj$response
  X <- getDesignmat(obj, full = T)  
  p <- ncol(X)
  selCourse <- selected(obj)
  
  ### get Upsilons if not already available
  if(is.null(Ups)) Ups <- getUpsilons(obj)
  # check for correctness:
  # sapply(2:length(Ups),function(i) sum(Ups[[i]]%*%Y - obj[i-1]$resid()))
  # should be approx 0 for all iterations
  
  ### get coefficient path
  sig <- getCoefPath(obj, what = "sign")
  
  ### function for comparison
  olsFun <- function(x) t(x)/sqrt(as.numeric(crossprod(x)))
  
  ### construct polyhedron characterisic matrix Gamma
  Gamma <- unlist(lapply(1:length(selCourse), function(i){
    # different indexing: Upsilons = Upsilon[0:(k-1)]
    k = selCourse[i]
    lapply(c(1:p)[-k], function(j){
      
      x <- rbind((sig[i]*olsFun(X[,k]) + olsFun(X[,j]))%*%Ups[[i]],
            (sig[i]*olsFun(X[,k]) - olsFun(X[,j]))%*%Ups[[i]]
      )
    
      rownames(x) <- c(paste(i,k,j,"+",sep="_"), paste(i,k,j,"-",sep="_"))
      return(x)
      
      })
  }), recursive=F)
  
  Gamma <- do.call("rbind", Gamma)
  
  base <- as.data.frame(do.call("rbind", lapply(rownames(Gamma), 
                                                function(x) strsplit(x,"_")[[1]])),
                        stringsAsFactors = FALSE)
  colnames(base) <- c("iteration", "selected", "comparison", "sign")
  base[,1:3] <- sapply(base[,1:3], as.numeric)
  
  attr(Gamma, "base") <- base

  # as long as covariance is diagonal, 
  # we dont need the variance
  # for the boundaries anyway
    
  if(returnComps){
    
    return(list(vT = vT, Gamma = Gamma))
    
  }else{
    
    return(lapply(1:length(vT), function(j) polyh_vlovup(vT=vT[[j]], Y = Y,
                                                         Gamma = Gamma)))
    
  }
    
}

polyh_vlovup <- function(vT, Y, Gamma, returnBounds = TRUE)
{
  
  z = as.numeric(vT %*% Y)
  vv = sum(vT^2)
  # sd = sqrt(var)*sqrt(vv) 
  rho = Gamma %*% t(vT) / vv
  vec = (- Gamma %*% Y + rho * z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))
  if(returnBounds) return(c(vlo,vup)) else
    return(c(which.max(vec[rho>0]), which.min(vec[rho<0])))
  
}