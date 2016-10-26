polyh_inf <- function(obj, vT, var, ncore, alpha, Ups = NULL, ...)
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
      
      rbind((sig[i]*olsFun(X[,k]) + olsFun(X[,j]))%*%Ups[[i]],
            (sig[i]*olsFun(X[,k]) - olsFun(X[,j]))%*%Ups[[i]]
      )
    
      })
  }), recursive=F)
  
  Gamma <- do.call("rbind", Gamma)
  # check for correctness:
  # all(Gamma%*%Y > 0)
  
  # attr(Gamma, "selected") <- rep(selCourse, each=p-1)
  # attr(Gamma, "notSelected") <- unlist(lapply(selCourse, function(k) c(1:p)[-k]))

  ### do inference
  ret <- lapply(1:length(vT), function(j) selectiveInf(v = vT[[j]], Y = Y, 
                                                       Gamma = Gamma,
                                                       sd = sqrt(var),
                                                       alpha = alpha))
  
  return(ret)
  
}