# extract the design matrix (matrices if split = TRUE) from mod
getDesignmat <- function(mod, split = FALSE, full = FALSE)
{
  
  if(full){
    
    ret <- do.call("cbind", lapply(mod$baselearner, extract, "design"))
    
  }else{
  
    if(class(mod)[1] == "glmboost"){
      
      # get design matrix
      ret <- extract(mod, "design")
      if(!1 %in% selected(mod) & colnames(ret)[1] == "(Intercept)")
        ret <- ret[,-1]
      
      if(split){ # recover single design matrices
      
        if(is.null(dim(ret))){
        
          ret <- list(ret)
          
        }else{
          
          ret <- lapply(1:ncol(ret), function(i) ret[,i])
        
        }
        
      }
      
    }else{
      
      # get single design matrices
      ret <- extract(mod, "design", expand = TRUE)
      
      if(!split) ret <- do.call("cbind", ret)
          
    }
    
  }
  
  return(ret)
  
}

#' Create test vector
#' 
#' @param obj mboost object
#' @param eps numeric; value to stabilize inversion
#' @return list of test vectors or matrices for each effect
#' @export 
#' 
# returns v^T for linear effects and \tilde{P}_g for group effects
getTestvector <- function(obj, eps = 1e-12)
{
  
  # get design matrix
  X <- getDesignmat(obj, split = T)
  nrcol <- sapply(X, "NCOL")
  n <- nrow(X[[1]])
  inds <- c(1, cumsum(nrcol)[-length(nrcol)]+1)
  inde <- cumsum(nrcol)

  # create list for result
  vT <- vector("list", length(nrcol))
  # get LS hat
  if(any(sapply(X, NCOL)==1))
  {
    
    if(n < sum(nrcol)) stop("p > n.")
    
    Xf <- do.call("cbind", X) 
    B <- crossprod(Xf) + diag(ncol(Xf))*eps
    Xplus <- try(solve(B) %*% t(Xf))
    
    if(class(Xplus)=="try-error"){  
      
      Rchol   <- try(chol(B))
      L1      <- backsolve(Rchol, t(Xf), transpose = TRUE)
      Xplus <- backsolve(Rchol, L1)
      
    }
    
  }
  
  for(j in 1:length(X)){
    
    s <- inds[j]
    # e <- inde[j]
    
    if(ncol(X[[j]])==1) # linear effect
    {  
      
      vT[[j]] <- t(Xplus[s, ]) 
    
    }else{
      
      
      ### handle smooth baselearner?
      
      # if(grepl("bbs", names(obj$baselearner)[j])){
      #   
      #   Xj <- X[[j]]
      #   # Pj <- P[[j]]
      #   # Vinv <- solve(crossprod(Xj) + lambda[2] * Pj)
      #   Vinvj <- Vinv[s:e,s:e]
      #   Hj <- Xj %*% Vinvj %*% t(Xj)
      #   vT[[j]] <- list(H = Hj, X = Xj, Vinv = Vinvj, 
      #                   r = sum(diag(Vinvj%*%crossprod(Xj))))
      #   
      # }else{
        
      ### code by barber for groups (needs to save matrix)
      # X0 = X[[j]]
      # XSh0 = do.call("cbind", X[-j])
      # X0perp = X0 - XSh0%*%solve(t(XSh0)%*%XSh0,t(XSh0)%*%X0) 
      # u = X0perp%*%solve(t(X0perp)%*%X0perp)%*%t(X0perp) 

      # from the selectiveInference package
      svdu_thresh <- function(x) {
        svdx <- svd(x)
        inds <- svdx$d > svdx$d[1] * sqrt(.Machine$double.eps)
        return(svdx$u[, inds, drop = FALSE])
      }
      
        Xj <- X[[j]]
        Xminusj <- svdu_thresh(do.call("cbind", X[-j]))
        vT[[j]] <- t(svdu_thresh(Xj - tcrossprod(Xminusj) %*% Xj))
        
      # }
      
    }
        
    
  }
  
  names(vT) <- names(X)
  
  return(vT)
  
}

getCoefPath <- function(obj, what = c("path", "incr", "sign"))
{
  
  what <- match.arg(what)
  
  signCourse <- if(inherits(obj, "glmboost"))
    sapply(1:mstop(obj), function(m) as.data.frame(obj[m]$coef())) else 
      sapply(1:mstop(obj),function(m) sapply(obj[m]$coef(),"[[",1))
  
  nams <- attr(signCourse[[length(signCourse)]], "names")
  signCourseS <- do.call("rbind",lapply(signCourse, function(sc){
    
    lenSc <- length(sc)
    
    if(lenSc<length(nams)){
      
      namSc <- names(sc)
      namsN <- nams[!nams%in%namSc]
      sc <- c(rep(0,length(namsN)),sc)
      names(sc) <- c(namsN, namSc)
      sc <- sc[nams]
      
    }
    
    unlist(sc)
    
  }))
  
  if(what == "path") return(signCourseS)
  
  signCoursePM <- apply(rbind(rep(0, ncol(signCourseS)), signCourseS), 2, diff)
  
  if(what == "incr") return(signCoursePM)
  
  # what == "sign"
  sig <- if(length(signCoursePM)==1) sign(signCoursePM) else rowSums(sign(signCoursePM))
  return(sig)
  
}

# make base-learners using character names
# blFun specifies the type of BL
makeBL <- function(charName, blFun, data, ...) {
  
  temp <- as.data.frame(data[, charName])
  colnames(temp) <- charName
  bl <- tryCatch(blFun(temp, ...), error=function(e) bols(temp))
  bl$set_names(charName)
  bl
  
}


getFmat <- function(mod, which)
{
  
  X <- extract(mod, which = which, "design")[[1]]
  K <- extract(mod, which = which, "penalty")[[1]]
  l <- extract(mod, which = which, "lambda")[[1]]
  solve(t(X)%*%X + l*K) %*% t(X)
  
}

poi <- function(x){ 
  
  voi <- ceiling(min(x)):floor(max(x))
  unique(sapply(voi, function(xoi) which.min((xoi-x)^2)))
  
}