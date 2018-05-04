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

# returns v^T for linear effects and \tilde{P}_g for group effects
getTestvector <- function(obj, eps, makeGroup = FALSE)
{
  
  # get design matrix
  X <- getDesignmat(obj, split = T)
  nrcol <- sapply(X, "ncol")
  n <- nrow(X[[1]])
  inds <- c(1, cumsum(nrcol)[-length(nrcol)]+1)
  inde <- cumsum(nrcol)
  
  # check for linear effects
  if(any(nrcol==1) & !makeGroup){
   
    Xf <- do.call("cbind", X) 
    
    if(n < sum(nrcol)){ # p > n
      
      stop("p > n.")
      # ups <- getUpsilons(obj)
      # Ztilde <- getZtilde(Ups = ups, mod = obj)
      # selmod <- selected(obj)
      # Ztilde_sel <- lapply(1:length(Ztilde), function(i) 10*Ztilde[[i]][[selmod[i]]])
      # Zdebiased <- lapply(Ztilde_sel, function(z) 2*z-z%*%z)
      # r <- sapply(Zdebiased, function(x) sum(diag(x)))
      # ri <- sapply(sort(unique(selected(obj))), function(i) r[max(which(selmod==i))])
      # XtX <- crossprod(Xf)
      # lambdas <- c(t(Xf)%*%resid(mod1) / unlist(coef(mod1)))
      # Xplus <- solve(crossprod(Xf) + diag(lambdas)) %*% t(Xf)
      
    }else{
    
      Xplus <- solve(crossprod(Xf) + diag(ncol(Xf))*eps) %*% t(Xf)
      
    }
  
  }
    
  vT <- vector("list", length(nrcol))
  
  # if(any(nrcol>1) & any(sapply(names(obj$baselearner), function(x) grepl("bbs", x)))){
  # 
  #   P <- extract(obj, "penalty")
  #   ups <- getUpsilons(obj)
  #   Ztilde <- getZtilde(Ups = ups, mod = obj)
  #   selmod <- selected(obj)
  #   Ztilde_sel <- lapply(1:length(Ztilde), function(i) 10*Ztilde[[i]][[selmod[i]]])
  #   Zdebiased <- lapply(Ztilde_sel, function(z) 2*z-z%*%z)
  #   r <- sapply(Zdebiased, function(x) sum(diag(x)))
  #   ri <- sapply(sort(unique(selected(obj))), function(i) r[max(which(selmod==i))])
  #   # options(mboost_lambdaMax = 1e30)
  #   # lambdas <- sapply(1:length(X), function(j) mboost:::df2lambda(X = X[[j]], df = ri[j], 
  #   #                              dmat = P[[j]], weights = rep(1, nrow(X[[j]])))[2])
  #   Xf <- do.call("cbind", X)
  #   # Pf <- bdiag(lapply(1:length(X), function(j) lambdas[j]*P[[j]]))
  #   # lambdaA <- mboost:::df2lambda(X = Xf, df = sum(ri), dmat = Pf, weights = rep(1, nrow(Xf)))[2]
  #   betas <- obj$coef()
  #   resid <- obj$resid()
  #   lambdas <- sapply(1:length(X), function(i) solve(t(betas[[i]])%*%crossprod(P[[i]])%*%betas[[i]]) %*% 
  #                       t(betas[[i]]) %*% t(P[[i]]) %*% t(X[[i]]) %*% resid)
  #   Pf <- bdiag(lapply(1:length(X), function(j) lambdas[j]*P[[j]]))
  #   # mboost:::df2lambda(X = Xf, df = NULL, lambda = 1, dmat = Pf, weights = rep(1, nrow(Xf)))[1]
  #   # lambdaA <- mboost:::df2lambda(X = Xf, df = sum(ri), dmat = Pf, weights = rep(1, nrow(Xf)))[2]
  #   ed <- eigen(crossprod(Xf) + Pf)
  #   ind <- ed$values > 1e-9
  #   vec <- ed$vectors
  #   vec <- t(t(vec[,ind])/sqrt(ed$val[ind]))
  #   Vinv <- tcrossprod(vec)
  #   # H <- Xf%*%Vinv%*%t(Xf)
  #   # sum(ri)
  #   # sum(diag(2*H-H%*%H))
  #   # Fi <- Vinv%*%crossprod(Xf)
  #   # sum(diag(2*Fi-Fi%*%Fi))
  #   # (ratio <- ri / c(sum(diag(2*Fi-Fi%*%Fi)[1:24]), sum(sum(diag(2*Fi-Fi%*%Fi)[25:31]))) )
  #   # lambdas[1] <- lambdas[1]*ratio[2]/ratio[1]
  #   # Pf <- bdiag(lapply(1:length(X), function(j) lambdas[j]*P[[j]]))
  #   # lambdaA <- mboost:::df2lambda(X = Xf, df = sum(ri), dmat = Pf, weights = rep(1, nrow(Xf)))[2]
  #   # Vinv <- solve(crossprod(Xf) + lambdaA*Pf + diag(ncol(Xf))*1e-10)
  #   
  # }
      
  for(j in 1:length(X)){
    
    s <- inds[j]
    e <- inde[j]
    
    if(ncol(X[[j]])==1 & !makeGroup) # linear effect
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

# make Baselearners using character names
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