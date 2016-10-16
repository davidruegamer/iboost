getPolyhedron <- function(mod, Sigma = NULL, testwith="ejTXAk", Ups, weights = NULL) {
  
  ##### general stuff
  
  Y <- mod$response
  n <- length(Y)
  nu <- mod$control$nu
  
  # Does this always work for L2-Boosting?
  X <- as.matrix(do.call("cbind", lapply(mod$baselearner, function(x)x$get_data())))
  
  if(!is.null(weights)) X[which(weights==0),] <- 0
  
  p <- ncol(X)
  selCourse <- selected(mod)
  
  if(missing(Ups)) Ups <- getUpsilons(mod)
  
  if(!is.null(weights)) Ups <- lapply(Ups, function(uu){ uu[which(weights==0),] <- 0; uu})
  
  if(FALSE){
    
    ##### create test vector
    
    testVec <- sapply(unique(selCourse),function(z) max(which(z==selCourse)))
    testVecSort <- testVec[order(selCourse[testVec])]
    testVecMin <- sapply(unique(selCourse),function(z) min(which(z==selCourse)))
    testVecSortMin <- testVecMin[order(selCourse[testVecMin])]
    
  }
  
  ##### create polyhedron
  
  #   W  = lapply(testVecSort, function(t){
  #     
  #     coef = selCourse[t]
  #     nu * hatMatsLeft[[coef]] %*% Lambda[[t]]
  #     
  #   })
  
  # Stuff for polyhedron paper
  
  # signs
  #   s = c(sign(as.numeric(t(X[,selCourse[1]])%*%Y)),
  #         sapply(selCourse[-1],function(k) sign(t(X[,k])%*%Upsilons[[k]]%*%Y)))
  
  signCourse <- sapply(1:mstop(mod),function(m) sapply(mod[m]$coef(),"[[","temp"))
  nams <- attr(signCourse[[length(signCourse)]],"names")
  signCourseS <- do.call("rbind",lapply(signCourse, function(sc){
    
    lenSc <- length(sc)
    
    if(lenSc<length(nams)){
      
      namSc <- names(sc)
      nams<- nams[!nams%in%namSc]
      sc <- c(rep(0,length(nams)),sc)
      names(sc) <- c(nams,namSc)
      sc <- sc[sort(names(sc))]
      
    }
    
    sc
    
  }))
  
  signCoursePM <- apply(rbind(rep(0,length(unique(selCourse))),signCourseS),2,diff)
  sig = if(length(signCoursePM)==1) sign(signCoursePM) else rowSums(sign(signCoursePM))
  
  Gamma <- unlist(lapply(1:length(selCourse), function(i){
    # different indexing: Upsilons = Upsilon[0:(k-1)]
    k = selCourse[i]
    lapply(c(1:p)[-k], function(j){
      rbind(((sig[i]*t(X[,k]))/as.numeric(crossprod(X[,k])) + X[,j]/
               as.numeric(crossprod(X[,j])))%*%Ups[[i]],
            ((sig[i]*t(X[,k]))/as.numeric(crossprod(X[,k])) - X[,j]/
               as.numeric(crossprod(X[,j])))%*%Ups[[i]])
    })
  }), recursive=F)
  
  attr(Gamma, "selected") <- rep(selCourse,each=p-1)
  attr(Gamma, "notSelected") <- unlist(lapply(selCourse,function(k)c(1:p)[-k]))
  names(Gamma) <- paste0("sel",rep(selCourse,each=p-1),"_notsel",unlist(lapply(selCourse,function(k)c(1:p)[-k])))
  
  #   Ups_List <- Upsilons[rep(1:length(Upsilons),each=2*(p-1))]
  #   
  #   Gamma = mapply(function(matInd,listEl) xFacKthStep[matInd,]%*%listEl, 
  #                  matInd=1:nrow(xFacKthStep), listEl=Ups_List, SIMPLIFY=FALSE)
  #   
  #  Gamma <- do.call("rbind",Gamma)
  
  return(Gamma)
  
  if(FALSE){ #### OLD STUFF TO CALCULATE BOUNDARIES
  
  ind <- selCourse[testVecSort]
  dMat <- solve(t(X[,ind])%*%X[,ind])%*%t(X[,ind])
  
  if(testwith=="Z"){ #etaT = W 
  }else{
    
    if(testwith=="ejTXAk") # test all effects at last iteration based on X with all previous selected cov
      # problem: testvector might have changed after the last selection of the ith cov
      etaT <- lapply(order(selCourse[testVecSort]),function(i)dMat[i,,drop=F]) else
        
        etaT <- lapply(testVecSort,function(iter){
          indI <- sort(unique(selCourse[1:iter]))
          i <- which(indI==selCourse[iter])
          (solve(t(X[,indI])%*%X[,indI])%*%t(X[,indI]))[i,,drop=F]
        })
      
  }
  
  if(testwith=="ejTXAk") times <- rep(max(testVecSort),length(testVecSort)) else times <- testVecSort 
  
  Gamma <- lapply(times,function(j)do.call("rbind",Gamma[1:(j*(p-1))]))
  
  #   etaSigEta <- lapply(etaT,function(e)e%*%t(e)*sigma)
  
  # Gamma <- lapply(Gamma,function(g)-1*g) # definition of polyhedron A by Lee & Taylor 2014
  
  # if(is.null(Sigma)) Sigma <- sd(residuals[length(residuals)])#*diag(n)
  
  GammaY <- lapply(Gamma,function(g)g%*%Y)
  vTv <- lapply(1:length(etaT),function(i) as.numeric(tcrossprod(etaT[[i]])))
  
  rho = lapply(1:length(etaT), function(i) as.numeric(Gamma[[i]]%*%t(etaT[[i]])) / vTv[[i]])
  # rho <- rho
  # rhojetaTY <- lapply(1:length(etaT),function(i)rho[[i]]*etaT[[i]]%*%Y)
  
  rhoBigger0 = lapply(rho,function(x)x>0)
  rhoSmaller0 = lapply(rho,function(x)x<0)
  rhoEqual0 = lapply(rho,function(x)x==0)
  
  Vup <- lapply(1:length(rho),function(i){
    
    ind <- rhoSmaller0[[i]]
    a <- rho[[i]][ind]
    aey <- a*etaT[[i]]%*%Y
    gy <- GammaY[[i]][ind]
    v <- (-as.numeric(gy)+aey)/a
    minj <- min(v[v>=as.numeric(etaT[[i]]%*%Y)], na.rm=T)
    j <- which(v==minj)
    return(list(minj=minj,j=j))
    
  })
  
  Vlo <- lapply(1:length(rho),function(i){
    
    ind <- rhoBigger0[[i]]
    a <- rho[[i]][ind]
    aey <- a*etaT[[i]]%*%Y
    gy <- GammaY[[i]][ind]
    v <- (-as.numeric(gy)+aey)/a
    maxj <- max(v[v<=as.numeric(etaT[[i]]%*%Y)], na.rm=T)
    j <- which(v==maxj)
    return(list(maxj=maxj,j=j))
    
  })
  
  V0 <- lapply(1:length(rho),function(i){
    
    ind <- rhoEqual0[[i]]
    gy <- -1*GammaY[[i]][ind]
    max(gy, na.rm=T)<=0
    
  })
  
  return(list(etaT=etaT,
              Vlo=sapply(Vlo,"[[","maxj"),
              Vup=sapply(Vup,"[[","minj"),
              jVlo=sapply(Vlo,"[[","j"),
              jVup=sapply(Vup,"[[","j"),
              V0=V0))
  
  
  # as.numeric(
  #     (
  #       #     pnorm(Vup/sigmaV2NormInv) - 
  #       #       pnorm(as.numeric(vT%*%Y)/sigmaV2NormInv)
  #       pnorm(as.numeric(vT%*%Y)/sigma) - pnorm(Vlo/sigma)
  #     ) / 
  #       (
  #         pnorm(Vup/sigma) - 
  #           pnorm(Vlo/sigma)
  #       )
  #   )
  
  #   return(list(desMats=desMats, lambdas=lambdas, pens=pens, hatMatsLeft=hatMatsLeft,
  #               hatMats=hatMats, Upsilons=Upsilons, residuals=residuals,
  #               Lambda=Lambda, testVec=testVec, testVecSort=testVecSort, W=W,
  #               WX=WX, betas=betas))
  # 
  #   return(W)
  
  }
  
}