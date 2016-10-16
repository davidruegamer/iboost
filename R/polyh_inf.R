polyh_inf <- function(obj, vT, is_congruent, B, var, ncore, alpha)
{
  
  Y <- obj$response
  n <- length(Y)
  nu <- obj$control$nu
  
  # Does this always work for L2-Boosting?
  X <- getDesignmat(obj)
  
  p <- ncol(X)
  selCourse <- selected(obj)
  
  Ups <- getUpsilons(obj)
  
  signCourse <- if(inherits(mod, "glmboost"))
    sapply(1:mstop(obj), function(m) obj[m]$coef()) else 
      sapply(1:mstop(obj),function(m) sapply(obj[m]$coef(),"[[","temp"))
    
  
  
  nams <- attr(signCourse[[length(signCourse)]], "names")
  signCourseS <- do.call("rbind",lapply(signCourse, function(sc){
    
    lenSc <- length(sc)
    
    if(lenSc<length(nams)){
      
      namSc <- names(sc)
      nams<- nams[!nams%in%namSc]
      sc <- c(rep(0,length(nams)),sc)
      names(sc) <- c(nams,namSc)
      sc <- sc[sort(names(sc))]
      
    }
    
    unlist(sc)
    
  }))
  
  signCoursePM <- apply(rbind(rep(0,length(unique(selCourse))), signCourseS), 2, diff)
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
  
  ret <- lapply(1:length(vT), function(vv){
    
    ret1 <- selectiveInference:::poly.pval(y = Y, G = Gamma[[length(Gamma)]], u = 0, v = t(vT[[j]]), sigma = sqrt(var))
    ret2 <- selectiveInference:::poly.int(y = Y, G = Gamma[[length(Gamma)]], u = 0, v = t(vT[[j]]), sigma = sqrt(var),
                                          alpha = alpha)
  
  }
  
  return(ret)
  
}