# creates v-testvector for OLS
vOLScreate <- function(iter, designmat, selCourse, n)
{
  
  selCov <- lapply(1:length(selCourse),function(i)unique(selCourse[1:i]))
  isFac <- FALSE
  
  # check for factor variables
  lenCheck <- sapply(designmat,nrow)!=n
  
  if(any(lenCheck)){
    
    for(j in which(lenCheck)){
      
      ind <- diff(c(as.numeric(rownames(designmat[[j]])),n+1))
      designmat[[j]] <- designmat[[j]][rep(1:nrow(designmat[[j]]),ind),]
      
    }
    
  }
  
  Xx <- do.call("cbind",designmat)
  colnames(Xx) <- rep(1:length(designmat),sapply(designmat,ncol))
  Xplus <- (Xx%*%solve(crossprod(Xx)))
  v <- Xplus[,colnames(Xplus)%in%which(selCourse[iter]==selCov[[iter]])] 
  #Xplus[,which(selCourse[i]==selCov[[i]])]
  
  listC <- !is.null(dim(v))
  if(listC){  
    #problV <- v[listC]
    #names(v) <- 1:length(v)
    #v[listC] <- lapply(problV,function(p)split(p,col(p)))
    v <- split(v,col(v))
    #v <- append(v[!listC],unlist(v[listC],recursive = F))
    #v <- v[order(names(v))]
    # isFac <- TRUE
  }
  
  return(v)
  
}

#### Splines

vBBScreate <- function(mod, which)
{
  
  X <- extract(mod, which = which, "design")[[1]]
  K <- extract(mod, which = which, "penalty")[[1]]
  l <- extract(mod, which = which, "lambda")[[1]]
  Xplus <- solve(t(X)%*%X + l*K)%*%t(X)
  
  return(lapply(1:nrow(Xplus),function(i)Xplus[i,]))
  
}

makeTestVecForSplines <- function(mod, which, isBBSC = FALSE)
{
  
  xvalDF <- mod$model.frame(which = which)[[1]]
  xvals <- seq(min(xvalDF[[1]]), max(xvalDF[[1]]), l=nrow(xvalDF))
  xvalDF[,1] <- xvals
  # yvals <- predict(mod, which = which, newdata = xvalDF)
  
  argsFromEnv <- environment(mod$baselearner[[which]]$dpp)$args
  
  Xnew <- mboost:::bsplines(xvals,
                            knots = argsFromEnv$knots[[1]]$knots,
                            boundary.knots =argsFromEnv$knots[[1]]$boundary.knots,
                            degree = argsFromEnv$degree,
                            Ts_constraint = argsFromEnv$Ts_constraint,
                            deriv = argsFromEnv$deriv, 
                            extrapolation = FALSE)
  
  if(isBBSC) Xnew <- Xnew %*% argsFromEnv$Z
    
  # # for other covariates, which have been selected in the model before which
  # if(selected(mod)[1]==which) Xcorr <- diag(nrow(xvalDF)) else
  #   Xcorr <- 
  
  X <- extract(mod, which = which, "design")[[1]]
  K <- extract(mod, which = which, "penalty")[[1]]
  l <- extract(mod, which = which, "lambda")[[1]]
  Xplus <- solve(t(X)%*%X + l*K)%*%t(X)
  
  testVec <- Xnew %*% Xplus
  
  return(list(xvalDF = xvalDF,
              testVec = testVec, 
              Xnew = Xnew,
              # Xcorr = Xcorr,
              Xplus = Xplus))
  
}