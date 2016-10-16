truncPnorm <- function(a,b,mu,sd,x,log.p=FALSE)
{
  #options(digits = 22)
  nom <- (pnorm(as.numeric(x-mu)/as.numeric(sd),log.p = log.p)-
            pnorm(as.numeric(a-mu)/as.numeric(sd),log.p = log.p))
  denom <- (pnorm(as.numeric(b-mu)/as.numeric(sd),log.p = log.p)-
              pnorm(as.numeric(a-mu)/as.numeric(sd),log.p = log.p))
  
  ind <- denom==0
  approxVals <- NULL
  if(any(ind) & !log.p) approxVals <- truncPnorm(a,b,mu[ind],sd,x[ind],log.p=TRUE)
  if(any(ind) & log.p){ 
    if(x[ind]>mu) approxVals <- 0 else approxVals <- 1
  }
  
  return( c((nom / denom)[!ind], approxVals ) )
}

truncRnorm <- function(a,b,mu,sd,x,log.p=FALSE)
{
  #options(digits = 22)
  nom <- (pnorm(as.numeric(x-mu)/as.numeric(sd),log.p = log.p)-
            pnorm(as.numeric(a-mu)/as.numeric(sd),log.p = log.p))
  denom <- (pnorm(as.numeric(b-mu)/as.numeric(sd),log.p = log.p)-
              pnorm(as.numeric(a-mu)/as.numeric(sd),log.p = log.p))
  
  ind <- denom==0
  approxVals <- NULL
  if(any(ind) & !log.p) approxVals <- truncPnorm(a,b,mu[ind],sd,x[ind],log.p=TRUE)
  if(any(ind) & log.p){ 
    if(x[ind]>mu) approxVals <- 0 else approxVals <- 1
  }
  
  return( c((nom / denom)[!ind], approxVals ) )
}


twoFractionPnorm <- function(a1,b1,a2,b2,mu,sd,x,log.p=FALSE)
{
  
  if(all(is.na(x))) return(NA)
  if(length(unique(x))==2 & any(is.na(x))) x <- unique(x)[!is.na(unique(x))]
  
  if(a2 <= x & x <= b2){
    
    a2t <- a1
    b2t <- b1
    
    a1 <- a2
    b1 <- b2
    a2 <- a2t
    b2 <- b2t
    

  }
  
  #options(digits = 22)
  nom <- ((pnorm(as.numeric(x-mu)/as.numeric(sd),log.p = log.p)-
            pnorm(as.numeric(a1-mu)/as.numeric(sd),log.p = log.p)) + 
            (pnorm(as.numeric(b2-mu)/as.numeric(sd),log.p = log.p)-
               pnorm(as.numeric(a2-mu)/as.numeric(sd),log.p = log.p)))
  denom <- ((pnorm(as.numeric(b1-mu)/as.numeric(sd),log.p = log.p)-
              pnorm(as.numeric(a1-mu)/as.numeric(sd),log.p = log.p)) + 
              (pnorm(as.numeric(b2-mu)/as.numeric(sd),log.p = log.p)-
                 pnorm(as.numeric(a2-mu)/as.numeric(sd),log.p = log.p)))
  
  ind <- denom==0
  approxVals <- NULL
  if(any(ind) & !log.p) approxVals <- twoFractionPnorm(a1,b1,a2,b2,mu[ind],sd,x[ind],log.p=TRUE)
  if(any(ind) & log.p){ 
    if(x[ind]>mu) approxVals <- 0 else approxVals <- 1
  }
  
  return( c((nom / denom)[!ind], approxVals ) )
}

inInterval <- function(value,int)int[[1]]<=value & int[[2]]>=value

moreFractionPnorm <- function(ints,mu,sd,x,log.p=FALSE)
{
  
  if(all(is.na(x))) return(NA)
  if(length(unique(x))==2 & any(is.na(x))) x <- unique(x)[!is.na(unique(x))]
  
  whichIntElemX <- apply(ints,1,function(oneIn)inInterval(x,oneIn))
  
  bs <- ints[,2]
  as <- ints[,1]
  
  bsNorm <- sapply(bs,function(b)pnorm(as.numeric(b-mu)/as.numeric(sd),log.p = log.p))
  asNorm <- sapply(as,function(a)pnorm(as.numeric(a-mu)/as.numeric(sd),log.p = log.p))
  
  denom <- sum(sapply(1:length(bs),function(i)bsNorm[i]-asNorm[i]))
  bsNorm[whichIntElemX] <- pnorm(as.numeric(x-mu)/as.numeric(sd),log.p = log.p)
  nom <- sum(sapply(1:length(bs),function(i)bsNorm[i]-asNorm[i]))

  ind <- denom==0
  approxVals <- NULL
  if(any(ind) & !log.p) approxVals <- moreFractionPnorm(ints,mu[ind],sd,x[ind],log.p=TRUE)
  if(any(ind) & log.p){ 
    if(x[ind]>mu) approxVals <- 0 else approxVals <- 1
  }
  
  return( c((nom / denom)[!ind], approxVals ) )
}


getCIbounds <- function(ints, sd, y, etaT, alpha=0.05, searchInPossDev=10, 
                        bySI=FALSE, gridpts=1000, griddepth=2, nrIterationsUniroot = 1000)
{

  
  getTk <- function(Rk)2*min(Rk,1-Rk)
  
  if(bySI) library("selectiveInference")

  muHat <- as.numeric(etaT%*%y)
  testSeq <- seq(muHat-searchInPossDev*sd,muHat+searchInPossDev*sd,length.out=1e6)
  
  
  if(nrow(ints)>1){
    
    if(nrow(ints)==2) warning("Got a combined interval") else stop("More than 2 intervals - not implemented yet.")
    Vlo1 <- ints[[1]]
    Vlo2 <- ints[[2]]
    Vup1 <- ints[[3]]
    Vup2 <- ints[[4]]
    
    
    funEval <- function(x) twoFractionPnorm(x=muHat, 
                                      mu=x,
                                      sd=as.numeric(sd),
                                      a1=Vlo1,
                                      b1=Vup1,
                                      a2=Vlo2,
                                      b2=Vup2)
    
    inverseF = function (f, lower = -10, upper = 10) {
      function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper,
                           extendInt = "yes", maxiter = nrIterationsUniroot)[1]
    }
    
    testVal <- funEval(testSeq)
    considerTestSeq <- testSeq[!is.nan(testVal) & testVal!=0 & testVal!=1]
    lowB <- min(considerTestSeq, na.rm=T)
    uppB <- max(considerTestSeq, na.rm=T)
    
    
    quantfun <- inverseF(f=funEval, lower = lowB, upper = uppB)
    
    upp = tryCatch(as.numeric(quantfun(alpha/2)),error=function(e)return(NA))
    low = tryCatch(as.numeric(quantfun(1-alpha/2)),error=function(e)return(NA))
  
  Rk = tryCatch(twoFractionPnorm(a1=Vlo1,b1=Vup1,a2=Vlo2,b2=Vup2,sd=sd,mu=0,x=muHat),
                error=function(e)return(NA))
  
  return(list(KI=c(low,upp),Tk=#Rk))# 
                getTk(Rk)))
    
  }else{
    
    Vlo <- ints[[1]]
    Vup <- ints[[2]]
  
  funEval <- function(x) truncPnorm(x=muHat, 
                                    mu=x,
                                    sd=as.numeric(sd),
                                    a=Vlo,
                                    b=Vup)

  if(bySI){
    
    z=muHat
    sd=as.numeric(sd)
    vlo=Vlo
    vup=Vup
    bits=NULL
    fun = function(x) { selectiveInference:::tnorm.surv(z,x,sd,vlo,vup,bits) }
    int = selectiveInference:::grid.search(testSeq,fun,alpha/2,1-alpha/2,gridpts,griddepth)
    low <- int[1]
    upp <- int[2]
    
    
  }else{
  
    inverseF = function (f, lower = -10, upper = 10) {
      function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper,
                           extendInt = "yes", maxiter = nrIterationsUniroot)[1]
    }
    
    testVal <- funEval(testSeq)
    considerTestSeq <- testSeq[!is.nan(testVal) & testVal!=0 & testVal!=1]
    lowB <- min(considerTestSeq, na.rm=T)
    uppB <- max(considerTestSeq, na.rm=T)


    quantfun <- inverseF(f=funEval, lower = lowB, upper = uppB)
  
    upp = tryCatch(as.numeric(quantfun(alpha/2)),error=function(e)return(NA))
    low = tryCatch(as.numeric(quantfun(1-alpha/2)),error=function(e)return(NA))

  }

  
  
  Rk = tryCatch(truncPnorm(a=Vlo,b=Vup,sd=sd,mu=0,x=muHat),
                error=function(e)return(NA))
  
  return(list(KI=c(low,upp),Tk=#Rk))# 
              getTk(Rk)))
  
  }
  
}

getPval <- function(ints, sd, y=NULL, testvecT = NULL, muHat = as.numeric(testvecT%*%y), noTrans = FALSE)
{
  
  getTk <- function(Rk)2*min(Rk,1-Rk)
  if(is.null(muHat)) muHat <- as.numeric(testvecT%*%y)
  
  if(nrow(ints)==1){
  
    Vlo <- ints[[1]]
    Vup <- ints[[2]]
  
    Rk = tryCatch(truncPnorm(a=Vlo,b=Vup,sd=sd,mu=0,x=muHat),
                error=function(e)return(NA))
  }else{
    
#     Vlo1 <- ints[[1]]
#     Vlo2 <- ints[[2]]
#     Vup1 <- ints[[3]]
#     Vup2 <- ints[[4]]
#     
#     Rk <- tryCatch(twoFractionPnorm(a1=Vlo1,b1=Vup1,a2=Vlo2,b2=Vup2,sd=sd,mu=0,x=muHat),
#                    error=function(e)return(NA))
#                    
#   }else{
    
    Rk <- tryCatch(moreFractionPnorm(ints=ints, sd=sd, mu=0, x=muHat),
                   error=function(e)return(NA))
    
  }
  
  if(noTrans) return(Rk) else return(getTk(Rk))
  
  
}
