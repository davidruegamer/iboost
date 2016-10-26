calcUpLow <- function(y, v, M, Sigma=NULL, PvOy=NULL, Pvy=NULL, vTy=NULL)
{
  
  vvT <- tcrossprod(v,v)
  if(is.null(vTy)) vTy <- as.numeric(crossprod(v,y))
  vTv <- as.numeric(crossprod(v,v))
  if(is.null(Pvy)){
    if(is.null(Sigma)) Pvy <- vvT%*%y / vTv else Pvy <- Sigma%*%vvT%*%y / as.numeric(t(v)%*%Sigma%*%v)
  }
  if(is.null(PvOy)) PvOy <- y - Pvy

  vTMPvOy <- sapply(M,function(m)as.numeric(t(Pvy)%*%m%*%PvOy))
  vTMv <- sapply(M,function(m)as.numeric(t(Pvy)%*%m%*%Pvy))
  PvOyTMPvOy <- sapply(M,function(m)as.numeric(t(PvOy)%*%m%*%PvOy))
  D <- sapply(1:length(M),function(i)suppressWarnings(sqrt(vTMPvOy[[i]]^2 - vTMv[[i]]*PvOyTMPvOy[[i]])))
  dev <- sapply(1:length(M),function(i)c((-1*vTMPvOy[[i]]-D[[i]])/(vTMv[[i]]), 
                                           (-1*vTMPvOy[[i]]+D[[i]])/(vTMv[[i]])))
  
  vTy2 <- sqrt(crossprod(Pvy))*sqrt(vTv)
  Vs <- sapply(dev,function(l)l*vTy2)
  
  devFr1 <- sapply(dev,function(d)ifelse(d>1,d-1,1-d))
  
  wUp <- which(devFr1==devFr1[dev>1][which.min(devFr1[dev>1])])
  wLo <- which(devFr1==devFr1[dev<1][which.min(devFr1[dev<1])])
  
  b=ifelse(length(wUp)==0,Inf,Vs[wUp])*sign(vTy)
  a=ifelse(length(wLo)==0,-Inf,Vs[wLo])*sign(vTy)
  
  Vup <- ifelse(sign(vTy)==1,b,a)
  Vlo <- ifelse(sign(vTy)==1,a,b)
  
  
  if(Vup<Vlo) stop("Something went wrong: Vup < Vlo")
  
  return(
    # list(PvOy=PvOy,
    #           Pvy=Pvy,
    #           #               devP=devP,
    #           #               devN=devN,
    #           vTy=vTy,
    #           vTv=vTv,
    #           v=v,
    #           Vs=Vs,
    #           dev=dev,
    #           devFr1=devFr1,
    #           #               Vlos=Vlos,
    #           #               Vups=Vups,
              #ints=
                Intervals(c(Vlo,Vup))
  )#)
  
}


# Version without M being a list -> TODO: merge both appropriately

getIntersection <- function(v,M,y,Sigma=NULL, Pvy=NULL, PvOy=NULL, vTy=NULL)
{
  
  vvT <- tcrossprod(v,v)
  if(is.null(vTy)) vTy <- as.numeric(crossprod(v,y))
  vTv <- as.numeric(crossprod(v,v))
  if(is.null(Pvy)){
    if(is.null(Sigma)) Pvy <- vvT%*%y / vTv else Pvy <- Sigma%*%vvT%*%y / as.numeric(t(v)%*%Sigma%*%v)
  }
  if(is.null(PvOy)) PvOy <- y - Pvy
  
  vTMPvOy <- as.numeric(t(Pvy)%*%M%*%PvOy)
  vTMv <- as.numeric(t(Pvy)%*%M%*%Pvy)
  PvOyTMPvOy <- as.numeric(t(PvOy)%*%M%*%PvOy)
  D <- suppressWarnings(sqrt(vTMPvOy^2 - vTMv*PvOyTMPvOy))
  # calculate tau_{1/2}
  dev <- c((-1*vTMPvOy-D)/(vTMv), (-1*vTMPvOy+D)/(vTMv))
  
  # should be equal to v^Ty but without sign (!)
  vTy2 <- as.numeric(sqrt(crossprod(Pvy))*sqrt(vTv))
  # mutliply tau_{1/2} with v^Ty to get the possible lengths
  # -> if tau == 1, Vs would be = v^Ty
  # -->CHECK: is this doing the right thing?
  Vs <- dev*vTy2
  
  # is it a deviation in the direction of vTy or not
  devFr1 <- sapply(dev,function(d)ifelse(d>1,d-1,1-d))
  
  # which is the smallest deviation bigger than 1 -> wUp
  # and which is the smallest deviation smaller than 1 -> wLo
  wUp <- which(devFr1==devFr1[dev>1][which.min(devFr1[dev>1])])
  wLo <- which(devFr1==devFr1[dev<1][which.min(devFr1[dev<1])])
  
  # boundaries are +- Inf or the lengths (Vs) corresponding to the
  # narrowest interval
  b=ifelse(length(wUp)==0,Inf,Vs[wUp])*sign(vTy)
  a=ifelse(length(wLo)==0,-Inf,Vs[wLo])*sign(vTy)
  
  # change Vup = b to Vup = a (and the same for Vlo) if vTy is negative
  # -> in this case, the maximal deviation is max(dev), which is the
  # factor multiplied with P_v y; if vTy < 0, P_v y is negative
  # and hence, P_v^\perp y + max(dev)*P_v y is a deviation < 0
  Vup <- ifelse(sign(vTy)==1,b,a)
  Vlo <- ifelse(sign(vTy)==1,a,b)
  
  
  if(Vup<Vlo) stop("Something went wrong - Vup < Vlo")
  
  return(list(ints=Intervals(c(Vlo,Vup)),vTy=vTy,v=v))
  
}