#### basic functions

createQ <- function(listOfIminusHatSq, indicesOfChoice, Upsilons, choiceIsBestBL = TRUE, notSelected=NULL)
{
  
  
  ss <-  1:length(listOfIminusHatSq)
  
  Q <- unlist(lapply(1:length(indicesOfChoice),function(i){
    
    selC <- indicesOfChoice[i]
    compPart <- if(choiceIsBestBL) ss[-selC] else notSelected
    
    M1minusM2 <-  lapply(compPart,function(ii) listOfIminusHatSq[[selC]]-listOfIminusHatSq[[ii]]) 
    lapply(M1minusM2, function(m) crossprod(crossprod(m,Upsilons[[i]]),Upsilons[[i]]))
    
  }), recursive = F)
  
  return(Q)
  
  
}

mCreate <- function(mod, conditionOn = c("last","notSelected","first","inBetween"), crit=c("none","GCV","AIC","gMDL"),
                    invertIneq = FALSE, Ups, Sigma=NULL, optimNotSel = FALSE, ...)
{
  
  # stuff for calculations
  hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
  selCourse <- selected(mod)
  # selCov <- lapply(1:mstop(mod),function(i)unique(selected(mod)[1:i]))
  Y <- mod$response
  n <- length(Y)
  p <- length(mod$basemodel)
  w <- mod$`(weights)`
  
  # define vector of iteration indices
  testVec <- if(!conditionOn[1]%in%c("first","inBetween")) 
    sapply(unique(selCourse),function(z) max(which(z==selCourse))) else 
      sapply(unique(selCourse),function(z) min(which(z==selCourse)))
  testVecSort <- testVec[order(selCourse[testVec])]
  
  notSelected <- (1:p)[!(1:p)%in%selCourse]
  
  # calculate ||(I-P)||^2
  IminusHatsq <- lapply(1:p,function(i)crossprod(diag(w)-hatMats[[i]]))
  
  M <- list()
  
  # create conditions resulting from selection
  if(conditionOn[1]!="notSelected"){
  
    M <- createQ(listOfIminusHatSq = IminusHatsq, indicesOfChoice = selCourse[1:max(testVecSort)],
               Upsilons = Ups, choiceIsBestBL = TRUE)
    
  }

  # create conditions resulting from not selected variables  
  if(conditionOn[1]%in%c("inBetween","notSelected") & length(notSelected) > 0){
    
    M <- append(M, 
                createQ(listOfIminusHatSq = IminusHatsq, indicesOfChoice = selCourse[1:max(testVecSort)],
                        Upsilons = Ups, choiceIsBestBL = FALSE, notSelected = notSelected))
    
  }
  
  Madd <- if(crit[1]!="none") mCreateCrit(mod=mod, crit=crit, Ups=Ups, Sigma=Sigma, ...) else NULL
               
  return(list(M=M,Madd=Madd,ind=testVecSort))
  
}

mCreateCrit <- function(mod, crit=c("none","GCV","AIC","gMDL"), Ups, eps=0.01, Sigma=NULL)
{
  
  hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
  # selCov <- lapply(1:mstop(mod),function(i)unique(selected(mod)[1:i]))
  Y <- mod$response
  n <- length(Y)
  p <- length(mod$basemodel)
  w <- mod$`(weights)`
  if(any(w!=1)) stop("Weights !=1 are given. Not implemented yet.")
  
  if(crit[1]!="none"){
    
    
    df <- sapply(1:length(Ups),function(i)sum(diag(diag(n)-Ups[[i]])))
    sigSqEst <- sapply(1:mstop(mod),function(i)as.numeric(crossprod(mod[i]$resid()))/n)
    
    if(crit[1]=="GCV"){
      
      kappa <- sapply(df,function(d)1-d/n)
      vals <- sapply(1:mstop(mod),function(i) sigSqEst[i]/kappa[i]^2)
      
    }else if(crit[1]=="AIC"){ # see Buehlmann & Hothorn
      
      k <- sapply(df,function(d)(1+d/n)/((1-d+2)/n))
      vals <- sapply(1:mstop(mod),function(i) log(sigSqEst[i]) + k[i])
      
    }else if(crit[1]=="gMDL"){  # see Buehlmann & Hothorn
      
      S = sapply(1:mstop(mod),function(i)n*sigSqEst[i]/(n-df[i]))
      F = sapply(1:mstop(mod),function(i)(crossprod(Y)-n*sigSqEst[i])/(df[i]*S[i]))
      vals <- sapply(1:mstop(mod),function(i)log(S[i]) + df[i]/n * log(F[i]))
      
    }else if(is.function(crit)){
      
      # ....
      
    }else{
      
      stop("crit != 'none', but does not correspond to a known criterion.")
      
    }
    
    stopIter <- ifelse(which.min(vals)==mstop(mod),min(which(vals/max(vals)<=eps),
                                                       mstop(mod)),which.min(vals))
    
    seqToC <- (1:mstop(mod))[-stopIter]
    if(crit[1]=="GCV"){
      G <- lapply(1:length(Ups),function(i)crossprod(Ups[[i]])/kappa[i]^2)
    }else{
      
      stop("Not implemented yet.")
      
    }
    Madd <- lapply(seqToC, function(i)G[[stopIter]]-G[[i]])
    warning("Selection Criterion might only be meaningful for a full simulation!")
    
  }
  
  return(list(Madd=Madd,stopIter=stopIter))
  
}