getLimitsFakeModel <- function(mod, listOfnewSelCourses = NULL, v, nrCores = 1, Sigma = NULL,
                               addNonSel = FALSE, sampleNr=0, Seed=1234, estSd=1, newYdis = 1,
                               conditionOn, diffCondForEachVar = FALSE, ...)
{
  
  if(addNonSel & diffCondForEachVar) stop("Combining addNonSel and diffCondForEachVar not meaningful.")
  if(!is.list(v)) v <- list(v)
  
  hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
  selCourse <- selected(mod)
  p <- length(mod$basemodel)
  Y <- mod$response
  n <- length(Y)
  Ups <- getUpsilons(mod)
  mst <- mstop(mod)
  
  set.seed(Seed)
  checkInputList <- !is.null(listOfnewSelCourses)
  if(!checkInputList & sampleNr > 0) 
    listOfnewSelCourses <- lapply(1:sampleNr,function(x)sample(selCourse))
  if(!checkInputList & sampleNr == 0) 
    listOfnewSelCourses <- lapply(v,function(vv)educatedGuess(mod = mod, v = vv, 
                                                              Sigma = Sigma, timesDistance = newYdis))
  
  if(length(v)>1 & sampleNr == 0 & !checkInputList) 
    listOfnewSelCourses <- unlist(listOfnewSelCourses, recursive = F)
  
  ss <- 1:p
  M <- NULL  
  res <- vector("list",length(v))
  notSelected <- ss[!ss%in%unique(selCourse)]
  
  testVec <- sapply(unique(selCourse),function(z) max(which(z==selCourse)))
  testVecSort <- testVec[order(selCourse[testVec])]
  
  IminusHatsq <- lapply(1:p,function(i)crossprod(diag(n)-hatMats[[i]]))
  if(!is.list(listOfnewSelCourses)) listOfnewSelCourses <- list(listOfnewSelCourses)
  fakeUps <- mclapply(listOfnewSelCourses, function(sc) getUpsilons(mod, otherSelCourse = sc), mc.cores = nrCores)
  
  if(addNonSel){
    
    if(is.list(listOfnewSelCourses)){
      
      if(length(listOfnewSelCourses)!=1) stop("Adding non-selected not implemented for lists") else
        listOfnewSelCourses <- unlist(listOfnewSelCourses)
      
    }
    
    if(checkInputList & length(listOfnewSelCourses)!=mst){
      
      warning("Using fillPath")
      selCourse <- fillPath(mod,listOfnewSelCourses,
                            fakeUps[[1]][[length(listOfnewSelCourses)]])
    
    }
      
    Mns <- createQ(listOfIminusHatSq = IminusHatsq, 
                   indicesOfChoice = selCourse, 
                   Upsilons = Ups, choiceIsBestBL = FALSE, notSelected=notSelected)
    
     
  }
    
  M <- lapply(1:length(fakeUps), function(j) 
              mCreate(mod, conditionOn = conditionOn, Ups=fakeUps[[j]], Sigma=Sigma, 
                      customSelCourse = listOfnewSelCourses[[j]]))
  
  if(!diffCondForEachVar) M <- lapply(M,"[[","M")
  if(addNonSel) M <- lapply(M,function(M)append(M,Mns))
  
  res <- mclapply(1:length(v),function(iv){
    
    if(diffCondForEachVar)
      combineIntervals(lapply(M,function(m)
        calcUpLow(v=v[[iv]],M=m$M[1:(m$ind[iv]*(p-1))],y=Y,Sigma=Sigma)[c("ints","vTy")]),
                       what="union") else
                         combineIntervals(lapply(M,function(m)
                           calcUpLow(v=v[[iv]],M=m,y=Y,Sigma=Sigma)[c("ints","vTy")]),
                                          what="union")
    
  }, mc.cores = nrCores)
  
  return(list(vals=res,testvec=v))
  
}

educatedGuess <- function(mod, v, timesDistance=1, Sigma=NULL)
{
  
  y <- mod$response
  orgSelC <- selected(mod)
  vvT <- tcrossprod(v,v)
  vTv <- as.numeric(crossprod(v,v))
  Pvy <- if(is.null(Sigma)) vvT%*%y / vTv else Sigma%*%vvT%*%y / t(v)%*%Sigma%*%v
  # PvOy <- y - Pvy
  newYs <- lapply(timesDistance, function(t) y + t * Pvy)
  
  newSelCs <- mclapply(newYs,function(yy)selected(mboost_fit(blg=mod$baselearner,
                                  response=as.numeric(yy),
                                  control = boost_control(mstop=mstop(mod),
                                                          nu=mod$control$nu))),
                       mc.cores=2)
  
  newSelCs <- lapply(newSelCs,function(nsc){
    
    allVars <- unique(selected(mod))
    elemsToInsert <- allVars[!allVars%in%nsc]
    if(length(elemsToInsert)>0) nsc[which(duplicated(nsc))[1:length(elemsToInsert)]] <- elemsToInsert
   
    nsc
     
  })
  
  return(newSelCs)
  
}


findBestPath <- function(mod, v, Sigma=NULL, estSd, addRestByBoosting = TRUE)
{
  
  selCourse <- selected(mod)
  selCov <- unique(selCourse)
  this.mstop <- mstop(mod)
  mstopIter <- mstop(mod)
  newPath <- c()
  notSelectedFake <- selCov
  n <- length(mod$ustart)
  nu <- mod$control$nu
  hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
  upsPart <- lapply(hatMats,function(x)(diag(n)-x*nu))
  UpsTillStep <- diag(n)
  ss <- 1:length(mod$basemodel)
  notSelected <- ss[!ss%in%selCov]
  
  remainingSelCov <- selCov
  
  for(i in 1:length(selCov)){
   
    chooseCov <- ifelse(i==length(selCov),1, 
      findBestCovInStep(mod = mod, v = v, estSd=estSd, Sigma=Sigma,
                        remainingCov = remainingSelCov, UpsTillStep = UpsTillStep)
    )
    
    newPath <- c(newPath,remainingSelCov[chooseCov])
    remainingSelCov <- remainingSelCov[remainingSelCov!=newPath[i]]
    UpsTillStep <- upsPart[[newPath[i]]] %*% UpsTillStep
    
  }
  
  if(addRestByBoosting){
   
    resids <- UpsTillStep%*%mod$response
    newPath <- c(newPath,
                 selected(mboost_fit(blg=mod$baselearner[selCov],
                                     response=as.numeric(resids),
                                     control = boost_control(mstop=mstop(mod),
                                                             nu=mod$control$nu)))
    )
                 
  }else{
  
#   for(i in (length(selCov)+1):this.mstop){
#     
#     chooseCov <- lapply(notSelected, function(ns) 
#       findBestComparison(mod, v, whichIsNotSelected = ns, estSd = estSd,
#                          Sigma = Sigma, UpsTillStep = UpsTillStep)
#     
#   }
  
  }
    
  return(newPath)
  
}

findBestCovInStep <- function(mod, v, remainingCov, UpsTillStep, Sigma, estSd)
{
  
  Y <- mod$response
  n <- length(Y)
  ss <- 1:length(mod$basemodel)
  hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
  IminusHatsq <- lapply(1:p,function(i)crossprod(diag(n)-hatMats[[i]]))
  
  M <- lapply(remainingCov, function(selC){
  
    M1minusM2 <- lapply(ss[-selC],function(ii) IminusHatsq[[selC]]-IminusHatsq[[ii]])
    lapply(M1minusM2, function(m) crossprod(crossprod(m,UpsTillStep),UpsTillStep))
  
  })

  strictestPerSelC <- lapply(M,function(m)calcUpLow(y=Y,v=v,M=m,Sigma=Sigma)[c("ints","vTy")])
  
  if(is.null(Sigma)) Sigma=as.numeric(sqrt(crossprod(v))*estSd)
  bestBoundary <- which.max(sapply(lapply(strictestPerSelC,function(iv)combineIntervals(iv[c("ints","vTy")], what="union")),
                         size))
  return(bestBoundary)
  
}

findBestComparison <- function(mod, v, whichIsNotSelected, UpsTillStep, Sigma, estSd)
{
  
  
  Y <- mod$response
  n <- length(Y)
  ss <- 1:length(mod$basemodel)
  hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
  IminusHatsq <- lapply(1:p,function(i)crossprod(diag(n)-hatMats[[i]]))

  M1minusM2 <- lapply(ss[-whichIsNotSelected],function(ii) IminusHatsq[[ii]]-IminusHatsq[[whichIsNotSelected]])
  
  M <- lapply(M1minusM2, function(m) crossprod(crossprod(m,UpsTillStep),UpsTillStep))
    
  bdr <- lapply(M,function(m)getIntersection(y=Y,v=v,M=m,Sigma=Sigma))
  
  if(is.null(Sigma)) Sigma=as.numeric(sqrt(crossprod(v))*estSd)
  bestBoundary <- which.max(sapply(bdr,function(iv)combineIntervals(iv, what="union"), size))
  
  return(bestBoundary)
  
}
  
  
# 
# 
# findMarginMaximizer <- function(mod, y=mod$response, v, lookAtSpecificSet = NULL,
#                                 chooseBL = 1, otherSelCourse=NULL, Sigma=NULL,
#                                 returnRparts = FALSE, Ups=NULL, fixedSecond = FALSE,
#                                 onlyConsiderBetterBLs = FALSE)
# {
#   
#   if(fixedSecond & !onlyConsiderBetterBLs) warning("Fixed the second BL with onlyConsiderBetterBLs = FALSE.")
#   
#   n <- length(mod$ustart)
#   p <- length(mod$basemodel)
#   ss=1:p
#   if(is.null(Ups)){
#     resid <- if(is.null(otherSelCourse)) mod$resid() else 
#       getUpsilons(mod, otherSelCourse=otherSelCourse)[[mstop(mod)]]%*%y
#   }else{
#     resid <- Ups%*%y
#   }
#   
#   vvT <- tcrossprod(v,v)
#   vTv <- as.numeric(crossprod(v,v))
#   rpar <- if(is.null(Sigma)) vvT%*%resid / vTv else 
#     Sigma%*%vvT%*%resid / t(v)%*%Sigma%*%v
#   rort <- resid - rpar
#   
#   hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
#   RortSq <- sapply(hatMats,function(h)crossprod(h%*%rort))
#   RparSq <- sapply(hatMats,function(h)crossprod(h%*%rpar))
#   
#   if(returnRparts) return(list(RortSq=RortSq, RparSq=RparSq))
#   
#   if(onlyConsiderBetterBLs){
#     
#     quadErr <- sapply(hatMats,function(h)crossprod(resid-h%*%resid))
#     ss <- ss[quadErr<=quadErr[chooseBL]]
#     
#   }
#   
#   if(!fixedSecond)
#     vals <- sapply(ss[-chooseBL],function(i)(RparSq[chooseBL]-RparSq[i])*(RortSq[i] - RortSq[chooseBL]))
#   else
#     vals <- sapply(ss[-chooseBL],function(i)(RparSq[i]-RparSq[chooseBL])*(RortSq[chooseBL] - RortSq[i]))
#   
#   names(vals) <- ss[-chooseBL]
#   return(vals)
#   
# }
# 
# 
# 
# 
# backtrackSearch <- function(mod, initialSelCourse = unique(selected(mod)), v, 
#                             nrCores = 25, Sigma = NULL, startPart=round(length(initialSelCourse)/1),
#                             maxNumberOfPathsToExtend = 50, estSd=1)
# {
# 
#     
#     initialBounds <- getLimitsFakeModel(mod = mod, v=v, listOfnewSelCourses = initialSelCourse[1:startPart],
#                                         addNonSel = FALSE, conditionOn = "first")$vals[[1]]
#     diffLow <- initialBounds["vTy"]-initialBounds["Vlo"]
#     diffUp <- initialBounds["Vup"]-initialBounds["vTy"]
#     betterBound <- ifelse(diffLow > diffUp, 1,2)
#     diff <- c(diffLow,diffUp)[betterBound]
#     
#     candidates = unique(selected(mod))
#     paths = as.list(candidates)
#     
#     while(all(sapply(paths,length)<length(candidates))){
#       
#       lim <- mclapply(paths, function(course)
#          getLimitsFakeModel(mod = mod, v=v, listOfnewSelCourses = course,
#                            addNonSel = FALSE, conditionOn = "first")$vals[[1]],
#          mc.cores=nrCores)
#       lims <- sapply(lim,"[[",betterBound)
#       if(length(paths[[1]])>=startPart){ 
#         
#         limsBest <- max(lims)
#         print(paste0(diff," vs ", limsBest))
#         diff <- max(diff,limsBest)
#         
#       }
#         
#       cond <- sapply(lim,function(l)abs(l[betterBound]-l[3])>diff)
#       if(sum(cond)==0 & length(paths)!=length(candidates)) return(initialBounds)
#       paths <- paths[cond]
#       print(sum(cond))
#       if(length(paths)>100) paths <- paths[order(lims[cond], decreasing = TRUE)][1:maxNumberOfPathsToExtend]
#       paths <- unlist(lapply(paths,function(courseBefore)
#         lapply(candidates[!candidates%in%courseBefore], function(ccc)
#           c(courseBefore,ccc))),recursive = F)
#       print(length(paths))
#       print(length(paths[[1]]))
#       
#     }
#     
#     valsFin <- combineIntervals(lapply(lim,as.list), wmin = which.max, wmax = which.min, Sigma=sqrt(crossprod(v))*estSd)
#     
#     return(valsFin)
#   
# }

backtrackSearch2 <- function(mod, initialSelCourse = unique(selected(mod)), vList, 
                             nrCores = 25, Sigma = NULL, startPart=round(length(initialSelCourse)/2),
                             maxNumberOfPathsToExtend = length(initialSelCourse), estSd=1, addNonSelAtEnd = FALSE,
                             returnResult = FALSE, resultIn = NULL, addCov = NULL)
{
  
  initialBounds <- getLimitsFakeModel(mod = mod, v=vList, 
                                      listOfnewSelCourses = initialSelCourse[1:startPart],
                                      addNonSel = FALSE, 
                                      conditionOn = "first", 
                                      nrCores=nrCores)$vals
  diffLow <- sapply(initialBounds,function(iB)iB["vTy"]-iB["Vlo"])
  diffUp <- sapply(initialBounds,function(iB)iB["Vup"]-iB["vTy"])
  betterBound <- sapply(1:length(diffLow),function(i)ifelse(diffLow[i] > diffUp[i], 1,2))
  diff <- sapply(1:length(betterBound),function(i)c(diffLow[i],diffUp[i])[betterBound[i]])
  
  if(is.null(resultIn)){
  
    candidates = unique(selected(mod))
    paths = as.list(candidates)
    
    while(length(paths[[1]])<length(candidates)){
      
      lim <- mclapply(paths, function(course)
        getLimitsFakeModel(mod = mod, v=vList, listOfnewSelCourses = course,
                           addNonSel = FALSE, conditionOn = "first")$vals,
        mc.cores=nrCores)
      
      cond <- checkCond(lim, betterBound, diff)# 
      tempCond <- apply(cond,2,function(cc)sum(cc>0)==(length(vList))^2)
      if(all(!tempCond) | length(tempCond)>maxNumberOfPathsToExtend){
        
        ocond <- order(apply(cond,2,function(cc)sum(cc)),decreasing = TRUE)
        nrElem <- min(maxNumberOfPathsToExtend,length(ocond))
        cond <- ocond[1:nrElem]
        print(nrElem)
        
        }else{
          
          cond <- tempCond
          print(sum(cond>0))
          
        }
      
      
      paths <- paths[cond]
      paths <- unlist(lapply(paths,function(courseBefore)
        lapply(candidates[!candidates%in%courseBefore], function(ccc)
          c(courseBefore,ccc))),recursive = F)
      print(length(paths))
      print(length(paths[[1]]))
      
    }
  
  }else paths <- resultIn
    
  if(returnResult) return(paths)
  
  
  if(!is.list(paths)) paths <- list(paths)
  
  diffLen <- mstop(mod)-length(paths[[1]])
  
#   paths <- unlist(lapply(paths,function(pp) lapply(unique(selected(mod)),function(cov)
#     c(pp,rep(cov,diffLen)))),recursive = F)

  
  
  if(!is.null(addCov)) paths <- lapply(paths,function(pp)c(pp,rep(selected(mod)[addCov],diffLen)))
  
  lim <- mclapply(paths, function(course)
                    getLimitsFakeModel(mod = mod, v=vList, listOfnewSelCourses = course,
                                       addNonSel = TRUE, conditionOn = "first")$vals,
                  mc.cores=nrCores)
  
  
  cond <- checkCond(lim, betterBound, diff)# 
  ocond <- order(apply(cond,2,function(cc)sum(cc)),decreasing = TRUE)
  paths <- paths[ocond]

  valsFin <- lapply(1:length(vList),function(i)
    combineIntervals(lapply(lim,function(l)as.list(l[[i]])), what="union"))
  
  return(valsFin)
  
}

checkCond <- function(listOfLists, indicatorListItem, comparison)
{
  
  sapply(listOfLists, function(lll){
    
    sapply(1:length(lll),function(i)abs(lll[[i]][indicatorListItem]-
                                              lll[[i]][3])-comparison[i])
    
  })
  
}

fillPath <- function(mod,path,UpsTillStep)
{
  
  selCov <- unique(selected(mod))
  
  resids <- UpsTillStep%*%mod$response
  path <- c(path, selected(mboost_fit(blg=mod$baselearner[selCov],
                                      response=as.numeric(resids),
                                      control = boost_control(mstop=mstop(mod)-length(path),
                                                              nu=mod$control$nu)))
  )
  
  return(path)
  
}

randomBoundaryExtension <- function(mod, numberOfPathsToAdd = 100, nrCores = 1)
{
  
  selCourse <- selected(mod)
  covarsMin <- min(sapply(unique(selCourse),function(z) min(which(z==selCourse))))
  
  #...
  
  selCourse <- fillPath(mod,listOfnewSelCourses,
                        fakeUps[[1]][[length(listOfnewSelCourses)]])
    
  #...
  
}


backtrackSearch3 <- function(mod, initialSelCourse = unique(selected(mod)), vList, 
                             nrCores = 25, Sigma = NULL, startPart=round(length(initialSelCourse)/2),
                             maxNumberOfPathsToExtend = length(initialSelCourse), estSd=1, addNonSelAtEnd = FALSE,
                             returnResult = FALSE, resultIn = NULL, addCov = NULL, condOn = "last")
{
  
  initialBounds <- getLimitsFakeModel(mod = mod, v=vList, 
                                      listOfnewSelCourses = initialSelCourse[1:startPart],
                                      addNonSel = FALSE, 
                                      conditionOn = condOn, 
                                      nrCores=nrCores)$vals
  sizeInit <- sapply(initialBounds,function(iB)size(iB$ints))
  
  if(is.null(resultIn)){
    
    candidates = unique(selected(mod))
    paths = as.list(candidates)
    
    while(length(paths[[1]])<length(candidates)){
      
      lim <- mclapply(paths, function(course)
        getLimitsFakeModel(mod = mod, v=vList, listOfnewSelCourses = course,
                           addNonSel = FALSE, conditionOn = condOn)$vals,
        mc.cores=nrCores)
      
      cond <- sapply(lim,function(ll)sum(sapply(lapply(ll,"[[","ints"),size)>sizeInit))
    
      ocond <- order(cond,decreasing = TRUE)
      nrElem <- min(maxNumberOfPathsToExtend,length(ocond))
      cond <- ocond[1:nrElem]
        
      paths <- paths[cond]
      paths <- unlist(lapply(paths,function(courseBefore)
        lapply(candidates[!candidates%in%courseBefore], function(ccc)
          c(courseBefore,ccc))),recursive = F)
      print(length(paths))
      print(length(paths[[1]]))
      
    }
    
  }else paths <- resultIn
  
  if(returnResult) return(paths)
  
  
  if(!is.list(paths)) paths <- list(paths)
  
  diffLen <- mstop(mod)-length(paths[[1]])
  
  #   paths <- unlist(lapply(paths,function(pp) lapply(unique(selected(mod)),function(cov)
  #     c(pp,rep(cov,diffLen)))),recursive = F)
  
  
  
  if(!is.null(addCov)) paths <- lapply(paths,function(pp)c(pp,rep(selected(mod)[addCov],diffLen)))
  
  lim <- mclapply(paths, function(course)
    getLimitsFakeModel(mod = mod, v=vList, listOfnewSelCourses = course,
                       addNonSel = addNonSelAtEnd, conditionOn = condOn)$vals,
    mc.cores=nrCores)
  
  stopifnot(length(lim)==1) # FIX: order limits and then combine all? of them
#   cond <- sapply(lim,function(ll)sum(sapply(lapply(ll,"[[","ints"),size)>sizeInit))
#   
#   ocond <- order(cond,decreasing = TRUE)
#   paths <- paths[ocond]
  
  valsFin <- lapply(1:length(vList),function(i)
    combineIntervals(lapply(lim,function(l)as.list(l[[i]])), what="union"))
  
  return(valsFin)
  
}

fitFakePathToTestVec <- function(path,orgPath)
{
  
  
  testVecMax <- sapply(unique(orgPath),function(z) max(which(z==orgPath)))
  
  for(i in 1:length(testVecMax)){

    act <- unique(path[1:max(which(unique(orgPath)[i]==path))])
    targ <- unique(orgPath[1:testVecMax[i]])
    path <- c(path,targ[!targ%in%act],unique(orgPath)[i])
    
    }
  
  return(path)
  
}