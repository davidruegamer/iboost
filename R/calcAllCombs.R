calcAllCombs <- function(mod, Ups, Sigma=NULL, estSd, nrCores = 10)
{
  
  # stuff for calculations
  dm <- extract(mod,"design")
  hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
  selCourse <- selected(mod)
  selCov <- lapply(1:mstop(mod),function(i)unique(selected(mod)[1:i]))
  Y <- mod$response
  n <- length(Y)
  p <- length(mod$basemodel)
  
  # define vector of iteration indices
  testVec <- sapply(unique(selCourse),function(z) max(which(z==selCourse)))
  testVecSort <- testVec[order(selCourse[testVec])]
  
  notSelected <- (1:p)[!(1:p)%in%selCourse]
  
  IminusHatsq <- lapply(1:p,function(i)crossprod(diag(n)-hatMats[[i]]))
  
  ss <- 1:p
  M <- NULL
  
  v = lapply(testVecSort,function(i) vOLScreate(iter = i,designmat = dm, selCourse = selCourse, n = n))
  
  
  M <- unlist(mclapply(1:max(testVecSort),function(i){
    
    selC <- selCourse[i]
    M1minusM2 <- lapply(notSelected,function(ii) IminusHatsq[[selC]]-IminusHatsq[[ii]])
    lapply(M1minusM2, function(m) crossprod(crossprod(m,Ups[[i]]),Ups[[i]]))
    
  }, mc.cores=nrCores), recursive = F)
  
  perStepConditions <- mclapply(1:max(testVecSort),function(i){
    
    selC <- selCourse[i]
    M1minusM2 <- lapply(ss[-selC],function(ii) IminusHatsq[[selC]]-IminusHatsq[[ii]]) 
    lapply(M1minusM2, function(m) crossprod(crossprod(m,Ups[[i]]),Ups[[i]]))
    
  }, mc.cores=nrCores)
  
  res <- mclapply(1:length(v),function(i){
    
    if(is.null(Sigma)) sd=sqrt(crossprod(v[[i]]))*estSd else sd=Sigma
    
    strictestForEachStep <- lapply(perStepConditions,function(setOfPminusOneCond)
      combineIntervals(lapply(setOfPminusOneCond,function(M)getIntersection(v=v[[i]],M=M,y=Y,Sigma=Sigma)),
                       what = "intersection"))
    unstrictestForEachCov <- lapply(testVecSort,function(cov){
      
      selCovI <- which(selCourse==selCourse[cov])
      combineIntervals(strictestForEachStep[selCovI], what = "union")
      
    })
    #print(paste0("==== ", i, " ===="))
#     print(paste0("lower: ", which.max(do.call("rbind",unstrictestForEachCov)[,1]),
#                  "; upper: ", which.min(do.call("rbind",unstrictestForEachCov)[,2])))
    strictestOverall <- combineIntervals(unstrictestForEachCov,
                                         what = "intersection")
    
    
    finList  <- if(length(notSelected)>0) list(strictestOverall,
              combineIntervals(lapply(M,function(m)
                getIntersection(v=v[[i]],M=m,y=Y,Sigma=Sigma)),
                what = "intersection")) else strictestOverall
    
    return(combineIntervals(finList, what = "intersection"))
    
  },mc.cores=nrCores)
  
  return(list(vals=res,testvec=v))
  
}