combineIntervals <- function(objFromSlice,
                 what="union",
                 onlyOneClosedInterval=FALSE,
                 ...)
{

  if(is.list(objFromSlice)){ if(is.list(objFromSlice[[1]])) 
    objFromSlice <- list(ints = lapply(objFromSlice,"[[","ints"), vTy = sapply(objFromSlice,"[[","vTy")[1])
  }
  
  vTy <- NULL
  
  if(!is.null(objFromSlice$vTy)) vTy <- objFromSlice$vTy else{
    
    if(onlyOneClosedInterval) stop("You have to supply vTy, if intervals should be checked.")
    
  }
  
  ints <- objFromSlice$ints
  
  if(!is.list(ints)) return(list(vTy=vTy,ints=ints))
  
  resInts <- if(what=="union") do.call("interval_union",ints) else do.call("interval_intersection",ints)
  
  if(onlyOneClosedInterval){
    
    incl <- interval_included(resInts,Intervals(c(vTy,vTy)))
    if(is.list(incl) & length(incl[[1]])==0) stop("No Interval includes vTy") else 
      return(list(vTy = vTy, ints = resInts[incl]))  
    
  }else 
    return(list(vTy = vTy, ints = resInts))
  
}

combineSliceObjects <- function(obj1, obj2, what="intersection", ...)
{
  
  stopifnot(is.list(obj1) & is.list(obj2) & 
              !is.null(obj1[[1]]$ints) & !is.null(obj2[[1]]$ints) & 
              length(obj1)==length(obj2))
  
  resObj <- vector("list",length(obj1))
  
  for(i in 1:length(obj1)){
    
    resObj[[i]] <- combineIntervals(list(obj1[[i]],obj2[[i]]), what=what, ...)
    
  }
  
  return(resObj)
  
  
}