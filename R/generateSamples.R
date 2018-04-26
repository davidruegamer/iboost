generateSamples <- function(
  refitFun,
  is_congruent,
  samplingFun,
  refPoint,
  dir,
  B = 1000,
  U = NULL
)
{

  # generate 'test statistics'
  rBs <- samplingFun(B)
  # generate responses
  if(is.null(U)) yBs <- sapply(rBs, function(rb) refPoint + dir*rb) else
    yBs <- sapply(rBs, function(rb) U%*%(refPoint + dir*rb))
    
  # check whether response yields congruent selection
  logvals <- apply(yBs, 2, function(yb){
    
    m <- try(refitFun(yb))
    
    if(class(m)=="try-error") 
      return(NA) else 
      return(is_congruent(m))
    
  })
  
  # get incidences, where model could not be fitted
  fitna <- is.na(logvals)
  logvals[fitna] <- FALSE
  # survivals
  survive <- rBs[logvals]
  if(length(survive)==0){
    
    search = FALSE
    warning("No sample survived.")
    
  }
  
  return(list(rBs = rBs, 
              logvals = logvals, 
              fitna = fitna, 
              survive = survive))
  
}