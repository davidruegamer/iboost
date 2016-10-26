getVloVup <- function(mod, Ups=NULL, v, Sigma=NULL,
                      # conditionOn = c("last", "notSelected", "first",
                      #                 "inBetween", "uniqSel"),
                      crit=c("none", "GCV", "AIC", "gMDL"), 
                      ncore, ...)
{
  
  # function that provides Vlo and Vup
  # checks that all M matrices are connected to v in the right way
  
  # get components
  y <- mod$response
  crit <- match.arg(crit)
  
  # compute Upsilons if neccessary
  if(is.null(Ups)) Ups <- getUpsilons(mod)
  
  # get M matrices
  MC <- mCreate(mod, #conditionOn = conditionOn, 
                Ups = Ups, 
                crit = crit, ...)
  M <- MC$M
  Madd <- MC$Madd
  ind <- MC$ind
  
  # define what to do
  
  forThisSeq <- 1:length(v) 
  
  if(!is.null(Madd)) M <- lapply(1:length(M), function(i) append(M[[i]], Madd)) 

  ul <- mclapply(forThisSeq, function(i) calcUpLow(y = y, v = v[[i]], M = M, Sigma = Sigma), 
    mc.cores = ncore)
    
  retList <- lapply(1:length(ul), function(i) list(ul = ul[[i]], ind = ind[i], v = v[[i]]))
  return(retList)
  
}