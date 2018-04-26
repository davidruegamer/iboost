getVloVup <- function(mod, Ups=NULL, v, Sigma=NULL,
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
  MC <- mCreate(mod, Ups = Ups, crit = crit, ...)
  M <- MC$M
  Madd <- MC$Madd
  ind <- MC$ind
  
  if(!is.null(Madd)) M <- lapply(1:length(M), function(i) append(M[[i]], Madd)) 

  ul <- calcUpLow(y = y, v = v, M = M, Sigma = Sigma)
    
  return(list(ul = ul, ind = ind, v = v))
  
}