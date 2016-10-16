getVloVup <- function(mod, Ups=NULL, v, Sigma=NULL,
                      conditionOn = c("last","notSelected","first","inBetween","uniqSel"),
                      crit=c("none","GCV","AIC","gMDL"), returnComps=FALSE, 
                      calcM = TRUE, Mextra=NULL, vert=FALSE,
                      ncore, ...)
{
  
  # function that provides Vlo and Vup
  # handles all M matrices to get connected to v in the right way
  
  # compute stuff for later
  dm <- extract(mod, "design")
  selCourse <- selected(mod)
  n <- length(mod$ustart)
  y <- mod$response
  
  # compute Upsilons if neccessary
  if(is.null(Ups) & calcM) Ups <- getUpsilons(mod)
  
  # get M matrices
  
  if(calcM){
    
    MC <- mCreate(mod, conditionOn = conditionOn, Ups = Ups, Sigma = Sigma, 
                  crit = crit[1], ...)
    M <- MC$M
    Madd <- MC$Madd
    ind <- MC$ind
    
  }else{
    
    if(is.null(Mextra) | !all(c("ind","v","M")%in%names(Mextra)))
      stop("If calcM is FALSE Mextra has to be supplied.")
    M <- Mextra$M
    Madd <- Mextra$Madd
    ind <- Mextra$ind
    
  }
  
  if(returnComps) return(list(ind=ind,v=v,M=M,Madd=Madd))
  
  # define what to do
  
  forThisSeq <- 1:length(v) # for factors no list of lists????
  
  if(!is.null(Madd)) M <- lapply(1:length(M), function(i) append(M[[i]], Madd)) 

  ul <- mclapply(forThisSeq,function(i)
    calcUpLow(y = y, v = v[[i]], M = M, Sigma = Sigma, vert = vert), mc.cores = ncore)
    
  
  return(list(ul=ul,ind=ind,v=v))
  
}

addZeros <- function(mat, dim)
{
  
  orgD <- nrow(mat)
  diff <- dim-orgD
  #nm <- 
  cbind(rbind(mat,matrix(rep(0,diff*orgD),ncol=orgD)),matrix(rep(0,diff*dim),nrow=dim))
  # diag(nm)[(orgD+1):dim] <- 1
  #return(nm)
  
}