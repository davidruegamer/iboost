# function to calculate the widest or narrowest interval

getFUL <- function(cc,wmin=which.min, wmax=which.max, 
                   pairwise = (wmin(1:2)==2), Sigma=NULL, returnIndex = FALSE,
                   combineIntervals = FALSE)
{
  
  
  
  Vlos <- sapply(cc,"[[","Vlo")
  Vups <- sapply(cc,"[[","Vup")
  vTys <- unique(sapply(cc,function(ccc) unique(ccc$vTy)))
  
  if(wmin(1:2)==2){ # if wmin==which.max
    
    includes <- Vlos <= vTys & vTys <= Vups
    
    if(sum(includes)!=length(Vlos)){
      
      VlosInclude <- Vlos[includes] 
      VlosNotInclude <- Vlos[!includes]
      VupsInclude <- Vups[includes]
      VupsNotInclude <- Vups[!includes]
      
      if(sum(includes)==0){
        
        Vlos <- Vups <- vTys
        
      }else{
        
        conditionOne <- VlosNotInclude<VupsNotInclude
        VlosNotInclude <- VlosNotInclude[conditionOne]
        VupsNotInclude <- VupsNotInclude[conditionOne]
        
        Vlos <- c(VlosInclude, VlosNotInclude[VlosNotInclude<max(VupsInclude)])
        Vups <- c(VupsInclude, VupsNotInclude[VupsNotInclude>min(VlosInclude)])
        
      }
      
      
    }
    
  }
  
  if(pairwise){
    
    wmin=which.min
    
    Rk <- sapply(1:length(Vlos),function(i)truncPnorm(a=Vlos[i],
                                                        b=Vups[i],
                                                        sd=Sigma,
                                                        mu=0,
                                                        x=vTys))
    wminX = wmaxX = wmin(2*pmin(Rk,1-Rk))
    
  }else{
  
    wminX <- wmin(Vups)
    wmaxX <- wmax(Vlos)
    
  }
  
  if(returnIndex) return(unique(c(wminX,wmaxX))) else
  
  return(list(Vlo=Vlos[wmaxX],
              Vup=Vups[wminX],
              vTy=vTys))
  
}