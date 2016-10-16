getDesignmat <- function(mod, split = FALSE)
{
  
  if(class(mod)[1] == "glmboost"){
    
    ret <- extract(mod,"design")
    asel <- new_order(mod$assign[mod$assign %in% unique(selected(mod))])
    if(split) ret <- lapply(unique(asel), function(i) ret[,which(asel == i)])
    
  }else{
    
    if(NOFACTORINMODEL){
      
      ret <- do.call("cbind", extract(mod, "design"))
      
    }else{
      
      
      
    }
  }
  
  return(ret)
  
}

new_order <- function(x) as.numeric(as.factor(x))
