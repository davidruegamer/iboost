# extract the design matrix (matrices if split = TRUE) from mod
getDesignmat <- function(mod, split = FALSE, full = FALSE)
{
  
  if(full){
    
    ret <- do.call("cbind", lapply(mod$baselearner, extract, "design"))
    
  }else{
  
    if(class(mod)[1] == "glmboost"){
      
      # get design matrix
      ret <- extract(mod,"design")
      if(!1 %in% selected(mod) & colnames(ret)[1] == "(Intercept)")
        ret <- ret[,-1]
      
      if(split){ # recover single design matrices
      
        if(is.null(dim(ret))){
        
          ret <- list(ret)
          
        }else{
          
          ret <- lapply(1:ncol(ret), function(i) ret[,i])
        
        }
        
      }
      
    }else{
      
      # get single design matrices
      ret <- extract(mod, "design")
      
      if(!all(sapply(ret,ncol)==1)){
        
        # factor bols baselearner must have at least 2 columns
        # therefore check # rows
        lenCheck <- sapply(ret, nrow) != length(mod$response)
        
        if(any(lenCheck)){
          
          # recover full matrix
          for(j in (1:length(ret))[lenCheck]){
            
            ind <- diff(c(as.numeric(rownames(ret[[j]])), 
                          length(mod$response) + 1))
            ret[[j]] <- ret[[j]][rep(1:nrow(ret[[j]]), ind),] 
            
          }
        }
      }
      
      if(!split) ret <- do.call("cbind", ret)
          
    }
    
  }
  
  return(ret)
  
}

# function to order values maintaining duplicates
# new_order <- function(x){ 
#   
#   minx <- min(x)
#   as.numeric(as.factor(x)) + minx - 1
#     
# }

# make Baselearners using character names
# blFun specifies the type of BL
makeBL <- function(charName, blFun, data, ...) {
  
  temp <- as.data.frame(data[, charName])
  colnames(temp) <- charName
  bl <- tryCatch(blFun(temp, ...), error=function(e) bols(temp))
  bl$set_names(charName)
  bl
  
}