#'
#' @param mod an object of class \code{mboost}.
#' @param otherSelCourse either \code{NULL} or \code{numeric}. If \code{NULL} Upsilon is computed for
#' \code{selected(mod)}, otherwise for otherSelCourse
#' @param last logical; whether to return only Upsilon for the last selection step. If FALSE, a list of
#' Upsilons for all iterations is returned
#' @param listFun function; used to combine design matrices, which are build with 
#' baselearner operators such as \code{%X%}.
#' @param returnUpsOnly logical; if FALSE, the hat matrices are also returned.
#'

getUpsilons <- function(mod, otherSelCourse = NULL, last = FALSE, 
                        listFun = kronecker, returnUpsOnly = TRUE)
{
  
  # computes Upsilons for all mstop(mod)-1 steps. 
  # The last Upsilon is not usually neccessary, but
  # can be calculated via upsPart[[selCourse[mstop(mod)]]]%*%Upsilons[[mstop(mod)]]%*%Y
  
  library(Matrix)
  # X <- as.matrix(do.call("cbind",lapply(mod$baselearner,function(x)x$get_data())))
  w <- mod$`(weights)`
  n <- length(w)
  # Y <- mod$response
  nu <- mod$control$nu
  
  selCourse <- if(is.null(otherSelCourse)) selected(mod) else otherSelCourse
  if(inherits(mod,"glmboost")) selCourse <- new_order(selCourse)
  
  ### part for hat matrix
  hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
  if(any(sapply(hatMats,is.null)) | !returnUpsOnly | (any(w!=1)) | inherits(mod, "glmboost")){
    
    if(any(w>1)) stop("Calculation for weights > 1 not implemented yet.")
    
    if(inherits(mod,"glmboost")){
      
      desMats <- getDesignmat(mod, split = T)
      
    }else{

      desMats <- lapply(1:length(mod$baselearner), function(i) 
        extract(mod, what="design", which=i)[[1]][as.logical(w),])
      desMats[sapply(desMats,class)=="list"] <- lapply(desMats[sapply(desMats,class)=="list"],
                                                       function(l)listFun(l[[1]],l[[2]]))
    }
    
    if(inherits(mod,"glmboost")){
      
      hatMatsLeft <- lapply(desMats,function(x)solve(t(x)%*%x)%*%t(x))
      
    }else{
      
      lambdas <- sapply(1:length(mod$baselearner), function(i) extract(mod, what="lambda", which=i)[[1]])
      
      if(any(lambdas!=0) & any(w!=1)) stop("Can't deal with given weights and penalized baselearners.")
      
      pens <- lapply(1:length(mod$baselearner), function(i) extract(mod, what="penalty", which=i)[[1]])
      hatMatsLeft <- lapply(1:length(pens),function(i)
        solve(Matrix::t(desMats[[i]])%*%desMats[[i]] + lambdas[i]*pens[[i]])%*%Matrix::t(desMats[[i]])
      )
      
    }
    
    hatMats <- mapply(function(x,h)x%*%h,x=desMats,h=hatMatsLeft,SIMPLIFY=FALSE)
    
  }
  
  
  upsPart <- lapply(hatMats,function(x)(diag(w)-x*nu))
  
  len <- length(selCourse)
  Upsilons <- vector("list",len)
  Upsilons[[1]] <- diag(w)
  
  if(len>1)
    for(i in 1:(len-1)) 
      Upsilons[[i+1]] <- upsPart[[selCourse[i]]] %*% Upsilons[[i]] 
  
  if(last) return(Upsilons[[length(Upsilons)]])
  
  if(returnUpsOnly) return(Upsilons) else return(list(Upsilons=Upsilons, hatMats=hatMats, hatMatsLeft=hatMatsLeft))
  
}

