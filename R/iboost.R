#'
#' @param obj mboost object
#' @param method character; if possible, choose custom method. See details.
#' @param var numeric; true or consistent estimator of variance
#' @param bootType character; type of bootstrap to employ if \code{method = 'boot'}
#' @param B numeric; number of bootstrap or informative samples.
#' @param alpha numeric; significance level for p-value / size of selective interval (\code{1 - alpha}).
#' @param ncore numeric; number of cores to use (via mclapply)
#' @param refit.mboost function; needed if obj was created by a direct call of \code{mboost_fit}.
#' 
#'
iboost <- function(obj, 
                   method = c("infsamp", "slice", "analytic", "boot"),
                   var,
                   bootType = c("param", "nonparam"),
                   B = 1000,
                   alpha = 0.05,
                   ncore = 1,
                   refit.mboost = NULL)
{
  
  ####### checks #######

  if(all(class(obj) != "mboost")) stop("obj must be of class mboost.")
  
  # condition <- match.arg(condition)
  method <- match.arg(method)
  if(alpha >= 1 | alpha <= 0.001) stop("Please provide an alpha value in between 0.001 and 1.")

  if(obj$family@name != "Squared Error (Regression)")
    stop("Inference for families other than Gaussian currently not available.")
  
  # check for not supported baselearners
  if(class(obj)[1] != "glmboost"){
    
    bls <- gsub("\\(.*\\)", "", names(obj$baselearner))
    if(any(bls != "bols")) stop("Inference currently restricted to models with linear baselearners only.")
    
  }else{
    
    bls <- (unique(obj$assign))
    
  }
  
  # check for group baselearners and restrict method
  nrcol <- if(class(obj)[1] == "glmboost") rle(obj$assign)$lengths else 
    unlist(lapply(extract(obj, "design"), "ncol"))
  if(any(nrcol > 1) & method == "analytic") stop("Method analytic not available for group baselearners.") 
  
  ####### prepare inference #######
  
  sel <- sort(unique(selected(obj)))
  
  # function definition for the check of congruent selection sets
  if(method %in% c("infsamp", "boot")) 
    is_congruent <- function(mod) setequal(selected(mod), sel)
  
  y <- obj$response
    
  # get design matrix and testvector(s)
  X <- getDesignmat(obj)
  Xplus <- solve(crossprod(X)) %*% t(X)
  vT <- lapply(1:length(sel), function(i) t(Xplus[i, ]))
  
  if(method == "boot"){
    
    # get testvector function
    extract_vT <- function(mod, which){
      
      # TODO - CHECK: is order correct? 
      X <- getDesignmat(mod)[, order(sel)]
      W <- mod$`(weights)`
      ( solve(crossprod(X * W, X)) %*% t(X) %*% diag(W) )[which, ]
      
    } 
  
  }
  
  if(is.null(refit.mboost)){
    
    refit.mboost <<- function(newY){
      
      call <- obj$call
      olddat <- eval(call$data, parent.frame())
      yname <- all.vars(as.formula(obj$call))[1]
      olddat[[yname]] <- newY
      assign(as.character(call$data), olddat)
      if(as.character(as.list(call)[[1]]) == "glmboost.formula")
        stop("Please provide a refit.mboost function.")
      eval(call)
      
    }
    
  }


  ####### call distribution helper functions #######
  
  res <- switch (method,
    boot = boot_inf(obj, evT = extract_vT, nrBL = length(sel), 
                    is_congruent, B, bootType, var, ncore),
    analytic = polyh_inf(obj, vT, is_congruent, B, var, ncore, alpha),
    slice = slice_inf(obj, vT, is_congruent, B, var, ncore),
    infsamp = infsamp(refit.mboost, y, vT, is_congruent, B, var, ncore)
  )
  
  ####### format results #######
  
  resDF <- list()
  
  for(j in 1:length(res)){
    
    vTv <- as.numeric(tcrossprod(vT[[j]]))
    mu <- as.numeric(vT[[j]] %*% y) 
    sqvTv <- sqrt(vTv)
    
    if(method == "infsamp"){
    
      ci <- qtnorm(c(alpha/2, 1-alpha/2), 
                   mean = mu,
                   sd = sqrt(var)*sqvTv,
                   lower = res[[j]][1],
                   upper = res[[j]][2])
      
    }else if(method == "boot"){
      
      ci <- quantile(res[[j]] * vTv, probs = c(alpha/2, 1-alpha/2))
      
    }else if(method == "analytic"){
      
      ci <- res[[j]][[2]]$int
      
    }else{ # slice
      
      
      
    }
    
    resDF[[j]] <- data.frame(lower = ci[1], mean = mu, upper = ci[2])
    
  }
  
  resDF <- do.call("rbind", resDF)
  
  ret <- list(dist = res,
              method = method,
              alpha = alpha,
              vT = vT,
              yorg = y,
              resDF = resDF)
  
  class(ret) <- "iboost"
  
  return(ret)
  
  
}



print.iboost <- function(ibo)
{
  
  if(method == "boot"){
    
    
  }else{
    
    
  }
  
  
}


