#'
#' @param obj mboost object
#' @param method character; if possible, choose custom method. See details.
#' @param var numeric; true or consistent estimator of variance
#' @param B numeric; number of informative samples.
#' @param alpha numeric; significance level for p-value / size of selective interval (\code{1 - alpha}).
#' @param ncore numeric; number of cores to use (via mclapply)
#' @param refit.mboost function; needed if obj was created by a direct call of \code{mboost_fit}.
#' @param eps factor to stabilize inverse of refit
#' @import parallel, Intervals
#' 
#' @examples 
#' 
#' library(mboost)
#' 
#' set.seed(0)
#' 
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n) + 0.25 * x1
#' x3 <- rnorm(n)
#' eta <- 3 * sin(x1) + x2^2
#' y <- eta + rnorm(n)
#' 
#' spline1 <- bbs(x1, knots = 20, df = 4)
#' knots.x2 <- quantile(x2, c(0.25, 0.5, 0.75))
#' spline2 <- bbs(x2, knots = knots.x2, df = 4)
#' spline3 <- bbs(x3, knots = 20, df = 4)
#' 
#' data <- data.frame(y=y, x1=x1, x2=x2, x3=x3)
#' 
#' mod1 <- mboost(y ~ spline1 + spline2 + spline3,
#' control=boost_control(mstop = 73), offset = 0, 
#' data = data)
#'
iboost <- function(obj, 
                   method = c("linesearch", 
                              "slice", 
                              "analytic",
                              "normalsamp",
                              "impsamp",
                              "unifsamp"),
                   vars,
                   varForSampling = vars[1],
                   alpha = 0.05,
                   ncore = 1,
                   B = 1000,
                   refit.mboost = NULL,
                   eps = 1e-12,
                   Ups = NULL,
                   checkBL = TRUE,
                   vT = NULL,
                   computeCI = TRUE,
                   fac = 2,
                   nInit = 20,
                   ...)
{
  
  # get method
  method <- match.arg(method)
  
  ####### checks #######

  # check model class
  if(all(class(obj) != "mboost")) stop("obj must be of class mboost.")
  
  # check that offset is (almost) zero
  if(abs(obj$offset) > 0.000001) stop(paste0("Please refit the model by substracting", 
                                  " the offset from the response before calling mboost / glmboost."))
    
  # check alpha level
  if(alpha >= 1 | alpha <= 0.001) stop("Please provide an alpha value in between 0.001 and 1.")

  # check that family == Gaussian
  if(obj$family@name != "Squared Error (Regression)")
    stop("Inference for families other than Gaussian currently not available.")
  
  # check for group baselearners and restrict method
  nrcol <- if(class(obj)[1] == "glmboost") rle(obj$assign)$lengths else 
    unlist(lapply(extract(obj, "design"), "ncol"))
  if(any(nrcol > 1) & method %in% c("analytic", 
                                    "normalsamp") & is.null(vT)) 
    stop("Method analytic not available for group base learner.") 

  # check for not supported base learners  
  learners <- gsub("(.*)\\(.*\\)","\\1",names(obj$baselearner))
  if(checkBL && any(sapply(learners, function(x) !x %in% c("bbs","bols"))))
    stop("Inference is currently restricted to model with linear and b-spline base learner only.")
  
  # check model call and data argument
  if(is.null(refit.mboost) && (is.null(obj$call) || class(obj$call$data)!="name") && method!="analytic")
    stop("Need data object as parameter in initial model call.")
  
  
  ####### prepare inference #######
  
  sel <- sort(unique(selected(obj)))
  
  # function definition for the check of congruent selection sets
  is_congruent <- function(mod) setequal(selected(mod), sel)
  
  y <- obj$response
  n <- length(y)  
  

  # testvector(s)
  if(is.null(vT)) vT <- getTestvector(obj, eps = eps)
  
  if(is.null(refit.mboost)){
    
    refit.mboost <- function(newY){
      
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

  var <- varForSampling
  
  ####### call distribution helper functions #######
  tstart <- Sys.time()
  res <- switch (method,
    analytic = polyh_inf(obj, vT, var, alpha, Ups = Ups),
    slice = mclapply(lapply(vT, t), function(v) getVloVup(mod = obj, v = v, Ups = Ups), 
                     mc.cores = ncore),
    linesearch = mclapply(vT, function(vt) 
      infsamp(refit.mboost, y, vT = vt, is_congruent, B, var, ...), 
                          mc.cores = ncore),
    normalsamp = mclapply(vT, function(vt) 
      normalsamp(refit.mboost, y, vT = vt, is_congruent, B, var, ...), 
                          mc.cores = ncore),
    impsamp = mclapply(vT, function(vt) 
      impsamp(refit.mboost, y, vT = vt, is_congruent, B, var, ...), 
                       mc.cores = ncore),
    unifsamp = mclapply(vT, function(vt) 
      unifsamp(refit.mboost, y, vT = vt, is_congruent, B, var, nInit = nInit, ...),
                        mc.cores = ncore)
  )
  dur <- as.numeric(Sys.time() - tstart)
  ####### format results #######
  
  ret <- vector("list", length(vars))
  
  for(i in 1:length(vars)){
    
    ret[[i]] <- format_iboost_res(res = res, alpha = alpha, 
                             method = method, vT = vT, 
                             this_var = vars[i],
                             computeCI = computeCI,
                             fac = fac, y = y)
    
    ret[[i]]$dur <- dur
    
  }

  
  return(ret)
  
  
}


format_iboost_res <- function(res, alpha, method, vT, this_var, computeCI, fac, y){
  
  resDF <- list()
  
  if(method=="unifsampTest"){
    
    method <- "unifsamp"
    computeCI <- FALSE
    
  }
  
  for(j in 1:length(res)){
    
    
    if(method %in% c("linesearch", "analytic")){
      
      res[[j]] <- sort(res[[j]])
      vlo <- res[[j]][1]
      vup <- res[[j]][2]
      
    }else if(method == "slice"){ # slice
      
      vlo <- res[[j]]$ul[[1]]
      vup <- res[[j]]$ul[[2]]
      
    }else if(method %in% c("normalsamp", "impsamp", "unifsamp", "impsamp")){
      
      if("extrcase" %in% names(attributes(res[[j]]))){
        resDF[[j]] <- data.frame(lower = NA, mean = res[[j]]$obsval, upper = NA, pval = 0,
                                 lowtrunc = NA, uptrunc = NA)
        next
      }
      
      if(sum(res[[j]]$logvals) == 0){
        
        warning("No sample survived. Most likely reason: truncated space very small.")
        resDF[[j]] <- data.frame(lower = NA, mean = res[[j]]$obsval, upper = NA, pval = NA,
                                 lowtrunc = NA, uptrunc = NA)
        next
        
      }
      
      vlo <- min(res[[j]]$rBs)
      vup <- max(res[[j]]$rBs)
      mu <- res[[j]]$obsval
      l1 <- res[[j]]$logvals
      
      survr <- res[[j]]$rBs[l1]
      l2 <- survr > mu
      survr_gr <- survr[l2]
      survr_le <- survr[!l2]
      
      pval <- sum(res[[j]]$weights(var = this_var)[l2]) / sum(res[[j]]$weights(var = this_var))
      
      vt <- vT[[j]]
      vTv <- as.numeric(tcrossprod(vT[[j]]))
      
      if(computeCI){
        
        # define binary search and get confidence intervaĺ
        if(method %in% c("normalsamp", "unifsamp"))
        {
          
          w <- res[[j]]$weights(var = this_var)
          if(all(w<1e-7)) w <- w / max(w)
          ftlo <- function(t) sum(exp(survr_gr * t / (this_var * vTv)) * w[l2]) / 
            sum(exp(survr * t / (this_var * vTv)) * w)
          ftup <- function(t) sum(exp(survr_le * t / (this_var * vTv)) * w[!l2]) / 
            sum(exp(survr * t / (this_var * vTv)) * w)
          
          if(length(survr_gr)==0){
            
            ci <- c(-Inf, Inf)
            
          }else{
            
            # check admissible values
            testvals <- seq(vlo - 4*sqrt(this_var * vTv), 
                            vup + 4*sqrt(this_var * vTv),
                            l = 1000)
            flovals <- sapply(testvals, ftlo)
            fupvals <- sapply(testvals, ftup)
            ll <- min(which(!is.na(flovals) & 
                              !is.nan(flovals)))
            lu <- max(which(!is.na(flovals) & 
                              !is.nan(flovals)))
            
            low <- try(uniroot(f = function(x) ftlo(x) - alpha/2, 
                               interval = testvals[c(ll,lu)],
                               extendInt = "upX")$root)
            
            up <- try(uniroot(f = function(x) ftup(x) - alpha/2, 
                              interval = testvals[c(ll,lu)],
                              extendInt = "downX")$root)
            
            if(class(low)=="try-error") low <- -Inf
            if(class(up)=="try-error") up <- Inf
            
            ci <- c(low, up)
            
          }
          
        }else{
          
          # impsamp
          w <- res[[j]]$weights(var = this_var)
          if(all(w<1e-7)) w <- w / max(w)
          ft <- function(t) sum(exp(survr_gr * t / this_var) * w[l2]) / 
            sum(exp(survr * t / this_var) * w)
          ibins <- function(d, fac){
            uniroot(f = function(x) ft(x) - d, 
                    interval = c(0, fac * mu),
                    extendInt = "upX")
          }
          
          low <- try(ibins(alpha, fac = fac)$root)
          if(class(low)=="try-error") low <- 0
          ci <- c(max(low, 0), Inf)
          
        }
        
      }else{
        
        ci <- c(NA, NA)
        
      }  
      
      # return two-sided p-value
      if(method != "imsamp") pval <- 2*min(1-pval,pval)
        
      resDF[[j]] <- data.frame(lower = ci[1], mean = mu, upper = ci[2], pval = pval,
                               lowtrunc = vlo, uptrunc = vup)
      
    }
    
    if(method %in% c("slice", "analytic", "linesearch")){
      
      df <- nrow(vT[[j]])
      lin <- df==1
      
      if(lin){
        
        vTv <- as.numeric(tcrossprod(vT[[j]]))
        mu <- as.numeric(vT[[j]] %*% y) 
        sqvTv <- sqrt(vTv)
        sd <- sqrt(this_var)*sqvTv
        gridrange = c(-fac,fac)*50
        
        pv <- selectiveInference:::tnorm.surv(mu, mean = 0, sd = sd, a = vlo, b = vup)
        if(mu<=0) pv <- 1 - pv
        
        xg = seq(gridrange[1]*sd, gridrange[2]*sd, length = 100)
        fun = function(x) { selectiveInference:::tnorm.surv(z = mu, mean = x, sd = sd, 
                                                            a = vlo, b = vup, bits = NULL) }
        
        if(computeCI) 
          ci = selectiveInference:::grid.search(xg, fun, alpha/2, 
                                                1-alpha/2, 100, 2) else
                                                  ci = c(NA, NA)
        
        
      }else{
        
        R <- vT[[j]] %*% y
        mu <- sqrt(sum(R^2))
        pv <- selectiveInference:::TC_surv(TC = mu, sigma = sqrt(this_var), 
                                           df = df, E = Intervals(c(vlo,vup)))
        ci <- c(NA, NA)
        
        
      }
      
      resDF[[j]] <- data.frame(lower = ci[1], mean = mu, upper = ci[2], pval = 2*min(pv, 1-pv),
                               lowtrunc = vlo, uptrunc = vup)   
      
    }

  }
  
  resDF <- do.call("rbind", resDF)
  if(is.null(vT)) rownames(resDF) <- sel
  
  ret <- list(dist = res,
              method = method, 
              alpha = alpha,
              vT = vT,
              yorg = y,
              resDF = resDF,
              var = this_var)
  
  class(ret) <- "iboost"
  
  
  return(ret)
  
}

