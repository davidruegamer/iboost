#' @title Valid Inference for Model-based Gradient Boosting Models
#' @param obj mboost object
#' @param method character; if possible, choose custom method. See details.
#' @param vars numeric vector; a single numeric value or vector of numeric values for the variance 
#' used in the linear model (preferably the true variance or an estimation from a consistent estimator). 
#' If NULL, the empirical response variance is used, which will result in rather conservative inference.
#' @param varForSampling variance used for generate new samples. Defaults to the first entry of \code{vars} 
#' if not given.
#' @param B numeric; number of samples drawn for inference.
#' @param alpha numeric; significance level for p-value / size of selective interval (\code{1 - alpha}).
#' @param ncore numeric; number of cores to use (via \code{mclapply})
#' @param refit.mboost function; this is needed if \code{obj} was created by a direct call to \code{mboost_fit}. 
#' In this case, \code{refit.mboost} should be a function of the response, refitting the exact model 
#' as given by \code{obj} with the response vector.
#' @param Ups list of residual matrix produces by \code{\link{getUpsilons}}.
#' @param checkBL logical; if \code{TRUE} checks whether base-learner only include linear, group and 
#' spline base-learner (which are currently supported)
#' @param vT list of test vectors as produced by \code{\link{getTestvector}}.
#' @param computeCI logical; whether or not to compute selective confidence intervals
#' @param returnSamples logical; whether or not (default = FALSE) to only return 
#' the samples produced. Per default, p-values (and intervals) are calculated using the 
#' samples using \code{\link{format_iboost_res}}.
#' @param which numeric; selects only certain base-learner, for which inference is conducted.
#' @param ... Further arguments passed to the sampling method.
#' 
#' @import parallel mboost intervals selectiveInference
#' @importFrom intervals Intervals
#' @importFrom methods as
#' @importFrom stats as.formula dnorm pnorm qnorm rnorm runif uniroot
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom msm ptnorm
#' 
#' @export
#' @details iboost provides inference for L_2-Boosting models fitted with \code{mboost} 
#' with linear, group or spline base-learner based on 
#' Ruegamer and Greven (2018) when \code{method = unifsamp}, 
#' Yang et al. (2016) when \code{method = impsamp}, 
#' Tibshirani et al. (2016) when \code{method = analytic},
#' Loftus and Taylor (2015) when \code{method = slice} and
#' two variations of the \code{unifsamp} approach when \code{method} is 
#' \code{normalsamp} or \code{normaladjsamp}.
#' Only the methods \code{impsamp} and \code{slice} can be used for testing 
#' group effects or whole spline functions. 
#' 
#' @references 
#' Ruegamer, D. and Greven, S. (2018), 
#' Valid Inferece for L2-Boosting, 
#' arXiv e-prints arXiv:1805.01852.
#' 
#' Yang, F., Barber, R. F., Jain, P. and Lafferty, J. (2016), 
#' Selective inference for group-sparse linear models, 
#' Advances in Neural Information Processing Systems, pp. 2469-2477. 
#' 
#' Tibshirani, R. J., Taylor, J., Lockhart, R. & Tibshirani, R. (2016), 
#' Exact post-selection inference for sequential regression procedures, 
#' Journal of the American Statistical Association 111(514), 600-620.
#' 
#' Loftus, J. R. & Taylor, J. E. (2015), 
#' Selective inference in regression models with groupsof variables, 
#' arXiv e-prints arXiv:1511.01478.
#' 
#' @return Returns an object of class \code{iboost} or, if \code{length(vars)>1}, a 
#' list of \code{iboost} objects for each variance.
#' An \code{iboost} object is a list containing the following items
#' \itemize{
#'   \item \code{dist}: a list obtained by the sampling procedure including \code{rB}, the 
#'   sampled values, \code{logvals}, logical values whether the corresponding \code{rB} yields 
#'   to a congruent model with the initial model fit, \code{obsval}, the actual observed value 
#'   in the initial model fit and corresponding \code{weights} of the importance sampling 
#'   procedure.
#'   \item \code{method}: name of the method used
#'   \item \code{alpha}: alpha level used for the confidence interval limits
#'   \item \code{vT}: the test vector(s)
#'   \item \code{yorg}: original response value
#'   \item \code{resDF}: a data.frame consisting of the \code{lower} and 
#'   \code{upper} confidence interval limits, the observed value \code{mean}, 
#'   the calculated p-value \code{pval} and the truncation limits of the effect 
#'   \code{lowtrunc} and \code{uptrunc}.
#'    \item \code{var}: the variance used for inference calculation
#'    \item \code{dur}: total duration of sampling in seconds
#'      
#' }
#' @description Function computes selective p-values (and confidence intervals) 
#' for \code{\link[mboost]{mboost}} objects. Currently \code{iboost} supports Gaussian family models 
#' (L2-Boosting) with linear, group and spline base-learners.
#' @examples 
#' 
#' if(require("mboost")){
#' 
#' set.seed(0)
#' 
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n) + 0.25 * x1
#' x3 <- rnorm(n)
#' eta <- 3 * sin(x1) + x2^2
#' y <- scale(eta + rnorm(n), scale = FALSE)
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
#' # calculate p-values and intervals for model with 
#' # fixed stopping iteration:
#' # this is done with only B = 100 samples for
#' # demonstrative purposes and should be increased 
#' # for actual research questions
#' res <- iboost(mod1, method = "impsamp", B = 100)
#' 
#' # do the same with crossvalidation
#' \dontrun{
#' 
#' fixFolds <- cv(weights = model.weights(mod1),
#' type = "kfold", B = 10)
#' cvr <- cvrisk(mod1, folds = fixFolds, papply = lapply)
#' modf <- mod1[mstop(cvr)]
#' 
#' # define corresponding refit function
#' modFun <- function(y){
#' 
#'  mod <- mboost_fit(response = y,                
#'                    blg = blList,
#'                    offset = 0, 
#'                    control = boost_control(mstop = 73))
#'  cvr <- cvrisk(mod, folds = fixFolds, papply = lapply)
#'  return(mod[mstop(cvr)])
#'  }
#'  
#'  # this will take a while
#' (res <- iboost(modf, refit.mboost = modFun, method = "impsamp", B = 1000))
#' 
#' }
#' }
#' 
iboost <- function(obj, 
                   method = c("unifsamp",
                              "impsamp",
                              "analytic",
                              "slice", 
                              "linesearch", 
                              "normalsamp",
                              "normaladjsamp"
                              ),
                   vars = NULL,
                   varForSampling = NULL,
                   B = 1000,
                   alpha = 0.05,
                   ncore = 1,
                   refit.mboost = NULL,
                   Ups = NULL,
                   checkBL = TRUE,
                   vT = NULL,
                   computeCI = TRUE,
                   returnSamples = FALSE,
                   which = NULL,
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
    stop("Inference for families other than Gaussian is currently not available.")
  
  # check for group base-learners and restrict method
  nrcol <- if(class(obj)[1] == "glmboost") rle(obj$assign)$lengths else 
    unlist(lapply(extract(obj, "design"), "ncol"))
  if(any(nrcol > 1) & method %in% c("analytic", 
                                    "normalsamp") & is.null(vT)) 
    stop("Method analytic not available for group base-learner.") 

  # check for not supported base-learners  
  learners <- gsub("(.*)\\(.*\\)","\\1",names(obj$baselearner))
  if(checkBL && any(sapply(learners, function(x) !x %in% c("bbs","bols"))))
    stop("Inference is currently restricted to model with linear and b-spline base-learner only.")
  
  # check model call and data argument
  if(is.null(refit.mboost) && (is.null(obj$call) || class(obj$call$data)!="name") && method!="analytic")
    stop("Need data object as parameter in initial model call.")
  
  # check variances
  if(is.null(vars)) vars <- var(obj$response) * (length(obj$response)-1)/length(obj$response)
  if(is.null(varForSampling)) varForSampling <- vars[1]
  
  ####### prepare inference #######
  
  sel <- sort(unique(selected(obj)))
  
  # function definition for the check of congruent selection sets
  is_congruent <- function(mod) setequal(selected(mod), sel)
  
  y <- obj$response
  n <- length(y)  
  

  # testvector(s)
  if(is.null(vT)) vT <- getTestvector(obj)
  if(!is.null(which)) vT <- vT[which]
  
  # check if method complies with vT
  if(any(sapply(vT, NROW)>1) & !method%in%c("impsamp", "slice"))
  {
    
    method <- "impsamp"
    warning(paste0("The supplied method does not allow for inference of grouped variable parameters.\n", 
                "method has changed to impsamp.\n", 
                "If linear effects are given and you wish",
                "to use a different method, please use the which argument to only select", 
                "effects of linear base-learners.")
    )
    
  }
  
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
    # linesearch = mclapply(vT, function(vt) 
    #   infsamp(refit.mboost, y, vT = vt, is_congruent, B, var, ...), 
    #                       mc.cores = ncore),
    normalsamp = mclapply(vT, function(vt) 
      normalsamp(refit.mboost, y, vT = vt, is_congruent, B, var, ...), 
                          mc.cores = ncore),
    normaladjsamp = mclapply(vT, function(vt) 
      normaladjsamp(refit.mboost, y, vT = vt, is_congruent, B, var), 
      mc.cores = ncore),
    impsamp = mclapply(vT, function(vt) 
      impsamp(refit.mboost, y, vT = vt, is_congruent, B, var, ...), 
                       mc.cores = ncore),
    unifsamp = mclapply(vT, function(vt) 
      unifsamp(refit.mboost, y, vT = vt, is_congruent, B, var, ...),
                        mc.cores = ncore)
  )
  dur <- as.numeric(difftime(Sys.time(), tstart))
  ####### format results #######
  
  if(returnSamples) return(list(res = res, 
                                alpha = alpha, 
                                method = method, 
                                vT = vT, 
                                this_var = vars[i],
                                computeCI = computeCI,
                                y = y))
  
  ret <- vector("list", length(vars))
  
  for(i in 1:length(vars)){
    
    ret[[i]] <- format_iboost_res(res = res, 
                                  alpha = alpha, 
                                  method = method, 
                                  vT = vT, 
                                  this_var = vars[i],
                                  computeCI = computeCI,
                                  y = y)
    
    ret[[i]]$dur <- dur
    
  }
  
  if(length(ret)==1) ret <- ret[[1]]

  return(ret)
  
  
}

#' @title Function to calculate selective inference based on sampled values
#' 
#' @description This function calculates p-values and confidence intervals 
#' based on sampled values and weights obtained by the \code{\link{iboost}} function.
#' 
#' @details The specific use of this function is to (internally) compute inference based 
#' on samples in the \code{\link{iboost}} function, to recalculate inference for 
#' given samples \code{dist} from an \code{iboost} object or to calculate inference in 
#' the first place when \code{\link{iboost}} has been ran with argument \code{returnSamples = TRUE}. 
#' In the first and third case, use \code{format_iboost_res}, otherwise \code{format_iboost}.
#' 
#' @examples 
#' if(require("mboost")){
#' 
#' set.seed(0)
#' 
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n) + 0.25 * x1
#' x3 <- rnorm(n)
#' eta <- 3 * sin(x1) + x2^2
#' y <- scale(eta + rnorm(n), scale = FALSE)
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
#' # calculate p-values and intervals for model with 
#' # fixed stopping iteration:
#' # this is done with only B = 100 samples for
#' # demonstrative purposes and should be increased 
#' # for actual research questions
#' res <- iboost(mod1, method = "impsamp", B = 100)
#' 
#' # recalculate inference for different variance or alpha level
#' format_iboost(res, alpha = 0.1, this_var = var(y)*5)
#' }
#' @param res list of sampled values either directly obtained in the \code{\link{iboost}} 
#' call or when using \code{\link{iboost}} with \code{returnSamples = TRUE}.
#' @param alpha numeric; level for confidence interval(s)
#' @param method character; method used for sampling
#' @param vT list of test vector(s)
#' @param this_var numeric; variance to be used for inference
#' @param computeCI logical; whether or not to compute confidence intervals
#' @param fac numeric; used in inverse search for confidence interval limits
#' @param y vector; response vector
#' @export
#' @rdname format_iboost
#' @aliases format_iboost
#' 
format_iboost_res <- function(
  res, 
  alpha, 
  method, 
  vT, 
  this_var, 
  computeCI, 
  fac = 2, 
  y
  ){
  
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
      
    }else if(method %in% c("normalsamp", "normaladjsamp", "impsamp", "unifsamp", "impsamp")){
      
      # if("extrcase" %in% names(attributes(res[[j]]))){
      #   resDF[[j]] <- data.frame(lower = NA, mean = res[[j]]$obsval, upper = NA, pval = 0,
      #                            lowtrunc = NA, uptrunc = NA)
      #   next
      # }
      
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
        if(method %in% c("normalsamp", "unifsamp", "normaladjsamp"))
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
      if(method != "impsamp") pval <- 2*min(1-pval,pval)
        
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
  if(is.null(vT)) rownames(resDF) <- names(res[[1]]$dist)
  
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

#' @export
#' @param obj an iboost object
#' @rdname format_iboost_res   
#' 
format_iboost <- function(
  obj, 
  alpha = 0.05,
  method = obj$method, 
  vT = obj$vT,
  this_var = obj$var,
  computeCI = TRUE,
  fac = 2,
  y = obj$yorg)
{
  
  format_iboost_res(obj$dist, alpha, method, vT, this_var, computeCI, fac, y)
  
}