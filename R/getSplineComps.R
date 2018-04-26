getSplineComps <- function(mod, mstop = mstop(mod), sigma, getLambda = c("refit_simple", "refit_mgcv"))
{

  getLambda <- match.arg(getLambda)
  if(getLambda == "refit_simple"){

    lambdas <- calcLambda(mod = mod, mstop = mstop)
    X <- extract(mod, "design")

    s1i <- 1:length(X)

    indsbbs <- grepl("bbs", names(X))
    indseff <- rep(1:length(X), sapply(X, ncol))
    Xbig <- do.call("cbind", X)
    D <- extract(mod, "penalty")
    Dbig <- bdiag(lapply(s1i, function(i) lambdas[[1]][i] * D[[i]]))
    V <- crossprod(Xbig) + Dbig 

    VbetajInv <- try(solve(as.matrix(V)) / sigma)
    
    if(class(VbetajInv)=="try-error"){
      
      VbetajInv <- solve(as.matrix(V + diag(ncol(Dbig))*options("mboost_eps")[[1]])) / sigma
      
    }

    vTs <- lapply(s1i[indsbbs], function(w) lapply(1:nrow(Xbig), function(x) Xbig[x, which(indseff==w), drop=F] %*%
                                                     (VbetajInv%*%t(Xbig))[which(indseff==w),]))

  }else{
    #
    modgam <- refitmboostwithgam(mod)
    lambdas <- list(modgam$sp)
    Xbig <- model.matrix(modgam)
    # D <- modgam$

  }

  return(list(
    Vi = VbetajInv,
    lambdas = lambdas,
    Dbig = Dbig,
    vTs = vTs)
  )

}

calcLambda <- function(mod, mstop)
{

  Dp <- extract(mod,"penalty")
  Xp <- extract(mod,"design")

  resList <- lapply(mstop, function(m){

    coefs <- coef(mod[m])
    nc <- names(coefs)
    Xs <- Xp[names(Xp)%in%nc]
    res <- mod[m]$resid()
    Ds <- Dp[names(Dp)%in%nc]

    Xtr <- lapply(Xs, function(x) t(x)%*%res)
    Dbeta <- lapply(1:length(coefs), function(i) Ds[[i]]%*%coefs[[i]])

    lambdas <- sapply(1:length(coefs), function(i) crossprod(Dbeta[[i]])^(-1)*t(Dbeta[[i]])%*%Xtr[[i]])
    names(lambdas) <- nc
    return(lambdas)

  })

  names(resList) <- mstop
  return(resList)

}

refitmboostwithgam <- function(mod)
{

  sel <- sort(unique(selected(mod)))
  bl <- mod$baselearner[sel]
  isSpline <- grepl("bbs", names(bl))
  inds <- (1:length(bl))[isSpline]
  envs <- lapply(inds, function(i) environment(bl[[i]]$dpp))
  degree <- sapply(envs, function(x) (x$args$degree))
  thisknots <- lapply(1:length(envs), function(i){

    x <- envs[[i]]
    boundary.knots <- x$args$knots[[1]][[2]]
    knots <- x$args$knots[[1]][[1]]
    dx <- diff(boundary.knots)/(length(knots) + 1)
    bk_lower <- seq(boundary.knots[1] - degree[i] * dx, boundary.knots[1],
                    length = degree[i] + 1)
    bk_upper <- seq(boundary.knots[2], boundary.knots[2] + degree[i] * dx,
                    length = degree[i] + 1)
    ## complete knot mesh
    c(bk_lower, knots, bk_upper)
  })


  diffs <- sapply(envs, function(x) (x$args$differences))
  nx <- sapply(envs, function(x) names(x$mf))
  names(thisknots) <- nx
  comps <- c()
  for(i in 1:length(sel)){

    comps[i] <- if(i %in% inds)
      paste0("s(",nx[i], ", bs='ps', k=", length(thisknots[[i]]) - degree[i] - 3,
             ", m=c(", degree[i] + 1, ",", diffs[i], "))") else
      names(bl[[i]]$model.frame())

  }

  y <- mod$response
  formula <- as.formula(paste0("y ~ -1 +", paste(comps, collapse = " + ")))
  eval(gam(formula, knots = thisknots), envir = environment(mod))

}
