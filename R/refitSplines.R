refitSplines <- function(mod, eps = 1e-12, 
                         lambdaInit = log(unlist(extract(mod, "lambda"))),
                         trace = 0, abstol = 1e-3, maxit = 500,
                         addConstraint = TRUE, extraRidge = TRUE)
{
  
  notPen <- lambdaInit==-Inf
  if(extraRidge) lambdaInit <- c(rep(0,length(lambdaInit)), lambdaInit)
  # extract smoothing parameter of unpenalized effects
  
  # penalties
  Pjs <- extract(mod, "penalty")
  # Pjs <- Pjs[!notPen]
  
  # penalty matrix depending on smoothing parameters
  ll <- length(lambdaInit)
  # which have actually not a zero penalization
  lambdaInit <- lambdaInit[c(rep(TRUE,ll/2),!notPen)]
  lln <- length(lambdaInit)
  
  # residual and hat matrix
  Ups <- getUpsilons(mod)
  X <- getDesignmat(mod, split = T)
  # add constraint to each effect?
  if(addConstraint){
    
    Z <- lapply(X, function(x){ 
      
      c <- colSums(x)
      qr_c <- qr(c)
      if( any(class(qr_c) == "sparseQR") ){
        rank_c <- qr_c@Dim[2]
      }else{
        rank_c <- qr_c$rank 
      } 
      Q <- qr.Q(qr_c, complete=TRUE)
      return(Q[  , (rank_c + 1) : ncol(Q) ])
      
    })
    X <- lapply(1:length(X), function(i) X[[i]] %*% Z[[i]])
    Pjs <- lapply(1:length(Pjs), function(i) t(Z[[i]]) %*% Pjs[[i]] %*% Z[[i]])
    
  }
  # Xc <- do.call("cbind", lapply(X, scale, scale=FALSE))
  ncols <- sapply(X, NCOL)
  XtXi <- lapply(X, crossprod)
  X <- do.call("cbind", X)
  n <- nrow(X)
  hatBoost <- diag(n)-Ups[[length(Ups)]]
  XtX <- crossprod(X)
  y <- mod$response
  
  # delta for helping inversion
  delta <- diag(rep(eps, ncol(X))) 

  if(extraRidge){
    
    P <- function(lambdas) 
    {
      
      lambdasTemp <- c(lambdas[1:(ll/2)],rep(-Inf,ll/2))
      lambdasTemp[ll/2+which(!notPen)] <- lambdas[ll/2+which(!notPen)]
      return(
        bdiag(lapply(1:(ll/2), function(j) exp(lambdasTemp[j]) * XtXi[[j]])) + 
          bdiag(lapply((ll/2+1):ll, function(j) exp(lambdasTemp[j]) * Pjs[[j-ll/2]]))
      )
    }
    
  }else{
  
    P <- function(lambdas) bdiag(lapply(1:length(lambdas), 
                                        function(j) exp(lambdas[j])*Pjs[[j]]))

  }

  V <- function(this_lambdas) XtX + delta + P(this_lambdas)
  
  invV <- function(this_lambdas, Vinput=NULL)
  {

    if(!is.null(Vinput)) this_V <- Vinput else
      this_V <- V(this_lambdas)
    eV <- eigen(this_V)
    (eV$vectors)%*%diag(1/eV$values)%*%t(eV$vectors)
    
  }
  
  # mse and gradient
  mse <- function(this_lambdas)
  {
    
    hatPLS <- X%*%invV(this_lambdas)%*%t(X)
    sum(diag(tcrossprod(hatPLS - hatBoost))) 
    
  }
  
  gradmse <- function(this_lambdas)
  {
    
    Vinv <- invV(this_lambdas)
    hatPLS <- X%*%Vinv%*%t(X)
    sapply(1:length(this_lambdas), function(i)
    {
      
      lambdas_dummy <- this_lambdas
      lambdas_dummy[-i] <- 0
      P_i <- P(lambdas_dummy)
      Htilde <- as.matrix(X%*%Vinv%*%Vinv%*%P_i%*%t(X))
      sum(diag(Htilde%*%hatBoost)) - sum(diag(hatPLS%*%Htilde))
      
    })
    
    
  }

  # find optimum
  if(length(Pjs)==1){
    om <- optim(par = 1, fn = mse)
  }else{
    om <- optim(par = lambdaInit, fn = mse, gr = gradmse, 
                control = list(trace = trace, abstol = abstol, 
                               maxit = maxit))
  }
  
  st <- c(1, cumsum(ncols[-length(ncols)])+1)
  en <- cumsum(ncols)
  Vinv <- invV(this_lambdas = NA, Vinput = crossprod(X) + delta + P(om$par))
  
  K <- Vinv%*%t(X)    
  
  fvals <- lapply(1:length(Pjs), function(i)
  {
    
    Xi <- X
    Xi[,-1*(st[i]:en[i])] <- 0
    hPLSi <- Xi%*%K
    return(as.numeric(hPLSi%*%y))
    
  })
  
  return(list(lambdas = exp(om$par), res = om, fitvals = fvals, Vinv = Vinv, X = X,
              st = st, en = en)) #, Z = Z))
  
}
