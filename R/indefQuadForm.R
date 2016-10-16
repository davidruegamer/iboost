indefQuadFrom <- function(X, A, a=rep(0,length(X)), d=0, mu, Sigma, B = NULL, roundEV=14)
{
  
#   # Example from Mohsenipour and Provost 2011
#   mu = c(100,0,-50,150,5)
#   a = c(-1,2,3,1,1)
#   d = 6
#   A = matrix(c(1,-0.9,-1,0,-5,-0.9,1,1,2,1,-1,1,2,3,1,0,2,3,-1,0,-5,1,1,0,1),ncol=5)
#   Sigma = matrix(c(rep(3,25)),ncol=5)
#   Sigma[3,3] <- 5
#   Sigma[4,] <- Sigma[,4] <- c(2,2,2,2,0)
#   Sigma[5,] <- Sigma[,5] <- c(0,0,0,0,1)
  
  
  library(expm)
  
  ind <- 1:length(X)
  
  if(is.null(B)){
    
    if(!is.matrix(Sigma)) B <- sqrt(Sigma)*diag(length(X)) else{
      
      B <- try(chol(Sigma))
      if(class(B)=="try-error"){
        
        E <- eigen(Sigma)
        rev <- round(E$values,roundEV)
        ind <- rev!=0
        B <- E$vectors[,ind]%*%diag(sqrt(rev[ind]))
        
      }
    
    }
  }
  
  c1 <- as.numeric(t(mu)%*%A%*%mu + t(a)%*%mu + d )
  SA <- Sigma%*%A
  
  # eigendecomposition and terms
  E <- eigen(t(B)%*%A%*%B)
  lambdaBiggerZero <- E$values[E$values>0]
  lambdaEqualZero <- E$values[E$values==0]
  lambdaSmallerZero <- E$values[E$values<0]
  P <- E$vectors # already normalized by eigen-function
  bStarT <- t(mu)%*%A%*%B%*%P # 
  mT <- t(a)%*%B%*%P
  nj <- 0.5*mT + bStarT
  theta <- length(lambdaEqualZero)
  nu1 <- nj[E$values>0]/lambdaBiggerZero
  nu2 <- nj[E$values<0]/lambdaSmallerZero
  
  Qfun <- function(X, mu, Sigma, s){
    
    sthCumulantFun <- function(s)
    {
      
      trASigma <- sum(diag(A%*%Sigma))
      if(s==1) return(trASigma + c1) else
      {
        
        SA2 <- SA%^%(s-2)
        SA1 <- SA2%*%SA      
        add <- if(!(sum(sapply(a,all.equal,0))==length(X)))
          0.25*t(a)%*%SA2%*%Sigma%*%a + t(a)%*%SA1%*%A%*%mu else 0
        
        return(
          2^(s-1)*factorial(s)%*%( trASigma^s / s +  t(mu)%*%SA1%*%A%*%mu + add) 
        )
        
      }
    }
    
    hthMomentFun <- function(h)
    {
      
      sum(sapply(0:(h-1),function(i){
        
        factorial(h-1) / (factorial(h-1-i)*factorial(i)) * sthCumulantFun(h-i) * hthMomentFun(i)
        
      }))
      
    }
    
  }
  
  
}