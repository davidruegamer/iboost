testStatWood <- function(X, V, rank, p)
{
  
  qrx <- qr(X,tol=0)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot,drop=FALSE]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)
  ## remove possible ambiguity from statistic...
  siv <- sign(ed$vectors[1,]);siv[siv==0] <- 1
  ed$vectors <- sweep(ed$vectors,2,siv,"*")
  
  k <- max(0,floor(rank)) 
  nu <- abs(rank - k)     ## fractional part of supplied edf
  
  if (nu>0) k1 <- k+1 else k1 <- k
  
  ## check that actual rank is not below supplied rank+1
  r.est <- sum(ed$values > max(ed$values)*.Machine$double.eps^.9)
  if (r.est<k1) {k1 <- k <- r.est;nu <- 0;rank <- r.est}
  
  ## Get the eigenvectors...
  # vec <- qr.qy(qrx,rbind(ed$vectors,matrix(0,nrow(X)-ncol(X),ncol(X))))
  vec <- ed$vectors
  if (k1<ncol(vec)) vec <- vec[,1:k1,drop=FALSE]
  
  ## deal with the fractional part of the pinv...
  if (nu>0&&k>0) {
    if (k>1) vec[,1:(k-1)] <- t(t(vec[,1:(k-1)])/sqrt(ed$val[1:(k-1)]))
    b12 <- .5*nu*(1-nu)
    if (b12<0) b12 <- 0
    b12 <- sqrt(b12)
    B <- matrix(c(1,b12,b12,nu),2,2)
    ev <- diag(ed$values[k:k1]^-.5,nrow=k1-k+1)
    B <- ev%*%B%*%ev
    eb <- eigen(B,symmetric=TRUE)
    rB <- eb$vectors%*%diag(sqrt(eb$values))%*%t(eb$vectors)
    # vec1 <- vec
    # vec1[,k:k1] <- t(rB%*%diag(c(-1,1))%*%t(vec[,k:k1]))
    vec[,k:k1] <- t(rB%*%t(vec[,k:k1]))
  } else {
    # vec1 <- 
      vec <- if (k==0) t(t(vec)*sqrt(1/ed$val[1])) else
      t(t(vec)/sqrt(ed$val[1:k]))
    if (k==1) rank <- 1
  }
  ## there is an ambiguity in the choise of test statistic, leading to slight
  ## differences in the p-value computation depending on which of 2 alternatives 
  ## is arbitrarily selected. Following allows both to be computed and p-values
  ## averaged (can't average test stat as dist then unknown) 
  d <- t(vec)%*%(R%*%p)
  # d1 <- t(vec1)%*%(R%*%p)
  # dorg = d1
  dstat <- sum(d^2) 
  # d1 <- sum(d1^2)
  ##d <- d1 ## uncomment to avoid averaging
  
  return(list(inn=t(vec)%*%R%*%t(X),dstat=dstat)) # list(d = d, d1 = d1))
  
}