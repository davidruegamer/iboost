# gcv.mboost <- function(mod, resMat=NULL){
#   
#   selCourse <- selected(mod)
#   selCovs <- sort(unique(selCourse))
#   n <- length(mod$ustart)
#   if(is.null(resMat)){
#   
#     hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
#     upsPart <- lapply(hatMats,function(x)(diag(n)-x*mod$control$nu))
#     resMat <- diag(n)
#     for(i in 1:(length(selCourse))) resMat <- upsPart[[which(selCovs==selCourse[i])]]%*%resMat
#     
#   }
#   
#   H = diag(n)-resMat
#   return(as.numeric(crossprod(mod$resid()))/(1-sum(diag(H))/n)^2)
#   
# }

library(Matrix)

gcv.mboost <- function(mod, Ups=NULL){
  
  H <- getHatMatrix(mod, Ups=Ups) 
  n <- length(mod$ustart)
  #return(as.numeric(crossprod(mod$resid()))/(1-sum(diag(H))/n)^2)
  return(log(as.numeric(crossprod(mod$resid()))/n) - 2*log(1-sum(Matrix::diag(H))/n))
}

aic.mboost <- function(mod, Ups=NULL){
  
  H <- getHatMatrix(mod, Ups=Ups)
  df <- sum(Matrix::diag(H))
  n <- length(mod$ustart)
  log(as.numeric(crossprod(mod$resid()))/n) + (1+df/n)/((1-df+2)/n)
  
}

aicc.mboost <- function(mod, Ups=NULL){
  
  H <- getHatMatrix(mod, Ups=Ups)
  df <- sum(Matrix::diag(H))
  n <- length(mod$ustart)
  log(as.numeric(crossprod(mod$resid()))/n) + 1 + (2+2*df)/(n-df-2)
  
}

bic.mboost <- function(mod, Ups=NULL){
  
  H <- getHatMatrix(mod, Ups=Ups)
  df <- sum(Matrix::diag(H))
  n <- length(mod$ustart)
  log(as.numeric(crossprod(mod$resid()))/n) + log(n*df)/n
  
}

gmdl.mboost <- function(mod, Ups=NULL){
  
  n <- length(mod$ustart)
  H <- getHatMatrix(mod, Ups=Ups)
  sigSqEst <- as.numeric(crossprod(mod$resid()))/n
  df <- sum(Matrix::diag(H))
  S = n*sigSqEst/(n-df)
  F = (crossprod(mod$response)-n*sigSqEst)/(df*S)
  return(log(S) + df/n * log(F))
  
}

getHatMatrix <- function(mod, Ups=NULL){
  
  n <- length(mod$ustart)
  if(is.null(Ups)) Ups <- getUpsilons(mod)
  if(is.list(Ups)) Ups <- Ups[[length(Ups)]] 
  Matrix::diag(n) - Ups
  
}

fpe.mboost <- function(mod, gamma=1, Ups=NULL){
  
  H <- getHatMatrix(mod, Ups=Ups)
  df <- sum(Matrix::diag(H))
  n <- length(mod$ustart)
  as.numeric(crossprod(mod$resid()))/n + gamma*df
  
  
}

crit.mboost <- function(mod, crits=c("AIC","AICC","BIC","GMDL","FPE"), gamma=1, Ups=NULL)
{
  
  vals <- c()

  # general stuff  
  n <- length(mod$ustart)
  H <- getHatMatrix(mod, Ups=Ups)
  sigSqEst <- as.numeric(crossprod(mod$resid()))/n
  df <- sum(Matrix::diag(H))

  # get criteria
  if("AIC"%in%crits) vals <- append(vals, log(sigSqEst) + (1+df/n)/((1-df+2)/n) )
  if("AICC"%in%crits) vals <- append(vals, log(sigSqEst) + 1 + (2+2*df)/(n-df-2) )
  if("BIC"%in%crits) vals <- append(vals, log(sigSqEst) + log(n*df)/n )
  if("GMDL"%in%crits){
    S = n*sigSqEst/(n-df)
    F = (crossprod(mod$response)-n*sigSqEst)/(df*S)
    vals <- append(vals, log(S) + df/n * log(F) )
  }
  if("FPE"%in%crits) vals <- append(vals, sigSqEst + gamma*df)
  
  names(vals) <- c("AIC","AICC","BIC","GMDL","FPE")[c("AIC","AICC","BIC","GMDL","FPE")%in%crits]
  vals
}

diagMinusMat <- function(mat)
{
  
  mat[diag(mat)] <- diag(mat)-1
  mat*(-1)
  
  
}

critEff.mboost <- function(mod, lastIterToConsider, gamma=1, listFun=kronecker, printPG=FALSE)
{
  
  if(lastIterToConsider > mstop(mod)) mod <- mod[lastIterToConsider]
  
  # general stuff  
  library(Matrix)
  crits <- c("AIC","AICC","BIC","GMDL","FPE")
  Y <- mod$response
  n <- length(mod$ustart)
  nu <- mod$control$nu
  selCourse <- selected(mod)
  hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
  if(any(sapply(hatMats,is.null))){
    
    desMats <- lapply(1:length(mod$baselearner), function(i) extract(mod, what="design", which=i)[[1]])
    desMats[sapply(desMats,class)=="list"] <- lapply(desMats[sapply(desMats,class)=="list"],
                                                     function(l)listFun(l[[1]],l[[2]]))
    
    if(inherits(mod,"glmboost")) hatMatsLeft <- lapply(desMats,function(x)solve(t(x)%*%x)%*%t(x)) else{
      
      lambdas <- sapply(1:length(mod$baselearner), function(i) extract(mod, what="lambda", which=i)[[1]])
      pens <- lapply(1:length(mod$baselearner), function(i) extract(mod, what="penalty", which=i)[[1]])
      hatMatsLeft <- lapply(1:length(pens),function(i)
        solve(Matrix::t(desMats[[i]])%*%desMats[[i]] + lambdas[i]*pens[[i]])%*%Matrix::t(desMats[[i]])
      )
      
    }
    hatMats <- mapply(function(x,h)x%*%h,x=desMats,h=hatMatsLeft,SIMPLIFY=FALSE)
  }
  upsPart <- lapply(hatMats,function(x)as(diagMinusMat(x*nu), "sparseMatrix"))
  len <- length(selCourse)
  # Upsilons <- vector("list",len)
  UpsilonsI <- as(Matrix::diag(n),"sparseMatrix")
  
  # matrix for results
  vals <- matrix(rep(NA,lastIterToConsider*length(crits)),ncol=length(crits))
  
  if(printPG) pb <- txtProgressBar(min = 1, max = lastIterToConsider)
  
  for(i in 1:lastIterToConsider){
    
    if(printPG) setTxtProgressBar(pb, i)
    
    dH <- 1-diag(UpsilonsI)
    sigSqEst <- as.numeric(crossprod(mod[i]$resid()))/n
    df <- sum(dH)
    
    # get criteria
    vals[i,1] <- log(sigSqEst) + (1+df/n)/((1-df+2)/n)
    vals[i,2] <- log(sigSqEst) + 1 + (2+2*df)/(n-df-2)
    vals[i,3] <- log(sigSqEst) + log(n*df)/n
    
    S = n*sigSqEst/(n-df)
    F = (crossprod(Y)-n*sigSqEst)/(df*S)
    
    vals[i,4] <- log(S) + df/n * log(F) 
    vals[i,5] <- sigSqEst + gamma*df
    
    # print(dim(upsPart[[selCourse[i]]]))
    # print(dim(UpsilonsI))
    if(i!=lastIterToConsider) UpsilonsI <- upsPart[[selCourse[i]]] %*% UpsilonsI else{
      if(printPG) close(pb)
    }
    
  }

  colnames(vals) <- crits
  rownames(vals) <- 1:lastIterToConsider
  vals
}

# 
# 
# #### Tests
# 
# ### correct residual matrix?
# 
# library(mboost)
# 
# n <- 100
# x1 <- scale(rnorm(n))
# x2 <- scale(rnorm(n))
# y <- 3 * x1 + x2 + rnorm(n)
# y <- as.numeric(y-mean(y))
# 
# mod <- mboost_fit(response=y, blg=list(bols(x1,intercept=FALSE),bols(x2,intercept=FALSE)))
# 
# selCourse <- selected(mod)
# selCovs <- sort(unique(selCourse))
# hatMats <- lapply(mod$basemodel,function(b)b$hatvalues())
# upsPart <- lapply(hatMats,function(x)(diag(n)-x*mod$control$nu))
# resMat <- diag(n)
# for(i in 1:(length(selCourse))) resMat <- upsPart[[which(selCovs==selCourse[i])]]%*%resMat
# all.equal(as.numeric(resMat%*%y),as.numeric(mod$resid()))
# 
# # -> Yes
# 
# ## Comparison GCV <-> CV
# 
# 
# set.seed(42)
# 
# cFr <- do.call("rbind", mclapply(1:100,function(sss){
#   
#   set.seed(sss)
#   
#   n <- 100
#   x1 <- rnorm(n)
#   x2 <- rnorm(n) 
#   y <- 3 * x1 + x2 + rnorm(n)
#   y <- y-mean(y)
#   
#   mod <- mboost(as.numeric(y)~bols(x1,intercept=FALSE)+bols(x2,intercept=FALSE),
#                 control = boost_control(mstop=200))
# 
#   cvStop <- mstop(cvrisk(mod))
#   gcvStop <- which.min(sapply(1:mstop(mod),function(i)gcv.mboost(mod[i])))
#   
#   return(c(cvStop,gcvStop))
#   
# },mc.cores=20))
# 
# boxplot(cFr) # gcv seems to be a better choice?
# at least in the context of the very small improvements in plot(cvrisk(mod[200]))
#   
#   