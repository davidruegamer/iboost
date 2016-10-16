getSplineBlDist <- function(mod, which=1)
{
  
  library(gtools)
  
  selCourse <- selected(mod)
  UpsAndHM <- getUpsilons(mod, returnUpsOnly = FALSE)
  indSum <- which(selCourse==which)
  nu <- mod$control$nu
  
  sumUps <- Reduce("+",lapply(indSum,function(l)UpsAndHM$Upsilons[[l]])) # not l-1 as Upsilon is already without last iter 
  E <- nu*UpsAndHM$hatMatsLeft[[which]]%*%sumUps
  sumVar <- nu^2*Reduce("+",lapply(indSum,function(l)
    UpsAndHM$hatMatsLeft[[which]]%*%UpsAndHM$Upsilons[[l]]%*%t(UpsAndHM$hatMatsLeft[[which]])%*%
      UpsAndHM$hatMatsLeft[[which]]%*%UpsAndHM$Upsilons[[l]]%*%t(UpsAndHM$hatMatsLeft[[which]])))
  combs <- combn(indSum,2)
  sumCov <- nu^2*Reduce("+",lapply(1:ncol(combs),function(i){
    l <- combs[1,i]
    k <- combs[2,i]
    UpsAndHM$hatMatsLeft[[which]]%*%UpsAndHM$Upsilons[[l]]%*%t(UpsAndHM$hatMatsLeft[[which]])%*%
      UpsAndHM$hatMatsLeft[[which]]%*%t(UpsAndHM$Upsilons[[k]])%*%t(UpsAndHM$hatMatsLeft[[which]])
  }))
  V <- sumVar + 2*sumCov
  return(list(E=E,V=V,sumVar=sumVar,sumCov=sumCov))
  
}




# 
# 
# set.seed(2010)
# library(mgcv)
# library(mboost)
# 
# res <- mclapply(1:1000,function(sss){
#   
#   set.seed(sss)
#   
#   n <- 100
#   x1 <- seq(0,1,length.out = n)
#   x2 <- rnorm(n)
#   sinX1 <- sin(13*x1)
#   
#   y <- sinX1 + rnorm(n)
#   
#   x1 <- scale(x1, center = T)
#   y <- as.numeric(scale(y, center = T))
#   
# #   modTrue <- gam(sinX1~-1+s(x1))
# #   trueBeta <- coef(modTrue)
# #   trueMean <- predict(modTrue)
# #   
#   ### set up base-learners
#   spline1 <- bbs(x1, knots = 5, df = 4)
#   spline2 <- bbs(x2, knots = 5, df = 4)
#   ### fit model, component-wise
#   mod1 <- mboost_fit(list(spline1,spline2), y, control = boost_control(mstop=500, nu=1))
#   
#   cvr <- cvrisk(mod1, fun=function(mod)return(list(risk=mod$risk()[1:mstop(mod)],
#                                                    sigEst=sapply(1:mstop(mod),function(i)
#                                                      sd(mod[i]$resid())))), mc.cores=5)
#   
#   ms <- which.min(rowSums(sapply(cvr,"[[","risk")))
#   mod1[ms]
#   if(length(unique(selected(mod1)))==1) return(NULL)
#   crossCVsigma <- rowMeans(sapply(cvr,"[[","sigEst"))[ms]
#   
#   
#   desMats <- lapply(1:length(mod1$baselearner), function(i) extract(mod1, what="design", which=i)[[1]])
#   lambdas <- sapply(1:length(mod1$baselearner), function(i) extract(mod1, what="lambda", which=i)[[1]])
#   pens <- lapply(1:length(mod1$baselearner), function(i) extract(mod1, what="penalty", which=i)[[1]])
#   hatMatsLeft <- lapply(1:length(pens),function(i)
#     solve(Matrix::t(desMats[[i]])%*%desMats[[i]] + lambdas[i]*pens[[i]])%*%Matrix::t(desMats[[i]])
#   )
#   
#   # Ve = lapply(1:length(hatMatsLeft),function(v)v%*%t(v))
#   Vp = lapply(1:length(pens),function(i)
#     solve(Matrix::t(desMats[[i]])%*%desMats[[i]] + lambdas[i]*pens[[i]]))
#   
#   res <- lapply(1:length(unique(selected(mod1))),function(i)getSplineBlDist(mod1, i))
#   
#   betas <- coef(mod1)
#   
#   chi9_Vp <- sapply(1:length(unique(selected(mod1))),function(i) t(betas[[i]])%*%Vp[[i]]%*%betas[[i]]*crossCVsigma)
#   chi9_Vboost <- sapply(1:length(unique(selected(mod1))),function(i) t(betas[[i]])%*%res[[i]]$V%*%betas[[i]]*crossCVsigma)
#   
#   return(list(chi9_Vp=chi9_Vp, chi9_Vboost=chi9_Vboost, betas=betas, res=res))
#   
# },mc.cores=5)
# 
# bound <- qchisq(0.95,df=9)
# chi9 <- do.call("rbind",lapply(res,function(x)x[c("chi9_Vp","chi9_Vboost")]))
# chi9p <- do.call("rbind",chi9) < bound

# -> besser: ueber Hat-Matrix EDF bestimmen und dann approx. Cov berechnen
