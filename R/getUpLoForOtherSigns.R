# getUpLoForOtherSigns <- function(y,v,M,Sigma=NULL,seqTau=seq(-10000,10000,length.out = 50),nrCores=1)
# {
#   
#   vvT <- tcrossprod(v,v)
#   vTy <- as.numeric(crossprod(v,y))
#   vTv <- as.numeric(crossprod(v,v))
#   if(is.null(Sigma)) Pvy <- vvT%*%y / vTv else Pvy <- Sigma%*%vvT%*%y / t(v)%*%Sigma%*%v
#   PvOy <- y - Pvy
#   
#   signsY <- sapply(seqTau, function(tau) sign(PvOy + tau*Pvy))
#   bounds <- mclapply(1:ncol(signsY),function(i){
#     
#     sy <- signsY[,i]
#     yy = sy*abs(y)
#     modmod <- mboost_fit(blg=mod$baselearner,
#                          response=as.numeric(yy),
#                          control = boost_control(mstop=mstop(mod),
#                                                  nu=mod$control$nu))
#     M <- mCreate(mod = modmod, Ups = getUpsilons(modmod),
#             Sigma = Sigma, conditionOn = conditionOn, crit = crit, ...)$M
# #     lapply(M,function(m) getIntersection(y = yy, v = v, M = m, Sigma = Sigma,
# #               Pvy=Pvy, PvOy=PvOy, vTy=vTy)[c("Vlo","Vup")])
#     calcUpLow(y = yy, v = v, M = M, Sigma = Sigma,
#               Pvy=Pvy, PvOy=PvOy, vTy=vTy)[c("Vlo","Vup")]
#   }, mc.cores=nrCores)
#   
#   Vlos <- unique(sapply(bounds,"[[","Vlo"))
#   Vups <- unique(sapply(bounds,"[[","Vup"))
#   
# }