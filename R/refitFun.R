refitFun <- function(mod, listFun=kronecker)
{
  
  desMats <- lapply(1:length(mod$baselearner), function(i) extract(mod, what="design", which=i)[[1]])
  desMats[sapply(desMats,class)=="list"] <- lapply(desMats[sapply(desMats,class)=="list"],
                                                   function(l)listFun(l[[1]],l[[2]]))
  Y <- mod$response - mod$offset
  
  Xbig <- as.data.frame(as.matrix(do.call("cbind",desMats)))
  colnames(Xbig) <- unlist(lapply(1:length(desMats),function(i)paste0("BL",i,"_coef",1:ncol(desMats[[i]]))))
  mod <- lm(Y~.-1,data=cbind(Y=Y,X=Xbig))
  return(mod)
  
}