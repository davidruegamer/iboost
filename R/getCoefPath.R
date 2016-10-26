getCoefPath <- function(obj, what = c("path", "incr", "sign"))
{
  
  what <- match.arg(what)
  
  signCourse <- if(inherits(obj, "glmboost"))
    sapply(1:mstop(obj), function(m) as.data.frame(obj[m]$coef())) else 
      sapply(1:mstop(obj),function(m) sapply(obj[m]$coef(),"[[",1))
  
  nams <- attr(signCourse[[length(signCourse)]], "names")
  signCourseS <- do.call("rbind",lapply(signCourse, function(sc){
    
    lenSc <- length(sc)
    
    if(lenSc<length(nams)){
      
      namSc <- names(sc)
      namsN <- nams[!nams%in%namSc]
      sc <- c(rep(0,length(namsN)),sc)
      names(sc) <- c(namsN, namSc)
      sc <- sc[nams]
      
    }
    
    unlist(sc)
    
  }))
  
  if(what == "path") return(signCourseS)
  
  signCoursePM <- apply(rbind(rep(0, ncol(signCourseS)), signCourseS), 2, diff)
  
  if(what == "incr") return(signCoursePM)
  
  # what == "sign"
  sig <- if(length(signCoursePM)==1) sign(signCoursePM) else rowSums(sign(signCoursePM))
  return(sig)
  
}