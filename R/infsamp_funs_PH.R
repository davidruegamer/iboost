infsamp_via_PH <- function(refitFun, y, vT, confun, var, ncore)
{
  
  Pv <- lapply(vT, function(x) crossprod(x) / as.numeric(tcrossprod(x)))
  vnorm <- lapply(vT, function(x) t(x) / as.numeric(tcrossprod(x)))
  yperp <- lapply(Pv, function(p) y - p %*% y) 
  r <- lapply(vT, function(x) as.numeric(x %*% y))
  
  lims <- list()
  
  for(j in 1:length(vT)){
    
    lims[[j]] <- findtrunclims_via_PH(refitFun, r = r[[j]], yperp = yperp[[j]],
                                      var = var, vnorm = vnorm[[j]], vT = vT[j], ncore = ncore, confun = confun)
    
    
  }
  
  return(lims)
  
}

findtrunclims_via_PH <- function(refitFun, r, yperp, var, vnorm, vT, ncore, confun, tol=1e-6)
{
  
  # initvalsLo <- r - 3^c(3:-8)
  # initvalsUp <- r + 3^c(3:-8)
  # loLog <- sapply(initvalsLo, function(val) confun(refitFun(as.numeric(yperp + val * Pv))))
  # upLog <- sapply(initvalsUp, function(val) confun(refitFun(as.numeric(yperp + val * Pv))))
  # initvalLo <- yperp + Pv * initvalsLo[min(which(loLog))]
  # initvalUp <- yperp + Pv * initvalsUp[max(which(upLog))]
  
  
  checkfun <- function(val) confun(refitFun(as.numeric(yperp + val * vnorm)))
  
  vlo <- binsearch_via_PH(ystart = yperp + vnorm * r, tol = tol, refitFun,
                          yperp = yperp, vnorm = vnorm, vT = vT, checkfun = checkfun, 
                          ncore = ncore, what = "lower", var = var)
  vup <- binsearch_via_PH(ystart = yperp + vnorm * r, tol = tol, refitFun, 
                          yperp = yperp, vnorm = vnorm, vT = vT, checkfun = checkfun, 
                          ncore = ncore, what = "upper", var = var)
  
  c(vlo,vup)
  
}



binsearch_via_PH <- function(ystart, tol, yperp, vnorm, vT, 
                             refitFun, checkfun, ncore, what, var)
{
  
  # count <- 1
  # this_tol <- 100
  # int_diff <- 0
  # alreadyInGap <- FALSE
  muInit <- as.numeric(vT[[1]]%*%ystart)
  
  getboundsfun <- function(val) polyh_inf(obj = refitFun(as.numeric(val * vnorm + yperp)), vT = vT, 
                                          var = var, ncore = ncore, 
                                          alpha = 0.05, # doesn't matter
                                          doInf=FALSE 
  )[[1]]
  
  newrval <- oldrval <-  polyh_inf(obj = refitFun(as.numeric(ystart)), vT = vT, var = var, ncore = ncore, 
                                   alpha = 0.05, # doesn't matter
                                   doInf=FALSE 
  )[[1]]
  
  range <- getSearchRange(mu = muInit, lower = newrval[1], 
                          upper = newrval[2], var = var, tol = tol)
  
  wu <- what=="upper"
  
  mincheck <- -Inf 
  maxcheck <- Inf
  len = 25
  
  while(mincheck==-Inf | maxcheck==Inf){
    
    # find minimum and maximum value
    
    # get some values to check their congruency
    startseq <- seq(range[1], range[2], l = len)
    # check values
    check_these_vals <- unlist(mclapply(startseq, function(x) checkfun(x), mc.cores = ncore))
    # use the smallest (mincheck) and largest (maxcheck) value, which is
    # closest to muInit without non-congruent values in between
    mincheck <- which(!check_these_vals & startseq < muInit) 
    if(length(mincheck)==0){
      mincheck <- -Inf
      range[1] <- startseq[min(which(startseq < muInit))]
    }else mincheck <- max(mincheck) + 1
    
    maxcheck <- which(!check_these_vals & startseq > muInit)
    
    if(length(maxcheck)==0){
      maxcheck <- Inf
      range[2] <- startseq[max(which(startseq > muInit))]
    }else maxcheck <- min(maxcheck) - 1
    
  }
  
  startseq <- startseq[c(mincheck, maxcheck)]
  
  stopatval <- ifelse(wu, startseq[2], startseq[1])
  
  endcheck <- getboundsfun(stopatval)[2 - wu]
  
  innercheck <- newrval[1 + wu]
  final_bound <- confirmedBound(int = sort(c(innercheck, endcheck)), checkfun = checkfun,
                                getboundsfun = getboundsfun, upperLog = wu, tol = tol,
                                sd = sqrt(var), otherBound = newrval[1 + !wu])
  return(final_bound[wu + 1])
  
}





getSearchRange <- function(mu, lower, upper, var, tol)
{
  
  atild <- pnorm(lower, mean = 0, sd = sqrt(var))
  btild <- pnorm(upper, mean = 0, sd = sqrt(var))
  qtild <- pnorm(mu, mean = 0, sd = sqrt(var))
  qma <- qtild - atild
  tola <- tol*(atild - 1)
  inA <- (tol*btild^2)/((1+tol)*btild + qtild)
  inB <- (qma + tola*atild) / (qma + tola)
  if(inA < tol) inA <- tol*10
  if(inB > 1-tol) inB <- 1 - tol*10
  a <- qnorm(inA, mean = 0, sd = sqrt(var))
  b <- qnorm(inB, mean = 0, sd = sqrt(var))
  return(c(a,b))
  
}




checkNegl <- function(int, sd, lower, upper, tol, print = FALSE)
{
  
  d <- c()
  
  if(is.null(nrow(int))){
    
    d <- ptnorm(q = int[2], mean = 0, sd = sd, lower = lower, upper = upper) - 
      ptnorm(q = int[1], mean = 0, sd = sd, lower = lower, upper = upper) 
    
  }else{
    
    for(i in 1:nrow(int)){
      
      d <- c(d,
             ptnorm(q = int[i,2], mean = 0, sd = sd, lower = lower, upper = upper) - 
               ptnorm(q = int[i,1], mean = 0, sd = sd, lower = lower, upper = upper) 
      )
      
    }
    
  }
  
  if(print) print(d)
  d < tol
  
}

confirmedBound <- function(int, checkfun, getboundsfun, nrInts = 100, probTol = 1e-5,
                           upperLog, tol, trace = TRUE, ncore = 25, sd, otherBound)
{
  
  ss <- seq(int[1], int[2], l = nrInts + 1)
  ss <- ss[-length(ss)] + diff(ss)/2
  
  checkvals <- unlist(mclapply(ss, function(x) checkfun(x), mc.cores = ncore))
  ints <- mclapply(ss[checkvals], function(x) getboundsfun(x), mc.cores = ncore)
  ints <- Intervals(do.call("rbind", ints))
  
  minconf <- min(ints)
  maxconf <- max(ints)
  
  lowB <- ifelse(upperLog, otherBound, minconf)
  upB <- ifelse(upperLog, maxconf, otherBound)
  
  ints <- interval_union(round(ints, -1*log(tol, base=10)))
  
  remints <- interval_complement(ints)
  remints <- remints[c(-1,-nrow(remints))]
  if(NROW(remints)==0) return(TRUE)
  cc <- checkNegl(remints, sd = sd, lower = lowB, upper = upB,
                  tol = probTol)
  remints <- remints[!cc,]
  if(NROW(remints)==0) return(TRUE)
  
  while(NROW(remints)>0){
    
    print(remints)
    
    int1 <- if(is.null(nrow(remints)) || nrow(remints) == 0) remints else as.numeric(remints[1,])
    ss <- seq(int1[1], int1[2], l = 11)
    ss <- ss[-length(ss)] + diff(ss)/2
    
    checkvals <- unlist(mclapply(ss, function(x) checkfun(x), mc.cores = ncore))
    ints <- mclapply(ss, function(x) getboundsfun(x), mc.cores = ncore)
    
    if(any(!checkvals)){
      
      if(upperLog){ 
        
        maxval <- ints[[min(which(!checkvals))]][1]
        delints <- sapply(1:nrow(remints), function(i) remints[i,1] > maxval)
        remints <- remints[!delints,]
        maxconf <- max(remints)
        
      }else{
        
        minval <- ints[[max(which(!checkvals))]][2]
        delints <- sapply(1:nrow(remints), function(i) remints[i,2] < minval)
        remints <- remints[!delints,]
        minconf <- min(remints)
        
      }
      
      lowB <- ifelse(upperLog, otherBound, maxconf)
      upB <- ifelse(upperLog, minconf, otherBound)
      
    }else{
      
      ints <- Intervals(do.call("rbind", ints))
      ints <- interval_union(round(ints, -1*log(tol, base=10)))
      
      remints1 <- interval_complement(ints)
      remints1 <- remints1[c(-1,-nrow(remints1))]
      
      if(NROW(remints1)==0){
        
        if(is.null(nrow(remints))) remints <- Intervals() else 
          remints <- remints[-1,]
        
      }else{
        
        cc <- checkNegl(remints1, sd = sd, lower = lowB, upper = upB,
                        tol = probTol)
        remints1 <- remints1[!cc,]
        
        if(is.null(nrow(remints))){
          
          remints <- remints1
          
        }else{
          
          remints <- remints[-1,]
          remints <- rbind(remints1, remints)
          
        }
        
      }
    }
  }
  
  return(c(minconf, maxconf))
  
}



