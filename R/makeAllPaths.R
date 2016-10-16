makeAllPaths <- function(selCourse, 
                         nuSmaller1 = any(duplicated(selCourse)), 
                         p = unique(selCourse), 
                         mstop = length(selCourse),
                         setEqualRleLength = TRUE,
                         selectedSetOnly = TRUE,
                         noVarOfNotSelectedSet = TRUE,
                         onlySamePos = TRUE
                         )
{
  
  library(gtools)
  
  if(p^mstop > 1e7) stop("Too many combinations!")
  if(sum(abs(diff(selCourse)))==0) return(as.data.frame(matrix(selCourse, nrow=1))) # there is only one path
  
  # generate all combinations
  ap <- as.data.frame(permutations(p, mstop, repeats.allowed = nuSmaller1))
  
  # delete duplicates or paths which dont change
  ap <- unique(ap[apply(ap,1,function(x)sum(abs(diff(x)))!=0),])
  
  # use only those vars, which were selected / not those, which were not selected
  if(selectedSetOnly) ap <- ap[apply(ap,1,function(x)all(selCov%in%x)),]
  if(noVarOfNotSelectedSet) ap <- ap[apply(ap,1,function(x)!any(notSel%in%x)),]
  
  # set equal RLE length and path
  if(setEqualRleLength){
    ap <- ap[apply(ap,1,function(x)length(rle(x)$values)==length(selRle) & 
                     sum(abs(rle(x)$values[1:length(selRle)]-selRle))==0),]
  }
    
  # use path only if the all unique variables are at the same position
  if(onlySamePos) ap <- ap[apply(ap,1,function(x)
    sum(abs(x[duplicated(selCourse)]-selCourse[duplicated(selCourse)]))==0),]  
  
  return(ap)
  
}