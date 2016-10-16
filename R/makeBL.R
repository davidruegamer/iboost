makeBL <- function(charName,blFun,data,...) {
  
  ########### make Baselearners using character names
  ########### blFun specifies the type of BL
  
  temp <- as.data.frame(data[,charName])
  colnames(temp) <- charName
  bl <- tryCatch(blFun(temp,...), error=function(e) bols(temp))
  bl$set_names(charName)
  bl
  
}