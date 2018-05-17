#' @title Generic methods for iboost objects
#' 
#' @description Generic methods which can be used for objects fitted with the \code{iboost} function
#' 
#' @param x iboost object
#'
#' @method print iboost
#' @rdname methodsIboost
#' 
print.iboost <- function(x, ...)
{

  print(x$resDF)
  
}