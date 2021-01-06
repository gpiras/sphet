#' @name listw2dgCMatrix
#' @aliases listw2dgCMatrix
#' @title Interface between Matrix class objects and weights list 
#' @description Interface between Matrix class objects and weights list 
#' @usage listw2dgCMatrix(listw, zero.policy = NULL) 
#' @param listw a \code{listw} object created for example by \code{nb2listw} 
#' @param zero.policy See \code{lagsarlm} for details
#' @return Matrix class object: a sparse Matrix
#' @author Gianfranco Piras \email{gpiras@mac.com}
#' @examples 
#' library(spdep)
#' data(columbus)
#' listw <- nb2listw(col.gal.nb)
#' spW <- listw2dgCMatrix(listw)
#' @keywords spatial
#' @export



listw2dgCMatrix<-function (listw, zero.policy=NULL) 
{
    if (!inherits(listw, "listw")) 
        stop("not a listw object")
    if (is.null(zero.policy))
        zero.policy <- get.ZeroPolicyOption()
    stopifnot(is.logical(zero.policy))
    n <- length(listw$neighbours)
    cardw <- card(listw$neighbours)
    p0 <- as.integer(c(0, cumsum(cardw)))
    scard <- sum(cardw)
    t<-unlist(listw$neighbours)
    if (zero.policy) t <- t[t > 0]
    t<-t-1
    res <- new("dgCMatrix", i = as.integer(t), p = p0,  Dim = as.integer(c(n,n)), x = unlist(listw$weights))
    res<-t(res)
}
