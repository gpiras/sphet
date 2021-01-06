#' @name read.gwt2dist
#' @aliases read.gwt2dist
#' @title Read distance ojbects
#' @description  The function reads "GWT" files (i.e. generated using \link{distance}. 
#' It will read also other more general formats (as for example .txt files).
#' @usage read.gwt2dist(file, region.id = NULL, skip = 1)
#' @param file name of file to be read
#' @param region.id variable that defines the ordering of the observations
#' @param skip number of lines to skip
#' @details The first line of a \code{'GWT'} file generally contains some information 
#' (e.g. the name of the shape file, the number of observations), 
#' in which case, \code{skip} should be equal to 1. 
#' When the \code{file} has a \code{'GWT'} extension, 
#' the number of observations is generally retrived from the first line. 
#' Alternatively, it is fixed to the length 
#' of the \link{unique} \code{region.id} variable.
#' @return An object of class \code{distance}
#' @author Gianfranco Piras \email{gpiras@mac.com}
#' @examples
#' \dontrun{
#' library(spdep)
#' data(columbus)
#' dist <- read.gwt2dist(file = system.file('extdata/knn10columbus.GWT',
#' package = "sphet"), region.id = columbus$POLYID)
#' }
#' @keywords spatial
#' @export


read.gwt2dist<-function(file, region.id=NULL, skip=1){

if (skip==1){
	con <- file(file, open = "r")
    firstline <- unlist(strsplit(readLines(con, 1), " "))
	  if (length(firstline) == 4) {
        n <- as.integer(firstline[2])
        shpfile <- firstline[3]
        ind <- firstline[4]
        if (ind != deparse(substitute(region.id))) 
            warning(paste("region.id not named", ind))
    }
    else if (length(firstline) == 1) {
        n <- as.integer(firstline[1])
        shpfile <- as.character(NA)
        ind <- as.character(NA)
        warning("Old-style GWT file")
    }
    else stop("Invalid header line format for GWT file")
    
close(con)
}
else 	n<-length(unique(region.id))

    if (n < 1) 
        stop("non-positive number of entities")
    nseq <- 1:n
    if (is.null(region.id)) 
        region.id <- nseq
    if (n != length(region.id)) 
        stop("Mismatch in dimensions of GWT file and region.id")
    if (length(unique(region.id)) != length(region.id)) 
        stop("non-unique region.id given")
    odij <- read.table(file, skip = skip)
	 regodij <- match(odij[, 1], region.id)
    regddij <- match(odij[, 2], region.id)
    odij <- cbind(regodij, regddij, odij[, 3])
    qorder <- order(odij[, 1], odij[, 2])
    odij <- odij[qorder, ]
    origvec <- unique(odij[, 1])
    weights<-odij[,3]
	    if (!all(nseq %in% origvec)) 
        warning(paste(paste(region.id[which(!(nseq %in% origvec))], 
            collapse = ", "), "are not origins"))
    destvec <- unique(odij[, 2])
   if (!all(nseq %in% destvec)) 
        warning(paste(paste(region.id[which(!(nseq %in% destvec))], 
            collapse = ", "), "are not destinations"))
    res <- vector(mode = "list", length = n)
    weights <- vector(mode = "list", length = n)
    rle.sn <- rle(odij[, 1]) ####computes the lengths and values
    cs1.sn <- cumsum(rle.sn$lengths)
    cs0.sn <- c(1, cs1.sn[1:(n - 1)] + 1)
    ii <- 1
    for (i in 1:n) {
        if (!is.na(rle.sn$value[ii]) && rle.sn$value[ii] == i) {
            res[[i]] <- as.integer(odij[cs0.sn[ii]:cs1.sn[ii], 
                2])
    #            print(res[[i]])
            weights[[i]] <- as.double(odij[cs0.sn[ii]:cs1.sn[ii], 
                3])
            ii <- ii + 1
        }
        else {
            res[[i]] <- as.integer(0)
        }
    }
#    result<-list(neigh=res, weights=weights)
    result<-res
    class(result)<-c("sphet","distance","nb", "GWT") ##modified 03/11/2010
    
    attr(result, "region.id") <- region.id
    attr(result, "neighbours.attrs") <- as.character(NA)
    attr(result, "weights.attrs") <- as.character(NA)
    attr(result, "GeoDa") <- list(dist = weights, shpfile = shpfile, 
        ind = ind)
    attr(result, "call") <- match.call()
    attr(result, "n") <- n
    return(result)
	}