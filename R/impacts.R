#' @title Generate impacts for spreg lag and sarar models
#'
#' @param obj An object of class sphet
#' @param ... Additional arguments to be passed
#'
#'
#' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#'
#' @examples
#' data(columbus, package="spdep")
#' listw <- spdep::nb2listw(col.gal.nb)
#' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
#'             het = TRUE, verbose = FALSE, model = "sarar")
#' summary(res)
#' effects <- impacts(res, listw = listw, R = 399)
#' summary(effects)
#' @export
impacts <- function(obj, ...){
  UseMethod("impacts", obj)
}

#' @title Generate impacts for objects of class sarar_gmm created in sphet
#'
#' @param obj A spreg spatial regression object created by \code{spreg} with model ="sarar"
#' @param tr A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression
#' @param R If given, simulations are used to compute distributions for the impact measures, returned as \code{mcmc} objects
#' @param listw a listw object
#' @param evalues vector of eigenvalues of spatial weights matrix for impacts calculations
#' @param tol Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix
#' @param empirical Argument passed to \code{mvrnorm} (default FALSE): if true, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix
#' @param Q default NULL, else an integer number of cumulative power series impacts to calculate if \code{tr} is given
#' @param ... Arguments passed through to methods in the \pkg{coda} package 
#'
#'
#'
#' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#'
#' @examples
#' data(columbus, package="spdep")
#' listw <- spdep::nb2listw(col.gal.nb)
#' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
#'             het = TRUE, verbose = FALSE, model = "sarar")
#' summary(res)
#' effects <- impacts(res, listw = listw,  R = 399)
#' summary(effects)
#' @export
#' @method impacts gstsls
impacts.gstsls <- function(obj, ..., tr=NULL, R=NULL, listw=NULL, evalues=NULL,
                           tol=1e-6, empirical=FALSE, Q=NULL) {

  
  
  object <- obj
  if (!is.null(object$endo)) stop("impacts for model with additional endogenous variables not yet available in sphet")
  
  
  if(isTRUE(object$Durbin) | class(object$Durbin) == "formula"){
    
    type <- "mixed"
  
    if (is.null(object$interval)) interval <- c(-1,0.999)
    
    coefs <- drop(object$coefficients)
    p2 <- length(coefs)
    lambda <- coefs[p2-1]
    rho <- coefs[p2]
    beta <- coefs[1:(p2-2)]

    p <- length(beta)
    p1 <- p + 1
    names(beta) <- rownames(object$coefficients)[1:(p2-2)]
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0L
    n <- length(object$residuals)
    
    irho <- length(beta) + 1
    drop2beta <- (length(beta) + 1) : p2
    
    Sigma <- object$var  
    rownames(Sigma) <- colnames(Sigma) <- names(coefs)
    
    #print(Sigma)
    if(isTRUE(object$Durbin)){
      zero_fill <- NULL
      dvars <- NULL
      
      if (iicept) b1 <- beta[-icept]
      else  b1 <- beta
      p <- length(b1)
      if (p %% 2 != 0) stop("non-matched coefficient pairs")
      P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
      bnames <- names(b1[1:(p/2)])
      
      if(!is.null(R)){
        mu <- c(beta, lambda, rho)
        #Sigma <- object$var
        
      }

    }  
    
    else{
      
      dvars <- NULL
      
      dn <- grep("lag_", names(beta)) #which of the names of beta has "lag_"
      dc <- beta[dn] # betas that are lagged
      beta1 <- beta[-dn] # all the betas (and)including the lagged) 
      xb <- beta1[which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") )]
      xb <- xb[order(names(xb))]
      l_xb <- length(xb)
      if(length(which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") ))>=1) xo <- beta1[ -which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") )]
      else xo <- beta1
      
      
      gamma <- dc[which( stringr::str_remove(names(dc),"lag_")  %in% names(beta1))]
      gamma <- gamma[order(names(gamma))]
      l_gamma <- length(gamma)
      
      if(length(which(stringr::str_remove(names(dc),"lag_") %in% names(beta1)))>=1) don <- dc[-which( stringr::str_remove(names(dc),"lag_")  %in% names(beta1))]
      else don <- dc
      l_don <- length(don)
      
      
      
      if (iicept) {
        xo <- xo[-icept]
        l_xo <- length(xo)  
        if(l_xb != 0){
          b2 <- c(xo, xb, rep(0, l_don), rep(0, l_xo), gamma, don)
          bnames <- c(names(xo), sort(names(xb)), names(don)) 
        } 
        else{
          b2 <- c(xo, rep(0, l_don), rep(0, l_xo),  don)
          bnames <- c(names(xo), names(don)) 
        }
        P <- matrix(b2, ncol = 2)
      } 
      else {
        xo <- xo
        l_xo <- length(xo)
        if(l_xb != 0){
          b2 <- c(xo, xb, rep(0, l_don), rep(0, l_xo), gamma, don)
          bnames <- c(names(xo), sort(names(xb)), names(don)) 
        } 
        else{
          b2 <- c(xo, rep(0, l_don), rep(0, l_xo),  don)
          bnames <- c(names(xo), names(don)) 
        }
        P <- matrix(b2, ncol = 2)
      }
      
      
      if(!is.null(R)){
        
        #first build the zero fills
        if (l_xo != 0 && l_don != 0 && l_xb !=0) {
          zero_fill <- list(1: l_xo, (l_xo+1):(l_xo+ l_xb),  l_don, 
                            l_xo, (l_xo + l_xb+1):  (l_xo + l_xb+l_gamma),
                            (l_xo + l_xb+l_gamma+1):(l_xo + l_xb+l_gamma+l_don))
          l_zero_fill <- 6
        }
        
        if(l_xo != 0 && l_don == 0 && l_xb != 0) {
          zero_fill <- list(1: l_xo, (l_xo+1):(l_xo+ l_xb),  
                            l_xo, (l_xo + l_xb+1):  (l_xo + l_xb+l_gamma))
          
          l_zero_fill <- 5
        }
        if(l_xo == 0 && l_don != 0 && l_xb !=0)  {
          zero_fill <- list(1: l_xb,  l_don,(l_xb+1):  (l_xb+l_gamma), (l_xb+l_gamma+1):( l_xb+l_gamma+l_don))
          l_zero_fill <- 4
        }
        
        if(l_xo != 0 && l_don != 0 && l_xb ==0)  {
          zero_fill <- list(1: l_xo, l_don, 
                            l_xo, (l_xo+1):(l_xo +l_don))
          l_zero_fill <- 3
        }
        
        # if(l_xo == 0 && l_don != 0 && l_xb ==0){
        #   zero_fill <- list(l_don, 1: l_don)
        #   l_zero_fill <- 2   
        #   
        # }
        
       # print(class(zero_fill))
        attr(zero_fill, "l_zero_fill") <- l_zero_fill
        #  print(l_zero_fill)
       # print(Sigma)
        if(iicept) mu <- c(beta[1], xo, xb, gamma, don, lambda, rho)
        else mu <- c(xo, xb, gamma, don, lambda, rho)
        #print(names(mu))
        #print(rownames(Sigma))
        Sigma <- Sigma[match(names(mu), rownames(Sigma)),match(names(mu), colnames(Sigma)) ]
        #
        #print(Sigma)
      }
      

      
    }
  }
  
  else{
    
    type <- "lag" 
    if (is.null(object$interval)) interval <- c(-1,0.999)
    
    coefs <- drop(object$coefficients)
    p2 <- length(coefs)
    lambda <- coefs[p2-1]
    rho <- coefs[p2]
    beta <- coefs[1:(p2-2)]
    names(beta) <- rownames(object$coefficients)[1:(p2-2)]
    p <- length(beta)
    p1 <- p + 1
    
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0L
    zero_fill <- NULL
    dvars <- NULL
    
    if (iicept) {
      P <- matrix(beta[-icept], ncol=1)
      bnames <- names(beta[-icept])
    } else {
      P <- matrix(beta, ncol=1)
      bnames <- names(beta)
    }
    
    n <- length(object$residuals)
    mu <- NULL
    Sigma <- NULL
    
    if (!is.null(R)) {
      mu <- c( beta, lambda, rho)
      #print(mu)
      Sigma <- object$var
      #print(Sigma)
    }      
    
    irho <- length(beta) + 1
    drop2beta <- (length(beta) + 1):p2
    
    # if (!requireNamespace("spatialreg", quietly=TRUE))
    #   stop("install spatialreg")
    
    
  }
  #print(beta)
  res <- spatialreg::intImpacts(rho=lambda, beta=beta, P=P, n=n, mu=mu,
                                Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
                                interval=interval, type = type, tr=tr, R=R, listw=listw, evalues=evalues,
                                tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p, zero_fill = zero_fill)
  
  attr(res, "iClass") <- class(object)
  
  res
  
  }

 

#' @title Generate impacts for objects of class lag_gmm created in sphet
#'
#' @param obj A spreg spatial regression object created by \code{spreg} with model ="lag"
#' @param tr A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression
#' @param R If given, simulations are used to compute distributions for the impact measures, returned as \code{mcmc} objects
#' @param listw a listw object
#' @param evalues vector of eigenvalues of spatial weights matrix for impacts calculations
#' @param tol Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix
#' @param empirical Argument passed to \code{mvrnorm} (default FALSE): if true, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix
#' @param Q default NULL, else an integer number of cumulative power series impacts to calculate if \code{tr} is given
#' @param ... Arguments passed through to methods in the \pkg{coda} package 

#'
#' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#'
#' @examples
#' require("sf", quietly=TRUE)
#' library(coda)
#' columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
#' col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])
#' listw <- spdep::nb2listw(col.gal.nb)
#' ev <- spatialreg::eigenw(listw)
#' W <- as(listw, "CsparseMatrix")
#' trMatc <- spatialreg::trW(W, type="mult")
#' trMC <- spatialreg::trW(W, type="MC")
#' #LAG
#' lobj_gm <- spreg(CRIME ~ INC + HOVAL, columbus, listw,
#'                 model = "lag")
#' summary(lobj_gm)
#' lobj_gmh <- spreg(CRIME ~ INC + HOVAL, columbus, listw,
#'                  model = "lag", het = TRUE)
#' summary(lobj_gmh)
#' set.seed(1)
#' impacts(lobj_gm, listw=listw)
#' impacts(lobj_gm, tr=trMatc)
#' impacts(lobj_gm, tr=trMC)
#' impacts(lobj_gm, evalues=ev)
#' impacts(lobj_gmh, listw=listw)
#' impacts(lobj_gmh, tr=trMatc)
#' impacts(lobj_gmh, tr=trMC)
#' impacts(lobj_gmh, evalues=ev)
#' #same impacts but different SD
#' summary(impacts(lobj_gm, evalues = ev, R = 1000))
#' summary(impacts(lobj_gmh, evalues = ev, R = 1000))
#' lobjIQ5_gm <- impacts(lobj_gm, tr=trMatc, R=1000, Q=5)
#' summary(lobjIQ5_gm, zstats=TRUE, short=TRUE)
#' summary(lobjIQ5_gm, zstats=TRUE, short=TRUE, reportQ=TRUE)
#' # LAG durbin TRUE
#' mobj_gm <- spreg(CRIME ~ INC + HOVAL, columbus, listw, Durbin=TRUE,
#'                 model = "lag")
#' summary(mobj_gm)
#' mobj_gmh <- spreg(CRIME ~ INC + HOVAL, columbus, listw, Durbin=TRUE,
#'                  model = "lag", het = TRUE)
#' summary(mobj_gmh)
#' impacts(mobj_gm, listw=listw)
#' impacts(mobj_gm, tr=trMatc)
#' impacts(mobj_gm, tr=trMC)
#' impacts(mobj_gm, evalues=ev)
#' summary(impacts(mobj_gm, evalues=ev, R=1000), short=TRUE, zstats=TRUE)
#' impacts(mobj_gmh, listw=listw)
#' impacts(mobj_gmh, tr=trMatc)
#' impacts(mobj_gmh, tr=trMC)
#' impacts(mobj_gmh, evalues=ev)
#' summary(impacts(mobj_gmh, tr=trMatc, R=1000), short=TRUE, zstats=TRUE)
#' #lag durbin = ~formula
#' mobj1_gm <- spreg(CRIME ~ INC + HOVAL, columbus, listw, Durbin= ~ INC,
#'                  model = "lag")
#' mobj1_gmh <- spreg(CRIME ~ INC + HOVAL, columbus, listw, Durbin= ~ INC,
#'                   model = "lag", het = TRUE)
#' impacts(mobj1_gm, tr=trMatc)
#' impacts(mobj1_gm, listw=listw)
#' summary(impacts(mobj_gm, evalues=ev, R=200), short=TRUE, zstats=TRUE)
#' summary(impacts(mobj1_gm, tr=trMatc, R=200), short=TRUE, zstats=TRUE)
#' mobj1_gm <- spreg(CRIME ~ HOVAL, columbus, listw, Durbin= ~ INC,
#'                  model = "lag")
#' summary(impacts(mobj1_gm, evalues=ev, R=200), short=TRUE, zstats=TRUE)
#' @export
#' @method impacts stsls_sphet
impacts.stsls_sphet <- function(obj, ..., tr=NULL, R=NULL, listw=NULL,
                                evalues=NULL, tol=1e-6, empirical=FALSE, Q=NULL) {
  
  
  object <- obj
  if (!is.null(object$endo)) stop("impacts for model with additional endogenous variables not yet available in sphet")
 if(isTRUE(object$Durbin) | class(object$Durbin) == "formula"){
  
   type <- "mixed"
   if (is.null(object$interval)) interval <- c(-1,0.999)
   
   coefs <- drop(object$coefficients)
   p2 <- length(coefs)
   lambda <- coefs[p2]
   beta <- coefs[1:(p2-1)]
   p <- length(beta)
   p1 <- p + 1
   names(beta) <- rownames(object$coefficients)[1:(p2-1)]
   icept <- grep("(Intercept)", names(beta))
   iicept <- length(icept) > 0L
   n <- length(object$residuals)
   
   irho <- length(beta) + 1
   drop2beta <- length(beta) + 1
   

  if(isTRUE(object$Durbin)){
    zero_fill <- NULL
    dvars <- NULL
 
       if (iicept) b1 <- beta[-icept]
     else  b1 <- beta
    p <- length(b1)
    if (p %% 2 != 0) stop("non-matched coefficient pairs")
    P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
    bnames <- names(b1[1:(p/2)])
    
    if(!is.null(R)){
    mu <- c(beta, lambda)
    Sigma <- object$var
    }
    # res <- spatialreg::intImpacts(rho=lambda, beta=beta, P=P, n=n, mu=mu,
    #                               Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
    #                               interval=interval, type = type, tr=tr, R=R, listw=listw, evalues=evalues,
    #                               tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p, 
    #                               zero_fill = zero_fill)
    # 
    # 
  }  
   
   else{
  
     dvars <- NULL
     
   dn <- grep("lag_", names(beta)) #which of the names of beta has "lag_"
   dc <- beta[dn] # betas that are lagged
   beta1 <- beta[-dn] # all the betas (and)including the lagged) 
   xb <- beta1[which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") )]
   xb <- xb[order(names(xb))]
   l_xb <- length(xb)
   if(length(which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") ))>=1) xo <- beta1[ -which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") )]
   else xo <- beta1
   
   
   gamma <- dc[which( stringr::str_remove(names(dc),"lag_")  %in% names(beta1))]
   gamma <- gamma[order(names(gamma))]
   l_gamma <- length(gamma)
   
 if(length(which(stringr::str_remove(names(dc),"lag_") %in% names(beta1)))>=1) don <- dc[-which( stringr::str_remove(names(dc),"lag_")  %in% names(beta1))]
 else don <- dc
   l_don <- length(don)
   
   
   
  if (iicept) {
      xo <- xo[-icept]
      l_xo <- length(xo)  
     if(l_xb != 0){
       b2 <- c(xo, xb, rep(0, l_don), rep(0, l_xo), gamma, don)
       bnames <- c(names(xo), sort(names(xb)), names(don)) 
     } 
      else{
        b2 <- c(xo, rep(0, l_don), rep(0, l_xo),  don)
        bnames <- c(names(xo), names(don)) 
      }
      P <- matrix(b2, ncol = 2)
    } 
  else {
    xo <- xo
    l_xo <- length(xo)
    if(l_xb != 0){
      b2 <- c(xo, xb, rep(0, l_don), rep(0, l_xo), gamma, don)
      bnames <- c(names(xo), sort(names(xb)), names(don)) 
    } 
    else{
      b2 <- c(xo, rep(0, l_don), rep(0, l_xo),  don)
      bnames <- c(names(xo), names(don)) 
    }
    P <- matrix(b2, ncol = 2)
    }
  
  
   if(!is.null(R)){
     
 #first build the zero fills
    if (l_xo != 0 && l_don != 0 && l_xb !=0) {
      zero_fill <- list(1: l_xo, (l_xo+1):(l_xo+ l_xb),  l_don, 
                     l_xo, (l_xo + l_xb+1):  (l_xo + l_xb+l_gamma),
                     (l_xo + l_xb+l_gamma+1):(l_xo + l_xb+l_gamma+l_don))
     l_zero_fill <- 6
    }
     
    if(l_xo != 0 && l_don == 0 && l_xb != 0) {
      zero_fill <- list(1: l_xo, (l_xo+1):(l_xo+ l_xb),  
                                                  l_xo, (l_xo + l_xb+1):  (l_xo + l_xb+l_gamma))
    
      l_zero_fill <- 5
    }
    if(l_xo == 0 && l_don != 0 && l_xb !=0)  {
      zero_fill <- list(1: l_xb,  l_don,(l_xb+1):  (l_xb+l_gamma), (l_xb+l_gamma+1):( l_xb+l_gamma+l_don))
    l_zero_fill <- 4
    }
     
     if(l_xo != 0 && l_don != 0 && l_xb ==0)  {
       zero_fill <- list(1: l_xo, l_don, 
                         l_xo, (l_xo+1):(l_xo +l_don))
       l_zero_fill <- 3
     }
     
     # if(l_xo == 0 && l_don != 0 && l_xb ==0){
     #   zero_fill <- list(l_don, 1: l_don)
     #   l_zero_fill <- 2   
     #   
     # }
 
# print(class(zero_fill))
   attr(zero_fill, "l_zero_fill") <- l_zero_fill
#   print(l_zero_fill)

   if(iicept) mu <- c(beta[1], xo, xb, gamma, don, lambda)
   else mu <- c(xo, xb, gamma, don, lambda)
  
      Sigma <- object$var[match(names(mu), rownames(object$var)),match(names(mu), rownames(object$var)) ]
  
   #print(Sigma)
   }
   
   #  if (!requireNamespace("spatialreg", quietly=TRUE))
   #  stop("install spatialreg")
   # res <- spatialreg::intImpacts(rho=lambda, beta=beta, P=P, n=n, mu=mu,
   #                               Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
   #                               interval=interval, type = type, tr=tr, R=R, listw=listw, evalues=evalues,
   #                               tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p, 
   #                               zero_fill = zero_fill)
   

   }
 }

  else{

    type <- "lag" 
    if (is.null(object$interval)) interval <- c(-1,0.999)
    
    coefs <- drop(object$coefficients)
    p2 <- length(coefs)
    lambda <- coefs[p2]
    beta <- coefs[1:(p2-1)]
    names(beta) <- rownames(object$coefficients)[1:(p2-1)]
    p <- length(beta)
    p1 <- p + 1
  
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0L
    zero_fill <- NULL
    dvars <- NULL

    if (iicept) {
      P <- matrix(beta[-icept], ncol=1)
      bnames <- names(beta[-icept])
    } else {
      P <- matrix(beta, ncol=1)
      bnames <- names(beta)
    }
    
    n <- length(object$residuals)
    mu <- NULL
    Sigma <- NULL
    
    if (!is.null(R)) {
        mu <- c( beta, lambda)
        Sigma <- object$var
    }      

    irho <- length(beta) + 1
    drop2beta <- length(beta) + 1
    
    #  if (!requireNamespace("spatialreg", quietly=TRUE))
    #  stop("install spatialreg")
    # res <- spatialreg::intImpacts(rho=lambda, beta=beta, P=P, n=n, mu=mu,
    #                               Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
    #                               interval=interval, type = type, tr=tr, R=R, listw=listw, evalues=evalues,
    #                               tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p)
    # 
    
  }

  # if (!requireNamespace("spatialreg", quietly=TRUE))
  #   stop("install spatialreg")
  res <- spatialreg::intImpacts(rho=lambda, beta=beta, P=P, n=n, mu=mu,
                                Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
                                interval=interval, type = type, tr=tr, R=R, listw=listw, evalues=evalues,
                                tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p, 
                                zero_fill = zero_fill)
    
  attr(res, "iClass") <- class(object)
  
  res
}




#' @title Generate impacts for objects of class ols_sphet created in sphet
#'
#' @param obj A spreg spatial regression object created by \code{spreg} with model ="lag"
#' @param tr A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression
#' @param R If given, simulations are used to compute distributions for the impact measures, returned as \code{mcmc} objects
#' @param listw a listw object
#' @param evalues vector of eigenvalues of spatial weights matrix for impacts calculations
#' @param tol Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix
#' @param empirical Argument passed to \code{mvrnorm} (default FALSE): if true, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix
#' @param Q default NULL, else an integer number of cumulative power series impacts to calculate if \code{tr} is given
#' @param ... Arguments passed through to methods in the \pkg{coda} package 

#'
#' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#'
#' @examples
#' data(columbus, package="spdep")
#' listw <- spdep::nb2listw(col.gal.nb)
#' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
#'             het = TRUE, verbose = FALSE, model = "ols", Durbin = TRUE)
#' summary(res)
#' effects <- impacts(res)
#' summary(effects)
#' @export
#' @method impacts ols_sphet

impacts.ols_sphet <- function(obj, ..., tr=NULL, R=NULL, listw=NULL,
                              evalues=NULL, tol=1e-6, empirical=FALSE, Q=NULL){
  
  
  if(isFALSE(obj$Durbin) && !is.formula(obj$Durbin)) 
    stop("Trying to evaluate impacts in linear models without spatially autocorrelated regressors")
  if (!is.null(obj$endo)) stop("impacts for model with additional endogenous variables not yet available in sphet")
  
  
     beta <- drop(obj$coefficients)
     Sigma <- obj$var
    p <- length(beta)
   names(beta) <- rownames(obj$coefficients)
    #print(beta)
    
   icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0L
    n <- length(obj$residuals)
    
     
    
    if(isTRUE(obj$Durbin)){

      zero_fill <- NULL
      dvars <- NULL
      
      if (iicept) {
        b1 <- beta[-icept]
        Sigma <- Sigma[-icept, -icept]
      }
      else{
        b1 <- beta
        Sigma <- Sigma
      } 
      p <- length(b1)

      if (p %% 2 != 0) stop("non-matched coefficient pairs")
      
      
     st.err <- sqrt(diag(Sigma))
     st.err.tot <- c()
     for(i in 1:(p/2))  st.err.tot[i] <-  sqrt(Sigma[i,i] + Sigma[(i+p/2),(i+p/2)] + 2*Sigma[i,(i+p/2)])
     #print(st.err.tot)
     #print(st.err)
     #print(names(b1)[1:(p/2)])
     #print(c(b1[((p/2)+1):p], st.err[((p/2)+1):p]))
     #print((length(p)/2))
     indirImps <- matrix(c(b1[((p/2)+1):p], st.err[((p/2)+1):p]), nrow = p/2, ncol = 2, byrow = F)
     rownames(indirImps) <- names(b1)[1:(p/2)]
     
     dirImps <- matrix(c(b1[1:(p/2)], st.err[1:(p/2)]), nrow = p/2, ncol = 2, byrow = F )
     rownames(dirImps) <- names(b1)[1:(p/2)]
     
     totImps <- cbind(indirImps[,1] + dirImps[,1], st.err.tot)
     
     colnames(indirImps) <- colnames(dirImps) <- colnames(totImps) <- c("Estimate", "Std. Error")
     
     mixedImps <- list(dirImps=dirImps, indirImps=indirImps,
                       totImps=totImps)
     res <- list()
     res$n <- n
     res$k <- obj$k
     attr(res, "mixedImps") <- mixedImps
     class(res) <- c("ols_sphet", "SLX")
     
     
     
    }  
    else{
      ##not yet implementeds_
      dvars <- NULL
      
      dn <- grep("lag_", names(beta)) #which of the names of beta has "lag_"
      dc <- beta[dn] # betas that are lagged
      beta1 <- beta[-dn] # all the betas (and)including the lagged) 
      xb <- beta1[which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") )]
      xb <- xb[order(names(xb))]
      l_xb <- length(xb)
      if(length(which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") ))>=1) xo <- beta1[ -which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") )]
      else xo <- beta1
      #print(l_xb)
      
      gamma <- dc[which( stringr::str_remove(names(dc),"lag_")  %in% names(beta1))]
      gamma <- gamma[order(names(gamma))]
      l_gamma <- length(gamma)
      
      if(length(which(stringr::str_remove(names(dc),"lag_") %in% names(beta1)))>=1) don <- dc[-which( stringr::str_remove(names(dc),"lag_")  %in% names(beta1))]
      else don <- dc
      l_don <- length(don)
      
      
      
      if (iicept) {
        xo <- xo[-icept]
        l_xo <- length(xo)  
        b1 <- c(xo, xb, rep(NA, l_don), rep(NA, l_xo), gamma, don)
        b2 <- c(xo, xb, rep(0, l_don), rep(0, l_xo), gamma, don)
        if(l_xb !=0) bnames <- c(names(xo), sort(names(xb)), paste("lag_", sort(names(xb)),sep=""), names(don))
        else bnames <- c(names(xo), names(don))
        Sigma <- Sigma[-icept, -icept]
        Sigma <- Sigma[match(bnames, rownames(Sigma)), match(bnames, rownames(Sigma))]
        #print(bnames)
        #print(rownames(Sigma))
        #print(Sigma)
        Sigma.x <- Sigma[match(names(xo), rownames(Sigma)), match(names(xo), rownames(Sigma))]
        Sigma.b <- Sigma[match(c(sort(names(xb)), paste("lag_", sort(names(xb)),sep="")), rownames(Sigma)),
                         match(c(sort(names(xb)), paste("lag_", sort(names(xb)),sep="")), rownames(Sigma))]
        Sigma.o <- Sigma[match(names(don), rownames(Sigma)), match(names(don), rownames(Sigma))]
      } 
      else {
        xo <- xo
        l_xo <- length(xo)
        b1 <- c(xo, xb, rep(NA, l_don), rep(NA, l_xo), gamma, don)
        b2 <- c(xo, xb, rep(0, l_don), rep(0, l_xo), gamma, don)
        if(l_xb !=0) bnames <- c(names(xo), sort(names(xb)), paste("lag_", sort(names(xb)),sep=""), names(don))
        else bnames <- c(names(xo), names(don))
        Sigma <- Sigma[match(bnames, rownames(Sigma)), match(bnames, rownames(Sigma))]
        Sigma.x <- Sigma[match(names(xo), rownames(Sigma)), match(names(xo), rownames(Sigma))]
        Sigma.b <- Sigma[match(c(sort(names(xb)), paste("lag_", sort(names(xb)),sep="")), rownames(Sigma)),
                         match(c(sort(names(xb)), paste("lag_", sort(names(xb)),sep="")), rownames(Sigma))]
        Sigma.o <- Sigma[match(names(don), rownames(Sigma)), match(names(don), rownames(Sigma))]
      }
      #  print(Sigma)
      #  print(sqrt(diag(Sigma)))
      
      st.err <- sqrt(diag(Sigma))
      p <- length(b1)
      if (p %% 2 != 0) stop("non-matched coefficient pairs")
      
      p1 <- ncol(Sigma.b)
      
      if(is.matrix(Sigma.x)) st.err.x <- sqrt(diag(Sigma.x))
      if(is.numeric(Sigma.x)) st.err.x <- sqrt(Sigma.x)
      if(!is.null(p1)) for (i in 1:(p1/2)) st.err.b <-  sqrt(Sigma.b[i,i] + Sigma.b[(i+p1/2),(i+p1/2)] + 2*Sigma.b[i,(i+p1/2)])
      if(is.matrix(Sigma.o)) st.err.o <-  sqrt(diag(Sigma.o))
      if(is.numeric(Sigma.o)) st.err.o <- sqrt(Sigma.o)
      

      if(!is.null(p1)) st.err.T <- c(st.err.x, st.err.b, st.err.o)
      else st.err.T <- c(st.err.x,  st.err.o)
      
      
      if(!is.null(p1)) st.err <- c(st.err[1:l_xo], st.err[(1+l_xo):(l_xo+l_xb)], rep(NA, l_don), 
                                   rep(NA, l_xo),st.err[(1+l_xo+l_xb):(l_xo + l_xb +l_gamma)], 
                                   st.err[(l_xo + l_xb +l_gamma+1):(l_xo + l_xb +l_gamma+l_don)])
      else st.err <- c(st.err[1:l_xo], rep(NA, l_xo), 
                       rep(NA, l_don), st.err[(1+l_xo):(l_xo + l_don)])

      
      indirImps <- matrix(c(b1[((p/2)+1):p], st.err[((p/2)+1):p]), nrow = p/2, ncol = 2, byrow = F)
      rownames(indirImps) <- names(b1)[1:(p/2)]
      
      dirImps <- matrix(c(b1[1:(p/2)], st.err[1:(p/2)]), nrow = p/2, ncol = 2, byrow = F )
      rownames(dirImps) <- names(b1)[1:(p/2)]
      
      totImps <- cbind(b2[1:(p/2)] + b2[((p/2)+1):p], st.err.T)
      
      colnames(indirImps) <- colnames(dirImps) <- colnames(totImps) <- c("Estimate", "Std. Error")
      
      mixedImps <- list(dirImps=dirImps, indirImps=indirImps,
                        totImps=totImps)
      res <- list()
      res$n <- n
      res$k <- obj$k
      attr(res, "mixedImps") <- mixedImps
      class(res) <- c("ols_sphet", "SLX")
      
      
      
      
      
      
    }    
  
    spatialreg::impacts.SLX(res)
}





#' @title Generate impacts for objects of class error_sphet created in sphet
#'
#' @param obj A spreg spatial regression object created by \code{spreg} with model ="lag"
#' @param tr A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression
#' @param R If given, simulations are used to compute distributions for the impact measures, returned as \code{mcmc} objects
#' @param listw a listw object
#' @param evalues vector of eigenvalues of spatial weights matrix for impacts calculations
#' @param tol Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix
#' @param empirical Argument passed to \code{mvrnorm} (default FALSE): if true, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix
#' @param Q default NULL, else an integer number of cumulative power series impacts to calculate if \code{tr} is given
#' @param ... Arguments passed through to methods in the \pkg{coda} package 

#'
#' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#'
#' @examples
#' require("sf", quietly=TRUE)
#' library(coda)
#' columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
#' col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])
#' listw <- spdep::nb2listw(col.gal.nb)
#' mobj_gme <- spreg(CRIME ~ INC + HOVAL, columbus, listw, Durbin=TRUE,
#'                   model = "error")
#' summary(mobj_gme)
#' mobj_gmhe <- spreg(CRIME ~ INC + HOVAL, columbus, listw, Durbin=TRUE,
#'                   model = "error", het = TRUE)
#' summary(mobj_gmhe)
#' summary(impacts(mobj_gme))
#' summary(impacts(mobj_gmhe))
#' mobj1_gme <- spreg(CRIME ~ INC + HOVAL, columbus, listw, Durbin= ~ INC,
#'                   model = "error")
#' mobj1_gmhe <- spreg(CRIME ~ INC + HOVAL, columbus, listw, Durbin= ~ INC,
#'                    model = "error", het = TRUE)
#' impacts(mobj1_gme)
#' impacts(mobj1_gmhe)
#' mobj1_gme <- spreg(CRIME ~ HOVAL, columbus, listw, Durbin= ~ INC,
#'                    model = "error")
#' summary(impacts(mobj1_gme))
#' @export
#' @method impacts error_sphet

impacts.error_sphet <- function(obj, ..., tr=NULL, R=NULL, listw=NULL,
                              evalues=NULL, tol=1e-6, empirical=FALSE, Q=NULL){
  
  
  if(isFALSE(obj$Durbin) && !is.formula(obj$Durbin)) 
    stop("Trying to evaluate impacts in error models without spatially autocorrelated regressors")
  if (!is.null(obj$endo)) stop("impacts for model with additional endogenous variables not yet available in sphet")
  
  
  cf <- drop(obj$coefficients)
  names(cf) <- rownames(obj$coefficients)
  rho <- length(cf)
  beta <- cf[-rho]
  Sigma <- obj$var[-rho,-rho]
  colnames(Sigma) <- rownames(Sigma) <- names(cf)[-rho]
  p <- length(beta)
  names(beta) <- rownames(obj$coefficients)[-rho]
  #print(Sigma)
  
  icept <- grep("(Intercept)", names(beta))
  iicept <- length(icept) > 0L
  n <- length(obj$residuals)
  
  
  
  if(isTRUE(obj$Durbin)){
    
    zero_fill <- NULL
    dvars <- NULL
    
    if (iicept) {
      b1 <- beta[-icept]
      Sigma <- Sigma[-icept, -icept]
    }
    else{
      b1 <- beta
      Sigma <- Sigma
    } 
    p <- length(b1)
    
    if (p %% 2 != 0) stop("non-matched coefficient pairs")
    
    
    st.err <- sqrt(diag(Sigma))
    st.err.tot <- c()
    for(i in 1:(p/2))  st.err.tot[i] <-  sqrt(Sigma[i,i] + Sigma[(i+p/2),(i+p/2)] + 2*Sigma[i,(i+p/2)])

    indirImps <- matrix(c(b1[((p/2)+1):p], st.err[((p/2)+1):p]), nrow = p/2, ncol = 2, byrow = F)
    rownames(indirImps) <- names(b1)[1:(p/2)]
    
    dirImps <- matrix(c(b1[1:(p/2)], st.err[1:(p/2)]), nrow = p/2, ncol = 2, byrow = F )
    rownames(dirImps) <- names(b1)[1:(p/2)]
    
    totImps <- cbind(indirImps[,1] + dirImps[,1], st.err.tot)
    
    colnames(indirImps) <- colnames(dirImps) <- colnames(totImps) <- c("Estimate", "Std. Error")
    
    mixedImps <- list(dirImps=dirImps, indirImps=indirImps,
                      totImps=totImps)
    res <- list()
    res$n <- n
    res$k <- obj$k
    attr(res, "mixedImps") <- mixedImps
    class(res) <- c("ols_sphet", "SLX")
    
    
    
  }  
  else{
    ##not yet implementeds_
    dvars <- NULL
    
    dn <- grep("lag_", names(beta)) #which of the names of beta has "lag_"
    dc <- beta[dn] # betas that are lagged
    beta1 <- beta[-dn] # all the betas (and)including the lagged) 
    xb <- beta1[which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") )]
    xb <- xb[order(names(xb))]
    l_xb <- length(xb)
    if(length(which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") ))>=1) xo <- beta1[ -which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") )]
    else xo <- beta1
    #print(l_xb)
    
    gamma <- dc[which( stringr::str_remove(names(dc),"lag_")  %in% names(beta1))]
    gamma <- gamma[order(names(gamma))]
    l_gamma <- length(gamma)
    
    if(length(which(stringr::str_remove(names(dc),"lag_") %in% names(beta1)))>=1) don <- dc[-which( stringr::str_remove(names(dc),"lag_")  %in% names(beta1))]
    else don <- dc
    l_don <- length(don)
    
    
    
    if (iicept) {
      xo <- xo[-icept]
      l_xo <- length(xo)  
      b1 <- c(xo, xb, rep(NA, l_don), rep(NA, l_xo), gamma, don)
      b2 <- c(xo, xb, rep(0, l_don), rep(0, l_xo), gamma, don)
      if(l_xb !=0) bnames <- c(names(xo), sort(names(xb)), paste("lag_", sort(names(xb)),sep=""), names(don))
      else bnames <- c(names(xo), names(don))
      Sigma <- Sigma[-icept, -icept]
      Sigma <- Sigma[match(bnames, rownames(Sigma)), match(bnames, rownames(Sigma))]
      #print(bnames)
      #print(rownames(Sigma))
      #print(Sigma)
      Sigma.x <- Sigma[match(names(xo), rownames(Sigma)), match(names(xo), rownames(Sigma))]
      Sigma.b <- Sigma[match(c(sort(names(xb)), paste("lag_", sort(names(xb)),sep="")), rownames(Sigma)),
                       match(c(sort(names(xb)), paste("lag_", sort(names(xb)),sep="")), rownames(Sigma))]
      Sigma.o <- Sigma[match(names(don), rownames(Sigma)), match(names(don), rownames(Sigma))]
    } 
    else {
      xo <- xo
      l_xo <- length(xo)
      b1 <- c(xo, xb, rep(NA, l_don), rep(NA, l_xo), gamma, don)
      b2 <- c(xo, xb, rep(0, l_don), rep(0, l_xo), gamma, don)
      if(l_xb !=0) bnames <- c(names(xo), sort(names(xb)), paste("lag_", sort(names(xb)),sep=""), names(don))
      else bnames <- c(names(xo), names(don))
      Sigma <- Sigma[match(bnames, rownames(Sigma)), match(bnames, rownames(Sigma))]
      Sigma.x <- Sigma[match(names(xo), rownames(Sigma)), match(names(xo), rownames(Sigma))]
      Sigma.b <- Sigma[match(c(sort(names(xb)), paste("lag_", sort(names(xb)),sep="")), rownames(Sigma)),
                       match(c(sort(names(xb)), paste("lag_", sort(names(xb)),sep="")), rownames(Sigma))]
      Sigma.o <- Sigma[match(names(don), rownames(Sigma)), match(names(don), rownames(Sigma))]
    }
  #  print(Sigma)
  #  print(sqrt(diag(Sigma)))
    
    st.err <- sqrt(diag(Sigma))
    p <- length(b1)
    if (p %% 2 != 0) stop("non-matched coefficient pairs")
    
    p1 <- ncol(Sigma.b)
    
    if(is.matrix(Sigma.x)) st.err.x <- sqrt(diag(Sigma.x))
    if(is.numeric(Sigma.x)) st.err.x <- sqrt(Sigma.x)
    if(!is.null(p1)) for (i in 1:(p1/2)) st.err.b <-  sqrt(Sigma.b[i,i] + Sigma.b[(i+p1/2),(i+p1/2)] + 2*Sigma.b[i,(i+p1/2)])
    if(is.matrix(Sigma.o)) st.err.o <-  sqrt(diag(Sigma.o))
    if(is.numeric(Sigma.o)) st.err.o <- sqrt(Sigma.o)

    #print(st.err)    
#print(st.err.x)
#print(st.err.b)
#print(st.err.o)


    if(!is.null(p1)) st.err.T <- c(st.err.x, st.err.b, st.err.o)
    else st.err.T <- c(st.err.x,  st.err.o)

    
if(!is.null(p1)) st.err <- c(st.err[1:l_xo], st.err[(1+l_xo):(l_xo+l_xb)], rep(NA, l_don), 
                rep(NA, l_xo),st.err[(1+l_xo+l_xb):(l_xo + l_xb +l_gamma)], 
                st.err[(l_xo + l_xb +l_gamma+1):(l_xo + l_xb +l_gamma+l_don)])
    else st.err <- c(st.err[1:l_xo], rep(NA, l_xo), 
                     rep(NA, l_don), st.err[(1+l_xo):(l_xo + l_don)])
#print(st.err[1:l_xo])    
#print(rep(NA, l_xo))
#print(rep(NA, l_don))
#print(st.err[(1+l_xo):(l_xo + l_don)])    
#print((1+l_xo):(l_xo + l_don))
    indirImps <- matrix(c(b1[((p/2)+1):p], st.err[((p/2)+1):p]), nrow = p/2, ncol = 2, byrow = F)
    rownames(indirImps) <- names(b1)[1:(p/2)]
    
    dirImps <- matrix(c(b1[1:(p/2)], st.err[1:(p/2)]), nrow = p/2, ncol = 2, byrow = F )
    rownames(dirImps) <- names(b1)[1:(p/2)]
    
    totImps <- cbind(b2[1:(p/2)] + b2[((p/2)+1):p], st.err.T)
    
    colnames(indirImps) <- colnames(dirImps) <- colnames(totImps) <- c("Estimate", "Std. Error")
    
    mixedImps <- list(dirImps=dirImps, indirImps=indirImps,
                      totImps=totImps)
    res <- list()
    res$n <- n
    res$k <- obj$k
    attr(res, "mixedImps") <- mixedImps
    class(res) <- c("ols_sphet", "SLX")
    
    
    
    
  }    
  
  spatialreg::impacts.SLX(res)
}

is.formula <- function(x){
  inherits(x,"formula")
}
# processXSample_2 <- function(x, drop2beta, type, iicept, icept, n, listw,
#                            irho, zero_fill, dvars) {
#   
#   rho <- x[irho]
#   SW <- spdep::invIrW(listw, rho)
#   beta <- x[-drop2beta]
#   
#   
#   if(is.list(zero_fill)){
#   
#         if (iicept) {
#       b1 <- beta[-icept]
#       
#     } else {
#       b1 <- beta
#       
#     }
#     #print(b1)
#     #print(b1[zero_fill[[1]]])
#     #print(b1[zero_fill[[2]]])
#     #print(zero_fill[[3]])
#     #print(zero_fill[[4]])
#     #print(b1[zero_fill[[5]]])
#     #print(b1[zero_fill[[6]]])
#     # b1 <- c(xo, xb, rep(0, l_don), rep(0, l_xo), gamma, don)
#     # zero_fill <- list(l_xo, l_xb, l_gamma, l_don)
#     b1_l <- c(b1[zero_fill[[1]]], b1[zero_fill[[2]]], rep(0, zero_fill[[3]]), rep(0, zero_fill[[4]]), b1[zero_fill[[5]]], b1[zero_fill[[6]]])
#     # else{
#     #   if (length(rep(0, zero_fill[[3]])!=0)) b1_l <- c(b1[zero_fill[[1]]], b1[zero_fill[[2]]], rep(0, zero_fill[[4]]), b1[zero_fill[[5]]], b1[zero_fill[[6]]])
#     #   else b1_l <- c(b1[zero_fill[[1]]], b1[zero_fill[[2]]], rep(0, zero_fill[[3]]),  b1[zero_fill[[5]]], b1[zero_fill[[6]]])
#     # }
#     #print(b1_l)
#     p <- length(b1_l)
#     if (p %% 2 != 0) stop("non-matched coefficient pairs")
#     P <- cbind(b1_l[1:(p/2)], b1_l[((p/2)+1):p])
#     
#     #print(P)
#   }
#   
# else{  
# 
#     if (type == "lag" || type == "sac") {
#     if (iicept) {
#       P <- matrix(beta[-icept], ncol=1)
#     } else {
#       P <- matrix(beta, ncol=1)
#     }
#     # print(P)
#     return(spdep::lagImpactsExact(SW, P, n))
#   } else if (type == "mixed" || type == "sacmixed") {
#     if (iicept) {
#       b1 <- beta[-icept]
#     } else {
#       b1 <- beta
#     }
#     #FIXME
#     if (!is.null(zero_fill)) {
#       if (length(zero_fill) > 0L) {
#         inds <- attr(dvars, "inds")
#         b1_long <- rep(0, 2*(dvars[1]-1))
#         b1_long[1:(dvars[1]-1L)] <- b1[1:(dvars[1]-1)]
#         for (i in seq(along=inds)) {
#           b1_long[(dvars[1]-1L)+(inds[i]-1L)] <- b1[(dvars[1]-1L)+i]
#         }
#         b1 <- b1_long
#         #            s_zero_fill <- sort(zero_fill, decreasing=TRUE)
#         #            for (i in s_zero_fill) {
#         #              b1 <- append(b1, values=0, after=i-1L)
#         #            }
#       }
#     }
#     p <- length(b1)
#     if (p %% 2 != 0) stop("non-matched coefficient pairs")
#     P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
#   }
#     return(mixedImpactsExact_2(SW, P, n, listw))
#   }
# }
# 
# intImpacts_2 <- function(rho, beta, P, n, mu, Sigma, irho, drop2beta, bnames,
#                        interval, type, tr, R, listw, evalues, tol, empirical, Q, icept, iicept, p,
#                        mess=FALSE, samples=NULL, zero_fill=NULL, dvars=NULL) {
#   # print(rho)
#   # print(beta)
#   # print(P)
#   # print(n)
#   # print(mu)
#   # print(Sigma)
#   #print(irho)
#   #print(drop2beta)
#   # print(bnames)
#   # print(interval)
#   # print(type)
#   # print(evalues)
#   # print(icept)
#   # print(iicept)
#   #print(p)
#   if (is.null(evalues)) {
#     if (is.null(listw) && is.null(tr))
#       stop("either tr or listw must be given")
#   } else {
#     if (!is.null(listw)) {
#       warning("evalues given: listw will be ignored")
#       listw <-NULL
#     }
#     if (!is.null(tr)) {
#       warning("evalues given: listw will be ignored")
#       tr <- NULL
#     }
#   }
#   timings <- list()
#   .ptime_start <- proc.time()
#   if (is.null(listw)) {
#     
#     q <- length(tr)-1L
#     g <- rho^(0:q)
#     T <- matrix(c(1, tr[-(q+1)]/n), nrow=1)
#     if (type == "mixed" || type == "sacmixed") {
#       T <- rbind(T, tr/n)
#     }
#     if (is.null(evalues)) {
#       res <- spatialreg::lagImpacts(T, g, P)
#       cmethod <- "trace"
#     } 
#     else {
#       if (type == "mixed" || type == "sacmixed")
#         stop("eigenvalue mixed impacts not available")
#       if (length(evalues) != n) stop("wrong eigenvalue vector length")
#       res <- lagImpacts_e(rho, P, n, evalues)
#       cmethod <- "evalues"        }
#     if (!is.null(Q)) {
#       if (!is.numeric(Q) || length(Q) > 1L) stop("Invalid Q argument")
#       if (Q > length(tr)) stop("Q larger than length of tr")
#       Qres <- lagDistrImpacts(T, g, P, q=as.integer(Q))
#       attr(res, "Qres") <- Qres
#     }
#     timings[["trace_impacts"]] <- proc.time() - .ptime_start
#     .ptime_start <- proc.time()
#     if (!is.null(R)) {
#       if (is.null(samples)) {
#         samples <- mvrnorm(n=R, mu=mu, Sigma=Sigma, tol=tol,
#                            empirical=empirical)
#         if (mess) samples[,irho] <- 1 - exp(samples[,irho])
#       }
#       if (!is.null(interval)) {
#         check <- ((samples[,irho] > interval[1]) & 
#                     (samples[,irho] < interval[2]))
#         if (any(!check)) samples <- samples[check,]
#       }
#       timings[["impacts_samples"]] <- proc.time() - .ptime_start
#       .ptime_start <- proc.time()
#       # type, iicept, icept, T, Q
#       sres <- apply(samples, 1, spreg::processSample, irho=irho,
#                     drop2beta=drop2beta, type=type, iicept=iicept,
#                     icept=icept, zero_fill=zero_fill, dvars=dvars, T=T, Q=Q, q=q,
#                     evalues=evalues)
#       timings[["process_samples"]] <- proc.time() - .ptime_start
#       .ptime_start <- proc.time()
#       # 100928 Eelke Folmer
#       if (length(bnames) == 1L) {
#         direct <- as.mcmc(t(matrix(sapply(sres, function(x) x$direct),
#                                    nrow=1)))
#         indirect <- as.mcmc(t(matrix(sapply(sres,
#                                             function(x) x$indirect), nrow=1)))
#         total <- as.mcmc(t(matrix(sapply(sres, function(x) x$total),
#                                   nrow=1)))
#       } else {
#         direct <- as.mcmc(t(sapply(sres, function(x) x$direct)))
#         indirect <- as.mcmc(t(sapply(sres, function(x) x$indirect)))
#         total <- as.mcmc(t(sapply(sres, function(x) x$total)))
#       }
#       colnames(direct) <- bnames
#       colnames(indirect) <- bnames
#       colnames(total) <- bnames
#       ssres <- list(direct=direct, indirect=indirect, total=total)
#       if (!is.null(Q)) {
#         Qdirect <- as.mcmc(t(sapply(sres, function(x)
#           attr(x, "Qres")$direct)))
#         Qindirect <- as.mcmc(t(sapply(sres, function(x) 
#           attr(x, "Qres")$indirect)))
#         Qtotal <- as.mcmc(t(sapply(sres, function(x) 
#           attr(x, "Qres")$total)))
#         Qnames <- c(sapply(bnames, function(x) 
#           paste(x, 1:Q, sep="__Q")))
#         if (length(Qnames) == 1L) {
#           Qdirect <- t(Qdirect)
#           Qindirect <- t(Qindirect)
#           Qtotal <- t(Qtotal)
#         }
#         colnames(Qdirect) <- Qnames
#         colnames(Qindirect) <- Qnames
#         colnames(Qtotal) <- Qnames
#         Qmcmc <- list(direct=Qdirect, indirect=Qindirect, total=Qtotal)
#         attr(ssres, "Qmcmc") <- Qmcmc
#       }
#       timings[["postprocess_samples"]] <- proc.time() - .ptime_start
#       res <- list(res=res, sres=ssres)
#     }
#     attr(res, "method") <- cmethod
#   } 
#   else {
#     # added checks 140304
#     stopifnot(length(listw$neighbours) == n)
#     V <- spdep::listw2mat(listw)
#     e <- eigen(V, only.values = TRUE)$values
#     if (is.complex(e)) interval <- 1/(range(Re(e)))
#     else interval <- 1/(range(e))
#     SW <- spdep::invIrW(listw, rho)
#     if (type == "lag" || type == "sac") res <- lagImpactsExact(SW, P, n)
#     else if (type == "mixed" || type == "sacmixed")
#       res <- mixedImpactsExact_2(SW, P, n, listw)
#     timings[["weights_impacts"]] <- proc.time() - .ptime_start
#     .ptime_start <- proc.time()
#     if (!is.null(R)) {
#       #print(Sigma)
#       samples <- mvrnorm(n=R, mu=mu, Sigma=Sigma, tol=tol,
#                          empirical=empirical)
#       check <- ((samples[,irho] > interval[1]) & 
#                   (samples[,irho] < interval[2]))
#       if (any(!check)) samples <- samples[check,]
#       #print(samples)
#       timings[["impacts_samples"]] <- proc.time() - .ptime_start
#       .ptime_start <- proc.time()
#       # type, iicept, icept, SW, n, listw
#       if(is.list(zero_fill)) sres <- apply(samples, 1, processXSample_2,
#                                              drop2beta=drop2beta, type=type, iicept=iicept,
#                                              icept=icept, n=n, listw=listw, irho=irho, zero_fill=zero_fill,
#                                              dvars=dvars)
#       
#       else sres <- apply(samples, 1, processXSample,
#                     drop2beta=drop2beta, type=type, iicept=iicept,
#                     icept=icept, n=n, listw=listw, irho=irho, zero_fill=zero_fill,
#                     dvars=dvars)
#       
#       timings[["process_samples"]] <- proc.time() - .ptime_start
#       .ptime_start <- proc.time()
#       if (length(bnames) == 1L) {
#         direct <- as.mcmc(t(matrix(sapply(sres, function(x) x$direct),
#                                    nrow=1)))
#         indirect <- as.mcmc(t(matrix(sapply(sres,
#                                             function(x) x$indirect), nrow=1)))
#         total <- as.mcmc(t(matrix(sapply(sres, function(x) x$total),
#                                   nrow=1)))
#       } else {
#         #print(sres)
#         direct <- as.mcmc(t(sapply(sres, function(x) x$direct)))
#         indirect <- as.mcmc(t(sapply(sres, function(x) x$indirect)))
#         total <- as.mcmc(t(sapply(sres, function(x) x$total)))
#       }
#       
#       colnames(direct) <- bnames
#       colnames(indirect) <- bnames
#       colnames(total) <- bnames
#       timings[["postprocess_samples"]] <- proc.time() - .ptime_start
#       res <- list(res=res, sres=list(direct=direct,
#                                      indirect=indirect, total=total))
#     }
#     attr(res, "method") <- "exact"
#   }
#   if (!is.null(R)) attr(res, "samples") <- list(samples=samples, irho=irho,
#                                                 drop2beta=drop2beta)
#   attr(res, "type") <- type
#   attr(res, "bnames") <- bnames
#   attr(res, "haveQ") <- !is.null(Q)
#   attr(res, "timings") <- do.call("rbind", timings)[, c(1,3)]
#   class(res) <- "lagImpact"
#   res
# }



