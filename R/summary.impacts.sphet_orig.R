#' #' Generate impacts for spreg lag and sarar models
#' #' 
#' #' 
#' #' @param object An object of class sphet
#' #' @param ... Additional arguments to be passed 
#' #' 
#' #' 
#' #' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#' #' 
#' #' @examples 
#' #' data(columbus, package="spdep")
#' #' listw <- spdep::nb2listw(col.gal.nb)
#' #' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
#' #'             het = TRUE, verbose = FALSE, model = "sarar")
#' #' summary(res)
#' #' effects <- impacts(res, n_mvn = 3000)
#' #' summary(effects)
#' impacts <- function(object, ...){
#'                       UseMethod("impacts", object)
#' }
#' #' describedIn impacts
#' #' @param object An object of class sphet
#' #' @param ... Additional arguments to be passed 
#' #' 
#' #' 
#' #' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#' #' 
#' #' @examples 
#' #' data(columbus, package="spdep")
#' #' listw <- spdep::nb2listw(col.gal.nb)
#' #' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
#' #'             het = TRUE, verbose = FALSE, model = "sarar")
#' #' summary(res)
#' #' effects <- impacts(res, n_mvn = 3000)
#' #' summary(effects)
#' impacts.sphet <- function(object, ...) {
#'    
#'   if(class(object$Durbin) == "formula") object$Durbin <- TRUE
#'  
#'   if(!is.null(object$endog)) stop("Impacts not yet implemented 
#'                                for models with additional endogenous variables")
#'  
#'   if(class(object)[[2]] %in% c("ols", "error_gmm")) stop("Impacts not needed
#'                                                  for error and ols models")
#'   
#'   
#'   if(!isTRUE(all.equal(as.numeric(Matrix::rowSums(object$listw)), rep(1,dim(object$listw)[2])))) stop("Impacts not yet implemented 
#'                                           for non row-standardized listw")
#'     
#'   ag <- class(object)[[2]]
#'  # print(ag)
#'   switch(match.arg(ag, c("lag_gmm", "sarar_gmm")),
#'          sarar_gmm = impacts.sarar(object, ...),
#'          lag_gmm = impacts.lag(object, ...)
#'          )
#' 
#' }
#'   # TODO in the future implement the VC matrix of the impacts 
#'   # TODO check two things in the impacts: 1) calculate the impacts only for the coefficients in parameter space 2) fix the Z test
#' 
#' #' Generate impacts
#' #' 
#' #' 
#' #' @param object An object of class sphet
#' #' @param mvn Inference based on MC simulation of mvn samples
#' #' @param n_mvn Number of samples
#' #' @param weg Eigenvalues can be provided by the user
#' #' @param inference should inference be performed
#' #' @param exact should inference be exact
#' #' @param trW either NULL or a trW from spatialreg
#' #' @param m number of traces of w - default is 30 
#' #' 
#' #' 
#' #' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#' #' 
#' #' @examples 
#' #' data(columbus, package="spdep")
#' #' listw <- spdep::nb2listw(col.gal.nb)
#' #' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
#' #'             het = TRUE, verbose = FALSE, model = "sarar")
#' #' summary(res)
#' #' effects <- impacts(res, n_mvn = 3000)
#' #' summary(effects)
#' impacts.sarar <-function(object, n_mvn = 3000, weg = NULL, 
#'                           inference = FALSE, exact = TRUE, m = 30, trW = NULL, mvn = TRUE){ 
#'    
#'    if(object$Durbin) {
#'       
#'       coefs <- drop(object$coefficients)
#'       p2 <- length(coefs)
#'       lambda <- coefs[(p2-1)]
#'       beta <- coefs[1:(p2-2)]
#'       rho <- coefs[p2]
#'       
#'       
#'       dn <- grep("lag_", names(beta))
#'       dc <- beta[dn]
#'       beta <- beta[-dn]
#'       
#'       xb <- beta[which(names(beta) %in% stringr::str_remove(names(dc),"lag_") )]
#'       xb <- xb[order(names(xb))]
#'       xo <- beta[ -which(names(beta) %in% stringr::str_remove(names(dc),"lag_") )]
#'       gamma <- dc[which( stringr::str_remove(names(dc),"lag_")  %in% names(beta))]
#'       gamma <- gamma[order(names(gamma))]
#'       don <- dc[-which( stringr::str_remove(names(dc),"lag_")  %in% names(beta))]
#'       
#'       icept <- grep("(Intercept)", names(xo))
#'       iicept <- length(icept) > 0
#'       Sigma <- object$var  
#' 
#'       if (iicept) {
#'          P <- xo[-icept]
#'          Sigma <- Sigma[-c(icept,p2),-c(icept,p2)]
#'          bnames <- c(names(xo[-icept]), stringr::str_remove(names(dc), "lag_"))
#'       } 
#'       else {
#'          P <- xo
#'          Sigma <- Sigma[-c(p2),-c(p2)]
#'          bnames <-c(names(xo), stringr::str_remove(names(dc), "lag_"))
#'       }
#'       
#'       n <- length(object$residuals)
#'       
#'       
#'       ATE <- c(AVTE_rn_Dur(c(P, lambda), c(0, lambda)), 
#'                AVTE_rn_Dur(c(xb, lambda), c(gamma, lambda)), 
#'                AVTE_rn_Dur(c(0, lambda), c(don, lambda)) )
#'       
#'       
#'       
#'       if(exact){
#'          
#'          if(is.null(weg)) weg <- eigen(object$listw)$values
#'  
#'       
#'          ADE <- c(ADRE_rn_Dur(c(P, lambda), weg, n, m, object$listw, c(0, lambda), exact = TRUE), 
#'                       ADRE_rn_Dur(c(xb, lambda), weg, n, m, object$listw, c(gamma, lambda), exact = TRUE), 
#'                       ADRE_rn_Dur(c(0, lambda), weg, n, m, object$listw, c(don, lambda), exact = TRUE) )
#'       }
#'       else{
#'          
#'          if(is.null(trW))   trW <-  smtrace(object$listw, m = m)
#' 
#'       
#'          ADE <- c(ADRE_rn_Dur(c(P, lambda), weg, n, m = m, object$listw, c(0, lambda), exact = F, trW = trW), 
#'                   ADRE_rn_Dur(c(xb, lambda), weg, n, m = m, object$listw, c(gamma, lambda), exact = F, trW = trW), 
#'                   ADRE_rn_Dur(c(0, lambda), weg, n, m = m, object$listw, c(don, lambda), exact = F, trW = trW) )
#' 
#'          
#'       }
#'       
#'       AIE <- ATE - ADE
#'       
#'       
#'       if (inference) { 
#'          
#'          mu <- c(P, xb, gamma, don, lambda)
#'          mc_splm <- mvtnorm::rmvnorm(n_mvn, mean = mu, sigma = Sigma)
#'          
#'          
#'          if (length(P) == 0 && length(don) == 0)       ATE_inf <- rbind(apply(mc_splm[ ,c((1 : length(xb)), ncol(mc_splm))], 1, AVTE_rn_Dur_fp) + 
#'                                                                            apply( mc_splm[ , c((length(xb) + 1) : ( ncol(mc_splm)-1), ncol(mc_splm))], 1, AVTE_rn_Dur_sp))
#'          
#'          else if(length(P) == 0 && length(don) != 0) ATE_inf <- rbind(apply(mc_splm[ ,c( (1 : length(xb)), ncol(mc_splm))], 1, AVTE_rn_Dur_fp) + 
#'                                                                          apply(mc_splm[ , c((length(xb) + 1) : (ncol(mc_splm) - 1), ncol(mc_splm))], 1, AVTE_rn_Dur_sp), 
#'                                                                       apply(mc_splm[ , c(length(xb) + length(gamma) +1) : (ncol(mc_splm) - 1), ncol(mc_splm)], 1, AVTE_rn_Dur_sp))
#'          
#'          else if(length(P) != 0 && length(don) == 0) ATE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, AVTE_rn_Dur_fp),
#'                                                                       apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, AVTE_rn_Dur_fp) + 
#'                                                                          apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( ncol(mc_splm) - 1), ncol(mc_splm))], 1, AVTE_rn_Dur_sp))  
#'          
#'          else ATE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, AVTE_rn_Dur_fp) 
#'                                + apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, AVTE_rn_Dur_sp),
#'                                apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, AVTE_rn_Dur_fp) + 
#'                                   apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( (length(P) + length(xb) + length(gamma))), ncol(mc_splm))], 1, AVTE_rn_Dur_sp),
#'                                apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, AVTE_rn_Dur_fp) + 
#'                                   apply( mc_splm[ , c( (length(P) +  length(xb) +length(gamma) +1) : (length(P) + length(xb) +length(gamma) +length(don)), ncol(mc_splm))], 1, AVTE_rn_Dur_sp))
#'          
#'          
#'          if(exact){ 
#'             
#'             
#'             if (length(P) == 0 && length(don) == 0)       ADE_inf <- rbind(apply(mc_splm[ ,c((1 : length(xb)), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE ) + 
#'                                                                               apply( mc_splm[ , c((length(xb) + 1) : ( ncol(mc_splm)-1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE))
#'             
#'             else if(length(P) == 0 && length(don) != 0) ADE_inf <- rbind(apply(mc_splm[ ,c( (1 : length(xb)), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE) + 
#'                                                                             apply(mc_splm[ , c((length(xb) + 1) : (ncol(mc_splm) - 1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE), 
#'                                                                          apply(mc_splm[ , c(length(xb) + length(gamma) +1) : (ncol(mc_splm) - 1), ncol(mc_splm)], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE))
#'             
#'             else if(length(P) != 0 && length(don) == 0) ADE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE),
#'                                                                          apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE) + 
#'                                                                             apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( ncol(mc_splm) - 1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE))  
#'             
#'             else ADE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE) 
#'                                   + apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE),
#'                                   apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE) + 
#'                                      apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( (length(P) + length(xb) + length(gamma))), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE),
#'                                   apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE) + 
#'                                      apply( mc_splm[ , c( (length(P) +  length(xb) +length(gamma) +1) : (length(P) + length(xb) +length(gamma) +length(don)), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE))
#'             
#'             
#'          }
#'          else{
#'             
#'             
#'             if (length(P) == 0 && length(don) == 0)       ADE_inf <- rbind(apply(mc_splm[ ,c((1 : length(xb)), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW ) + 
#'                                                                               apply( mc_splm[ , c((length(xb) + 1) : ( ncol(mc_splm)-1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW))
#'             
#'             else if(length(P) == 0 && length(don) != 0) ADE_inf <- rbind(apply(mc_splm[ ,c( (1 : length(xb)), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW) + 
#'                                                                             apply(mc_splm[ , c((length(xb) + 1) : (ncol(mc_splm) - 1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW), 
#'                                                                          apply(mc_splm[ , c(length(xb) + length(gamma) +1) : (ncol(mc_splm) - 1), ncol(mc_splm)], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW))
#'             
#'             else if(length(P) != 0 && length(don) == 0) ADE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW),
#'                                                                          apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW) + 
#'                                                                             apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( ncol(mc_splm) - 1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW))  
#'             
#'             else ADE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW) 
#'                                   + apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW),
#'                                   apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW) + 
#'                                      apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( (length(P) + length(xb) + length(gamma))), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW),
#'                                   apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW) + 
#'                                      apply( mc_splm[ , c( (length(P) +  length(xb) +length(gamma) +1) : (length(P) + length(xb) +length(gamma) +length(don)), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW))
#'             
#'             
#'             
#'          }
#'          
#'          AIE_inf <- ATE_inf - ADE_inf
#'          
#'          res <- list(ATE, ADE, AIE, ATE_inf, ADE_inf, AIE_inf) 
#'       }
#'       else res <- list(ATE, ADE, AIE)
#'       
#'       
#'    }
#'     
#'   else{
#'     
#'   coefs <- drop(object$coefficients)
#'   p2 <- length(coefs)
#'   lambda <- coefs[(p2-1)]
#'   beta <- coefs[1:(p2-2)]
#'   p <- length(beta)
#'   p1 <- p + 1
#'   rho <- coefs[p2]
#'   icept <- grep("(Intercept)", names(beta))
#'   iicept <- length(icept) > 0
#'   Sigma <- object$var  
#'   n <- length(object$residuals)
#'   
#'   if (iicept) {
#'     P <- beta[-icept]
#'     Sigma <- Sigma[-c(icept,p2),-c(icept,p2)]
#'     bnames <- names(beta[-icept])
#'   } 
#'   else {
#'     P <- beta
#'     Sigma <- Sigma[-c(p2),-c(p2)]
#'     bnames <- names(beta)
#'   }
#'   
#'   ATE <- AVTE_rn(c(P, lambda))
#'   
#'    if(exact){
#'     
#'       if(is.null(weg)) weg <- eigen(object$listw)$values
#'   ADE <- ADRE_rn(c(P, lambda), weg, n, exact = TRUE)
#'  }
#'   else{
#'     
#'   if(is.null(trW))   trW <-  smtrace(object$listw, m = m)
#' 
#'      ADE <- ADRE_rn(c(P, lambda), exact = FALSE, trW = trW, m = m, n = n)
#'   
#'      }
#'   AIE <- ATE - ADE
#'   
#'    if (inference) { 
#' 
#'         mu <- c(P, lambda)
#'   mc_splm <- mvtnorm::rmvnorm(n_mvn, mean = mu, sigma = Sigma)
#' 
#'   
#'   ATE_inf <- apply(mc_splm, 1, AVTE_rn)
#'   
#'   if(exact) ADE_inf <- apply(mc_splm, 1, ADRE_rn, weg, n, exact = TRUE) 
#'   else     ADE_inf <- apply(mc_splm, 1, ADRE_rn, exact = FALSE, trW = trW, m = m, n = n)
#'      
#'   AIE_inf <- ATE_inf - ADE_inf
#'      
#'   res <- list(ATE, ADE, AIE, ATE_inf, ADE_inf, AIE_inf) 
#'   }
#'   else res <- list(ATE, ADE, AIE)
#'   
#'   }
#'   
#'    attr(res, "inference") <- inference
#'    class(res) <- c("impacts_sphet") 
#' 
#'   invisible(res)
#'  }
#'  
#' #' Generate impacts
#' #' 
#' #' 
#' #' @param object An object of class sphet
#' #' @param mvn Inference based on MC simulation of mvn samples
#' #' @param n_mvn Number of samples
#' #' @param weg Eigenvalues can be provided by the user
#' #' @param inference should inference be performed
#' #' @param exact should inference be exact
#' #' @param trW either NULL or a trW from spatialreg
#' #' @param m number of traces of w - default is 30 
#' #' 
#' #' 
#' #' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#' #' 
#' #' @examples 
#' #' data(columbus, package="spdep")
#' #' listw <- spdep::nb2listw(col.gal.nb)
#' #' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
#' #'             het = TRUE, verbose = FALSE, model = "sarar")
#' #' summary(res)
#' #' effects <- sphet::impacts(res, n_mvn = 3000)
#' #' summary(effects)
#' impacts.lag <-function(object, n_mvn = 3000, weg = NULL, 
#'                        inference = FALSE, exact = TRUE, m = 30, trW = NULL, mvn = TRUE){ 
#'    
#'    if(object$Durbin){ 
#'       
#'       coefs <- drop(object$coefficients)
#'       p2 <- length(coefs)
#'       lambda <- coefs[p2]
#'       beta <- coefs[1:(p2-1)]
#'       p <- length(beta)
#'       p1 <- p + 1
#'       
#'       
#'       dn <- grep("lag_", names(beta))
#'       dc <- beta[dn]
#'       beta <- beta[-dn]
#'       
#'       xb <- beta[which(names(beta) %in% stringr::str_remove(names(dc),"lag_") )]
#'       xb <- xb[order(names(xb))]
#'       xo <- beta[ -which(names(beta) %in% stringr::str_remove(names(dc),"lag_") )]
#'       gamma <- dc[which( stringr::str_remove(names(dc),"lag_")  %in% names(beta))]
#'       gamma <- gamma[order(names(gamma))]
#'       don <- dc[-which( stringr::str_remove(names(dc),"lag_")  %in% names(beta))]
#'       
#'       icept <- grep("(Intercept)", names(xo))
#'       iicept <- length(icept) > 0
#'       Sigma <- as.matrix(object$var)
#'       #print(Sigma)
#'       
#'       if (iicept) {
#'          P <- xo[-icept]
#'          Sigma <- Sigma[-c(icept),-c(icept)]
#'          bnames <- c(names(xo[-icept]), stringr::str_remove(names(dc), "lag_"), names(don))
#'       } 
#'       else {
#'          P <- xo
#'          Sigma <- Sigma
#'          bnames <-c(names(xo), stringr::str_remove(names(dc), "lag_"), names(don))
#'       }
#'       
#'       n <- length(object$residuals)
#'       
#'       
#'       ATE <- c(AVTE_rn_Dur(c(P, lambda), c(0, lambda)), 
#'                AVTE_rn_Dur(c(xb, lambda), c(gamma, lambda)), 
#'                AVTE_rn_Dur(c(0, lambda), c(don, lambda)) )
#'       
#'       
#'       
#'       if(exact){
#'          
#'          if(is.null(weg)) weg <- eigen(object$listw)$values
#'          
#'          
#'          ADE <- c(ADRE_rn_Dur(c(P, lambda), weg, n, m, object$listw, c(0, lambda), exact = TRUE), 
#'                   ADRE_rn_Dur(c(xb, lambda), weg, n, m, object$listw, c(gamma, lambda), exact = TRUE), 
#'                   ADRE_rn_Dur(c(0, lambda), weg, n, m, object$listw, c(don, lambda), exact = TRUE) )
#'       }
#'       else{
#'          
#'          if(is.null(trW))   trW <-  smtrace(object$listw, m = m)
#'          
#'          
#'          ADE <- c(ADRE_rn_Dur(c(P, lambda), weg, n, m = m, object$listw, c(0, lambda), exact = F, trW = trW), 
#'                   ADRE_rn_Dur(c(xb, lambda), weg, n, m = m, object$listw, c(gamma, lambda), exact = F, trW = trW), 
#'                   ADRE_rn_Dur(c(0, lambda), weg, n, m = m, object$listw, c(don, lambda), exact = F, trW = trW) )
#'          
#'       }
#'       
#'       AIE <- ATE - ADE
#'       
#'       
#'       if (inference) { 
#'          
#'          mu <- c(P, xb, gamma, don, lambda)
#'          mc_splm <- mvtnorm::rmvnorm(n_mvn, mean = mu, sigma = Sigma)
#'          
#' 
#' if (length(P) == 0 && length(don) == 0)       ATE_inf <- rbind(apply(mc_splm[ ,c((1 : length(xb)), ncol(mc_splm))], 1, AVTE_rn_Dur_fp) + 
#'                                                                   apply( mc_splm[ , c((length(xb) + 1) : ( ncol(mc_splm)-1), ncol(mc_splm))], 1, AVTE_rn_Dur_sp))
#'          
#'          else if(length(P) == 0 && length(don) != 0) ATE_inf <- rbind(apply(mc_splm[ ,c( (1 : length(xb)), ncol(mc_splm))], 1, AVTE_rn_Dur_fp) + 
#'                                                                     apply(mc_splm[ , c((length(xb) + 1) : (ncol(mc_splm) - 1), ncol(mc_splm))], 1, AVTE_rn_Dur_sp), 
#'                                                                     apply(mc_splm[ , c(length(xb) + length(gamma) +1) : (ncol(mc_splm) - 1), ncol(mc_splm)], 1, AVTE_rn_Dur_sp))
#'          
#'          else if(length(P) != 0 && length(don) == 0) ATE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, AVTE_rn_Dur_fp),
#'                                                                      apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, AVTE_rn_Dur_fp) + 
#'                                                                         apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( ncol(mc_splm) - 1), ncol(mc_splm))], 1, AVTE_rn_Dur_sp))  
#'    
#'    else ATE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, AVTE_rn_Dur_fp) 
#'                                            + apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, AVTE_rn_Dur_sp),
#'                                            apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, AVTE_rn_Dur_fp) + 
#'                                              apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( (length(P) + length(xb) + length(gamma))), ncol(mc_splm))], 1, AVTE_rn_Dur_sp),
#'                                              apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, AVTE_rn_Dur_fp) + 
#'                                              apply( mc_splm[ , c( (length(P) +  length(xb) +length(gamma) +1) : (length(P) + length(xb) +length(gamma) +length(don)), ncol(mc_splm))], 1, AVTE_rn_Dur_sp))
#'                          
#' 
#' #         print(head(ATE_inf)[,1:5])
#'          
#'          if(exact){ 
#'             
#'             
#'             if (length(P) == 0 && length(don) == 0)       ADE_inf <- rbind(apply(mc_splm[ ,c((1 : length(xb)), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE ) + 
#'                                                                               apply( mc_splm[ , c((length(xb) + 1) : ( ncol(mc_splm)-1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE))
#'             
#'             else if(length(P) == 0 && length(don) != 0) ADE_inf <- rbind(apply(mc_splm[ ,c( (1 : length(xb)), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE) + 
#'                                                                             apply(mc_splm[ , c((length(xb) + 1) : (ncol(mc_splm) - 1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE), 
#'                                                                          apply(mc_splm[ , c(length(xb) + length(gamma) +1) : (ncol(mc_splm) - 1), ncol(mc_splm)], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE))
#'             
#'             else if(length(P) != 0 && length(don) == 0) ADE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE),
#'                                                                          apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE) + 
#'                                                                             apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( ncol(mc_splm) - 1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE))  
#'             
#'             else ADE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE) 
#'                                   + apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE),
#'                                   apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE) + 
#'                                      apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( (length(P) + length(xb) + length(gamma))), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE),
#'                                   apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= TRUE) + 
#'                                      apply( mc_splm[ , c( (length(P) +  length(xb) +length(gamma) +1) : (length(P) + length(xb) +length(gamma) +length(don)), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= TRUE))
#'             
#'             
#'          }
#'          else{
#'             
#' 
#'             if (length(P) == 0 && length(don) == 0)       ADE_inf <- rbind(apply(mc_splm[ ,c((1 : length(xb)), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW ) + 
#'                                                                               apply( mc_splm[ , c((length(xb) + 1) : ( ncol(mc_splm)-1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW))
#'             
#'             else if(length(P) == 0 && length(don) != 0) ADE_inf <- rbind(apply(mc_splm[ ,c( (1 : length(xb)), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW) + 
#'                                                                             apply(mc_splm[ , c((length(xb) + 1) : (ncol(mc_splm) - 1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW), 
#'                                                                          apply(mc_splm[ , c(length(xb) + length(gamma) +1) : (ncol(mc_splm) - 1), ncol(mc_splm)], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW))
#'             
#'             else if(length(P) != 0 && length(don) == 0) ADE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW),
#'                                                                          apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW) + 
#'                                                                             apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( ncol(mc_splm) - 1), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW))  
#'             
#'             else ADE_inf <- rbind(apply(mc_splm[ ,c(1:length(P), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW) 
#'                                   + apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW),
#'                                   apply( mc_splm[ ,c( ((length(P) + 1) : (length(P) + length(xb))), ncol(mc_splm))], 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW) + 
#'                                      apply( mc_splm[ , c((length(P) + length(xb) + 1) : ( (length(P) + length(xb) + length(gamma))), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW),
#'                                   apply(cbind(0, mc_splm[ ,c(ncol(mc_splm))]), 1, ADRE_rn_Dur_fp, weg, n, m, object$listw, exact= F, trW) + 
#'                                      apply( mc_splm[ , c( (length(P) +  length(xb) +length(gamma) +1) : (length(P) + length(xb) +length(gamma) +length(don)), ncol(mc_splm))], 1, ADRE_rn_Dur_sp, weg, n, m, object$listw, exact= F, trW))
#'             
#'             
#'             
#'          }
#'          
#'          AIE_inf <- ATE_inf - ADE_inf
#'          
#'          res <- list(ATE, ADE, AIE, ATE_inf, ADE_inf, AIE_inf) 
#'       }
#'       else res <- list(ATE, ADE, AIE)
#'       
#'       
#'       
#'       
#'       
#'       }
#'    
#'    else{
#'      
#'      coefs <- drop(object$coefficients)
#'      p2 <- length(coefs)
#'      lambda <- coefs[p2]
#'      beta <- coefs[1:(p2-1)]
#'      p <- length(beta)
#'      p1 <- p + 1
#'      icept <- grep("(Intercept)", names(beta))
#'      iicept <- length(icept) > 0
#'      Sigma <- as.matrix(object$var)  
#' 
#'      if (iicept) {
#'        P <- beta[-icept]
#'        Sigma <- Sigma[-c(icept),-c(icept)]
#'        bnames <- names(beta[-icept])
#'      } 
#'      else {
#'        P <- beta
#'        Sigma <- Sigma
#'        bnames <- names(beta)
#'      }
#'      
#'      n <- length(object$residuals)
#'      
#'      if(exact){
#'         
#'         if(is.null(weg)) weg <- eigen(object$listw)$values
#'         
#'         
#'         ATE <- AVTE_rn(c(P, lambda))
#'         ADE <- ADRE_rn(c(P, lambda), weg, n, exact = TRUE)
#'         AIE <- ATE - ADE
#'      }
#'      else{
#'         
#'         if(is.null(trW))   trW <-  smtrace(object$listw, m = m)
#'         
#'        
#'         ATE <- AVTE_rn(c(P, lambda))
#'         ADE <- ADRE_rn(c(P, lambda), exact = FALSE, trW = trW, m = m, n = n)
#'         AIE <- ATE - ADE
#'         
#'      }
#'      
#'      if (inference) { 
#'         
#'         mu <- c(P, lambda)
#'         mc_splm <- mvtnorm::rmvnorm(n_mvn, mean = mu, sigma = Sigma)
#'         
#'         if(exact){ 
#'            
#'            ATE_inf <- apply(mc_splm, 1, AVTE_rn)
#'            ADE_inf <- apply(mc_splm, 1, ADRE_rn, weg, n, exact = TRUE)
#'            AIE_inf <- ATE_inf - ADE_inf
#'            
#'         }
#'         else{
#'            
#'            ATE_inf <- apply(mc_splm, 1, AVTE_rn)
#'            ADE_inf <- apply(mc_splm, 1, ADRE_rn, exact = FALSE, trW = trW, m = m, n = n)
#'            AIE_inf <- ATE_inf - ADE_inf
#'            
#'         }
#'         res <- list(ATE, ADE, AIE, ATE_inf, ADE_inf, AIE_inf) 
#'      }
#'      else res <- list(ATE, ADE, AIE)
#' 
#'    }
#'    attr(res, "inference") <- inference
#'    class(res) <- c("impacts_sphet") 
#'    
#'    invisible(res)
#'  }
#'  
#'  #' Summary for spreg impacts
#'  #' 
#'  #' 
#'  #' @param object An object of class impacts.sphet
#'  #' @param ... Additional arguments that can be passed to the function
#'  #' @return Summary of the Average Total, Average Direct, and Average Indirect Effects
#'  #' @method summary impacts_sphet
#'  #' @examples 
#'  #' data(columbus, package = "spdep")
#'  #' listw <- spdep::nb2listw(col.gal.nb)
#'  #' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
#'  #'             het = TRUE, verbose = FALSE, model = "sarar")
#'  #' summary(res)
#'  #' effects <- impacts(res, n_mvn = 3000)
#'  #' summary(effects)
#'  summary.impacts_sphet <- function(object, ...){
#'   
#'    cat("\n ---------------------- \n")
#'    cat("\n ---------------------- \n")
#'    cat("\n Impacts estimates:\n")
#'    cat("\n ---------------------- \n")
#'    cat("\n ---------------------- \n")
#'    
#'   # print(cbind(ATE,apply(ATE_inf, 2, sd), (ATE - apply(ATE_inf, 2, mean) )/apply(ATE_inf, 2, sd)))
#'    #print(cbind(obj[[1]], sd(obj[[4]]), obj[[1]]/sd(obj[[4]]) ))
#'    eff <- cbind(object[[1]], object[[2]], object[[3]])
#'    colnames(eff) <- c("ATE", "ADE", "AIE")
#'    print(eff)
#'    
#' if(attributes(object)$inference){   
#'    cat("\n ---------------------- \n")
#'    cat("\n Inference for ATE:\n")
#'    cat("\n ---------------------- \n")
#'    
#'    mean_ate <- apply(object[[4]], 1, mean)
#'   # print(head(mean_ate))
#'    sd_ate <- apply(object[[4]], 1, sd)
#'    z_ate <- (object[[1]])/sd_ate
#'    p_ate <- 2*pnorm(abs(z_ate), lower.tail = F)
#'    total <- cbind(object[[1]], sd_ate, z_ate, p_ate)
#'    colnames(total) <- c("ATE", "sd", "z-stat","p-val")
#'    print(total)
#'    cat("\n Summary stats for the simulated ATE samples:\n")
#'    a <- matrix(nrow = nrow(total), ncol = 6)
#'   for (i in 1:nrow(total)) a[i,] <- c(quantile(object[[4]][i,], c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(object[[4]][i,]))
#'    colnames(a) <- c("0.025", "Q1", "Median", "Q3", "0.975","Mean")
#'    rownames(a) <- rownames(eff)
#'    print(a)
#'    
#'    cat("\n ---------------------- \n")
#'    cat("\n Inference for ADE:\n")
#'    cat("\n ---------------------- \n")
#'    
#'    mean_ade <- apply(object[[5]], 1, mean)
#'    # print(head(mean_ate))
#'    sd_ade <- apply(object[[5]], 1, sd)
#'    z_ade <- (object[[2]])/sd_ade
#'    p_ade <- 2*pnorm(abs(z_ade), lower.tail = F)
#'    direct <- cbind(object[[2]], sd_ade, z_ade, p_ade)
#'    colnames(direct) <- c("ADE", "sd", "z-stat","p-val")
#'    print(direct)
#'    cat("\n Summary stats for the simulated ADE samples:\n")
#'    
#'    for (i in 1:nrow(total)) a[i,] <- c(quantile(object[[5]][i,], c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(object[[5]][i,]))
#'    colnames(a) <- c("0.025", "Q1", "Median", "Q3", "0.975","Mean")
#'    rownames(a) <- rownames(eff)
#'    print(a)
#'    
#'    cat("\n ---------------------- \n")
#'    cat("\n Inference for AIE:\n")
#'    cat("\n ---------------------- \n")
#'    
#'    mean_aie <- apply(object[[6]], 1, mean)
#'    sd_aie <- apply(object[[6]], 1, sd)
#'    z_aie <- (object[[3]])/sd_aie
#'       p_aie <- 2*pnorm(abs(z_aie), lower.tail = F)
#'    indirect <- cbind(object[[3]], sd_aie, z_aie, p_aie)
#'    dim(indirect)
#'    colnames(indirect) <- c("AIE", "sd", "z-stat","p-val")
#'    print(indirect)
#'    cat("\n Summary stats for the simulated AIE samples:\n")
#'    
#'    for (i in 1:nrow(total)) a[i,] <- c(quantile(object[[6]][i,], c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(object[[6]][i,]))
#'    colnames(a) <- c("0.025", "Q1", "Median", "Q3", "0.975","Mean")
#'    rownames(a) <- rownames(eff)
#'    print(a)
#' }
#'    
#'    }
#'  
#'  
#' AVTE_rn <- function(x){
#'   x[-length(x)]/(1-x[length(x)])
#' }
#' AVTE_rn_Dur <- function(x,y){
#'    x[-length(x)]/(1-x[length(x)]) + y[-length(y)]/(1-y[length(y)]) 
#' }
#' AVTE_rn_Dur_fp <- function(x){
#'    x[-length(x)]/(1-x[length(x)]) 
#' }
#' AVTE_rn_Dur_sp <- function(y){
#'    y[-length(y)]/(1-y[length(y)]) 
#' }
#' 
#' 
#' ADRE_rn <- function(x, egv, n, m, exact, trW){
#'  if (exact) Re(x[-length(x)] *(1/n)*sum(1/(1-x[length(x)]*egv)))
#'    else x[-length(x)] *(1/n) *(sum(as.numeric(sapply(x[length(x)], function(s) s^seq(1,m)), m )/n *trW) + n)
#' }
#' ADRE_rn_Dur <- function(x, egv, n, m, listw, y, exact, trW){
#'    if(exact) Re(x[-length(x)] *(1/n)*sum(1/(1-x[length(x)]*egv))) +
#'    (1/n)*y[-length(y)] *sum(diag(solve(Diagonal(n)- y[length(y)]* listw) %*% listw))
#'    else x[-length(x)] *(1/n) *(sum(as.numeric(sapply(x[length(x)], function(s) s^seq(1,m)), m )/n *trW) + n) +
#'       y[-length(y)] *(1/n) * (sum(as.numeric(sapply(y[length(y)], function(s) s^seq(1,(m-1))), m )/n *trW[-1]) )
#' }
#' ADRE_rn_Dur_fp <- function(x, egv, n, m, listw, exact, trW){
#'    if(exact) Re(x[-length(x)] *(1/n)*sum(1/(1-x[length(x)]*egv))) 
#'    else x[-length(x)] *(1/n) *(sum(as.numeric(sapply(x[length(x)], function(s) s^seq(1,m)), m )/n *trW) + n) 
#'       
#' }
#' ADRE_rn_Dur_sp <- function(y, egv, n, m, listw, exact, trW){
#'    if(exact)  (1/n)*y[-length(y)] *sum(diag(solve(Diagonal(n)- y[length(y)]* listw) %*% listw))
#'    else y[-length(y)] *(1/n) * (sum(as.numeric(sapply(y[length(y)], function(s) s^seq(1,(m-1))), m )/n *trW[-1]) )
#' }
#' smtrace <- function(listw, m){
#'    listwp <- listw
#'    tr <- numeric(m)
#'    for (i in 1:m) {
#'       tr[i] <- sum(diag(listwp))
#'       listwp <- listw %*% listwp
#'    }
#'    return(tr)
#' }
