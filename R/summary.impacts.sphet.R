#' Generate impacts for spreg lag and sarar models
#' 
#' 
#' @param obj An object of class sphet
#' @param mvn Inference based on MC simulation of mvn samples
#' @param n_mvn Number of samples
#' @param weg Eigenvalues can be provided by the user
#' @return Estimate of the Average Total, Average Direct, and Average Indirect Effects
#' 
#' @examples 
#' data(columbus, package="spdep")
#' listw <- spdep::nb2listw(col.gal.nb)
#' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
#'             het = TRUE, verbose = FALSE, model = "sarar")
#' summary(res)
#' effects <- impacts(res, n_mvn = 3000)
#' summary(effects)
impacts <- function(obj, mvn = TRUE, n_mvn = 3000, weg = NULL) {
  UseMethod("impacts", obj)
}

impacts.sphet <- function(obj, mvn = TRUE, n_mvn = 3000, weg = NULL) {

  if(!is.null(obj$endog)) stop("Impacts not yet implemented 
                               for models with additional endogenous variables")
 
  if(class(obj)[[2]] %in% c("ols", "error_gmm")) stop("Impacts not implemented
                                                 for error and ols models")
  
  if(!isTRUE(all.equal(Matrix::rowSums(obj$listw), rep(1,dim(obj$listw)[2])))) stop("Impacts not yet implemented 
                                          for non row-standardized listw")
    
  ag <- class(obj)[[2]]
  switch(match.arg(ag, c("lag_gmm", "sarar_gmm")),
         sarar_gmm = impacts.sarar(obj, mvn, n_mvn),
         lag_gmm = impacts.lag(obj, mvn, n_mvn),
         )

}
  
 impacts.sarar <-function(obj, mvn = TRUE, n_mvn = n_mvn, weg = NULL ){ 
   
   if(obj$Durbin) stop("Impacts not yet implemented
                      for Durbin models")
   
  else{
    
  coefs <- drop(obj$coefficients)
  p2 <- length(coefs)
  lambda <- coefs[(p2-1)]
  beta <- coefs[1:(p2-2)]
  p <- length(beta)
  p1 <- p + 1
  rho <- coefs[p2]
  icept <- grep("(Intercept)", names(beta))
  iicept <- length(icept) > 0
  Sigma <- obj$var  

  if (iicept) {
    P <- beta[-icept]
    Sigma <- Sigma[-c(icept,p2),-c(icept,p2)]
    bnames <- names(beta[-icept])
  } 
  else {
    P <- beta
    Sigma <- Sigma[-c(p2),-c(p2)]
    bnames <- names(beta)
  }
  
 if(is.null(weg)) weg <- eigen(obj$listw)$values
  n <- length(obj$residuals)

  ATE <- AVTE_rn(c(P, lambda))
  ADE <- ADRE_rn(c(P, lambda), weg, n)
  AIE <- ATE - ADE
  
  mu <- c(P, lambda)
  mc_splm <- mvtnorm::rmvnorm(n_mvn, mean = mu, sigma = Sigma)
  
  ATE_inf <- apply(mc_splm, 1, AVTE_rn)
  ADE_inf <- apply(mc_splm, 1, ADRE_rn, weg, n)
  AIE_inf <- ATE_inf - ADE_inf
  
  res <- list(ATE, ADE, AIE, ATE_inf, ADE_inf, AIE_inf) 
  
  
  }
  #print(ATE)
   class(res) <- c("impacts_sphet") 

  # z_ate <- (ATE - apply(ATE_inf, 2, mean) )/apply(ATE_inf, 2, sd)
  # z_ade <- (ADE - apply(ADE_inf, 2, mean) )/apply(ADE_inf, 2, sd)
  # z_aie <- (AIE - apply(AIE_inf, 2, mean) )/apply(AIE_inf, 2, sd)
 
  #invisible(x)
# print(res)
  invisible(res)
 }
 
 impacts.lag <-function(obj, mvn = TRUE, n_mvn = n_mvn, weg = NULL ){ 
   
   if(obj$Durbin) stop("Impacts not yet implemented
                      for Durbin models")
   
   else{
     
     coefs <- drop(obj$coefficients)
    # print(coefs)
     p2 <- length(coefs)
     lambda <- coefs[p2]
    # print(lambda)
     beta <- coefs[1:(p2-1)]
    # print(beta)
     p <- length(beta)
     p1 <- p + 1
     icept <- grep("(Intercept)", names(beta))
     iicept <- length(icept) > 0
     Sigma <- obj$var  
     #print(Sigma)
     if (iicept) {
       P <- beta[-icept]
       Sigma <- Sigma[-c(icept),-c(icept)]
       bnames <- names(beta[-icept])
     } 
     else {
       P <- beta
       Sigma <- Sigma
       bnames <- names(beta)
     }
     
     if(is.null(weg)) weg <- eigen(obj$listw)$values
     n <- length(obj$residuals)
     
     ATE <- AVTE_rn(c(P, lambda))
     ADE <- ADRE_rn(c(P, lambda), weg, n)
     AIE <- ATE - ADE
     
     mu <- c(P, lambda)
     mc_splm <- mvtnorm::rmvnorm(n_mvn, mean = mu, sigma = as.matrix(Sigma))
     
     ATE_inf <- apply(mc_splm, 1, AVTE_rn)
     ADE_inf <- apply(mc_splm, 1, ADRE_rn, weg, n)
     AIE_inf <- ATE_inf - ADE_inf
     
     res <- list(Av_tot = ATE, Av_dir = ADE, av_ind = AIE, 
                 sim_ate = ATE_inf, sim_ade =  ADE_inf, sim_aie =  AIE_inf) 
     
     
   }
   #print(ATE)
   class(res) <- c("impacts_sphet") 
   
   # z_ate <- (ATE - apply(ATE_inf, 2, mean) )/apply(ATE_inf, 2, sd)
   # z_ade <- (ADE - apply(ADE_inf, 2, mean) )/apply(ADE_inf, 2, sd)
   # z_aie <- (AIE - apply(AIE_inf, 2, mean) )/apply(AIE_inf, 2, sd)
   
   #invisible(x)
   # print(res)
   invisible(res)
 }
 
 #' Summary for spreg impacts
 #' 
 #' 
 #' @param object An object of class impacts.sphet
 #' @param ... Additional arguments that can be passed to the function
 #' @return Summary of the Average Total, Average Direct, and Average Indirect Effects
 #' @method summary impacts_sphet
 #' @examples 
 #' data(columbus, package = "spdep")
 #' listw <- spdep::nb2listw(col.gal.nb)
 #' res <- spreg(CRIME~HOVAL + INC, data=columbus , listw= listw,
 #'             het = TRUE, verbose = FALSE, model = "sarar")
 #' summary(res)
 #' effects <- impacts(res, n_mvn = 3000)
 #' summary(effects)
 summary.impacts_sphet <- function(object, ...){
   cat("\n ---------------------- \n")
   cat("\n ---------------------- \n")
   cat("\n Impacts estimates:\n")
   cat("\n ---------------------- \n")
   cat("\n ---------------------- \n")
   
  # print(cbind(ATE,apply(ATE_inf, 2, sd), (ATE - apply(ATE_inf, 2, mean) )/apply(ATE_inf, 2, sd)))
   #print(cbind(obj[[1]], sd(obj[[4]]), obj[[1]]/sd(obj[[4]]) ))
   eff <- cbind(object[[1]], object[[2]], object[[3]])
   colnames(eff) <- c("ATE", "ADE", "AIE")
   print(eff)
   
   
   cat("\n ---------------------- \n")
   cat("\n Inference for ATE:\n")
   cat("\n ---------------------- \n")
   
   mean_ate <- apply(object[[4]], 1, mean)
  # print(head(mean_ate))
   sd_ate <- apply(object[[4]], 1, sd)
   z_ate <- (object[[1]])/sd_ate
   p_ate <- 2*pnorm(abs(z_ate), lower.tail = F)
   total <- cbind(object[[1]], sd_ate, z_ate, p_ate)
   colnames(total) <- c("ATE", "sd", "z-stat","p-val")
   print(total)
   cat("\n Summary stats for the simulated ATE samples:\n")
   a <- matrix(nrow = nrow(total), ncol = 6)
  for (i in 1:nrow(total)) a[i,] <- c(quantile(object[[4]][i,], c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(object[[4]][i,]))
   colnames(a) <- c("0.025", "Q1", "Median", "Q3", "0.975","Mean")
   rownames(a) <- rownames(eff)
   print(a)
   
   cat("\n ---------------------- \n")
   cat("\n Inference for ADE:\n")
   cat("\n ---------------------- \n")
   
   mean_ade <- apply(object[[5]], 1, mean)
   # print(head(mean_ate))
   sd_ade <- apply(object[[5]], 1, sd)
   z_ade <- (object[[2]])/sd_ade
   p_ade <- 2*pnorm(abs(z_ade), lower.tail = F)
   direct <- cbind(object[[2]], sd_ade, z_ade, p_ade)
   colnames(direct) <- c("ADE", "sd", "z-stat","p-val")
   print(direct)
   cat("\n Summary stats for the simulated ADE samples:\n")
   
   for (i in 1:nrow(total)) a[i,] <- c(quantile(object[[5]][i,], c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(object[[5]][i,]))
   colnames(a) <- c("0.025", "Q1", "Median", "Q3", "0.975","Mean")
   rownames(a) <- rownames(eff)
   print(a)
   
   cat("\n ---------------------- \n")
   cat("\n Inference for AIE:\n")
   cat("\n ---------------------- \n")
   
   mean_aie <- apply(object[[6]], 1, mean)
   sd_aie <- apply(object[[6]], 1, sd)
   z_aie <- (object[[3]])/sd_aie
      p_aie <- 2*pnorm(abs(z_aie), lower.tail = F)
   indirect <- cbind(object[[3]], sd_aie, z_aie, p_aie)
   dim(indirect)
   colnames(indirect) <- c("AIE", "sd", "z-stat","p-val")
   print(indirect)
   cat("\n Summary stats for the simulated AIE samples:\n")
   
   for (i in 1:nrow(total)) a[i,] <- c(quantile(object[[6]][i,], c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(object[[6]][i,]))
   colnames(a) <- c("0.025", "Q1", "Median", "Q3", "0.975","Mean")
   rownames(a) <- rownames(eff)
   print(a)
   
   }
 
 
AVTE_rn <- function(x){
  x[-length(x)]/(1-x[length(x)])
}
ADRE_rn <- function(x, egv, n){
 Re(x[-length(x)] *(1/n)*sum(1/(1-x[length(x)]*egv)))
}



