vc_impacts_formula_lag <- function(obj, ev, tr = NULL, prt = T){
  
  type <- "lag" 
  if (is.null(obj$interval)){
  if(!is.null(tr)) interval <- 1/c(-1,1)
    else  interval <- 1/range(Re(ev))
  }
  
  betas <- coefficients(obj)
  p <- length(betas)
  beta <- betas[-p]
  names(beta) <- rownames(obj$coefficients)[1:(p-1)]
  lambda <- betas[p]
 
  if((lambda > interval[2] ) | (lambda < interval[1])) warning("Value of the spatial parameter outside of parameter space")     
  
  Sigma <- obj$var
  n <- length(obj$residuals)
  
  icept <- grep("(Intercept)", names(beta))
  iicept <- length(icept) > 0L
  
  if (iicept) {
    beta <- beta[-icept]
    Sigma <- Sigma[-icept, -icept]
    bnames <- names(beta)
  } 
  
  else  bnames <- names(beta)

  p1 <- length(beta)  
  pl <- ncol(Sigma)
  ############ ATE #######
  ATE <- c()
  for(i in 1:p1) ATE <- c(ATE, beta[i]/(1-lambda)) 
  
  p_l_ATE <- c()
  for(i in 1:p1) p_l_ATE <- c(p_l_ATE, beta[i]/(1-lambda)^2)
  p_b_ATE <- 1/(1-lambda)
  der_ATE <- matrix(cbind(rep(p_b_ATE, p1), p_l_ATE), nrow = p1, ncol = 2)
  se_ATE <- c()
  for(i in 1:p1) se_ATE <- c(se_ATE, sqrt(as.numeric(matrix(der_ATE[i,], nrow = 1, ncol = 2) %*% 
                                                       Sigma[c(i,pl),c(i,pl)] %*% t(matrix(der_ATE[i,], nrow = 1, ncol = 2)))))
 
  
  ############ ADE######
  if(!is.null(tr)){
    powl <- vector("numeric", length = length(tr))
    for(i in 1:length(tr)) powl[i] <-  lambda^i
    tr_G <- n +sum(tr * powl)
    scl <- 2:length(tr)
    powl <- vector("numeric", length = (length(tr)-1))
    for(i in 1:length(powl)) powl[i] <-  lambda^i
    dv_l <- sum(tr * c(0,powl)*c(0,scl))
  }
  else{
  tr_G <- sum(1/(1-lambda*ev))
  dv_l <- sum(ev/(1-lambda*ev)^2)
  }
  
  ADE <- c()
  for(i in 1:p1) ADE <- c(ADE, Re((1/n)*beta[i]*tr_G)) 
  ## Partial Derivatives
  p_l_ADE <- c()
  for(i in 1:p1) p_l_ADE <- c(p_l_ADE, Re((1/n)*beta[i]*dv_l))
  p_b_ADE <- Re((1/n)*tr_G)
  der_ADE <- matrix(cbind(rep(p_b_ADE, p1), p_l_ADE), nrow = p1, ncol = 2)
  se_ADE <- c()
  for(i in 1:p1) se_ADE <- c(se_ADE, sqrt(as.numeric(matrix(der_ADE[i,], nrow = 1, ncol = 2) %*% 
                                                    Sigma[c(i,pl),c(i,pl)] %*% t(matrix(der_ADE[i,], nrow = 1, ncol = 2)))))


  
  ######### AIE ######
  AIE <- ATE - ADE
  der_AIE <- der_ATE - der_ADE
  se_AIE <- c()
  for(i in 1:p1) se_AIE <- c(se_AIE, sqrt(as.numeric(matrix(der_AIE[i,], nrow = 1, ncol = 2) %*% 
                                        Sigma[c(i,pl),c(i,pl)] %*% t(matrix(der_AIE[i,], nrow = 1, ncol = 2)))))
  
  #############################################
  
  
  if(prt){
  
  cat("Impact Measures (lag, KP_formula):\n")
  tb <- cbind(ADE, AIE, ATE)
  colnames(tb) <- c("Direct", "Indirect", "Total")
  rownames(tb) <- bnames
  print(tb)
  cat("========================================================\n")
  cat("Results based on Kelejian and Piras formula:\n")
  cat("========================================================\n")
  cat("Analytical standard errors\n")
  se <- cbind(se_ADE, se_AIE, se_ATE)
  colnames(se) <- c("Direct", "Indirect", "Total")
  rownames(se) <- bnames
  print(se)
  cat("\nAnalytical z-values:\n")
  mat <- cbind(ADE, AIE, ATE)/se
  colnames(mat) <- c("Direct", "Indirect", "Total")
  print(mat)
  cat("\nAnalytical p-values:\n")
  xx <- apply(2*(1-pnorm(abs(mat))), 2, format.pval)
  if (length(bnames) == 1L) {
    xx <- matrix(xx, ncol=3)
    colnames(xx) <- c("Direct", "Indirect", "Total")
  }
  rownames(xx) <- bnames
 print(xx, quote=FALSE)
  }
  else{
    tb <- cbind(ADE, AIE, ATE)
    colnames(tb) <- c("Direct", "Indirect", "Total")
    rownames(tb) <- bnames
    se <- cbind(se_ADE, se_AIE, se_ATE)
    colnames(se) <- c("Direct", "Indirect", "Total")
    rownames(se) <- bnames
    mat <- cbind(ADE, AIE, ATE)/se
    colnames(mat) <- c("Direct", "Indirect", "Total")
    rownames(mat) <- bnames
    xx <- apply(2*(1-pnorm(abs(mat))), 2, format.pval)
    if (length(bnames) == 1L) {
      xx <- matrix(xx, ncol=3)
      colnames(xx) <- c("Direct", "Indirect", "Total")
    }
    rownames(xx) <- bnames
  }
  results <- list(tb, se, mat, xx)
 invisible(results)
}
vc_impacts_formula_sarar <- function(obj, ev, tr = NULL, prt = T){
  
  type <- "sarar" 
  
  if (is.null(obj$interval)){
    if(!is.null(tr)) interval <- 1/c(-1,1)
    else  interval <- 1/range(Re(ev))
  }
  
  betas <- coefficients(obj)
  p <- length(betas)
  beta <- betas[1:(p-2)]
  names(beta) <- rownames(obj$coefficients)[1:(p-2)]
  lambda <- betas[p-1]
  
  if((lambda > interval[2] ) | (lambda < interval[1])) warning("Value of the spatial parameter outside of parameter space")     
  
  Sigma <- obj$var
  n <- length(obj$residuals)
  
  icept <- grep("(Intercept)", names(beta))
  iicept <- length(icept) > 0L
  
  if (iicept) {
    beta <- beta[-icept]
    Sigma <- Sigma[-icept, -icept]
    bnames <- names(beta)
  } 
  
  else  bnames <- names(beta)
  
  p1 <- length(beta)  
  pl <- ncol(Sigma)-1
  
  ############ ATE #######
  ATE <- c()
  for(i in 1:p1) ATE <- c(ATE, beta[i]/(1-lambda)) 
  
  p_l_ATE <- c()
  for(i in 1:p1) p_l_ATE <- c(p_l_ATE, beta[i]/(1-lambda)^2)
  p_b_ATE <- 1/(1-lambda)
  
  der_ATE <- matrix(cbind(rep(p_b_ATE, p1), p_l_ATE), nrow = p1, ncol = 2)
  se_ATE <- c()
  for(i in 1:p1) se_ATE <- c(se_ATE, sqrt(as.numeric(matrix(der_ATE[i,], nrow = 1, ncol = 2) %*% 
                                                       Sigma[c(i,pl),c(i,pl)] %*% t(matrix(der_ATE[i,], nrow = 1, ncol = 2)))))
  ############ ADE######
  if(!is.null(tr)){
    powl <- vector("numeric", length = length(tr))
    for(i in 1:length(tr)) powl[i] <-  lambda^i
    tr_G <- n +sum(tr * powl)
    scl <- 2:length(tr)
    powl <- vector("numeric", length = (length(tr)-1))
    for(i in 1:length(powl)) powl[i] <-  lambda^i
    dv_l <- sum(tr * c(0,powl)*c(0,scl))
  }
  else{
    tr_G <- sum(1/(1-lambda*ev))
    dv_l <- sum(ev/(1-lambda*ev)^2)
  }
  
  ADE <- c()
  for(i in 1:p1) ADE <- c(ADE, Re((1/n)*beta[i]*tr_G)) 
  
  ## Partial Derivatives
  p_l_ADE <- c()
  for(i in 1:p1) p_l_ADE <- c(p_l_ADE, Re((1/n)*beta[i]*dv_l))
  p_b_ADE <- Re((1/n)*tr_G)
  der_ADE <- matrix(cbind(rep(p_b_ADE, p1), p_l_ADE), nrow = p1, ncol = 2)
  se_ADE <- c()
  for(i in 1:p1) se_ADE <- c(se_ADE, sqrt(as.numeric(matrix(der_ADE[i,], nrow = 1, ncol = 2) %*% 
                                                       Sigma[c(i,pl),c(i,pl)] %*% t(matrix(der_ADE[i,], nrow = 1, ncol = 2)))))
  
  
  
  ######### AIE ######
  AIE <- ATE - ADE
  der_AIE <- der_ATE - der_ADE
  se_AIE <- c()
  for(i in 1:p1) se_AIE <- c(se_AIE, sqrt(as.numeric(matrix(der_AIE[i,], nrow = 1, ncol = 2) %*% 
                                        Sigma[c(i,pl),c(i,pl)] %*% t(matrix(der_AIE[i,], nrow = 1, ncol = 2)))))
  
  #############################################
  
  if(prt){
    
  cat("Impact Measures (lag, KP_formula):\n")
  tb <- cbind(ADE, AIE, ATE)
  colnames(tb) <- c("Direct", "Indirect", "Total")
  rownames(tb) <- bnames
  print(tb)
  cat("========================================================\n")
  cat("Results based on Kelejian and Piras formula:\n")
  cat("========================================================\n")
  cat("Analytical standard errors\n")
  se <- cbind(se_ADE, se_AIE, se_ATE)
  colnames(se) <- c("Direct", "Indirect", "Total")
  rownames(se) <- bnames
  print(se)
  cat("\nAnalytical z-values:\n")
  mat <- cbind(ADE, AIE, ATE)/se
  colnames(mat) <- c("Direct", "Indirect", "Total")
  print(mat)
  cat("\nAnalytical p-values:\n")
  xx <- apply(2*(1-pnorm(abs(mat))), 2, format.pval)
  if (length(bnames) == 1L) {
    xx <- matrix(xx, ncol=3)
    colnames(xx) <- c("Direct", "Indirect", "Total")
  }
  rownames(xx) <- bnames
  print(xx, quote=FALSE)
  }
  else{
    tb <- cbind(ADE, AIE, ATE)
    colnames(tb) <- c("Direct", "Indirect", "Total")
    rownames(tb) <- bnames
    se <- cbind(se_ADE, se_AIE, se_ATE)
    colnames(se) <- c("Direct", "Indirect", "Total")
    rownames(se) <- bnames
    mat <- cbind(ADE, AIE, ATE)/se
    colnames(mat) <- c("Direct", "Indirect", "Total")
    rownames(mat) <- bnames
    xx <- apply(2*(1-pnorm(abs(mat))), 2, format.pval)
    if (length(bnames) == 1L) {
      xx <- matrix(xx, ncol=3)
      colnames(xx) <- c("Direct", "Indirect", "Total")
    }
    rownames(xx) <- bnames
  }
  results <- list(tb, se, mat, xx)
  invisible(results)
}


vc_impacts_formula_lag_mixed <- function(obj, ev, tr = NULL, prt = T){
  
  type <- "mixed" 
  
  if (is.null(obj$interval)){
    if(!is.null(tr)) interval <- 1/c(-1,1)
    else  interval <- 1/range(Re(ev))
  }
  
if(isTRUE(obj$Durbin)){  
  
  betas <- coefficients(obj)
  p <- length(betas)
  beta <- betas[-p]
  names(beta) <- rownames(obj$coefficients)[1:(p-1)]
  lambda <- betas[p]
  
  if((lambda > interval[2] ) | (lambda < interval[1])) warning("Value of the spatial parameter outside of parameter space")     
  
  Sigma <- obj$var
  n <- length(obj$residuals)
  
  icept <- grep("(Intercept)", names(beta))
  iicept <- length(icept) > 0L
  
  if (iicept) {
    beta <- beta[-icept]
    Sigma <- Sigma[-icept, -icept]
    bnames <- names(beta)[1:(length(beta)/2)]
  } 
  
  else  bnames <- names(beta)[1:(length(beta)/2)]
  
  p1 <- length(beta)/2  
  pl <- ncol(Sigma)
  P <- matrix(beta, ncol = 2)
 
   ############ ATE #######
  ATE <- c()
  for(i in 1:p1) ATE <- c(ATE, P[i,1]/(1-lambda) + P[i,2]/(1-lambda)) 
  
  p_l_ATE <- c()
  for(i in 1:p1) p_l_ATE <- c(p_l_ATE, P[i,1]/(1-lambda)^2+ P[i,2]/(1-lambda)^2)
  p_b_ATE <- 1/(1-lambda)
  p_g_ATE <- 1/(1-lambda) 
  der_ATE <- matrix(cbind(rep(p_b_ATE, p1), rep(p_g_ATE, p1), p_l_ATE) , nrow = p1, ncol = 3)
  
  se_ATE <- c()
  for(i in 1:p1) se_ATE <- c(se_ATE, sqrt(as.numeric(matrix(der_ATE[i,], nrow = 1, ncol = 3) %*% 
                                                       Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                       t(matrix(der_ATE[i,], nrow = 1, ncol = 3)))))
  
  
  ############ ADE######
  
  if(!is.null(tr)){
    powl <- vector("numeric", length = length(tr))
    for(i in 1:length(tr)) powl[i] <-  lambda^i
    tr_G <- n +sum(tr * powl)
    powl <- vector("numeric", length = (length(tr)-1))
    for(i in 0:(length(tr)-1)) powl[i] <-  lambda^i
    powl <- c(0,powl)
    tr_H <- sum(tr * powl)
    powl <- vector("numeric", length = (length(tr)-2))
    for(i in 1:(length(tr)-2)) powl[i] <-  lambda^i
    powl <- c(0,1,powl)
    scl <- 0:(length(tr)-1)
    tr_H2 <- sum(tr * scl*powl)
    scl <- 2:length(tr)
    powl <- vector("numeric", length = (length(tr)-1))
    for(i in 1:length(powl)) powl[i] <-  lambda^i
    dv_l <- sum(tr * c(0,powl)*c(0,scl))
    
  }
  else{
    tr_G <- sum(1/(1-lambda*ev))
    tr_H <- sum(ev/(1-lambda*ev))
    tr_H2 <- sum(ev^2/(1-lambda*ev)^2)
    dv_l <- sum(ev/(1-lambda*ev)^2)
  }
  
  
  ADE <- c()
  for(i in 1:p1) ADE <- c(ADE, Re((1/n)*P[i,1]*tr_G) + Re((1/n)*P[i,2]*tr_H)) 
  
  p_l_ADE <- c()
  for(i in 1:p1) p_l_ADE <- c(p_l_ADE, Re((1/n)*P[i,1] * dv_l) + Re((1/n)*P[i,2] * tr_H2))
  p_b_ADE <- Re((1/n)*tr_G)
  p_g_ADE <- Re((1/n)*tr_H)
  der_ADE <- matrix(cbind(rep(p_b_ADE, p1), rep(p_g_ADE, p1), p_l_ADE), nrow = p1, ncol = 3)
  
  se_ADE <- c()
  for(i in 1:p1) se_ADE <- c(se_ADE, sqrt(as.numeric(matrix(der_ADE[i,], nrow = 1, ncol = 3) %*% 
                                                       Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                       t(matrix(der_ADE[i,], nrow = 1, ncol = 3)))))

  ######### AIE ######
  AIE <- ATE - ADE
  der_AIE <- der_ATE - der_ADE
  se_AIE <- c()
  for(i in 1:p1) se_AIE <- c(se_AIE, sqrt(as.numeric(matrix(der_AIE[i,], nrow = 1, ncol = 3) %*% 
                                                       Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                       t(matrix(der_AIE[i,], nrow = 1, ncol = 3)))))
  
  #############################################
}
  else{
   
  
    coefs <- drop(obj$coefficients)
    p2 <- length(coefs)
    lambda <- coefs[p2]
    beta <- coefs[1:(p2-1)]
    p <- length(beta)
    p1 <- p + 1
    names(beta) <- rownames(obj$coefficients)[1:(p2-1)]
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0L
    n <- length(obj$residuals)
    Sigma <- obj$var

    if((lambda > interval[2] ) | (lambda < interval[1])) warning("Value of the spatial parameter outside of parameter space")    
    
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
      Sigma <- Sigma[-icept, -icept]
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
    
    p1 <- length(b2)/2  
    pl <- ncol(Sigma)
    
    adS <- length(which(b2 == 0))
    cnmS <- colnames(Sigma)
    Sigma <- rbind(cbind(Sigma , matrix(0, nrow = nrow(Sigma), ncol = adS)), matrix(0, nrow = adS, ncol =ncol(Sigma) + adS))
    colnames(Sigma) <- rownames(Sigma) <- c(cnmS, rep("", adS)) 
    ## Modifying sigma
    b2nam <- c(names(b2), "lambda") 
    #####
    Sigma <- Sigma[match(b2nam, rownames(Sigma)), match(b2nam, rownames(Sigma))]
   
    ############ ATE #######
    ATE <- c()
    for(i in 1:p1) ATE <- c(ATE, P[i,1]/(1-lambda) + P[i,2]/(1-lambda)) 
    
    p_l_ATE <- c()
    for(i in 1:p1) p_l_ATE <- c(p_l_ATE, P[i,1]/(1-lambda)^2+ P[i,2]/(1-lambda)^2)
    p_b_ATE <- 1/(1-lambda)
    p_g_ATE <- 1/(1-lambda) 
    der_ATE <- matrix(cbind(rep(p_b_ATE, p1), rep(p_g_ATE, p1), p_l_ATE) , nrow = p1, ncol = 3)
    
    se_ATE <- c()
    for(i in 1:p1) se_ATE <- c(se_ATE, sqrt(as.numeric(matrix(der_ATE[i,], nrow = 1, ncol = 3) %*% 
                                                         Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                         t(matrix(der_ATE[i,], nrow = 1, ncol = 3)))))
    
    
    ############ ADE######
    if(!is.null(tr)){
      powl <- vector("numeric", length = length(tr))
      for(i in 1:length(tr)) powl[i] <-  lambda^i
      tr_G <- n +sum(tr * powl)
      powl <- vector("numeric", length = (length(tr)-1))
      for(i in 0:(length(tr)-1)) powl[i] <-  lambda^i
      powl <- c(0,powl)
      tr_H <- sum(tr * powl)
      powl <- vector("numeric", length = (length(tr)-2))
      for(i in 1:(length(tr)-2)) powl[i] <-  lambda^i
      powl <- c(0,1,powl)
      scl <- 0:(length(tr)-1)
      tr_H2 <- sum(tr * scl*powl)
      scl <- 2:length(tr)
      powl <- vector("numeric", length = (length(tr)-1))
      for(i in 1:length(powl)) powl[i] <-  lambda^i
      dv_l <- sum(tr * c(0,powl)*c(0,scl))
      
    }
    else{
      tr_G <- sum(1/(1-lambda*ev))
      tr_H <- sum(ev/(1-lambda*ev))
      tr_H2 <- sum(ev^2/(1-lambda*ev)^2)
      dv_l <- sum(ev/(1-lambda*ev)^2)
    }
    
    ADE <- c()
    for(i in 1:p1) ADE <- c(ADE, Re((1/n)*P[i,1]*tr_G) + Re((1/n)*P[i,2]*tr_H)) 
    
    p_l_ADE <- c()
    for(i in 1:p1) p_l_ADE <- c(p_l_ADE, Re((1/n)*P[i,1] * dv_l) + Re((1/n)*P[i,2] * tr_H2))
    p_b_ADE <- Re((1/n)*tr_G)
    p_g_ADE <- Re((1/n)*tr_H)
    der_ADE <- matrix(cbind(rep(p_b_ADE, p1), rep(p_g_ADE, p1), p_l_ADE), nrow = p1, ncol = 3)
    
    se_ADE <- c()
    for(i in 1:p1) se_ADE <- c(se_ADE, sqrt(as.numeric(matrix(der_ADE[i,], nrow = 1, ncol = 3) %*% 
                                                         Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                         t(matrix(der_ADE[i,], nrow = 1, ncol = 3)))))
    
    ######### AIE ######
    AIE <- ATE - ADE
    der_AIE <- der_ATE - der_ADE
    se_AIE <- c()
    for(i in 1:p1) se_AIE <- c(se_AIE, sqrt(as.numeric(matrix(der_AIE[i,], nrow = 1, ncol = 3) %*% 
                                                         Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                         t(matrix(der_AIE[i,], nrow = 1, ncol = 3)))))
    
    #############################################    
     
  }
  
  if(prt){
  
  cat("Impact Measures (lag, KP_formula):\n")
  tb <- cbind(ADE, AIE, ATE)
  colnames(tb) <- c("Direct", "Indirect", "Total")
  rownames(tb) <- bnames
  print(tb)
  cat("========================================================\n")
  cat("Results based on Kelejian and Piras formula:\n")
  cat("========================================================\n")
  cat("Analytical standard errors\n")
  se <- cbind(se_ADE, se_AIE, se_ATE)
  colnames(se) <- c("Direct", "Indirect", "Total")
  rownames(se) <- bnames
  print(se)
  cat("\nAnalytical z-values:\n")
  mat <- cbind(ADE, AIE, ATE)/se
  colnames(mat) <- c("Direct", "Indirect", "Total")
  rownames(mat) <- bnames
  print(mat)
  cat("\nAnalytical p-values:\n")
  xx <- apply(2*(1-pnorm(abs(mat))), 2, format.pval)
  if (length(bnames) == 1L) {
    xx <- matrix(xx, ncol=3)
    colnames(xx) <- c("Direct", "Indirect", "Total")
  }
  rownames(xx) <- bnames
  print(xx, quote=FALSE)
  }
  else{
    tb <- cbind(ADE, AIE, ATE)
    colnames(tb) <- c("Direct", "Indirect", "Total")
    rownames(tb) <- bnames
    se <- cbind(se_ADE, se_AIE, se_ATE)
    colnames(se) <- c("Direct", "Indirect", "Total")
    rownames(se) <- bnames
    mat <- cbind(ADE, AIE, ATE)/se
    colnames(mat) <- c("Direct", "Indirect", "Total")
    rownames(mat) <- bnames
    xx <- apply(2*(1-pnorm(abs(mat))), 2, format.pval)
    if (length(bnames) == 1L) {
      xx <- matrix(xx, ncol=3)
      colnames(xx) <- c("Direct", "Indirect", "Total")
    }
    rownames(xx) <- bnames
  }
  results <- list(tb, se, mat, xx)
  invisible(results)
}
vc_impacts_formula_sarar_mixed <- function(obj, ev, tr = NULL, prt = T){
  
  type <- "mixed" 
  
  if (is.null(obj$interval)){
    if(!is.null(tr)) interval <- 1/c(-1,1)
    else  interval <- 1/range(Re(ev))
  }
  
  if(isTRUE(obj$Durbin)){  
    
    betas <- coefficients(obj)
    p <- length(betas)
    beta <- betas[1:(p-2)]
    names(beta) <- rownames(obj$coefficients)[1:(p-2)]
    lambda <- betas[p-1]
    
    
    
    if((lambda > interval[2] ) | (lambda < interval[1])) warning("Value of the spatial parameter outside of parameter space")     
    
    Sigma <- obj$var[-nrow(obj$var), -nrow(obj$var)]
    n <- length(obj$residuals)
    
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0L
    
    if (iicept) {
      beta <- beta[-icept]
      Sigma <- Sigma[-icept, -icept]
      bnames <- names(beta)[1:(length(beta)/2)]
    } 
    
    else  bnames <- names(beta)[1:(length(beta)/2)]
    
    p1 <- length(beta)/2  
    pl <- ncol(Sigma)
    P <- matrix(beta, ncol = 2)
    ############ ATE #######
    ATE <- c()
    for(i in 1:p1) ATE <- c(ATE, P[i,1]/(1-lambda) + P[i,2]/(1-lambda)) 
    
    p_l_ATE <- c()
    for(i in 1:p1) p_l_ATE <- c(p_l_ATE, P[i,1]/(1-lambda)^2+ P[i,2]/(1-lambda)^2)
    p_b_ATE <- 1/(1-lambda)
    p_g_ATE <- 1/(1-lambda) 
    der_ATE <- matrix(cbind(rep(p_b_ATE, p1), rep(p_g_ATE, p1), p_l_ATE) , nrow = p1, ncol = 3)
    
    se_ATE <- c()
    for(i in 1:p1) se_ATE <- c(se_ATE, sqrt(as.numeric(matrix(der_ATE[i,], nrow = 1, ncol = 3) %*% 
                                                         Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                         t(matrix(der_ATE[i,], nrow = 1, ncol = 3)))))
    
    
    ############ ADE######
    if(!is.null(tr)){
      powl <- vector("numeric", length = length(tr))
      for(i in 1:length(tr)) powl[i] <-  lambda^i
      tr_G <- n +sum(tr * powl)
      powl <- vector("numeric", length = (length(tr)-1))
      for(i in 0:(length(tr)-1)) powl[i] <-  lambda^i
      powl <- c(0,powl)
      tr_H <- sum(tr * powl)
      powl <- vector("numeric", length = (length(tr)-2))
      for(i in 1:(length(tr)-2)) powl[i] <-  lambda^i
      powl <- c(0,1,powl)
      scl <- 0:(length(tr)-1)
      tr_H2 <- sum(tr * scl*powl)
      scl <- 2:length(tr)
      powl <- vector("numeric", length = (length(tr)-1))
      for(i in 1:length(powl)) powl[i] <-  lambda^i
      dv_l <- sum(tr * c(0,powl)*c(0,scl))
      
    }
    else{
      tr_G <- sum(1/(1-lambda*ev))
      tr_H <- sum(ev/(1-lambda*ev))
      tr_H2 <- sum(ev^2/(1-lambda*ev)^2)
      dv_l <- sum(ev/(1-lambda*ev)^2)
    }
    
    ADE <- c()
    for(i in 1:p1) ADE <- c(ADE, Re((1/n)*P[i,1]*tr_G) + Re((1/n)*P[i,2]*tr_H)) 
    
    p_l_ADE <- c()
    for(i in 1:p1) p_l_ADE <- c(p_l_ADE, Re((1/n)*P[i,1] * dv_l) + Re((1/n)*P[i,2] * tr_H2))
    p_b_ADE <- Re((1/n)*tr_G)
    p_g_ADE <- Re((1/n)*tr_H)
    der_ADE <- matrix(cbind(rep(p_b_ADE, p1), rep(p_g_ADE, p1), p_l_ADE), nrow = p1, ncol = 3)
    
    se_ADE <- c()
    for(i in 1:p1) se_ADE <- c(se_ADE, sqrt(as.numeric(matrix(der_ADE[i,], nrow = 1, ncol = 3) %*% 
                                                         Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                         t(matrix(der_ADE[i,], nrow = 1, ncol = 3)))))
    
    ######### AIE ######
    AIE <- ATE - ADE
    der_AIE <- der_ATE - der_ADE
    se_AIE <- c()
    for(i in 1:p1) se_AIE <- c(se_AIE, sqrt(as.numeric(matrix(der_AIE[i,], nrow = 1, ncol = 3) %*% 
                                                         Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                         t(matrix(der_AIE[i,], nrow = 1, ncol = 3)))))
    
    #############################################
  }
  else{
    
    coefs <- drop(obj$coefficients)
    cnames <- names(coefs)
    p2 <- length(coefs)
    lambda <- coefs[p2-1]
    beta <- coefs[1:(p2-2)]
    p <- length(beta)
    p1 <- p + 1
    names(beta) <- rownames(obj$coefficients)[1:(p2-2)]
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0L
    n <- length(obj$residuals)
    
    Sigma <- obj$var[-p2,-p2]
    #print(Sigma)
    if((lambda > interval[2] ) | (lambda < interval[1])) warning("Value of the spatial parameter outside of parameter space")    
    
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
      Sigma <- Sigma[-icept, -icept]
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
    
    p1 <- length(b2)/2  
    pl <- ncol(Sigma)
    
    adS <- length(which(b2 == 0))
    cnmS <- colnames(Sigma)
    Sigma <- rbind(cbind(Sigma , matrix(0, nrow = nrow(Sigma), ncol = adS)), matrix(0, nrow = adS, ncol =ncol(Sigma) + adS))
    colnames(Sigma) <- rownames(Sigma) <- c(cnmS, rep("", adS)) 
    ## Modifying sigma
    b2nam <- c(names(b2), "lambda") 
    #####
    Sigma <- Sigma[match(b2nam, rownames(Sigma)), match(b2nam, rownames(Sigma))]
    
    ############ ATE #######
    ATE <- c()
    for(i in 1:p1) ATE <- c(ATE, P[i,1]/(1-lambda) + P[i,2]/(1-lambda)) 
    
    p_l_ATE <- c()
    for(i in 1:p1) p_l_ATE <- c(p_l_ATE, P[i,1]/(1-lambda)^2+ P[i,2]/(1-lambda)^2)
    p_b_ATE <- 1/(1-lambda)
    p_g_ATE <- 1/(1-lambda) 
    der_ATE <- matrix(cbind(rep(p_b_ATE, p1), rep(p_g_ATE, p1), p_l_ATE) , nrow = p1, ncol = 3)
    
    se_ATE <- c()
    for(i in 1:p1) se_ATE <- c(se_ATE, sqrt(as.numeric(matrix(der_ATE[i,], nrow = 1, ncol = 3) %*% 
                                                         Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                         t(matrix(der_ATE[i,], nrow = 1, ncol = 3)))))
    
    
    ############ ADE######
    if(!is.null(tr)){
      powl <- vector("numeric", length = length(tr))
      for(i in 1:length(tr)) powl[i] <-  lambda^i
      tr_G <- n +sum(tr * powl)
      powl <- vector("numeric", length = (length(tr)-1))
      for(i in 0:(length(tr)-1)) powl[i] <-  lambda^i
      powl <- c(0,powl)
      tr_H <- sum(tr * powl)
      powl <- vector("numeric", length = (length(tr)-2))
      for(i in 1:(length(tr)-2)) powl[i] <-  lambda^i
      powl <- c(0,1,powl)
      scl <- 0:(length(tr)-1)
      tr_H2 <- sum(tr * scl*powl)
      scl <- 2:length(tr)
      powl <- vector("numeric", length = (length(tr)-1))
      for(i in 1:length(powl)) powl[i] <-  lambda^i
      dv_l <- sum(tr * c(0,powl)*c(0,scl))
      
    }
    else{
      tr_G <- sum(1/(1-lambda*ev))
      tr_H <- sum(ev/(1-lambda*ev))
      tr_H2 <- sum(ev^2/(1-lambda*ev)^2)
      dv_l <- sum(ev/(1-lambda*ev)^2)
    }
    
    ADE <- c()
    for(i in 1:p1) ADE <- c(ADE, Re((1/n)*P[i,1]*tr_G) + Re((1/n)*P[i,2]*tr_H)) 
    
    p_l_ADE <- c()
    for(i in 1:p1) p_l_ADE <- c(p_l_ADE, Re((1/n)*P[i,1] * dv_l) + Re((1/n)*P[i,2] * tr_H2))
    p_b_ADE <- Re((1/n)*tr_G)
    p_g_ADE <- Re((1/n)*tr_H)
    der_ADE <- matrix(cbind(rep(p_b_ADE, p1), rep(p_g_ADE, p1), p_l_ADE), nrow = p1, ncol = 3)
    
    se_ADE <- c()
    for(i in 1:p1) se_ADE <- c(se_ADE, sqrt(as.numeric(matrix(der_ADE[i,], nrow = 1, ncol = 3) %*% 
                                                         Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                         t(matrix(der_ADE[i,], nrow = 1, ncol = 3)))))
    
    ######### AIE ######
    AIE <- ATE - ADE
    der_AIE <- der_ATE - der_ADE
    se_AIE <- c()
    for(i in 1:p1) se_AIE <- c(se_AIE, sqrt(as.numeric(matrix(der_AIE[i,], nrow = 1, ncol = 3) %*% 
                                                         Sigma[c(i,i+p1,pl),c(i,i+p1,pl)] %*% 
                                                         t(matrix(der_AIE[i,], nrow = 1, ncol = 3)))))
    
    #############################################
    
  }
  
  if(prt){
  
  cat("Impact Measures (lag, KP_formula):\n")
  tb <- cbind(ADE, AIE, ATE)
  colnames(tb) <- c("Direct", "Indirect", "Total")
  rownames(tb) <- bnames
   print(tb)
  cat("========================================================\n")
  cat("Results based on Kelejian and Piras formula:\n")
  cat("========================================================\n")
  cat("Analytical standard errors\n")
  se <- cbind(se_ADE, se_AIE, se_ATE)
  colnames(se) <- c("Direct", "Indirect", "Total")
  rownames(se) <- bnames
  print(se)
  cat("\nAnalytical z-values:\n")
  mat <- cbind(ADE, AIE, ATE)/se
  colnames(mat) <- c("Direct", "Indirect", "Total")
  rownames(mat) <- bnames
  print(mat)
  cat("\nAnalytical p-values:\n")
  xx <- apply(2*(1-pnorm(abs(mat))), 2, format.pval)
  if (length(bnames) == 1L) {
    xx <- matrix(xx, ncol=3)
    colnames(xx) <- c("Direct", "Indirect", "Total")
  }
  rownames(xx) <- bnames
  print(xx, quote=FALSE)
  }
  else{
    tb <- cbind(ADE, AIE, ATE)
    colnames(tb) <- c("Direct", "Indirect", "Total")
    rownames(tb) <- bnames
    se <- cbind(se_ADE, se_AIE, se_ATE)
    colnames(se) <- c("Direct", "Indirect", "Total")
    rownames(se) <- bnames
    mat <- cbind(ADE, AIE, ATE)/se
    colnames(mat) <- c("Direct", "Indirect", "Total")
    rownames(mat) <- bnames
    xx <- apply(2*(1-pnorm(abs(mat))), 2, format.pval)
    if (length(bnames) == 1L) {
      xx <- matrix(xx, ncol=3)
      colnames(xx) <- c("Direct", "Indirect", "Total")
    }
    rownames(xx) <- bnames
  }
  results <- list(tb, se, mat, xx)
  invisible(results)
  }
