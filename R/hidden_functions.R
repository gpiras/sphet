sarargmm <- function(formula, data, listw, listw2, endog, 
                     instruments, lag.instr, initial.value, 
                     het, verbose, na.action,
                     step1.c, control, HAC, cl, Durbin = NULL){
  
  mt <- terms(formula, data = data)
  mf <- lm(formula, data, na.action = na.action, method = "model.frame")
  na.act <- attr(mf, 'na.action')
  
  y <- c(model.extract(mf, "response"))
  x <- model.matrix(mt, mf)
  
  if (length(y) != nrow(x)) 
    stop("x and y have different length")
  
  if (any(is.na(y))) 
    stop("NAs in dependent variable")
  if (any(is.na(x))) 
    stop("NAs in independent variable")
  
  n <- nrow(x)
  k <- ncol(x)	
  xcolnames <- colnames(x)
  #print(xcolnames)
  K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1] == 1), 2, 1)
  
  if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
  if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
  if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
  if(inherits(listw,"Matrix"))  Ws <- listw	
  
  
  if (nrow(x) != nrow(Ws)) stop("Input data and weights have different dimension")
  
  if(!is.null(listw2)) {
     
    if(!inherits(listw2,c("listw", "Matrix", "matrix"))) stop("listw2 format unknown")
    if(inherits(listw2,"listw"))  Ws2<-listw2dgCMatrix(listw2)	
    if(inherits(listw2,"matrix"))  Ws2<-Matrix(listw2)	
    
    if(identical(listw, listw2)) Ws2 <- Ws
    twow <- TRUE    
  }
  else {
    Ws2 <- Ws
    twow <- FALSE
  }
  
  
  if(Durbin == TRUE | class(Durbin) == "formula"  ){
    if(class(Durbin) == "formula"){
      ## For each value of K (i.e, 1 or 2) there are four cases:
      # 1) formula = y ~ x1 + x2 + x3, Durbin = ~ x2 + x3 (only a subset of the x's are lagged and all of them are also in the formula)
      # 2) formula = y ~ x1 + x2 + x3, Durbin = ~ x1 + x2 + x3  (same as Durbin = TRUE)
      # 3) formula = y ~ x1 + x2 + x3, Durbin = ~ x4 (only variables that show up only lagged)
      # 4) formula = y ~ x1 + x2 + x3, Durbin = ~ x3 + x4 (one variable is in formula and one not)
      xdur <- as.matrix(lm(Durbin, data, na.action=na.action, method="model.frame"))
      if(sum(match(xcolnames, colnames(xdur), nomatch = 0))==0  && K != 2 && k != 1){
        
        wxdur <- as.matrix(Ws %*% xdur)
        wwxdur <- as.matrix(Ws %*% wxdur)
        wwwxdur <- as.matrix(Ws %*% wwxdur)
        
        if (k > 1 || (k == 1 && K == 1)) {
          
          wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
          for (i in K:k) {
            Wx <- Ws %*% x[, i]
            wx[, (i - (K - 1))] <- as.matrix(Wx)
            
          }
          wwx <- as.matrix(Ws %*% wx)  
        } 
        
        if(K==2){
          Hin <- cbind(x, wx, wwx, wxdur, wwxdur, wwwxdur)
          x <- cbind(x, wxdur)  
          colnames(x) <-  c(xcolnames, paste("lag_", colnames(xdur), sep=""))
          
        } 
        else {
          Hin <- cbind(1, x, wx, wwx, wxdur, wwxdur, wwwxdur)
          x <- cbind(x, wxdur)
          colnames(x) <-  c(xcolnames, paste("lag_", colnames(xdur), sep=""))
          
        }
        
        if(!is.null(listw2)) {
          w2H <- as.matrix(Ws2 %*% Hin[,-1])            
          Hin <- cbind(Hin, w2H)
        }
        
      }  
      else{
        #no intercept
        if(K==1){
          pos.xdur <- which(colnames(xdur) %in% xcolnames)
          pos.x <- which(xcolnames %in% colnames(xdur))
          
          if(all(is.na(match(colnames(xdur), xcolnames)))) onx <- as.matrix(x)
          else onx <- as.matrix(x[,-pos.x])
          
          if(dim(onx)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(onx) <- nmonx <-  xcolnames
            else names(onx) <- nmonx <-  xcolnames[-pos.x] 
          } 
          else nmonx <- NULL
          
          if(all(is.na(match(colnames(xdur), xcolnames)))) onxl <- as.matrix(xdur)
          else onxl <- as.matrix(xdur[, -pos.xdur])
          
          if(dim(onxl)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(onxl) <- nmonxl <- colnames(xdur)
            else names(onxl) <- nmonxl <- colnames(xdur)[-pos.xdur]  
          } 
          else nmonxl <- NULL
          
          onboth <- as.matrix(xdur[,pos.xdur])
          if(dim(onboth)[2] != 0){
            if(all(is.na(match(colnames(xdur), xcolnames))))  names(onboth) <- nmonb <- colnames(xdur)
            else  names(onboth) <- nmonb <- colnames(xdur)[pos.xdur] 
          }
          else nmonb <- NULL
          
          wonx <- as.matrix(Ws %*% onx)
          wwonx <- as.matrix(Ws %*% wonx)
          
          wonxl <- as.matrix(Ws %*% onxl)
          if(dim(wonxl)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur),sep = "") 
            else names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur)[-pos.xdur],sep = "") 
          } 
          else nmwonxl <- NULL
          wwonxl <- as.matrix(Ws %*% wonxl)
          wwwonxl <- as.matrix(Ws %*% wwonxl)
          
          wonboth <- as.matrix(Ws %*% onboth)
          if(dim(wonboth)[2]!=0) {
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur),sep = "") 
            else names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur)[pos.xdur],sep = "") 
          }
          else nmwonb <- NULL
          wwonboth <- as.matrix(Ws %*% wonboth)
          wwwonboth <- as.matrix(Ws %*% wwonboth)
          
          
          x <- cbind(onx, onboth, wonboth, wonxl)
          colnames(x) <- c(nmonx, nmonb, nmwonb, nmwonxl)
          Hin <- cbind(1, x, wonx, wwonx, wwonboth, wwwonboth,wwonxl,wwwonxl )
          
          
          if(!is.null(listw2)){ 
            
            w2H <- Ws2 %*% Hin[,-1]          
            Hin <- cbind(Hin, w2H)                  	          
            
          }
          
        }
        
        else{
          
          if(k !=1){
            pos.xdur <- which(colnames(xdur) %in% xcolnames)
            pos.x <- which(xcolnames %in% colnames(xdur))
            
            if(all(is.na(match(colnames(xdur), xcolnames)))) onx <- as.matrix(x)
            else onx <- as.matrix(x[,-pos.x]) 
            if(dim(onx)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onx) <- nmonx <-  xcolnames
              else names(onx) <- nmonx <-  xcolnames[-pos.x]
            }  
            else nmonx <- NULL
            
            
            if(all(is.na(match(colnames(xdur), xcolnames)))) onxl <- as.matrix(xdur)
            else  onxl <- as.matrix(xdur[, -pos.xdur])
            if(dim(onxl)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onxl) <- nmonxl <- colnames(xdur)
              else names(onxl) <- nmonxl <- colnames(xdur)[-pos.xdur]
            }   
            else nmonxl <- NULL
            
            onboth <- as.matrix(xdur[,pos.xdur])
            if(dim(onboth)[2] != 0) {
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onboth) <- nmonb <- colnames(xdur)
              else names(onboth) <- nmonb <- colnames(xdur)[pos.xdur] 
            } 
            else nmonb <- NULL
            wonx <- as.matrix(Ws %*% onx[,-1])
            wwonx <- as.matrix(Ws %*% wonx)
            
            wonxl <- as.matrix(Ws %*% onxl)
            if(dim(wonxl)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur),sep = "") 
              else names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur)[-pos.xdur],sep = "") 
            } 
            else nmwonxl <- NULL
            wwonxl <- as.matrix(Ws %*% wonxl)
            wwwonxl <- as.matrix(Ws %*% wwonxl)
            
            wonboth <- as.matrix(Ws %*% onboth)
            if(dim(wonboth)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur),sep = "") 
              else names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur)[pos.xdur],sep = "") 
            } 
            else nmwonb <- NULL
            wwonboth <- as.matrix(Ws %*% wonboth)
            wwwonboth <- as.matrix(Ws %*% wwonboth)
            
            
            x <- cbind(onx, onboth, wonboth, wonxl)
            colnames(x) <- c(nmonx,nmonb,nmwonb,nmwonxl)
            Hin <- cbind(x, wonx, wwonx, wwonboth, wwwonboth,wwonxl,wwwonxl)
            
            if(!is.null(listw2)){ 
              
              w2H <- Ws2 %*% Hin[,-1]          
              Hin <- cbind(Hin, w2H)                  	          
              
            }

          }
          else{
            
            wxdur <-  as.matrix(Ws %*% xdur)
            wwxdur <- as.matrix(Ws %*% wxdur)
            wwwxdur <- as.matrix(Ws %*% wwxdur)
            x <- cbind(x, wxdur)
            colnames(x) <- c(xcolnames, paste("lag_",colnames(xdur),sep =""))
            Hin <- cbind(x, wwxdur, wwwxdur)
            
            if(!is.null(listw2)){ 
              
              w2H <- Ws2 %*% Hin[,-1]          
              Hin <- cbind(Hin, w2H)                  	          
              
            }
          }
        }
      }
    }
    
    else{
      
      if (k > 1 || (k == 1 && K == 1)) {
        
        wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
        for (i in K:k) {
          Wx <- Ws %*% x[, i]
          wx[, (i - (K - 1))] <- as.matrix(Wx)
        }
        wwx <- as.matrix(Ws %*% wx)  
        wwwx <- as.matrix(Ws  %*%  wwx)
      } 
      
      if(K==2){
        
        Hin <- cbind(x, wx, wwx, wwwx)
        x <- cbind(x, wx)  
        colnames(x) <-  c(xcolnames, paste("lag_", xcolnames[-1], sep=""))
        
      } 
      else {
        Hin <- cbind(1, x, wx, wwx, wwwx)
        x <- cbind(x, wx)  
        colnames(x) <-  c(xcolnames, paste("lag_", xcolnames, sep=""))
        
      }
      
      if(!is.null(listw2)) {
        
        w2H <- as.matrix(Ws2 %*% Hin[,-1])
        Hin <- cbind(Hin, w2H)
      }
    
    }
  }
  else{  
    
    if (k > 1 || (k == 1 && K == 1)) {
      wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
      for (i in K:k) {
        Wx <- Ws %*% x[, i]
        wx[, (i - (K - 1))] <- as.matrix(Wx)
      }
      wwx <- Ws %*% wx                    					         
    }
    
    if(K==2)    Hin <- cbind(x, wx, wwx)
    else        Hin <- cbind(1, x, wx, wwx)
    x <- x
    
    if(!is.null(listw2)) {
      
      w2H <- as.matrix(Ws2 %*% Hin[,-1])
      Hin <- cbind(Hin, w2H)
      
    }
    
  }
  
  
  wy <- Ws %*% y	
  colnames(wy) <- "lambda"
  
  if (!is.null(endog) && is.null(instruments)) stop("No instruments specified for the endogenous variable in the model")
  
  if (!is.null(endog)) {
    endog <- as.matrix(lm(endog, data, na.action = na.action, method = "model.frame"))
    instruments <- as.matrix(lm(instruments, data, na.action = na.action, method = "model.frame"))
    
    if(lag.instr) {
      winst <- Ws %*% instruments
      wwinst <- Ws %*% winst	
      if(twow){
        w2i <- Ws2 %*% instruments 
        w2wi <- Ws2 %*% winst 
        w2wwi <- Ws2 %*% wwinst 	
        AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst), as.matrix(w2i), as.matrix(w2wi),as.matrix(w2wwi))        
      }
      else  AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst))        
    }
    else  AddH <- instruments        
    
    Hmat <- cbind(Hin, AddH)
    Zmat <- cbind(x, endog, as.matrix(wy))            
    colnames(Zmat) <- c(colnames(x), colnames(endog), colnames(wy))               
  } 
  else {
    Zmat <- cbind(x, as.matrix(wy))
    colnames(Zmat) <- c(colnames(x), colnames(wy))
    Hmat <- Hin 
  }
  
  firststep <- spatial.ivreg(y = y , Zmat = Zmat, Hmat = Hmat, het = het, HAC = HAC)
  ubase <- residuals(firststep)
  
  if (initial.value == "SAR"){
    Wubase <- Ws %*% ubase
    pars <- coefficients(lm(as.numeric(ubase) ~ as.numeric(Wubase)-1))
  }
  else pars <-initial.value
  
  
  
  if(het){
    Ggmat <- gg_het(Ws2, ubase, n)
    optres <- nlminb(pars, optimfunct, lower = -0.9 + .Machine$double.eps , 
                     upper = 0.9 -  .Machine$double.eps, control = control, 
                     v = Ggmat, verbose = verbose)
    rhotilde <- optres$par
    
    if(step1.c){
      gmm.weghts1.c <- psirhorho_het(rhotilde, ubase, Hmat, Zmat, Ws2, step1.c = TRUE)
      optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, 
                       vcmat = gmm.weghts1.c$Phiinv, verbose = verbose, 
                       lower= -0.9 + .Machine$double.eps , 
                       upper= 0.9 -  .Machine$double.eps, control = control)	
      
      rhotilde <- optres$par
      gmm.weghts1.c <- psirhorho_het(rhotilde, ubase, Hmat, Zmat, Ws2, step1.c = TRUE)
      vcmat_2sls <- Omega_het(rhotilde, gmm.weghts1.c$Pmat, gmm.weghts1.c$A1, 
                              gmm.weghts1.c$A2, gmm.weghts1.c$a.vec1, 
                              gmm.weghts1.c$a.vec2, Hmat, Ggmat$bigG, 
                              gmm.weghts1.c$Phiinv, gmm.weghts1.c$epsilon, 
                              gmm.weghts1.c$Zstar, Ws2, step1.c = TRUE)
      coeff_2sls <- as.matrix(c(coefficients(firststep), rhotilde))
      rownames(coeff_2sls) <- c(colnames(Zmat), 'rho')
      s2_2sls <- crossprod(ubase)/(n-k)
      model.data <- data.frame(cbind(y,x[,-1]))
      vc_mat <- vcmat_2sls$Omega
      rownames(vc_mat) <- colnames(vc_mat) <- rownames(coeff_2sls)
      method <- "gmm sarar"
      
      k <- nrow(coeff_2sls)
      R <- matrix(0,1,k)
      R[,((k-1):k)] <- 1
      Rbeta <- R %*% coeff_2sls
      Rvar <- R %*% vcmat_2sls$Omega %*% t(R)
      stat <- as.numeric(t(Rbeta) %*% Rbeta/Rvar)
      pval <- pchisq(stat, df = 1, lower.tail = FALSE)
      W <- list(stat = stat, pval = pval)
      results_2sls <- list(coefficients = coeff_2sls, var = vc_mat, s2 = s2_2sls, 
                           call = cl, residuals = as.numeric(ubase), model = model.data, 
                           method = method, W = W, firststep = firststep$coefficients, 
                           init.rho = rhotilde)
      class(results_2sls) <- c("sphet", "sarar", "gstsls")
    }
    
  }
  else{
    
    Ggmat <- gg_hom(Ws2, ubase, n)
    optres <- nlminb(pars, optimfunct, lower= -0.9 + .Machine$double.eps , 
                     upper= 0.9 -  .Machine$double.eps, control = control, 
                     v = Ggmat, verbose = verbose)
    rhotilde <- optres$par
    
  }
  
  yt  <- y - rhotilde * Ws2 %*% y
  wZmat <- Ws2 %*% Zmat
  Zt <- Zmat - rhotilde * wZmat
  
  secondstep <- spatial.ivreg(y =yt , Zmat = Zt, Hmat = Hmat, het = het, HAC = HAC)
  delta <- coefficients(secondstep)
  utildeb <- y - Zmat %*% delta
  
  if(het){
    
    Ggmat <- gg_het(Ws2, utildeb, n)
    gmm.weghts <- psirhorho_het(rhotilde, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)
    optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat= gmm.weghts$Phiinv, 
                     verbose = verbose, lower= -0.9 + .Machine$double.eps , 
                     upper= 0.9 -  .Machine$double.eps, control = control)	
    rhofin<-optres$par
    gmm.weghts <- psirhorho_het(rhofin, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)
    
    vcmat <- Omega_het(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, 
                       gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, 
                       Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, 
                       gmm.weghts$Zstar, Ws2, step1.c = FALSE)
    
  }
  else{
    
    Ggmat<-gg_hom(Ws2, utildeb, n)
    
    gmm.weghts <- psirhorho_hom(rhotilde, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )
    optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat = gmm.weghts$Phiinv, 
                     verbose = verbose, lower= -0.9 + .Machine$double.eps , 
                     upper= 0.9 -  .Machine$double.eps, control = control)	
    
    rhofin <- optres$par
    gmm.weghts <- psirhorho_hom(rhofin, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )
  
    vcmat <- Omega_hom(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, 
                       gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, 
                       gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar)
  }
  
  vc_Omega <- vcmat$Omega
  
  coeff <- as.matrix(c(as.numeric(delta), rhofin))
  rownames(coeff) <- c(colnames(Zmat), 'rho')
  rownames(vc_Omega) <- colnames(vc_Omega) <- rownames(coeff)
  s2 <- crossprod(utildeb)/(n-k)
  
  model.data<-data.frame(cbind(y,x[,-1]))
  
  method<-"gmm sarar"
  
  k <- nrow(coeff)
  R <- matrix(0,1,k)
  R[,((k-1):k)] <- 1
  Rbeta <- R %*% coeff
  Rvar <- R %*% vcmat$Omega %*% t(R)
  stat <- as.numeric(t(Rbeta) %*% Rbeta/Rvar)
  pval <- pchisq(stat, df = 1, lower.tail = FALSE)
  W <- list(stat = stat, pval = pval)
  
  
  
  if(het && step1.c) results <- list(coefficients = coeff, var = vc_Omega, 
                                     s2 = s2, call = cl, residuals = as.numeric(utildeb), 
                                     model = model.data, method = method, W = W, firststep = firststep$coefficients, 
                                     init.rho = rhotilde,  twosls = results_2sls, Durbin = Durbin, endog = endog)
  else  results <- list(coefficients = coeff, var = vc_Omega, s2 = s2, call = cl, 
                        residuals = as.numeric(utildeb), model = model.data, method = method, W = W, 
                        firststep = firststep$coefficients, init.rho = rhotilde, Durbin = Durbin, endog = endog)
  
  results$listw <- Ws
  
  class(results)<-c("sphet", "sarar_gmm", "gstsls") #remember to change to sarar when impacts will be developed
  return(results)
}

laggmm <- function(formula, data, listw, listw2, endog, 
                   instruments, lag.instr, 
                   het, verbose, na.action, HAC, cl, Durbin = NULL){
  
  
  mt <- terms(formula,data = data)
  mf <- lm(formula, data, na.action = na.action, method = "model.frame")
  na.act <- attr(mf, 'na.action')
  
  y <- c(model.extract(mf, "response"))
  x <- model.matrix(mt,mf)

  if (length(y)!=nrow(x)) 
    stop("x and y have different length")
  
  if (any(is.na(y))) 
    stop("NAs in dependent variable")
  if (any(is.na(x))) 
    stop("NAs in independent variable")
  
  n <- nrow(x)
  k <- ncol(x)	
  xcolnames <- colnames(x)
  K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)

  if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
  if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
  if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
  if(inherits(listw,"Matrix"))  Ws <- listw	
  
  if (nrow(x) != nrow(Ws))
    stop("Input data and weights have different dimension")
  
  if(Durbin == TRUE | class(Durbin) == "formula"  ){
    if(class(Durbin) == "formula"){
      ## For each value of K (i.e, 1 or 2) there are four cases:
      # 1) formula = y ~ x1 + x2 + x3, Durbin = ~ x2 + x3 (only a subset of the x's are lagged and all of them are also in the formula)
      # 2) formula = y ~ x1 + x2 + x3, Durbin = ~ x1 + x2 + x3  (same as Durbin = TRUE)
      # 3) formula = y ~ x1 + x2 + x3, Durbin = ~ x4 (only variables that show up only lagged)
      # 4) formula = y ~ x1 + x2 + x3, Durbin = ~ x3 + x4 (one variable is in formula and one not)
      xdur <- as.matrix(lm(Durbin, data, na.action=na.action, method="model.frame"))
      
      if(sum(match(xcolnames, colnames(xdur), nomatch = 0))==0  && K != 2 && k != 1){
       
        wxdur <- as.matrix(Ws %*% xdur)
        wwxdur <- as.matrix(Ws %*% wxdur)
        wwwxdur <- as.matrix(Ws %*% wwxdur)
        
        if (k > 1 || (k == 1 && K == 1)) {
          
          wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
          for (i in K:k) {
            Wx <- Ws %*% x[, i]
            wx[, (i - (K - 1))] <- as.matrix(Wx)
            
          }
          wwx <- as.matrix(Ws %*% wx)  
        } 
      
        if(K==2){
         # Hin <- cbind(x, wx, wwx, wwwx, wxdur, wwxdur, wwwxdur)
          Hin <- cbind(x, wx, wwx, wxdur, wwxdur, wwwxdur)
          x <- cbind(x, wxdur)  
          colnames(x) <-  c(xcolnames, paste("lag_", colnames(xdur), sep=""))
          # print(head(x))
          # print(head(Hin))
        } 
        else {
          #Hin <- cbind(1, x, wx, wwx, wwwx, wxdur, wwxdur, wwwxdur)
          Hin <- cbind(1, x, wx, wwx, wxdur, wwxdur, wwwxdur)
          x <- cbind(x, wxdur)  
          colnames(x) <-  c(xcolnames, paste("lag_", colnames(xdur), sep=""))
          # print(head(x))
          # print(head(Hin))
        }
      }  
      else{
       
        if(K==1 ){
       
           pos.xdur <- which(colnames(xdur) %in% xcolnames)
           pos.x <- which(xcolnames %in% colnames(xdur))
        
           if(all(is.na(match(colnames(xdur), xcolnames)))) onx <- as.matrix(x)
           else onx <- as.matrix(x[,-pos.x])
          
           if(dim(onx)[2]!=0){
             if(all(is.na(match(colnames(xdur), xcolnames)))) names(onx) <- nmonx <-  xcolnames
               else names(onx) <- nmonx <-  xcolnames[-pos.x] 
           } 
           else nmonx <- NULL
           
         if(all(is.na(match(colnames(xdur), xcolnames)))) onxl <- as.matrix(xdur)
           else onxl <- as.matrix(xdur[, -pos.xdur])
           
           if(dim(onxl)[2]!=0){
             if(all(is.na(match(colnames(xdur), xcolnames)))) names(onxl) <- nmonxl <- colnames(xdur)
               else names(onxl) <- nmonxl <- colnames(xdur)[-pos.xdur]  
           } 
           else nmonxl <- NULL
           
           onboth <- as.matrix(xdur[,pos.xdur])
           if(dim(onboth)[2] != 0){
             if(all(is.na(match(colnames(xdur), xcolnames))))  names(onboth) <- nmonb <- colnames(xdur)
                else  names(onboth) <- nmonb <- colnames(xdur)[pos.xdur] 
           }
           else nmonb <- NULL
            
           wonx <- as.matrix(Ws %*% onx)
           wwonx <- as.matrix(Ws %*% wonx)
           
           wonxl <- as.matrix(Ws %*% onxl)
           if(dim(wonxl)[2]!=0){
             if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur),sep = "") 
               else names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur)[-pos.xdur],sep = "") 
           } 
           else nmwonxl <- NULL
           wwonxl <- as.matrix(Ws %*% wonxl)
           wwwonxl <- as.matrix(Ws %*% wwonxl)
           
           wonboth <- as.matrix(Ws %*% onboth)
           
           if(dim(wonboth)[2]!=0) {
             if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur),sep = "") 
               else names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur)[pos.xdur],sep = "") 
             }
           else nmwonb <- NULL
           wwonboth <- as.matrix(Ws %*% wonboth)
           wwwonboth <- as.matrix(Ws %*% wwonboth)
             
           
         x <- cbind(onx, onboth, wonboth, wonxl)
         colnames(x) <- c(nmonx,nmonb,nmwonb,nmwonxl)
         Hin <- cbind(1, x, wonx, wwonx, wwonboth, wwwonboth,wwonxl,wwwonxl )
         # print(head(x))
         # print(head(Hin))
        } 
        else{
          
          if(k !=1){
            
            pos.xdur <- which(colnames(xdur) %in% xcolnames)
            pos.x <- which(xcolnames %in% colnames(xdur))
            
            if(all(is.na(match(colnames(xdur), xcolnames)))) onx <- as.matrix(x)
            else onx <- as.matrix(x[,-pos.x]) 
          
            if(dim(onx)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onx) <- nmonx <-  xcolnames
              else names(onx) <- nmonx <-  xcolnames[-pos.x]
            }  
            else nmonx <- NULL
          
            if(all(is.na(match(colnames(xdur), xcolnames)))) onxl <- as.matrix(xdur)
            else  onxl <- as.matrix(xdur[, -pos.xdur])
            
            if(dim(onxl)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onxl) <- nmonxl <- colnames(xdur)
              else names(onxl) <- nmonxl <- colnames(xdur)[-pos.xdur]
            }   
            else nmonxl <- NULL
          
              onboth <- as.matrix(xdur[,pos.xdur])
            if(dim(onboth)[2] != 0) {
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onboth) <- nmonb <- colnames(xdur)
                else names(onboth) <- nmonb <- colnames(xdur)[pos.xdur] 
            } 
            else nmonb <- NULL
           
            wonx <- as.matrix(Ws %*% onx[,-1])
            wwonx <- as.matrix(Ws %*% wonx)
            
            wonxl <- as.matrix(Ws %*% onxl)
            if(dim(wonxl)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur),sep = "") 
                else names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur)[-pos.xdur],sep = "") 
            } 
            else nmwonxl <- NULL
            wwonxl <- as.matrix(Ws %*% wonxl)
            wwwonxl <- as.matrix(Ws %*% wwonxl)
            
            wonboth <- as.matrix(Ws %*% onboth)
            
            if(dim(wonboth)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur),sep = "") 
                else names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur)[pos.xdur],sep = "") 
            } 
            else nmwonb <- NULL
            wwonboth <- as.matrix(Ws %*% wonboth)
            wwwonboth <- as.matrix(Ws %*% wwonboth)
            
            
            x <- cbind(onx, onboth, wonboth, wonxl)
            colnames(x) <- c(nmonx,nmonb,nmwonb,nmwonxl)
            Hin <- cbind(x, wonx, wwonx, wwonboth, wwwonboth,wwonxl,wwwonxl )
            # print(head(x))
            # print(head(Hin))
                      }
        else{
      
          wxdur <-  as.matrix(Ws %*% xdur)
          wwxdur <- as.matrix(Ws %*% wxdur)
          wwwxdur <- as.matrix(Ws %*% wwxdur)
          x <- cbind(x, wxdur)
          colnames(x) <- c(xcolnames, paste("lag_",colnames(xdur),sep =""))
          Hin <- cbind(x, wwxdur, wwwxdur)
          # print(head(x))
          # print(head(Hin))
        }
        }
      }
    
    }
    
    else{
      
      if (k > 1 || (k == 1 && K == 1)) {
        
        wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
        for (i in K:k) {
          Wx <- Ws %*% x[, i]
          wx[, (i - (K - 1))] <- as.matrix(Wx)
        }
        wwx <- as.matrix(Ws %*% wx)  
        wwwx <- as.matrix(Ws  %*%  wwx)
        #wwwwx <- as.matrix
         
      } 
      
      if(K==2){
        Hin <- cbind(x,wx,wwx,wwwx)
        x <- cbind(x, wx)  
        colnames(x) <-  c(xcolnames, paste("lag_", xcolnames[-1], sep=""))
        # print(head(x))
         #print(head(Hin))
      } 
      else {
        Hin <- cbind(1,x,wx,wwx,wwwx)
        x <- cbind(x, wx)  
        colnames(x) <-  c(xcolnames, paste("lag_", xcolnames, sep=""))
        # print(head(x))
         #print(head(Hin))
      }
    }
  }
  else{  
   
     if (k > 1 || (k == 1 && K == 1)) {
      wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
      for (i in K:k) {
        Wx <- Ws %*% x[, i]
        wx[, (i - (K - 1))] <- as.matrix(Wx)
      }
      wwx <- Ws %*% wx                    					         
    }
    
    if(K==2)    Hin <- cbind(x, wx, wwx)
    else        Hin <- cbind(1, x, wx, wwx)
    x <- x
    # print(head(x))
    # print(head(Hin))
  }
  
  
  wy <- Ws %*% y	
  colnames(wy)<-"lambda"
  if (!is.null(endog) && is.null(instruments)) stop("No instruments specified for the endogenous variable in the model")
  
  if (!is.null(endog)) {
    
    endog <- as.matrix(lm(endog, data, na.action=na.action, method="model.frame"))
    instruments <- as.matrix(lm(instruments, data, na.action=na.action, method="model.frame"))
    
    if(lag.instr) {
      winst <- Ws %*% instruments
      wwinst<- Ws %*% winst	
      AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst))        
    }
    else  AddH <- instruments        
    
    Hmat <- cbind(Hin, AddH)
    Zmat<- cbind(x, endog, as.matrix(wy))            
    colnames(Zmat) <- c(colnames(x), colnames(endog), colnames(wy))               
  }
  else {
    Zmat<- cbind(x, as.matrix(wy))                    
    Hmat <- Hin
  }
  
  results <-spatial.ivreg(y, Zmat, Hmat, het, HAC)
  #print(results$coefficients)
  #print(results$var)
  model.data <- data.frame(cbind(y, x[, -1]))
  results$call <- cl
  results$listw <- Ws
  results$model <- model.data
  results$type <- NULL
  results$bandwidth <- NULL
  results$method <- "gmm lag"
  results$HAC <- FALSE
  results$Durbin <- Durbin
  results$endog <- endog
  class(results) <- c("sphet", "lag_gmm", "stsls_sphet") #change to lag gmm
  
  return(results)
  
}

errorgmm <- function(formula, data, listw, listw2, endog, 
                     instruments, lag.instr, initial.value, 
                     het, verbose, na.action,
                     step1.c, control, HAC, cl, Durbin = NULL){
  
  mt <- terms(formula,data=data)
  mf <- lm(formula, data, na.action = na.action, method = "model.frame")
  na.act <- attr(mf,'na.action')
  
  y <- c(model.extract(mf,"response"))
  x <- model.matrix(mt,mf)
  
  if (length(y)!=nrow(x)) 
    stop("x and y have different length")
  
  if (any(is.na(y))) 
    stop("NAs in dependent variable")
  if (any(is.na(x))) 
    stop("NAs in independent variable")
  
  n <- nrow(x)
  k <- ncol(x)	
  xcolnames <- colnames(x)
  
  K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)
  
  if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
  if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
  if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
  if(inherits(listw,"Matrix"))  Ws <- listw	
  
  if (nrow(x) != nrow(Ws))
    stop("Input data and weights have different dimension")
  
  if(Durbin == TRUE | class(Durbin) == "formula"  ){
    if(class(Durbin) == "formula"){
      ## For each value of K (i.e, 1 or 2) there are four cases:
      # 1) formula = y ~ x1 + x2 + x3, Durbin = ~ x2 + x3 (only a subset of the x's are lagged and all of them are also in the formula)
      # 2) formula = y ~ x1 + x2 + x3, Durbin = ~ x1 + x2 + x3  (same as Durbin = TRUE)
      # 3) formula = y ~ x1 + x2 + x3, Durbin = ~ x4 (only variables that show up only lagged)
      # 4) formula = y ~ x1 + x2 + x3, Durbin = ~ x3 + x4 (one variable is in formula and one not)
      xdur <- as.matrix(lm(Durbin, data, na.action=na.action, method="model.frame"))

            if(sum(match(xcolnames, colnames(xdur), nomatch = 0))==0  && K != 2 && k != 1){
      
                wxdur <- as.matrix(Ws %*% xdur)
                wwxdur <- as.matrix(Ws %*% wxdur)
        
        if (k > 1 || (k == 1 && K == 1)) {
          
          wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
          for (i in K:k) {
            Wx <- Ws %*% x[, i]
            wx[, (i - (K - 1))] <- as.matrix(Wx)
          }
          
        } 
        
        if(K==2){
          Hin <- cbind(x, wx, wxdur, wwxdur)
          x <- cbind(x, wxdur)  
          colnames(x) <-  c(xcolnames, paste("lag_", colnames(xdur), sep=""))
          # print(head(x))
          # print(head(Hin))
        } 
        else {
          Hin <- cbind(1, x, wx, wxdur, wwxdur)
          x <- cbind(x, wxdur)  
          colnames(x) <-  c(xcolnames, paste("lag_", colnames(xdur), sep=""))
          # print(head(x))
          # print(head(Hin))
        }
      }  
      
      else{
        
        if(K==1 ){
          
          pos.xdur <- which(colnames(xdur) %in% xcolnames)
          pos.x <- which(xcolnames %in% colnames(xdur))
          
          if(all(is.na(match(colnames(xdur), xcolnames)))) onx <- as.matrix(x)
          else onx <- as.matrix(x[,-pos.x])
          
          if(dim(onx)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(onx) <- nmonx <-  xcolnames
            else names(onx) <- nmonx <-  xcolnames[-pos.x] 
          } 
          else nmonx <- NULL
         
          if(all(is.na(match(colnames(xdur), xcolnames)))) onxl <- as.matrix(xdur)
          else onxl <- as.matrix(xdur[, -pos.xdur])
          
          if(dim(onxl)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(onxl) <- nmonxl <- colnames(xdur)
            else names(onxl) <- nmonxl <- colnames(xdur)[-pos.xdur]  
          } 
          else nmonxl <- NULL
          
          onboth <- as.matrix(xdur[,pos.xdur])
          if(dim(onboth)[2] != 0){
            if(all(is.na(match(colnames(xdur), xcolnames))))  names(onboth) <- nmonb <- colnames(xdur)
            else  names(onboth) <- nmonb <- colnames(xdur)[pos.xdur] 
          }
          else nmonb <- NULL
         
          wonx <- as.matrix(Ws %*% onx)
          
          wonxl <- as.matrix(Ws %*% onxl)
          if(dim(wonxl)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur),sep = "") 
            else names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur)[-pos.xdur],sep = "") 
          } 
          else nmwonxl <- NULL
          wwonxl <- as.matrix(Ws %*% wonxl)
          
          wonboth <- as.matrix(Ws %*% onboth)
        
          if(dim(wonboth)[2]!=0) {
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur),sep = "") 
            else names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur)[pos.xdur],sep = "") 
          }
          else nmwonb <- NULL
          wwonboth <- as.matrix(Ws %*% wonboth)
          
          
          x <- cbind(onx, onboth, wonboth, wonxl)
          colnames(x) <- c(nmonx,nmonb,nmwonb,nmwonxl)
          Hin <- cbind(1, x, wonx, wwonboth, wwonxl )
          # print(head(x))
          # print(head(Hin))
          # 
          
        } 
        else{
         
           if(k !=1){
            
            pos.xdur <- which(colnames(xdur) %in% xcolnames)
            pos.x <- which(xcolnames %in% colnames(xdur))
            
            if(all(is.na(match(colnames(xdur), xcolnames)))) onx <- as.matrix(x)
            else onx <- as.matrix(x[,-pos.x]) 
            
            if(dim(onx)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onx) <- nmonx <-  xcolnames
              else names(onx) <- nmonx <-  xcolnames[-pos.x]
            }  
            else nmonx <- NULL
            
            if(all(is.na(match(colnames(xdur), xcolnames)))) onxl <- as.matrix(xdur)
            else  onxl <- as.matrix(xdur[, -pos.xdur])
            
            if(dim(onxl)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onxl) <- nmonxl <- colnames(xdur)
              else names(onxl) <- nmonxl <- colnames(xdur)[-pos.xdur]
            }   
            else nmonxl <- NULL
           
             onboth <- as.matrix(xdur[,pos.xdur])
            if(dim(onboth)[2] != 0) {
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onboth) <- nmonb <- colnames(xdur)
              else names(onboth) <- nmonb <- colnames(xdur)[pos.xdur] 
            } 
            else nmonb <- NULL
             
            wonx <- as.matrix(Ws %*% onx[,-1])
            
            wonxl <- as.matrix(Ws %*% onxl)
            if(dim(wonxl)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur),sep = "") 
              else names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur)[-pos.xdur],sep = "") 
            } 
            else nmwonxl <- NULL
            wwonxl <- as.matrix(Ws %*% wonxl)
            
            wonboth <- as.matrix(Ws %*% onboth)
           
            if(dim(wonboth)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur),sep = "") 
              else names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur)[pos.xdur],sep = "") 
            } 
            else nmwonb <- NULL
            wwonboth <- as.matrix(Ws %*% wonboth)
            
            
            x <- cbind(onx, onboth, wonboth, wonxl)
            colnames(x) <- c(nmonx,nmonb,nmwonb,nmwonxl)
            Hin <- cbind(x, wonx, wwonboth, wwonxl )
            # print(head(x))
            # print(head(Hin))
          }
          else{
           
            wxdur <-  as.matrix(Ws %*% xdur)
            wwxdur <- as.matrix(Ws %*% wxdur)
           
            x <- cbind(x, wxdur)
            colnames(x) <- c(xcolnames, paste("lag_",colnames(xdur),sep =""))
            Hin <- cbind(x, wwxdur)
            # print(head(x))
            # print(head(Hin))
          }
        }
      }
      
    }
    
    else{
      
      if (k > 1 || (k == 1 && K == 1)) {
        
        wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
        for (i in K:k) {
          Wx <- Ws %*% x[, i]
          wx[, (i - (K - 1))] <- as.matrix(Wx)
        }
        wwx <- as.matrix(Ws %*% wx)  
        
      } 
      
      if(K==2){
        Hin <- cbind(x,wx,wwx)
        x <- cbind(x, wx)  
        colnames(x) <-  c(xcolnames, paste("lag_", xcolnames[-1], sep=""))
        # print(head(x))
        # print(head(Hin))
      } 
      else {
        Hin <- cbind(1,x,wx,wwx)
        x <- cbind(x, wx)  
        colnames(x) <-  c(xcolnames, paste("lag_", xcolnames, sep=""))
        # print(head(x))
        # print(head(Hin))
      }
    }
  }
  else{  
    
    if (k > 1 || (k == 1 && K == 1)) {
      wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
      for (i in K:k) {
        Wx <- Ws %*% x[, i]
        wx[, (i - (K - 1))] <- as.matrix(Wx)
      }
     
    }
    
    if(K==2)    Hin <- cbind(x, wx)
    else        Hin <- cbind(1, x, wx)
    x <- x
    # print(head(x))
    # print(head(Hin))
  }
  
  
  if (!is.null(endog) && is.null(instruments)) stop("No instruments specified for the endogenous variable in the model")
  
  if (!is.null(endog)) {
    endog <- as.matrix(lm(endog, data, na.action=na.action, method="model.frame"))
    instruments <- as.matrix(lm(instruments, data, na.action=na.action, method="model.frame"))
    if(lag.instr) { 
      winst <- Ws %*% instruments           
      AddH<- cbind(instruments, as.matrix(winst))
    }
    else AddH<- cbind(instruments)
    
    Zmat<- cbind(x, endog)            
    colnames(Zmat) <- c(colnames(x), colnames(endog)) 
    Hmat<-cbind(Hin, AddH)
    
  }
  else {
    
    
    Hmat <- Hin
    Zmat<- as.matrix(x)
    
  }
  
  firststep <- spatial.ivreg(y = y , Zmat = Zmat, Hmat = Hmat, HAC = HAC, het = het)
  ubase <- residuals(firststep)
  
  if (initial.value=="SAR"){
    Wubase<-Ws %*% ubase
    pars<-coefficients(lm(as.numeric(ubase) ~ as.numeric(Wubase)-1))
  }
  else pars <- initial.value
  
 
  
  if(het){
    Ggmat <- gg_het(Ws, ubase, n)
  
    
    optres <- nlminb(pars, optimfunct, lower= -0.9 + .Machine$double.eps , 
                    upper= 0.9 -  .Machine$double.eps, control= control, 
                    v = Ggmat, verbose = verbose)
    rhotilde<-optres$par
    
    
    if(step1.c){
      gmm.weghts1.c<-psirhorho_het(rhotilde, ubase, Hmat, Zmat, Ws, step1.c = TRUE)
    
      optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat = gmm.weghts1.c$Phiinv, 
                       verbose = verbose, lower = -0.9 + .Machine$double.eps, 
                       upper = 0.9 -  .Machine$double.eps, control = control)	
      rhotilde <- optres$par
      gmm.weghts1.c <- psirhorho_het(rhotilde, ubase, Hmat, Zmat, Ws, step1.c = TRUE)
      vcmat_2sls <- Omega_het(rhotilde, gmm.weghts1.c$Pmat, gmm.weghts1.c$A1, 
                              gmm.weghts1.c$A2, gmm.weghts1.c$a.vec1, 
                              gmm.weghts1.c$a.vec2, Hmat, Ggmat$bigG, 
                              gmm.weghts1.c$Phiinv, gmm.weghts1.c$epsilon, 
                              gmm.weghts1.c$Zstar, Ws, step1.c = TRUE)
      
      
      coeff_2sls <- as.matrix(c(coefficients(firststep), rhotilde))
      rownames(coeff_2sls)<-c(colnames(Zmat), 'rho')
      s2_2sls<-crossprod(ubase)/(n-k)
      
      
      model.data<-data.frame(cbind(y,x[,-1]))
      
      method<-"gm spatial"
      
      k <- nrow(coeff_2sls)
      R <- matrix(0,1,k)
      R[,((k-1):k)] <- 1
      Rbeta <- R%*%coeff_2sls
      Rvar <- R %*% vcmat_2sls$Omega %*% t(R)
      stat <- as.numeric(t(Rbeta) %*% Rbeta/Rvar)
      pval <- pchisq(stat,df=1,lower.tail=FALSE)
      W <- list(stat = stat, pval = pval)
      
      
      
      results_2sls <- list(coefficients = coeff_2sls, var = vcmat_2sls$Omega, 
                           s2 = s2_2sls, call = cl, residuals = as.numeric(ubase), 
                           model = model.data, method = method, W = W, 
                           firststep = firststep$coefficients, init.rho = rhotilde, 
                           Durbin = Durbin)
      
      class(results_2sls)<-c("sphet", "error_gmm", "error_sphet")
      
    }
    
    
  }
  else{
    
    Ggmat<-gg_hom(Ws, ubase, n)
    optres <- nlminb(pars, optimfunct, control = control, 
                     v = Ggmat, verbose = verbose, lower = -0.9 + .Machine$double.eps , 
                     upper = 0.9 -  .Machine$double.eps)
    rhotilde <- optres$par
     
  }
  
  yt  <- y - rhotilde * Ws %*% y
  wZmat <- Ws %*% Zmat
  Zt <- Zmat - rhotilde * wZmat
  
  secondstep <- spatial.ivreg(y = yt , Zmat = Zt, Hmat = Hmat, het = het, HAC = HAC)
  delta <- coefficients(secondstep)
  utildeb <- y - Zmat %*% delta
 
  if(het){
    
    
    Ggmat <- gg_het(Ws, utildeb, n)
    gmm.weghts <- psirhorho_het(rhotilde, utildeb, Hmat, Zmat, Ws, step1.c = FALSE)
    
    optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, 
                     vcmat = gmm.weghts$Phiinv, verbose = verbose, 
                     lower = -0.9 + .Machine$double.eps , 
                     upper = 0.9 -  .Machine$double.eps, control = control)	
    
    rhofin <- optres$par
    gmm.weghts <- psirhorho_het(rhofin, utildeb, Hmat, Zmat, Ws, step1.c = FALSE)
    
    vcmat <- Omega_het(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, 
                       gmm.weghts$A2, gmm.weghts$a.vec1, 
                       gmm.weghts$a.vec2, Hmat, 
                       Ggmat$bigG, gmm.weghts$Phiinv, 
                       gmm.weghts$epsilon, gmm.weghts$Zstar, Ws, step1.c = FALSE)
  }
  else{
    
    Ggmat <- gg_hom(Ws, utildeb, n)
    gmm.weghts <- psirhorho_hom(rhotilde, utildeb, Hmat, Zmat, Ws, Ggmat$d, Ggmat$v.vec )
    optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, 
                     vcmat = gmm.weghts$Phiinv, verbose = verbose, 
                      control = control, lower = -1 + .Machine$double.eps , 
                     upper = 1 -  .Machine$double.eps)	
    
    rhofin <- optres$par
    gmm.weghts<-psirhorho_hom(rhofin, utildeb, Hmat, Zmat, Ws, Ggmat$d, Ggmat$v.vec)
    
    vcmat <- Omega_hom(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, 
                       gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, 
                       Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar)
    
    
  }
  
  coeff <- as.matrix(c(as.numeric(delta), rhofin))
  rownames(coeff) <- c(colnames(Zmat), 'rho')
  s2 <- crossprod(utildeb)/(n-k)
  
  
  model.data <- data.frame(cbind(y,x[,-1]))
  
  method <- "gmm error"
  
  k <- nrow(coeff)
  R <- matrix(0,1,k)
  R[,((k-1):k)] <- 1
  Rbeta <- R %*% coeff
  Rvar <- R %*% vcmat$Omega %*% t(R)
  stat <- as.numeric(t(Rbeta) %*% Rbeta/Rvar)
  pval <- pchisq(stat,df = 1, lower.tail = FALSE)
  W <- list(stat = stat, pval = pval)
  
  
  
  if(het && step1.c) results <- list(coefficients = coeff,var = vcmat$Omega, 
                                   s2 = s2, call = cl, residuals = as.numeric(utildeb), 
                                   model = model.data, method = method, W = W, 
                                   firststep = firststep$coefficients, 
                                   init.rho = rhotilde,  twosls = results_2sls, 
                                   Durbin = Durbin, endog = endog)
  
  else  results <- list(coefficients = coeff, var = vcmat$Omega, s2 = s2, call = cl, 
                      residuals = as.numeric(utildeb), model = model.data,
                      method = method, W = W, firststep = firststep$coefficients, 
                      init.rho = rhotilde, Durbin = Durbin, endog = endog)
  
  
  class(results) <- c("sphet", "error_gmm", "error_sphet")# gmm error
  
  return(results)
  
  
}

laghac <- function(formula, data, listw, listw2, endog, 
                   instruments, lag.instr,  verbose, 
                   na.action,  het, HAC, distance, 
                   type, bandwidth, cl, Durbin = NULL){
  
  mt <- terms(formula, data = data)
  mf <- lm(formula, data, na.action = na.action, method="model.frame")
  na.act <- attr(mf,'na.action')
  
  y <- c(model.extract(mf,"response"))
  x <- model.matrix(mt,mf)
  
  if (length(y)!=nrow(x)) 
    stop("x and y have different length")
  
  if (any(is.na(y))) 
    stop("NAs in dependent variable")
  if (any(is.na(x))) 
    stop("NAs in independent variable")
  
  if(HAC){ 
    if(is.null(distance)) stop("No distance measure specified")
    if(!inherits(distance,"distance")) 
      stop("The distance measure is not a distance object")
    if(!(type %in% c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"))) stop("Unknown kernel")
  }	
  
  n <- nrow(x)
  k <- ncol(x)	
  xcolnames <- colnames(x)
  
  K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)
  
  if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
  if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
  if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
  if(inherits(listw,"Matrix"))  Ws <- listw	
  
  if (nrow(x) != nrow(Ws))
    stop("Input data and weights have different dimension")
  
  
  if(Durbin == TRUE | class(Durbin) == "formula"  ){
    if(class(Durbin) == "formula"){
      ## For each value of K (i.e, 1 or 2) there are four cases:
      # 1) formula = y ~ x1 + x2 + x3, Durbin = ~ x2 + x3 (only a subset of the x's are lagged and all of them are also in the formula)
      # 2) formula = y ~ x1 + x2 + x3, Durbin = ~ x1 + x2 + x3  (same as Durbin = TRUE)
      # 3) formula = y ~ x1 + x2 + x3, Durbin = ~ x4 (only variables that show up only lagged)
      # 4) formula = y ~ x1 + x2 + x3, Durbin = ~ x3 + x4 (one variable is in formula and one not)
      xdur <- as.matrix(lm(Durbin, data, na.action=na.action, method="model.frame"))
      
      if(sum(match(xcolnames, colnames(xdur), nomatch = 0))==0  && K != 2 && k != 1){
        
        wxdur <- as.matrix(Ws %*% xdur)
        wwxdur <- as.matrix(Ws %*% wxdur)
        wwwxdur <- as.matrix(Ws %*% wwxdur)
        
        if (k > 1 || (k == 1 && K == 1)) {
          
          wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
          for (i in K:k) {
            Wx <- Ws %*% x[, i]
            wx[, (i - (K - 1))] <- as.matrix(Wx)
            
          }
          wwx <- as.matrix(Ws %*% wx)  
        } 
        
        if(K==2){
          Hin <- cbind(x, wx, wwx, wxdur, wwxdur, wwwxdur)
          x <- cbind(x, wxdur)  
          colnames(x) <-  c(xcolnames, paste("lag_", colnames(xdur), sep=""))
        #   print(head(x))
        #   print(head(Hin))
         } 
        else {
          Hin <- cbind(1, x, wx, wwx, wxdur, wwxdur, wwwxdur)
          x <- cbind(x, wxdur)  
          colnames(x) <-  c(xcolnames, paste("lag_", colnames(xdur), sep=""))
          # print(head(x))
          # print(head(Hin))
        }
      }  
      
      else{
        
        if(K==1 ){
          
          pos.xdur <- which(colnames(xdur) %in% xcolnames)
          pos.x <- which(xcolnames %in% colnames(xdur))
          
          if(all(is.na(match(colnames(xdur), xcolnames)))) onx <- as.matrix(x)
          else onx <- as.matrix(x[,-pos.x])
          
          if(dim(onx)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(onx) <- nmonx <-  xcolnames
            else names(onx) <- nmonx <-  xcolnames[-pos.x] 
          } 
          else nmonx <- NULL
          
          if(all(is.na(match(colnames(xdur), xcolnames)))) onxl <- as.matrix(xdur)
          else onxl <- as.matrix(xdur[, -pos.xdur])
          
          if(dim(onxl)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(onxl) <- nmonxl <- colnames(xdur)
            else names(onxl) <- nmonxl <- colnames(xdur)[-pos.xdur]  
          } 
          else nmonxl <- NULL
          
          onboth <- as.matrix(xdur[,pos.xdur])
          if(dim(onboth)[2] != 0){
            if(all(is.na(match(colnames(xdur), xcolnames))))  names(onboth) <- nmonb <- colnames(xdur)
            else  names(onboth) <- nmonb <- colnames(xdur)[pos.xdur] 
          }
          else nmonb <- NULL
         
          wonx <- as.matrix(Ws %*% onx)
          wwonx <- as.matrix(Ws %*% wonx)
          
          wonxl <- as.matrix(Ws %*% onxl)
          if(dim(wonxl)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur),sep = "") 
            else names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur)[-pos.xdur],sep = "") 
          } 
          else nmwonxl <- NULL
          wwonxl <- as.matrix(Ws %*% wonxl)
          wwwonxl <- as.matrix(Ws %*% wwonxl)
          
          wonboth <- as.matrix(Ws %*% onboth)
          
          if(dim(wonboth)[2]!=0) {
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur),sep = "") 
            else names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur)[pos.xdur],sep = "") 
          }
          else nmwonb <- NULL
          wwonboth <- as.matrix(Ws %*% wonboth)
          wwwonboth <- as.matrix(Ws %*% wwonboth)
          
          
          x <- cbind(onx, onboth, wonboth, wonxl)
          colnames(x) <- c(nmonx,nmonb,nmwonb,nmwonxl)
          Hin <- cbind(1, x, wonx, wwonx, wwonboth, wwwonboth,wwonxl,wwwonxl )
          # print(head(x))
          # print(head(Hin))
          # 
          
        } 
        else{
          
          if(k !=1){
            
            pos.xdur <- which(colnames(xdur) %in% xcolnames)
            pos.x <- which(xcolnames %in% colnames(xdur))
            
            if(all(is.na(match(colnames(xdur), xcolnames)))) onx <- as.matrix(x)
            else onx <- as.matrix(x[,-pos.x]) 
            
            if(dim(onx)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onx) <- nmonx <-  xcolnames
              else names(onx) <- nmonx <-  xcolnames[-pos.x]
            }  
            else nmonx <- NULL
            
            if(all(is.na(match(colnames(xdur), xcolnames)))) onxl <- as.matrix(xdur)
            else  onxl <- as.matrix(xdur[, -pos.xdur])
            
            if(dim(onxl)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onxl) <- nmonxl <- colnames(xdur)
              else names(onxl) <- nmonxl <- colnames(xdur)[-pos.xdur]
            }   
            else nmonxl <- NULL
             
            onboth <- as.matrix(xdur[,pos.xdur])
            if(dim(onboth)[2] != 0) {
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onboth) <- nmonb <- colnames(xdur)
              else names(onboth) <- nmonb <- colnames(xdur)[pos.xdur] 
            } 
            else nmonb <- NULL
            
            wonx <- as.matrix(Ws %*% onx[,-1])
            wwonx <- as.matrix(Ws %*% wonx)
            
            wonxl <- as.matrix(Ws %*% onxl)
            if(dim(wonxl)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur),sep = "") 
              else names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur)[-pos.xdur],sep = "") 
            } 
            else nmwonxl <- NULL
            wwonxl <- as.matrix(Ws %*% wonxl)
            wwwonxl <- as.matrix(Ws %*% wwonxl)
            
            wonboth <- as.matrix(Ws %*% onboth)
            
            if(dim(wonboth)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur),sep = "") 
              else names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur)[pos.xdur],sep = "") 
            } 
            else nmwonb <- NULL
            wwonboth <- as.matrix(Ws %*% wonboth)
            wwwonboth <- as.matrix(Ws %*% wwonboth)
            
            
            x <- cbind(onx, onboth, wonboth, wonxl)
            colnames(x) <- c(nmonx,nmonb,nmwonb,nmwonxl)
            Hin <- cbind(x, wonx, wwonx, wwonboth, wwwonboth,wwonxl,wwwonxl )
            # print(head(x))
            # print(head(Hin))
          }
          else{
            
            wxdur <-  as.matrix(Ws %*% xdur)
            wwxdur <- as.matrix(Ws %*% wxdur)
            wwwxdur <- as.matrix(Ws %*% wwxdur)
            x <- cbind(x, wxdur)
            colnames(x) <- c(xcolnames, paste("lag_",colnames(xdur),sep =""))
            Hin <- cbind(x, wwxdur, wwwxdur )
            # print(head(x))
            # print(head(Hin))
          }
        }
      }
    
    }
    
    else{
      
      if (k > 1 || (k == 1 && K == 1)) {
        
        wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
        for (i in K:k) {
          Wx <- Ws %*% x[, i]
          wx[, (i - (K - 1))] <- as.matrix(Wx)
        }
        wwx <- as.matrix(Ws %*% wx)  
        wwwx <- as.matrix(Ws  %*%  wwx)
      } 
      
      if(K==2){
        Hin <- cbind(x,wx,wwx,wwwx)
        x <- cbind(x, wx)  
        colnames(x) <-  c(xcolnames, paste("lag_", xcolnames[-1], sep=""))
        # print(head(x))
        # print(head(Hin))
      } 
      else {
        Hin <- cbind(1,x,wx,wwx,wwwx)
        x <- cbind(x, wx)  
        colnames(x) <-  c(xcolnames, paste("lag_", xcolnames, sep=""))
        # print(head(x))
        # print(head(Hin))
      }
    }
  }
  else{  
    
    if (k > 1 || (k == 1 && K == 1)) {
      wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
      for (i in K:k) {
        Wx <- Ws %*% x[, i]
        wx[, (i - (K - 1))] <- as.matrix(Wx)
      }
      wwx <- Ws %*% wx                    					         
    }
    
    if(K==2)    Hin <- cbind(x, wx, wwx)
    else        Hin <- cbind(1, x, wx, wwx)
    x <- x
    # print(head(x))
    # print(head(Hin))
  }
  
  wy<-Ws %*% y	
  colnames(wy)<-"lambda"
  
  
  if (!is.null(endog) && is.null(instruments)) stop("No instruments specified for the endogenous variable in the model")
  
  
  if (!is.null(endog)) {
    endog <- as.matrix(lm(endog, data, na.action=na.action, method="model.frame"))
    instruments <- as.matrix(lm(instruments, data, na.action=na.action, method="model.frame"))
    if(lag.instr) {
      winst <- Ws %*% instruments
      wwinst<- Ws %*% winst	
      AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst))        
    }
    else  AddH <- instruments        
    Hmat <- cbind(Hin, AddH)
    Zmat<- cbind(x, endog, as.matrix(wy))            
    colnames(Zmat) <- c(colnames(x), colnames(endog), colnames(wy))               
  }
  else {
    Zmat<- cbind(x, as.matrix(wy))                    
    Hmat <- Hin
  }
  
  
  results <- spatial.ivreg(y =y , Zmat = Zmat, Hmat = Hmat, het = het,  HAC=HAC, type=type, bandwidth=bandwidth, distance=distance)	
  model.data <- data.frame(cbind(y, x[, -1]))
  results$call <- cl
  results$listw <- Ws
  results$model <- model.data
  results$type <- type
  results$bandwidth <- bandwidth
  results$method <- "s2slshac"
  results$HAC <- HAC
  results$endog <- endog
  class(results) <- c("sphet", "lag_gmm", "stsls_sphet")
  return(results)
  
}

olshac <- function(formula, data, endog, instruments, listw, 
                   na.action, het, HAC, distance, type, bandwidth, cl, Durbin = NULL){
  
  #if(!isTRUE(HAC) || !isTRUE(het) || !is.null(Durbin))  
  #extract model objects	
  mt <- terms(formula,data=data)
  mf <- lm(formula, data, na.action=na.action, method="model.frame")
  na.act <- attr(mf, 'na.action')
  
  #generates x and y 
  y <- c(model.extract(mf, "response"))
  x <- model.matrix(mt, mf)
  
  #checks on teh dimensions of x and y 	
  if (length(y)!=nrow(x)) 
    stop("x and y have different length")
  
  #check that X and y does not have missing values	
  if (any(is.na(y))) 
    stop("NAs in dependent variable")
  if (any(is.na(x))) 
    stop("NAs in independent variable")
 
  Ws <- NULL
#
  if (!is.null(Durbin)){ 
  if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
  if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
  if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
  if(inherits(listw,"Matrix"))  Ws <- listw	
  if (nrow(x) != nrow(Ws)) stop("Input data and weights have different dimension")
  }
  
  if(!is.null(endog)) model <- 'ols.end'
  else model <- "ols"
  
  if(HAC){
    if(is.null(distance)) stop("No distance measure specified")
    if(!inherits(distance,"distance")) 
      stop("The distance measure is not a distance object")
    if(!(type %in% c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"))) stop("Unknown kernel")
  }	
  
  
  
  #fix the dimensions of the problem
  n <- nrow(x)
  k <- ncol(x)	
  xcolnames <- colnames(x)
  
  K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)
  
  if(Durbin == TRUE | class(Durbin) == "formula"  ){
    if(class(Durbin) == "formula"){
      ## For each value of K (i.e, 1 or 2) there are four cases:
      # 1) formula = y ~ x1 + x2 + x3, Durbin = ~ x2 + x3 (only a subset of the x's are lagged and all of them are also in the formula)
      # 2) formula = y ~ x1 + x2 + x3, Durbin = ~ x1 + x2 + x3  (same as Durbin = TRUE)
      # 3) formula = y ~ x1 + x2 + x3, Durbin = ~ x4 (only variables that show up only lagged)
      # 4) formula = y ~ x1 + x2 + x3, Durbin = ~ x3 + x4 (one variable is in formula and one not)
      xdur <- as.matrix(lm(Durbin, data, na.action=na.action, method="model.frame"))
      
      if(sum(match(xcolnames, colnames(xdur), nomatch = 0))==0  && K != 2 && k != 1){
      
        wxdur <- as.matrix(Ws %*% xdur)
        wwxdur <- as.matrix(Ws %*% wxdur)
        
        if (k > 1 || (k == 1 && K == 1)) {
          
          wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
          for (i in K:k) {
            Wx <- Ws %*% x[, i]
            wx[, (i - (K - 1))] <- as.matrix(Wx)
          }
          
        } 
        
        if(K==2){
          Hin <- cbind(x, wx, wxdur, wwxdur)
          x <- cbind(x, wxdur)  
          colnames(x) <-  c(xcolnames, paste("lag_", colnames(xdur), sep=""))
          # print(head(x))
          # print(head(Hin))
        } 
        else {
          
          Hin <- cbind(1, x, wx, wxdur, wwxdur)
          x <- cbind(x, wxdur)  
          
          colnames(x) <-  c(xcolnames, paste("lag_", colnames(xdur), sep=""))
          # print(head(x))
          # print(head(Hin))
        }
      }  
      
      else{
        
        if(K==1 ){
      
          pos.xdur <- which(colnames(xdur) %in% xcolnames)
          pos.x <- which(xcolnames %in% colnames(xdur))
          
          
          if(all(is.na(match(colnames(xdur), xcolnames)))) onx <- as.matrix(x)
          else onx <- as.matrix(x[,-pos.x])
          
          if(dim(onx)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(onx) <- nmonx <-  xcolnames
            else names(onx) <- nmonx <-  xcolnames[-pos.x] 
          } 
          else nmonx <- NULL
          
          if(all(is.na(match(colnames(xdur), xcolnames)))) onxl <- as.matrix(xdur)
          else onxl <- as.matrix(xdur[, -pos.xdur])
          
          if(dim(onxl)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(onxl) <- nmonxl <- colnames(xdur)
            else names(onxl) <- nmonxl <- colnames(xdur)[-pos.xdur]  
          } 
          else nmonxl <- NULL
          
          onboth <- as.matrix(xdur[,pos.xdur])
          if(dim(onboth)[2] != 0){
            if(all(is.na(match(colnames(xdur), xcolnames))))  names(onboth) <- nmonb <- colnames(xdur)
            else  names(onboth) <- nmonb <- colnames(xdur)[pos.xdur] 
          }
          else nmonb <- NULL
          
          wonx <- as.matrix(Ws %*% onx)
  
          
          wonxl <- as.matrix(Ws %*% onxl)
          if(dim(wonxl)[2]!=0){
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur),sep = "") 
            else names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur)[-pos.xdur],sep = "") 
          } 
          else nmwonxl <- NULL
          wwonxl <- as.matrix(Ws %*% wonxl)
          
          wonboth <- as.matrix(Ws %*% onboth)
         
          if(dim(wonboth)[2]!=0) {
            if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur),sep = "") 
            else names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur)[pos.xdur],sep = "") 
          }
          else nmwonb <- NULL
          wwonboth <- as.matrix(Ws %*% wonboth)
          
          
          x <- cbind(onx, onboth, wonboth, wonxl)
          colnames(x) <- c(nmonx,nmonb,nmwonb,nmwonxl)
          Hin <- cbind(1, x, wonx, wwonboth, wwonxl )
          # print(head(x))
          # print(head(Hin))
          # 
          
        } 
        else{
          
          if(k !=1){
          
            pos.xdur <- which(colnames(xdur) %in% xcolnames)
            pos.x <- which(xcolnames %in% colnames(xdur))
            
            if(all(is.na(match(colnames(xdur), xcolnames)))) onx <- as.matrix(x)
            else onx <- as.matrix(x[,-pos.x]) 
          
            if(dim(onx)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onx) <- nmonx <-  xcolnames
              else names(onx) <- nmonx <-  xcolnames[-pos.x]
            }  
            else nmonx <- NULL
            
            if(all(is.na(match(colnames(xdur), xcolnames)))) onxl <- as.matrix(xdur)
            else  onxl <- as.matrix(xdur[, -pos.xdur])
            
            if(dim(onxl)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onxl) <- nmonxl <- colnames(xdur)
              else names(onxl) <- nmonxl <- colnames(xdur)[-pos.xdur]
            }   
            else nmonxl <- NULL
               
            onboth <- as.matrix(xdur[,pos.xdur])
            if(dim(onboth)[2] != 0) {
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(onboth) <- nmonb <- colnames(xdur)
              else names(onboth) <- nmonb <- colnames(xdur)[pos.xdur] 
            } 
            else nmonb <- NULL
             
            wonx <- as.matrix(Ws %*% onx[,-1])
            wonxl <- as.matrix(Ws %*% onxl)
            if(dim(wonxl)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur),sep = "") 
              else names(wonxl) <- nmwonxl <- paste("lag_", colnames(xdur)[-pos.xdur],sep = "") 
            } 
            else nmwonxl <- NULL
            wwonxl <- as.matrix(Ws %*% wonxl)
            
            wonboth <- as.matrix(Ws %*% onboth)
            if(dim(wonboth)[2]!=0){
              if(all(is.na(match(colnames(xdur), xcolnames)))) names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur),sep = "") 
              else names(wonboth) <- nmwonb <- paste("lag_",colnames(xdur)[pos.xdur],sep = "") 
            } 
            else nmwonb <- NULL
            wwonboth <- as.matrix(Ws %*% wonboth)
            
            
            x <- cbind(onx, onboth, wonboth, wonxl)
            colnames(x) <- c(nmonx,nmonb,nmwonb,nmwonxl)
            Hin <- cbind(x, wonx, wwonboth, wwonxl )
            # print(head(x))
            # print(head(Hin))
          }
          else{
            
            wxdur <-  as.matrix(Ws %*% xdur)
            wwxdur <- as.matrix(Ws %*% wxdur)
            
            x <- cbind(x, wxdur)
            colnames(x) <- c(xcolnames, paste("lag_",colnames(xdur),sep =""))
            Hin <- cbind(x, wwxdur)
            # print(head(x))
            # print(head(Hin))
          }
        }
      }

    }
    
    else{
      
      if (k > 1 || (k == 1 && K == 1)) {
        
        wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
        for (i in K:k) {
          Wx <- Ws %*% x[, i]
          wx[, (i - (K - 1))] <- as.matrix(Wx)
        }
        wwx <- as.matrix(Ws %*% wx)  
        
      } 
      
      if(K==2){
        Hin <- cbind(x,wx,wwx)
        x <- cbind(x, wx)  
        colnames(x) <-  c(xcolnames, paste("lag_", xcolnames[-1], sep=""))
        # print(head(x))
        # print(head(Hin))
      } 
      else {
        Hin <- cbind(1,x,wx,wwx)
        x <- cbind(x, wx)  
        colnames(x) <-  c(xcolnames, paste("lag_", xcolnames, sep=""))
        # print(head(x))
        # print(head(Hin))
      }
    }
  }
  else{  
    
    if (k > 1 || (k == 1 && K == 1)) {
      wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
      for (i in K:k) {
        Wx <- Ws %*% x[, i]
        wx[, (i - (K - 1))] <- as.matrix(Wx)
      }
     
    }
    
    if(K==2)    Hin <- cbind(x, wx)
    else        Hin <- cbind(1, x, wx)
    x <- x
    # print(head(x))
    # print(head(Hin))
  }
  

  if(model == "ols.end" ){
   
    endog <- as.matrix(lm(endog, data, na.action=na.action, method="model.frame"))	
    instruments <- as.matrix(lm(instruments, data, na.action=na.action, method="model.frame"))	
    AddH<- cbind(instruments)
    Zmat<- cbind(x, endog)            
    Hmat <- cbind(Hin, AddH)
    
    
    if(HAC) results <- spatial.ivreg(y =y , Zmat = Zmat, Hmat = Hmat, 
                            HAC = HAC, type = type, bandwidth = bandwidth, distance = distance)	
    else results <-spatial.ivreg(y =y , Zmat = Zmat, Hmat = Hmat, 
                                 het = het)	
     
    model.data <- data.frame(cbind(y, x[, -1]))
    results$call <- cl
    results$model <- model.data
    results$type <- type
    results$bandwidth <- bandwidth
    results$method <- "s2slshac"
    results$HAC <- HAC
    results$endog <- TRUE
    class(results) <- c("sphet", "ols_sphet")
    
  }
  
  if(model == "ols"){
  if(HAC)  results <- hac.ols(y =y , x = x, HAC=HAC, type=type, bandwidth=bandwidth, distance=distance, het = FALSE)	
   else  results <-  hac.ols(y =y , x = x, het = het) 
    model.data <- data.frame(cbind(y, x[, -1]))
    results$call <- cl
    results$listw <- Ws
    results$model <- model.data
    results$type <- type
    results$bandwidth <- bandwidth
    results$method <- "olshac"
    results$HAC <- HAC
    results$Durbin <- Durbin
    results$endog <- FALSE
    class(results) <- c("sphet", "ols_sphet")
  }
  
  return(results)
}

