

# sarargmm <- function(formula, data, listw, listw2, endog, 
#                        instruments, lag.instr, initial.value, 
#                        het, verbose, na.action,
#                        step1.c, control, HAC, cl, Durbin = NULL){
#   
#   mt <- terms(formula, data = data)
#   mf <- lm(formula, data, na.action = na.action, method = "model.frame")
#   na.act <- attr(mf, 'na.action')
#   
#   y <- c(model.extract(mf, "response"))
#   x <- model.matrix(mt, mf)
#   
#   if (length(y) != nrow(x)) 
#     stop("x and y have different length")
#   
#   if (any(is.na(y))) 
#     stop("NAs in dependent variable")
#   if (any(is.na(x))) 
#     stop("NAs in independent variable")
#   
#   n <- nrow(x)
#   k <- ncol(x)	
#   xcolnames <- colnames(x)
#   #print(xcolnames)
#   K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1] == 1), 2, 1)
#   
#   if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
#   if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
#   if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
#   if(inherits(listw,"Matrix"))  Ws <- listw	
#   
#   
#   if (nrow(x) != nrow(Ws)) stop("Input data and weights have different dimension")
#   
#   if (k > 1) {
#     wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
#     for (i in K:k) {
#       Wx <- Ws %*% x[, i]
#       wx[, (i - (K - 1))] <- as.matrix(Wx)
#     }
#     wwx <- Ws %*% wx                    					         
#   }
#   
#   if(!is.null(Durbin)){
#     x.dur <- as.matrix(lm(Durbin, data, na.action = na.action, method = "model.frame"))
#    # print(colnames(x.dur))
#     # print(match(colnames(x.dur), xcolnames, nomatch = 0))
# if(sum(match(colnames(x.dur), xcolnames, nomatch = 0))!= 0) stop("Lagged explanatory variables cannot be endogenous")    
#     wx.dur <- Ws %*% x.dur
#     wwx.dur <- Ws %*% wx.dur
#     wwwx.dur <- Ws %*% wwx.dur
#     addx <- cbind(x.dur, wx.dur)
#     x.lag.names <- c(colnames(x.dur), paste("lag_", colnames(x.dur), sep=""))
#     Inx.dur <- cbind(as.numeric(wwx.dur),as.numeric(wwwx.dur))
#   }
#   
#   
#   if(is.null(listw2)) {
#     twow <- FALSE		
#     Ws2 <- Ws
#   } 
#   else{ 
#     twow <- TRUE	
#     if(!inherits(listw2,c("listw", "Matrix", "matrix"))) stop("listw2 format unknown")
#     if(inherits(listw2,"listw"))  Ws2<-listw2dgCMatrix(listw2)	
#     if(inherits(listw2,"matrix"))  Ws2<-Matrix(listw2)	
#     
#     if(identical(listw, listw2)){ 
#       twow <- FALSE		
#       Ws2 <- Ws
#     }
#     
#     
#     if (k > 1) {
#       w2x <- matrix(nrow = n, ncol = (k  - (K - 1)))
#       for (i in K:k) {
#         W2x <- Ws2 %*% x[, i]          
#         w2x[, (i - (K - 1))] <- as.matrix(W2x)
#       }
#       w2wx <- Ws2 %*% wx                   	
#       w2wwx <- Ws2 %*% wwx                    	          
#       
#     }
#     
#     if(!is.null(Durbin)){
#       w2x.dur <- Ws2 %*% x.dur
#       print(x.dur)
#       print(w2x.dur)
#       w2wx.dur <- Ws2 %*% wx.dur
#       w2wwx.dur <- Ws2 %*% wwx.dur
#       w2wwwx.dur <- Ws2 %*% wwwx.dur
#       Inw2x.dur <- cbind(as.numeric(w2x.dur), as.numeric(w2wx.dur), as.numeric(w2wwx.dur), as.numeric(w2wwwx.dur))
#     }
#     
#   }
#   
#   
#   
#   wy <- Ws %*% y	
#   colnames(wy) <- "lambda"
#   
#   if (!is.null(endog) && is.null(instruments)) stop("No instruments specified for the endogenous variable in the model")
#   
#   if (!is.null(endog)) {
#     endog <- as.matrix(lm(endog, data, na.action = na.action, method = "model.frame"))
#     if(!is.null(Durbin) && (sum(match(colnames(x.dur), colnames(endog), nomatch = 0))!= 0)) stop("Lagged explanatory variables cannot be endogenous")
#     instruments <- as.matrix(lm(instruments, data, na.action = na.action, method = "model.frame"))
#     if(lag.instr) {
#       winst <- Ws %*% instruments
#       wwinst <- Ws %*% winst	
#       if(twow){
#         w2i <- Ws2 %*% instruments 
#         w2wi <- Ws2 %*% winst 
#         w2wwi <- Ws2 %*% wwinst 	
#         AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst), as.matrix(w2i), as.matrix(w2wi),as.matrix(w2wwi))        
#       }
#       else  AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst))        
#     }
#     else  AddH <- instruments        
#     if (K==2) {
#       if(twow){ 
#         if(is.null(Durbin)) Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx), AddH)
#         else Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx), addx, Inx.dur, Inw2x.dur ,AddH)
#       }
#       else{
#         if(is.null(Durbin))  Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx),  AddH)
#         else  Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur, AddH)
#       } 
#     }
#     else {
#       if(twow){
#         if (is.null(Durbin)) Hmat <- cbind(1,x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx), AddH)
#         else Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx), addx, Inx.dur, Inw2x.dur ,AddH)
#       } 
#       else {
#         if (is.null(Durbin)) Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), AddH)
#         else Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur, AddH)
#       } 
#     }
#     
#     if (is.null(Durbin)){
#       Zmat<- cbind(x, endog, as.matrix(wy))            
#       colnames(Zmat) <- c(colnames(x), colnames(endog), colnames(wy))               
#     }
#     else {
#       Zmat<- cbind(x, addx, endog, as.matrix(wy))
#       colnames(Zmat) <- c(colnames(x), x.lag.names , colnames(endog), colnames(wy))
#     }
#   } 
#   else {
#     if (is.null(Durbin)){
#       Zmat<- cbind(x, as.matrix(wy))
#       colnames(Zmat) <- c(colnames(x), colnames(wy))
#     }
#     else {
#       Zmat<- cbind(x, addx, as.matrix(wy))
#       colnames(Zmat) <- c(colnames(x), x.lag.names, colnames(wy))
#     }
#     if (K==2){
#       if(twow) {
#         if(is.null(Durbin)) Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx))
#         else Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx), addx, Inx.dur, Inw2x.dur)
#       }
#       else {
#         if(is.null(Durbin)) Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx)) 
#         else Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur) 
#       } 
#     }
#     else {
#       if(twow) {
#         if(is.null(Durbin))  Hmat <- cbind(1,x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx))
#         else Hmat <- cbind(1,x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx), addx, Inx.dur, Inw2x.dur)
#       }
#       else  {
#         if(is.null(Durbin))   Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx)) 
#         else Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur) 
#       }
#     }
#   }
#   
#   firststep<-spatial.ivreg(y = y , Zmat = Zmat, Hmat = Hmat, het = het, HAC = HAC)
#   ubase<-residuals(firststep)
#   
#   if (initial.value=="SAR"){
#     Wubase<- Ws %*% ubase
#     pars<-coefficients(lm(as.numeric(ubase) ~ as.numeric(Wubase)-1))
#   }
#   else pars<-initial.value
#   
#   
#   
#   if(het){
#     Ggmat <- gg_het(Ws2, ubase, n)
#     optres <-nlminb(pars, optimfunct, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps, control= control, v = Ggmat, verbose = verbose)
#     rhotilde <- optres$par
#     
#     
#     if(step1.c){
#       gmm.weghts1.c <- psirhorho_het(rhotilde, ubase, Hmat, Zmat, Ws2, step1.c = TRUE)
#       optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat= gmm.weghts1.c$Phiinv, verbose = verbose, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps, control = control)	
#       
#       rhotilde<-optres$par
#       gmm.weghts1.c<-psirhorho_het(rhotilde, ubase, Hmat, Zmat, Ws2, step1.c = TRUE)
#       vcmat_2sls <- Omega_het(rhotilde, gmm.weghts1.c$Pmat, gmm.weghts1.c$A1, gmm.weghts1.c$A2, gmm.weghts1.c$a.vec1, gmm.weghts1.c$a.vec2, Hmat, Ggmat$bigG, gmm.weghts1.c$Phiinv, gmm.weghts1.c$epsilon, gmm.weghts1.c$Zstar, Ws2, step1.c = TRUE)
#       coeff_2sls <- as.matrix(c(coefficients(firststep), rhotilde))
#       rownames(coeff_2sls)<-c(colnames(Zmat), 'rho')
#       s2_2sls<-crossprod(ubase)/(n-k)
#       model.data<-data.frame(cbind(y,x[,-1]))
#       method<-"gm sarar spatial"
#       
#       k<-nrow(coeff_2sls)
#       R<-matrix(0,1,k)
#       R[,((k-1):k)]<-1
#       Rbeta<-R%*%coeff_2sls
#       Rvar<-R%*% vcmat_2sls$Omega %*%t(R)
#       stat<-as.numeric(t(Rbeta)%*% Rbeta/Rvar)
#       pval <- pchisq(stat,df=1,lower.tail=FALSE)
#       W<-list(stat=stat,pval=pval)
#       results_2sls <- list(coefficients=coeff_2sls,var=vcmat_2sls$Omega, s2=s2_2sls, call=cl, residuals=as.numeric(ubase), model=model.data,method=method,W=W, firststep=firststep$coefficients, init.rho = rhotilde)
#       class(results_2sls)<-c("sphet", "gstsls")
#     }
#     
#     
#   }
#   else{
#     
#     Ggmat<-gg_hom(Ws2, ubase, n)
#     optres <- nlminb(pars, optimfunct, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control = control, v = Ggmat, verbose = verbose)
#     rhotilde<-optres$par
#     # print(rhotilde)	
#   }
#   
#   yt  <- y - rhotilde * Ws2 %*% y
#   wZmat <- Ws2 %*% Zmat
#   Zt <- Zmat - rhotilde * wZmat
#   
#   secondstep<-spatial.ivreg(y =yt , Zmat = Zt, Hmat = Hmat, het = het, HAC = HAC)
#   delta <- coefficients(secondstep)
#   utildeb <- y - Zmat %*% delta
#   
#   if(het){
#     
#     
#     Ggmat<-gg_het(Ws2, utildeb, n)
#     
#     gmm.weghts<-psirhorho_het(rhotilde, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)
#     
#     # gmm.weghts<-psirhorho_het_mod(rhotilde, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)
#     
#     # print(gmm.weghts$Phiinv)
#     optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat= gmm.weghts$Phiinv, verbose = verbose, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps, control = control)	
#     
#     rhofin<-optres$par
#     gmm.weghts<-psirhorho_het(rhofin, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)
#     
#     # gmm.weghts<-psirhorho_het_mod(rhofin, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)
#     
#     vcmat <- Omega_het(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar, Ws2, step1.c = FALSE)
#     # vcmat <- Omega_het_mod(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar, Ws2, step1.c = FALSE)
#   }
#   else{
#     
#     Ggmat<-gg_hom(Ws2, utildeb, n)
#     
#     gmm.weghts<-psirhorho_hom(rhotilde, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )
#     
#     # gmm.weghts<-psirhorho_hom_mod(rhotilde, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )
#     
#     
#     optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat= gmm.weghts$Phiinv, verbose = verbose, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps, control = control )	
#     
#     rhofin<-optres$par
#     # print(rhofin)
#     gmm.weghts<-psirhorho_hom(rhofin, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )
#     # gmm.weghts<-psirhorho_hom_mod(rhofin, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )
#     
#     vcmat <- Omega_hom(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar)
#     
#     # vcmat <- Omega_hom_mod(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2,gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar)
#     
#   }
#   
#   
#   coeff <- as.matrix(c(as.numeric(delta), rhofin))
#   rownames(coeff)<-c(colnames(Zmat), 'rho')
#   s2<-crossprod(utildeb)/(n-k)
#   
#   model.data<-data.frame(cbind(y,x[,-1]))
#   
#   method<-"gmm spatial"
#   
#   k<-nrow(coeff)
#   R<-matrix(0,1,k)
#   R[,((k-1):k)]<-1
#   Rbeta<-R%*%coeff
#   Rvar<-R%*% vcmat$Omega %*%t(R)
#   stat<-as.numeric(t(Rbeta)%*% Rbeta/Rvar)
#   pval <- pchisq(stat,df=1,lower.tail=FALSE)
#   W<-list(stat=stat,pval=pval)
#   
#   
#   
#   if(het && step1.c) results<-list(coefficients=coeff,var=vcmat$Omega, s2=s2, call=cl, residuals=as.numeric(utildeb), model=model.data,method=method,W=W, firststep=firststep$coefficients, init.rho = rhotilde,  twosls = results_2sls)
#   else  results<-list(coefficients=coeff,var=vcmat$Omega, s2=s2, call=cl, residuals=as.numeric(utildeb), model=model.data,method=method,W=W, firststep=firststep$coefficients, init.rho = rhotilde)
#   
#   
#   class(results)<-c("sphet", "gstsls") #remember to change to sarar when impacts will be developed
#   return(results)
# }
# 


# laghac <- function(formula, data, listw, listw2, endog, 
#                    instruments, lag.instr,  verbose, 
#                    na.action,  het, HAC, distance, 
#                    type, bandwidth, cl, Durbin = NULL){
#   
#   mt <- terms(formula, data = data)
#   mf <- lm(formula, data, na.action = na.action, method="model.frame")
#   na.act <- attr(mf,'na.action')
#   
#   y <- c(model.extract(mf,"response"))
#   x <- model.matrix(mt,mf)
#   
#   if (length(y)!=nrow(x)) 
#     stop("x and y have different length")
#   
#   if (any(is.na(y))) 
#     stop("NAs in dependent variable")
#   if (any(is.na(x))) 
#     stop("NAs in independent variable")
#   
#   if(HAC){ 
#     if(is.null(distance)) stop("No distance measure specified")
#     if(!inherits(distance,"distance")) 
#       stop("The distance measure is not a distance object")
#     if(!(type %in% c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"))) stop("Unknown kernel")
#   }	
#   
#   n <- nrow(x)
#   k <- ncol(x)	
#   xcolnames <- colnames(x)
#   
#   K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)
#   
#   if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
#   if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
#   if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
#   if(inherits(listw,"Matrix"))  Ws <- listw	
#   
#   if (nrow(x) != nrow(Ws))
#     stop("Input data and weights have different dimension")
#   
#   if (k > 1) {
#     wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
#     for (i in K:k) {
#       Wx <- Ws %*% x[, i]
#       wx[, (i - (K - 1))] <- as.matrix(Wx)
#     }
#     wwx <- Ws %*% wx                    					         
#   }
#   
#   if(!is.null(Durbin)){
#     x.dur <- as.matrix(lm(Durbin, data, na.action = na.action, method = "model.frame"))
#     if(sum(match(colnames(x.dur), xcolnames, nomatch = 0))!= 0) stop("Explanatory variables to be lagged cannot be specified in the main formula")    
#     wx.dur <- Ws %*% x.dur
#     wwx.dur <- Ws %*% wx.dur
#     wwwx.dur <- Ws %*% wwx.dur
#     addx <- cbind(x.dur, wx.dur)
#     x.lag.names <- c(colnames(x.dur), paste("lag_", colnames(x.dur), sep=""))
#     Inx.dur <- cbind(as.numeric(wwx.dur),as.numeric(wwwx.dur))
#   } 
#   
#   wy<-Ws %*% y	
#   colnames(wy)<-"lambda"
#   
#   
#   if (!is.null(endog) && is.null(instruments)) stop("No instruments specified for the endogenous variable in the model")
#   
#   
#   if (!is.null(endog)) {
#     endog <- as.matrix(lm(endog, data, na.action=na.action, method="model.frame"))
#     if(!is.null(Durbin) && (sum(match(colnames(x.dur), colnames(endog), nomatch = 0))!= 0)) stop("Lagged explanatory variables cannot be endogenous") 
#     instruments <- as.matrix(lm(instruments, data, na.action=na.action, method="model.frame"))
#     if(lag.instr) {
#       winst <- Ws %*% instruments
#       wwinst<- Ws %*% winst	
#       AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst))        
#     }
#     else  AddH <- instruments        
#     if (K==2){
#       if (is.null(Durbin)) Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), AddH)
#       else Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur, AddH)
#       
#     } 
#     else {
#       if (is.null(Durbin)) Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), AddH)
#       else Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur, AddH)
#     }
#     if (is.null(Durbin)){  
#       Zmat<- cbind(x, endog, as.matrix(wy))            
#       colnames(Zmat) <- c(colnames(x), colnames(endog), colnames(wy))               
#     }
#     else{
#       Zmat<- cbind(x, addx, endog, as.matrix(wy))            
#       colnames(Zmat) <- c(colnames(x), x.lag.names, colnames(endog), colnames(wy))
#     }
#   }
#   else {
#     
#     if(is.null((Durbin)))  Zmat<- cbind(x, as.matrix(wy))                    
#     else Zmat<- cbind(x, addx, as.matrix(wy))
#     
#     if (K==2) {
#       if(is.null((Durbin))) Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx)) 
#       else Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur)
#     }
#     else {
#       if(is.null((Durbin)))  Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx)) 
#       else Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur)
#     }
#     
#     
#   }
#   
#   
#   results <- spatial.ivreg(y =y , Zmat = Zmat, Hmat = Hmat, het = het,  HAC=HAC, type=type, bandwidth=bandwidth, distance=distance)	
#   model.data <- data.frame(cbind(y, x[, -1]))
#   results$call <- cl
#   results$model <- model.data
#   results$type <- type
#   results$bandwidth <- bandwidth
#   results$method <- "gmm spatial"
#   results$HAC <- HAC
#   class(results) <- c("sphet", "stsls_sphet")#change to lag hac
#   
#   return(results)
#   
# }
# 
# 





















# errorgmm <- function(formula, data, listw, listw2, endog, 
#                      instruments, lag.instr, initial.value, 
#                      het, verbose, na.action,
#                      step1.c, control, HAC, cl, Durbin = NULL){
#   
#   mt <- terms(formula,data=data)
#   mf <- lm(formula, data, na.action = na.action, method = "model.frame")
#   na.act <- attr(mf,'na.action')
#   
#   y <- c(model.extract(mf,"response"))
#   x <- model.matrix(mt,mf)
#   
#   if (length(y)!=nrow(x)) 
#     stop("x and y have different length")
#   
#   if (any(is.na(y))) 
#     stop("NAs in dependent variable")
#   if (any(is.na(x))) 
#     stop("NAs in independent variable")
#   
#   n <- nrow(x)
#   k <- ncol(x)	
#   xcolnames <- colnames(x)
#   
#   K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)
#   
#   if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
#   if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
#   if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
#   if(inherits(listw,"Matrix"))  Ws <- listw	
#   
#   if (nrow(x) != nrow(Ws))
#     stop("Input data and weights have different dimension")
#   
#   if (k > 1) {
#     wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
#     for (i in K:k) {
#       Wx <- Ws %*% x[, i]
#       wx[, (i - (K - 1))] <- as.matrix(Wx)
#     }
#     wwx <- Ws %*% wx                    					         
#   }
#   
#   if(!is.null(Durbin)){
#     x.dur <- as.matrix(lm(Durbin, data, na.action = na.action, method = "model.frame"))
#     if(sum(match(colnames(x.dur), xcolnames, nomatch = 0))!= 0) stop("Explanatory variables to be lagged cannot be specified in the main formula")    
#     wx.dur <- Ws %*% x.dur
#     wwx.dur <- Ws %*% wx.dur
#     # wwwx.dur <- Ws %*% wwx.dur
#     addx <- cbind(x.dur, wx.dur)
#     x.lag.names <- c(colnames(x.dur), paste("lag_", colnames(x.dur), sep=""))
#     Inx.dur <- as.numeric(wwx.dur)
#   }
#   
#   
#   if (!is.null(endog) && is.null(instruments)) stop("No instruments specified for the endogenous variable in the model")
#   
#   if (!is.null(endog)) {
#     endog <- as.matrix(lm(endog, data, na.action=na.action, method="model.frame"))
#     instruments <- as.matrix(lm(instruments, data, na.action=na.action, method="model.frame"))
#     if(lag.instr) { 
#       winst <- Ws %*% instruments           
#       AddH<- cbind(instruments, as.matrix(winst))
#     }
#     else AddH<- cbind(instruments)
#     
#     if(is.null(Durbin)){   
#       Zmat<- cbind(x, endog)            
#       colnames(Zmat) <- c(colnames(x), colnames(endog)) 
#     }
#     else {
#       Zmat<- cbind(x, addx, endog)            
#       colnames(Zmat) <- c(colnames(x), x.lag.names, colnames(endog)) 
#     }
#     if (K==2) {
#       if(is.null(Durbin))  Hmat<-cbind(x, wx, AddH) 
#       else Hmat<-cbind(x, wx, addx, Inx.dur, AddH)
#     }
#     else {
#       if(is.null(Durbin))  Hmat<-cbind(1, x, wx, AddH)
#       else Hmat<-cbind(1, x, wx, addx, Inx.dur, AddH)
#     }
#   }
#   else {
#     
#     if (K==2) {
#       if(is.null(Durbin)) Hmat<-cbind(x,wx)
#       else Hmat<-cbind(x, wx, addx, Inx.dur)
#     }
#     else {
#       if (is.null(Durbin)) Hmat<-cbind(1, x,wx)	
#       else Hmat<-cbind(1, x, wx, addx, Inx.dur) 
#     }
#     if(is.null((Durbin)))  Zmat<- x
#     else Zmat <- cbind(x, addx)
#   }
#   
#   firststep<-spatial.ivreg(y = y , Zmat = Zmat, Hmat = Hmat, HAC = HAC, het = het)
#   ubase<-residuals(firststep)
#   
#   if (initial.value=="SAR"){
#     Wubase<-Ws %*% ubase
#     pars<-coefficients(lm(as.numeric(ubase) ~ as.numeric(Wubase)-1))
#   }
#   else pars<-initial.value
#   
#   
#   if(het){
#     
#     Ggmat<-gg_het(Ws, ubase, n)
#     
#     optres <-nlminb(pars, optimfunct, lower= -0.9 + .Machine$double.eps , 
#                     upper= 0.9 -  .Machine$double.eps, control= control, 
#                     v = Ggmat, verbose = verbose)
#     
#     #list(abs.tol = abs.tol, rel.tol = rel.tol)
#     
#     rhotilde<-optres$par
#     # print(rhotilde)
#     
#     if(step1.c){
#       gmm.weghts1.c<-psirhorho_het(rhotilde, ubase, Hmat, Zmat, Ws, step1.c = TRUE)
#       # print(gmm.weghts1.c$Phiinv)
#       # gmm.weghts1.c<-psirhorho_het_mod(rhotilde, ubase, Hmat, Zmat, Ws2, step1.c = TRUE)
#       
#       optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat= gmm.weghts1.c$Phiinv, 
#                        verbose = verbose, lower= -0.9 + .Machine$double.eps, 
#                        upper= 0.9 -  .Machine$double.eps, control = control)	
#       
#       rhotilde<-optres$par
#       gmm.weghts1.c<-psirhorho_het(rhotilde, ubase, Hmat, Zmat, Ws, step1.c = TRUE)
#       
#       # gmm.weghts1.c<-psirhorho_het_mod(rhotilde, ubase, Hmat, Zmat, Ws2, step1.c = TRUE)
#       
#       vcmat_2sls <- Omega_het(rhotilde, gmm.weghts1.c$Pmat, gmm.weghts1.c$A1, 
#                               gmm.weghts1.c$A2, gmm.weghts1.c$a.vec1, 
#                               gmm.weghts1.c$a.vec2, Hmat, Ggmat$bigG, 
#                               gmm.weghts1.c$Phiinv, gmm.weghts1.c$epsilon, 
#                               gmm.weghts1.c$Zstar, Ws, step1.c = TRUE)
#       
#       
#       coeff_2sls <- as.matrix(c(coefficients(firststep), rhotilde))
#       rownames(coeff_2sls)<-c(colnames(Zmat), 'rho')
#       s2_2sls<-crossprod(ubase)/(n-k)
#       
#       
#       model.data<-data.frame(cbind(y,x[,-1]))
#       
#       method<-"gm spatial"
#       
#       k<-nrow(coeff_2sls)
#       R<-matrix(0,1,k)
#       R[,((k-1):k)]<-1
#       Rbeta<-R%*%coeff_2sls
#       Rvar<-R%*% vcmat_2sls$Omega %*%t(R)
#       stat<-as.numeric(t(Rbeta)%*% Rbeta/Rvar)
#       pval <- pchisq(stat,df=1,lower.tail=FALSE)
#       W<-list(stat=stat,pval=pval)
#       
#       
#       
#       results_2sls <- list(coefficients=coeff_2sls,var=vcmat_2sls$Omega, 
#                            s2=s2_2sls, call=cl, residuals=as.numeric(ubase), 
#                            model=model.data,method=method, W=W, 
#                            firststep=firststep$coefficients, init.rho = rhotilde)
#       
#       class(results_2sls)<-c("sphet", "gstsls")
#       
#     }
#     
#     
#   }
#   else{
#     
#     Ggmat<-gg_hom(Ws, ubase, n)
#     optres <- nlminb(pars, optimfunct, lower= -0.9 + .Machine$double.eps , 
#                      upper= 0.9 -  .Machine$double.eps,control = control, 
#                      v = Ggmat, verbose = verbose)
#     rhotilde<-optres$par
#     # print(rhotilde)	
#   }
#   
#   yt  <- y - rhotilde * Ws %*% y
#   wZmat <- Ws %*% Zmat
#   Zt <- Zmat - rhotilde * wZmat
#   
#   secondstep<-spatial.ivreg(y =yt , Zmat = Zt, Hmat = Hmat, het = het, HAC = HAC)
#   delta <- coefficients(secondstep)
#   utildeb <- y - Zmat %*% delta
#   
#   if(het){
#     
#     
#     Ggmat<-gg_het(Ws, utildeb, n)
#     
#     gmm.weghts<-psirhorho_het(rhotilde, utildeb, Hmat, Zmat, Ws, step1.c = FALSE)
#     
#     # gmm.weghts<-psirhorho_het_mod(rhotilde, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)
#     
#     # print(gmm.weghts$Phiinv)
#     optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, 
#                      vcmat= gmm.weghts$Phiinv, verbose = verbose, 
#                      lower= -0.9 + .Machine$double.eps , 
#                      upper= 0.9 -  .Machine$double.eps, control = control)	
#     
#     rhofin<-optres$par
#     gmm.weghts<-psirhorho_het(rhofin, utildeb, Hmat, Zmat, Ws, step1.c = FALSE)
#     
#     # gmm.weghts<-psirhorho_het_mod(rhofin, utildeb, Hmat, Zmat, Ws2, step1.c = FALSE)
#     
#     vcmat <- Omega_het(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, 
#                        gmm.weghts$A2, gmm.weghts$a.vec1, 
#                        gmm.weghts$a.vec2, Hmat, 
#                        Ggmat$bigG, gmm.weghts$Phiinv, 
#                        gmm.weghts$epsilon, gmm.weghts$Zstar, Ws, step1.c = FALSE)
#     # vcmat <- Omega_het_mod(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar, Ws2, step1.c = FALSE)
#   }
#   else{
#     
#     Ggmat<-gg_hom(Ws, utildeb, n)
#     
#     gmm.weghts<-psirhorho_hom(rhotilde, utildeb, Hmat, Zmat, Ws, Ggmat$d, Ggmat$v.vec )
#     
#     # gmm.weghts<-psirhorho_hom_mod(rhotilde, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )
#     
#     
#     optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, 
#                      vcmat= gmm.weghts$Phiinv, verbose = verbose, 
#                      lower= -0.9 + .Machine$double.eps , 
#                      upper= 0.9 -  .Machine$double.eps, control = control )	
#     
#     rhofin<-optres$par
#     # print(rhofin)
#     gmm.weghts<-psirhorho_hom(rhofin, utildeb, Hmat, Zmat, Ws, Ggmat$d, Ggmat$v.vec )
#     # gmm.weghts<-psirhorho_hom_mod(rhofin, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )
#     
#     vcmat <- Omega_hom(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, 
#                        gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, 
#                        Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar)
#     
#     # vcmat <- Omega_hom_mod(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2,gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar)
#     
#   }
#   
#   coeff <- as.matrix(c(as.numeric(delta), rhofin))
#   rownames(coeff)<-c(colnames(Zmat), 'rho')
#   s2<-crossprod(utildeb)/(n-k)
#   
#   
#   model.data<-data.frame(cbind(y,x[,-1]))
#   
#   method<-"gmm spatial"
#   
#   k<-nrow(coeff)
#   R<-matrix(0,1,k)
#   R[,((k-1):k)]<-1
#   Rbeta<-R%*%coeff
#   Rvar<-R%*% vcmat$Omega %*%t(R)
#   stat<-as.numeric(t(Rbeta)%*% Rbeta/Rvar)
#   pval <- pchisq(stat,df=1,lower.tail=FALSE)
#   W<-list(stat=stat,pval=pval)
#   
#   
#   
#   if(het && step1.c) results<-list(coefficients=coeff,var=vcmat$Omega, 
#                                    s2=s2, call=cl, residuals=as.numeric(utildeb), 
#                                    model=model.data,method=method,W=W, 
#                                    firststep=firststep$coefficients, 
#                                    init.rho = rhotilde,  twosls = results_2sls)
#   
#   else  results<-list(coefficients=coeff,var=vcmat$Omega, s2=s2, call=cl, 
#                       residuals=as.numeric(utildeb), model=model.data,
#                       method=method,W=W, firststep=firststep$coefficients, 
#                       init.rho = rhotilde)
#   
#   
#   class(results)<-c("sphet", "gstsls")# gmm error
#   
#   return(results)
#   
#   
# }
# 
# 
# 




# laggmm <- function(formula, data, listw, listw2, endog, 
#                    instruments, lag.instr, 
#                    het, verbose, na.action, HAC, cl, Durbin){
#   
#   
#   mt <- terms(formula,data = data)
#   mf <- lm(formula, data, na.action = na.action, method = "model.frame")
#   na.act <- attr(mf, 'na.action')
#   
#   y <- c(model.extract(mf, "response"))
#   x <- model.matrix(mt,mf)
#   
#   # if((Durbin == TRUE  | class(Durbin) == "formula" ) && !is.null(endog)) 
#   #   stop("Lagged explanatory variables cannot be endogenous")
#   # 
#   
#   if (length(y)!=nrow(x)) 
#     stop("x and y have different length")
#   
#   if (any(is.na(y))) 
#     stop("NAs in dependent variable")
#   if (any(is.na(x))) 
#     stop("NAs in independent variable")
#   
#   n <- nrow(x)
#   k <- ncol(x)	
#   xcolnames <- colnames(x)
#   
#   K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)
#   
#   if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
#   if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
#   if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
#   if(inherits(listw,"Matrix"))  Ws <- listw	
#   
#   if (nrow(x) != nrow(Ws))
#     stop("Input data and weights have different dimension")
#   
# if(Durbin == TRUE | class(Durbin) == "formula"  ){
#   if(class(Durbin) == "formula"){
#     xdur <- as.matrix(lm(Durbin, data, na.action=na.action, method="model.frame"))
#     inxdur <- ifelse(K == 1, which(xcolnames %in%colnames(xdur)), which(colnames(xdur)  %in% xcolnames)-1)
#     if (k > 1) {
#       wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
#       for (i in K:k) {
#         Wx <- Ws %*% x[, i]
#         wx[, (i - (K - 1))] <- as.matrix(Wx)
#       }
#       wwx <- as.matrix(Ws %*% wx)  
#       wwwx <- as.matrix(Ws  %*%  wwx[,inxdur])
#     } 
#     
#     
#     if(K==2)   Hin <- cbind(x,wx,wwx,wwwx)
#     else Hin <- cbind(1,x,wx,wwx,wwwx)
#     x <- cbind(x, wx[,indxdur])      
#   }
# 
#   else{
#     if (k > 1) {
#       wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
#       for (i in K:k) {
#         Wx <- Ws %*% x[, i]
#         wx[, (i - (K - 1))] <- as.matrix(Wx)
#       }
#       wwx <- as.matrix(Ws %*% wx)  
#       wwwx <- as.matrix(Ws  %*%  wwx)
#     } 
#     
#     if(K==2) Hin <- cbind(x,wx,wwx,wwwx)
#     else Hin <- cbind(1,x,wx,wwx,wwwx)
#     x <- cbind(x, wx)
#   }
# }
#   else{  
#   if (k > 1) {
#     wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
#     for (i in K:k) {
#       Wx <- Ws %*% x[, i]
#       wx[, (i - (K - 1))] <- as.matrix(Wx)
#     }
#     wwx <- Ws %*% wx                    					         
#   }
# if(K==2)    Hin <- cbind(x, wx, wwx)
# else        Hin <- cbind(1, x, wx, wwx)
# x <- x
# }
#   
#   
#   print(wx)
#   # if(!is.null(Durbin)){
#   #   # print(Durbin)
#   #   x.dur <- as.matrix(lm(Durbin, data, na.action = na.action, method = "model.frame"))
#   #   #print(colnames(x.dur))
#   #   #print(xcolnames)
#   #   #print(match(colnames(x.dur), xcolnamprint(sum(match(colnames(x.dur), xcolnames, nomatch = 0)))
#   #   if(sum(match(colnames(x.dur), xcolnames, nomatch = 0))!= 0) stop("Explanatory variables to be lagged cannot be specified in the main formula")    
#   #   wx.dur <- Ws %*% x.dur
#   #   wwx.dur <- Ws %*% wx.dur
#   #   wwwx.dur <- Ws %*% wwx.dur
#   #   addx <- cbind(x.dur, wx.dur)
#   #   x.lag.names <- c(colnames(x.dur), paste("lag_", colnames(x.dur), sep=""))
#   #   Inx.dur <- cbind(as.numeric(wwx.dur),as.numeric(wwwx.dur))
#   # }
#   # 
#   
#   wy<-Ws %*% y	
#   colnames(wy)<-"lambda"
#   if (!is.null(endog) && is.null(instruments)) stop("No instruments specified for the endogenous variable in the model")
#   
#   if (!is.null(endog)) {
#     
#     endog <- as.matrix(lm(endog, data, na.action=na.action, method="model.frame"))
#     instruments <- as.matrix(lm(instruments, data, na.action=na.action, method="model.frame"))
#     
#     if(lag.instr) {
#       winst <- Ws %*% instruments
#       wwinst<- Ws %*% winst	
#       AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst))        
#     }
#     else  AddH <- instruments        
# 
#     Hmat <- cbind(Hin, AddH)
#     Zmat<- cbind(x, endog, as.matrix(wy))            
#     colnames(Zmat) <- c(colnames(x), colnames(endog), colnames(wy))               
#   }
#   else {
#     Zmat<- cbind(x, as.matrix(wy))                    
#     Hmat <- Hin
#   }
#   
#   results <-spatial.ivreg(y, Zmat, Hmat, het, HAC)
#   #print(results$coefficients)
#   #print(results$var)
#   model.data <- data.frame(cbind(y, x[, -1]))
#   results$call <- cl
#   results$model <- model.data
#   results$type <- NULL
#   results$bandwidth <- NULL
#   results$method <- "gmm spatial"
#   results$HAC <- FALSE
#   class(results) <- c("sphet", "stsls_sphet") #change to lag gmm
#   
#   return(results)
#   
# }

# laggmm <- function(formula, data, listw, listw2, endog, 
#                    instruments, lag.instr, 
#                    het, verbose, na.action, HAC, cl, Durbin = FALSE){
#   
#   mt <- terms(formula,data = data)
#   mf <- lm(formula, data, na.action = na.action, method = "model.frame")
#   na.act <- attr(mf, 'na.action')
#   
#   y <- c(model.extract(mf, "response"))
#   x <- model.matrix(mt,mf)
#   
#     
#   
#   if (length(y)!=nrow(x)) 
#     stop("x and y have different length")
#   
#   if (any(is.na(y))) 
#     stop("NAs in dependent variable")
#   if (any(is.na(x))) 
#     stop("NAs in independent variable")
#   
#   n <- nrow(x)
#   k <- ncol(x)	
#   xcolnames <- colnames(x)
#   
#   K <- ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)
#   
#   if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
#   if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
#   if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
#   if(inherits(listw,"Matrix"))  Ws <- listw	
#   
#   if (nrow(x) != nrow(Ws))
#     stop("Input data and weights have different dimension")
#   
#   if (k > 1) {
#     
#     wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
#     for (i in K:k) {
#       Wx <- Ws %*% x[, i]
#       wx[, (i - (K - 1))] <- as.matrix(Wx)
#     }
#     wwx <- Ws %*% wx                    					         
#   }
#   print(wx)
#   if(!is.null(Durbin)){
#    # print(Durbin)
#     x.dur <- as.matrix(lm(Durbin, data, na.action = na.action, method = "model.frame"))
#     #print(colnames(x.dur))
#     #print(xcolnames)
#     #print(match(colnames(x.dur), xcolnamprint(sum(match(colnames(x.dur), xcolnames, nomatch = 0)))
#     if(sum(match(colnames(x.dur), xcolnames, nomatch = 0))!= 0) stop("Explanatory variables to be lagged cannot be specified in the main formula")    
#     wx.dur <- Ws %*% x.dur
#     wwx.dur <- Ws %*% wx.dur
#     wwwx.dur <- Ws %*% wwx.dur
#     addx <- cbind(x.dur, wx.dur)
#     x.lag.names <- c(colnames(x.dur), paste("lag_", colnames(x.dur), sep=""))
#     Inx.dur <- cbind(as.numeric(wwx.dur),as.numeric(wwwx.dur))
#   }
#   
#   
#   wy<-Ws %*% y	
#   colnames(wy)<-"lambda"
#   if (!is.null(endog) && is.null(instruments)) stop("No instruments specified for the endogenous variable in the model")
#   
#   if (!is.null(endog)) {
#     endog <- as.matrix(lm(endog, data, na.action=na.action, method="model.frame"))
# if(!is.null(Durbin) && (sum(match(colnames(x.dur), colnames(endog), nomatch = 0))!= 0)) stop("Lagged explanatory variables cannot be endogenous")
#     instruments <- as.matrix(lm(instruments, data, na.action=na.action, method="model.frame"))
#     if(lag.instr) {
#       winst <- Ws %*% instruments
#       wwinst<- Ws %*% winst	
#       AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst))        
#     }
#     else  AddH <- instruments        
#     if (K==2){
#       print(wx)
#       if (is.null(Durbin)) Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), AddH)
#       else Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur, AddH)
#       
#     } 
#     else {
#       if (is.null(Durbin)) Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), AddH)
#       else Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur, AddH)
#     }
#     if (is.null(Durbin)){  
#       Zmat<- cbind(x, endog, as.matrix(wy))            
#       colnames(Zmat) <- c(colnames(x), colnames(endog), colnames(wy))               
#     }
#     else{
#       Zmat<- cbind(x, addx, endog, as.matrix(wy))            
#       colnames(Zmat) <- c(colnames(x), x.lag.names, colnames(endog), colnames(wy))
#     }
#   }
#   else {
#     if(is.null((Durbin)))  Zmat<- cbind(x, as.matrix(wy))                    
#     else Zmat<- cbind(x, addx, as.matrix(wy))
#     
#     if (K==2) {
#       if(is.null((Durbin))) Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx)) 
#       else Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur)
#     }
#     else {
#       if(is.null((Durbin)))  Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx)) 
#       else Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), addx, Inx.dur)
#     }
#   }
#   
#   results <-spatial.ivreg(y, Zmat, Hmat, het, HAC)
#   #print(results$coefficients)
#   #print(results$var)
#   model.data <- data.frame(cbind(y, x[, -1]))
#   results$call <- cl
#   results$model <- model.data
#   results$type <- NULL
#   results$bandwidth <- NULL
#   results$method <- "gmm spatial"
#   results$HAC <- FALSE
#   class(results) <- c("sphet", "stsls_sphet") #change to lag gmm
#   
#   return(results)
#   
# }
