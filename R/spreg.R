spreg<-function(formula, data=list(), listw, listw2=NULL, 
                endog = NULL, instruments= NULL, lag.instr = FALSE, 
                initial.value=0.2, 
                model = c("sarar", "lag", "error", "ivhac", "ols"), 
                het = FALSE, verbose=FALSE, na.action = na.fail,  
                HAC = FALSE, distance = NULL, 
                type =  c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"), bandwidth="variable", 
                step1.c = FALSE, control = list()){
 
  
switch(match.arg(model),
       sarar = sarargmm(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                      instruments = instruments, lag.instr = lag.instr, initial.value = initial.value, 
                      het = het, verbose = verbose, na.action = na.action,
                      step1.c = step1.c, control = control, HAC = HAC),
       lag = laggmm(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                    instruments = instruments, lag.instr = lag.instr, 
                     het = het, verbose = verbose, na.action = na.action, HAC = HAC),
        error = errorgmm(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                        instruments = instruments, lag.instr = lag.instr, initial.value = initial.value, 
                        het = het, verbose = verbose, na.action = na.action,
                        step1.c = step1.c, control = control, HAC = HAC),
        ivhac = laghac(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                       instruments = instruments, lag.instr = lag.instr,  verbose = verbose, 
                       na.action = na.action, het = het, HAC = HAC, distance = distance, 
                       type = type, bandwidth = bandwidth),
        ols = olshac(formula = formula, data = data, endog = endog, 
                     instruments= instruments,  
                     na.action = na.action, het = het, HAC = HAC, distance = distance, 
                     type = type, bandwidth = bandwidth),
        stop("Argument model incorrectly specified")
  )

  }


impacts.gstsls <- function(obj, ..., tr=NULL, R=NULL, listw=NULL, evalues=NULL,
    tol=1e-6, empirical=FALSE, Q=NULL) {
    if (!is.null(obj$listw_style) && obj$listw_style != "W") 
        stop("Only row-standardised weights supported")
coefs <- drop(obj$coefficients)
    p2 <- length(coefs)
    rho <- coefs[(p2-1)]
    beta <- coefs[1:(p2-2)]
    p <- length(beta)
	p1 <- p + 1
    lambda <- coefs[p2]
# rho is lag coef., lambda is error coef (reversed from function use)
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0
    if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
    } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
    }
    n <- length(obj$residuals)
    mu <- c(rho, beta)
    Sigmawor <- obj$var[-p2,-p2]
	Sigma <- Sigmawor[c(p1, (1:(p1-1))), c(p1, (1:(p1-1)))]    
	irho <- 1
    drop2beta <- 1
    if (!requireNamespace("spatialreg", quietly=TRUE))
        stop("install spatialreg")
    res <- spatialreg::intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu,
        Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
        interval=NULL, type="lag", tr=tr, R=R, listw=listw, evalues=evalues,
        tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "iClass") <- class(obj)
    res
}






impacts.stsls_sphet <- function(obj, ..., tr=NULL, R=NULL, listw=NULL,
    evalues=NULL, tol=1e-6, empirical=FALSE, Q=NULL) {
    if (!is.null(obj$listw_style) && obj$listw_style != "W") 
        stop("Only row-standardised weights supported")
    coefs <- drop(obj$coefficients)
    p2 <- length(coefs)
    rho <- coefs[(p2)]
    beta <- coefs[1:(p2-1)]
# rho is lag coef., lambda is error coef (reversed from function use)
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0
    if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
    } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
    }
    p <- length(beta)
    n <- length(obj$residuals)
    mu <- c(rho, beta)
    Sigma <- obj$var[c(p2, (1:(p2-1))), c(p2, (1:(p2-1)))]
    irho <- 1
    drop2beta <- 1
    if (!requireNamespace("spatialreg", quietly=TRUE))
        stop("install spatialreg")
    res <- spatialreg::intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu,
        Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
        interval=NULL, type="lag", tr=tr, R=R, listw=listw, evalues=evalues,
        tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "iClass") <- class(obj)
    res
}
