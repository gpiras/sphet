spreg<-function(formula, data=list(), listw, listw2=NULL, 
                endog = NULL, instruments= NULL, lag.instr = FALSE, 
                initial.value=0.2, 
                model = c("sarar", "lag", "error", "ivhac", "ols"), 
                het = FALSE, verbose=FALSE, na.action = na.fail,  
                HAC = FALSE, distance = NULL, 
                type =  c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"), bandwidth="variable", 
                step1.c = FALSE, control = list(), Durbin = FALSE){
 
 
  cl = match.call()
switch(match.arg(model),
       sarar = sarargmm(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                      instruments = instruments, lag.instr = lag.instr, initial.value = initial.value, 
                      het = het, verbose = verbose, na.action = na.action,
                      step1.c = step1.c, control = control, HAC = HAC, cl = cl, Durbin = Durbin),
       lag = laggmm(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                    instruments = instruments, lag.instr = lag.instr, 
                     het = het, verbose = verbose, na.action = na.action, HAC = HAC, cl = cl, Durbin = Durbin),
        error = errorgmm(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                        instruments = instruments, lag.instr = lag.instr, initial.value = initial.value, 
                        het = het, verbose = verbose, na.action = na.action,
                        step1.c = step1.c, control = control, HAC = HAC, cl = cl, Durbin = Durbin),
        ivhac = laghac(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                       instruments = instruments, lag.instr = lag.instr,  verbose = verbose, 
                       na.action = na.action, het = het, HAC = HAC, distance = distance, 
                       type = type, bandwidth = bandwidth, cl = cl, Durbin = Durbin),
        ols = olshac(formula = formula, data = data, listw = listw, 
                     endog = endog, instruments= instruments,  
                     na.action = na.action, het = het, HAC = HAC, distance = distance, 
                     type = type, bandwidth = bandwidth, cl = cl, Durbin = Durbin),
        stop("Argument model incorrectly specified")
  )

  }


