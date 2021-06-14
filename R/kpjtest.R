#' @name kpjtest
#' @aliases kpjtest 
#' @title Kelejian and Piras J-test
#' @description  The function calculate the Kelejian and Piras J-test for spatial models.
#'  Both models (under the null and under the alternative) 
#'  can be specified with additional endogenous variables, and additional instruments. 
#'  The model under the null allows for heteroskedasticity as well as spatial autocorrelation:
#'    
#'    \deqn{y=\lambda W y + X \beta + u  }
#'  \deqn{u=Re}  with \deqn{e ~ N(0,\sigma^2_i) }
#'  
#'  Note that when R reduces to an identity matrix, the error term, while still heteroskedastic, is not spatially autocorrelated. 
#'  On the other hand, when the \eqn{\sigma^2_i} are all the same (and R is an identity matrix) than the error term is neither heteroskedastic nor autocorrelated.  
#'
#' @usage kpjtest(H0model, H1model, data = list(), listw0 = NULL, listw1 = NULL, 
#'          endogH0 = NULL, endogH1 = NULL, instrumentsH0 = NULL, instrumentsH1 = NULL, 
#'          lag.instr = FALSE, model = "lag", het = FALSE, HAC = F, 
#'          distance = NULL, type = "Epanechnikov",
#'          bandwidth = "variable",  na.action = na.fail)
#'
#' @param H0model Formula object for the specification of the model under the null
#' @param H1model Formula object for the specification of the model under the alternative
#' @param data an object of class \link{data.frame}. An optional data frame containing the variables
#'    in the model
#' @param listw0 an object of class \code{listw}, \code{matrix}, or \code{Matrix}. The spatial weighting matrix under the null model
#' @param listw1 an object of class \code{listw}, \code{matrix}, or \code{Matrix}. The spatial weighting matrix under the alternative model
#' @param endogH0 additional endogenous variables under the null model. Default \code{NULL}. If not \code{NULL} should be specified as a formula with no dependent variable (endog = ~ x1 + x2). Note the ~ before the expression
#' @param endogH1 additional endogenous variables under the alternative model. Default \code{NULL}. If not \code{NULL} should be specified as a formula with no dependent variable (endog = ~ x1 + x2). Note the ~ before the expression
#' @param instrumentsH0 external instruments for the null model. Default \code{NULL}. If not \code{NULL} should be specified 
#'    as a formula with no dependent variable (instruments = ~ x1 + x2). Note the ~ before the expression
#' @param instrumentsH1 external instruments for the alternative model. Default \code{NULL}. If not \code{NULL} should be specified 
#'    as a formula with no dependent variable (instruments = ~ x1 + x2). Note the ~ before the expression
#' @param lag.instr should the external instruments be spatially lagged?
#' @param  model one of \code{lag}, or \code{sarar}. The current version of the function only implements the lag model.  
#' @param het default FALSE: if TRUE uses the methods developed for heteroskedasticity
#' @param  na.action a function which indicates what should happen when the data contains missing values. See \link{lm} for details
#' @param HAC perform the HAC estimator of Kelejian and Prucha, 2007 on the null (and augmented) model.
#' @param distance an object of class \code{distance} created for example by \link{read.gwt2dist} 
#'    The object contains the specification of the distance measure 
#'    to be employed in the estimation of the VC matrix. See Details. 
#' @param type One of \code{c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular")}. 
#'    The type of Kernel to be used. See Details. 
#' @param  bandwidth "variable" (default) - or numeric when a fixed bandwidth is specified by the user.
#' 
#' 
#' @details   In order to calculate the J-test, the function follows a few steps:
#'    \itemize{
#'      \item The alternative model is estimated by S2SLS. 
#'      \item Based on the estimated parameters in the previous step,  obtain a prediction based on the alternative models of the dependent vector in the null model. The predictor is based on the right hand side of the model. 
#'      \item Use these predicted values of the dependent variable based on the alternative models into the null model to obtain the augmented model. 
#'      \item Estimate the augmented model by 2SLS using all of the instruments relating to the null model as well as all of the instruments relating to the alternative models. 
#'      \item Test for the statistical significance of the predicted value. If it is not significant, accept the null model. If it is  significant, reject the null and conclude that the true model is the alternative models. 
#'    }
#'  The output is an object of class \code{sphet} where the last row of the table of coefficients is the prediction. 
#'  
#'  When the model is heteroskedastic as well as spatially autocorrelated, an HAC procedure is employed. 
#'  The default sets the bandwith for each observation to the maximum distance for that observation (i.e.
#'                                                                                                   the max of each element of the list of distances). 
#'  
#'  Six different kernel functions are implemented:
#'    \itemize{
#'      \item \code{'Epanechnikov'}: \eqn{K(z) = 1-z^2}
#'      \item \code{'Triangular'}: \eqn{K(z) = 1-z} 
#'      \item \code{'Bisquare'}: \eqn{K(z) = (1-z^2)^2}
#'      \item \code{'Parzen'}: \eqn{K(z) = 1-6z^2+6 |z|^3} if  \eqn{z \leq 0.5} and  
#'      \eqn{ K(z) = 2(1-|z|)^3} if \eqn{0.5 < z \leq 1} 
#'      \item \code{'TH'} (Tukey - Hanning):  \eqn{ K(z) = \frac{1+ \cos(\pi z)}{2}}
#'      \item \code{'Rectangular'}:  \eqn{ K(z) = 1}
#'      \item \code{'QS'} (Quadratic Spectral): \eqn{K(z) = \frac{25}{12\pi^2z^2} 
#'      (\frac{\sin(6\pi z)/5)}{6\pi z/5} - \cos(6\pi z)/5)}). 
#'    }
#'  
#'  If the kernel type is not one of the six implemented, the function will terminate with an error message. 
#'  The spatial two stage least square estimator is based on the matrix of instruments \eqn{H=[X,WX,W^2X^2]}.
#'
#'
#' @return A list object of class \code{sphet}
#'  \item{coefficients}{Generalized Spatial two stage least squares coefficient estimates of \eqn{\delta} and GM estimator for \eqn{\rho}. }
#'  \item{var}{variance-covariance matrix of the estimated coefficients}
#'  \item{s2}{GS2SLS residuals variance}
#'  \item{residuals}{GS2SLS residuals}
#'  \item{yhat}{difference between GS2SLS residuals and response variable}
#'  \item{call}{the call used to create this object}
#'  \item{model}{the model matrix of data}
#'  \item{method}{\code{'s2slshac'}}
#'  
#' @references
#'  
#'  Kelejian and Piras (2017). \emph{Spatial Econometrics}. Academic Press. ISBN: 978-0-12-813387-3 
#'  
#'  
#'  Gianfranco Piras (2010). sphet: Spatial Models with Heteroskedastic Innovations in R. \emph{Journal of Statistical Software}, 35(1), 1-21. \url{https://www.jstatsoft.org/v35/i01/}.
#'  
#'  Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of Estimation Methods for Spatial Econometrics. \emph{Journal of Statistical Software}, 63(18), 1-36. \url{https://www.jstatsoft.org/v63/i18/}.
#'  
#' @author  Gianfranco Piras \email{gpiras@mac.com}
#' 
#' @examples 
#'  library(spdep)
#'  library(sphet)
#'  data(boston)
#'  boslw <- nb2listw(boston.soi)
#'  
#'  Bos.Knn <- knearneigh(boston.utm, k = 5)
#'  bos.nb <- knn2nb(Bos.Knn)
#'  boslw2 <- nb2listw(bos.nb)
#'  
#'  fm <- log(MEDV) ~ CRIM + ZN + INDUS + CHAS
#'  fm2 <- log(MEDV) ~ CRIM + ZN + INDUS + RM + AGE 
#'  
#'  test <- kpjtest(fm, fm2, data = boston.c, 
#'  listw0 = boslw, listw1 = boslw2, model = "lag")
#'  
#' @keywords spatial
#' @export








# The assumption is that listw1 and listw2 are different, otherwise there would be no reason to use the spatial J-test

kpjtest <- function(H0model, H1model, data = list(), listw0 = NULL, 
                    listw1 = NULL, endogH0 = NULL, endogH1 = NULL, 
                    instrumentsH0 = NULL, instrumentsH1 = NULL, 
                    lag.instr = FALSE, model = "lag", 
                    het = FALSE, HAC = F, distance = NULL, 
                    type = "Epanechnikov", bandwidth = "variable",  
                    na.action = na.fail){
	 cl <- match.call()
	if(is.null(listw0)) 
		stop("listw for the null model was not specified")
	if(is.null(listw1))
		stop("listw for the alternative model was not specified")
	
	
	 if (!inherits(listw0, c("listw", "Matrix", "matrix"))) 
            stop("listw format unknown")
        if (inherits(listw0, "listw"))  Ws0 <- listw2dgCMatrix(listw0)
        if (inherits(listw0, "matrix")) Ws0 <- Matrix(listw0)
        if (inherits(listw0, "Matrix")) Ws0 <- listw0

	 if (!inherits(listw1, c("listw", "Matrix", "matrix"))) 
            stop("listw format unknown")
        if (inherits(listw1, "listw"))  Ws1 <- listw2dgCMatrix(listw1)
        if (inherits(listw1, "matrix")) Ws1 <- Matrix(listw1)
        if (inherits(listw1, "Matrix")) Ws1 <- listw1
	
	if((dim(Ws0)[1] != dim(Ws1)[1]) &&  (dim(Ws0)[2] != dim(Ws1)[2])) 
			stop("listw0 and listw1 must have the same dimension")
			
	if( all(Ws0 == Ws1)) 
			stop("listw0 and listw1 cannot be the same")		
	
	
if(model == "lag") res <- jtestlag(H0model, H1model, data = data, 
                                   listw0 = Ws0, listw1 = Ws1, 
                                   endogH0 = endogH0, endogH1 = endogH1, 
                                   instrumentsH0 = instrumentsH0, 
                                   instrumentsH1 = instrumentsH1, 
                                   lag.instr = lag.instr, model = model,  
                                   HAC = HAC, distance = distance, 
                                   type = type, bandwidth = bandwidth,  
                                   na.action = na.action, het = het, cl = cl)


if(model == "sarar") stop ("Method for sarar not yet implemented")
	
	
	return(res)
}



jtestlag <- function(H0model, H1model, data, listw0, listw1, endogH0, endogH1, instrumentsH0, instrumentsH1, lag.instr, model,  HAC, distance, type, bandwidth,  na.action, het, cl){
	
Alt <- spreg(H1model, data = data, listw = listw1, 
             endog = endogH1, instruments = instrumentsH1, 
             lag.instr = lag.instr, model = "lag", het = F, Durbin = FALSE)	

mt <- terms(H1model, data = data)
mf <- lm(H1model, data, na.action = na.action, method = "model.frame")
y1 <- c(model.extract(mf, "response"))
x1 <- model.matrix(mt, mf)
wy1 <- as.matrix(listw1 %*% y1)
reg <- as.matrix(cbind(x1, wy1))
yp <- reg %*% coefficients(Alt)
data$yp <- yp
# H0fm <- as.formula(paste(names(mf)[1], paste(names(mf)[-1], collapse = " + "), sep = " ~" ))

H0 <- augmented(H0model, data = data, listw0 = listw0, listw1 = listw1, yp = yp, x1 = x1, endogH0 = endogH0, endogH1 = endogH1, instrumentsH0 = instrumentsH0, instrumentsH1 = instrumentsH1, lag.instr = lag.instr, model = model,  HAC = HAC, distance = distance, type = type, bandwidth = bandwidth,  na.action = na.fail, het = het, cl = cl)	
# print(H0)
return(H0)
	
}

# # jtestsarar <- function(){
	
	
# }



augmented <- function(H0model, data, listw0, listw1, yp, x1, endogH0, endogH1, instrumentsH0, instrumentsH1, lag.instr, model,  HAC, distance, type, bandwidth,  na.action, het, cl){
	
	
	    if (!is.null(endogH0) && is.null(instrumentsH0)) 
        stop("No instruments specified for the endogenous variable in the null model")

	
				if(colnames(x1)[1] == "(Intercept)")	 x1 <- x1[,-1]
				else x1 <- x1

						mt <- terms(H0model, data = data)
						mf <- lm(H0model, data, na.action = na.action, method = "model.frame")

						y0 <- c(model.extract(mf, "response"))
						x0 <- model.matrix(mt, mf)
						w0y0 <- listw0 %*% y0
						
				if(colnames(x0)[1] == "(Intercept)")	 x0 <- x0[,-1]
				else x0 <- x0
						
						xaug <- cbind(x0 , x1[,-which(colnames(x1) %in% colnames(x0))])
						w0x0 <- listw0 %*% x0
						w02x0 <- listw0 %*% w0x0
						w1x1 <- listw1 %*% x1
						w12x1 <- listw1 %*% w1x1

	

		if(is.null(instrumentsH0) && is.null(instrumentsH1) ){
		
						Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1)	
		}


		if(!is.null(instrumentsH0) && !is.null(instrumentsH1) ){
			
			instrumentsH0 <- as.matrix(lm(instrumentsH0, data, na.action = na.action, method = "model.frame"))
			instrumentsH1 <- as.matrix(lm(instrumentsH1, data, na.action = na.action, method = "model.frame"))
			instaug <- cbind(instrumentsH0 , instrumentsH1[,-which(colnames(instrumentsH1) %in% colnames(instrumentsH0))])
			
			if(lag.instr){
				w0inH0 <- listw0 %*% instrumentsH0
				ww0inH0 <- listw0 %*% w0inH0
				w1inH1 <- listw1 %*% instrumentsH1
				ww1inH1 <- listw1 %*% w1inH1
				Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instaug, w0inH0, ww0inH0, w1inH1, ww1inH1)			
			}	
							Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instaug )			
		}


		if(!is.null(instrumentsH0) && is.null(instrumentsH1) ){
						
						instrumentsH0 <- as.matrix(lm(instrumentsH0, data, na.action = na.action, method = "model.frame"))
						
						if(lag.instr){
				w0inH0 <- listw0 %*% instrumentsH0
				ww0inH0 <- listw0 %*% w0inH0
				Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instrumentsH0, w0inH0, ww0inH0)			
			}
		Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instrumentsH0)				
		}
	
	
			if(is.null(instrumentsH0) && !is.null(instrumentsH1) ){
						
						instrumentsH1 <- as.matrix(lm(instrumentsH1, data, na.action = na.action, method = "model.frame"))
						
						if(lag.instr){
				w1inH1 <- listw1 %*% instrumentsH1
				ww1inH1 <- listw1 %*% w1inH1
				Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instrumentsH1, w1inH1, ww1inH1)			
			}
		Haug <- cbind(1, xaug, w0x0, w02x0, w1x1, w12x1, instrumentsH1)				
					
		}


Haug <- as.matrix(Haug)

if(!is.null(endogH0)){
	endogH0 <- as.matrix(lm(endogH0, data, na.action = na.action, method = "model.frame"))
	Zaug <- cbind(1, x0, as.matrix(w0y0), endogH0, as.matrix(yp))
}
else  Zaug <- cbind(1, x0, as.matrix(w0y0), as.matrix(yp))
colnames(Zaug) <- c("(Intercept)", colnames(x0), "w0y", "yp")

results <- spatial.ivreg(y0, Zaug, Haug, het = het, HAC = HAC, distance = distance, type = type, bandwidth = bandwidth)
    
 model.data <- data.frame(cbind(as.matrix(y0), as.matrix(Zaug)))
        results$call <- cl
        results$model <- model.data
        results$type <- type
        results$bandwidth <- bandwidth
        results$method <- "s2slshac"
        results$HAC <- HAC
         class(results) <- c("sphet", "stsls_sphet")
# print(summary(results))

# df <- nrow(Zaug) - ncol(Zaug)
# HHaug <- crossprod(Haug)
# HHaugi <- solve(HHaug)
# fp <- Haug %*% HHaugi
# sp <- crossprod(Haug,Zaug)
# Zaugp <- fp %*% sp
# daug <- solve(crossprod(Zaugp), crossprod(Zaugp, y0))

# raug <- y0 - Zaug %*% daug

# s2 <- as.numeric(crossprod(raug)/df)
# vc <- s2 * solve(crossprod(Zaugp))
# print(ret)

# class(results) <- c("KP_jtest", "sphet", "stsls_sphet")
return(results)
}



# summary.KP_jtest <- function(x, digits= max(3, getOption("digits") - 2),...){
# cp <- coefficients(x)
# vc <- x$var
# test <- cp[length(cp)]^2/vc[length(cp),length(cp)]

# STAT <- qchisq(0.05, 1, lower.tail=FALSE)

# # ret <- list(test, STAT)	
	
# cat("\n Kelejian and Piras J-test for spatial models:\n")
# cat("\n H0 vs H1:\n")
# cat("\n ---------------------------------------------\n")
# cat("\n Stat:\n")
# print(test)
# cat("\n Chi-sq (df = 1):\n")
# print(STAT)

# }