#' @name spreg
#' @aliases spreg
#' @title GM estimation of a Cliff-Ord type model with Heteroskedastic Innovations
#' @description Multi step GM/IV estimation of a linear Cliff and Ord -type of model 
#' of the form:
#' \deqn{y=\lambda W y + X \beta + u  }
#' \deqn{u=\rho W u + e}  with \deqn{e ~ N(0,\sigma^2_i) }
#' 
#' The model allows for spatial lag in the dependent variable
#' and disturbances. The innovations in the disturbance process are assumed 
#' heteroskedastic of an unknown form.
#' @usage
#' spreg(formula, data = list(), listw, listw2 = NULL, 
#'         endog = NULL, instruments = NULL, 
#'         lag.instr = FALSE, initial.value = 0.2, q = 2,
#'         model = c("sarar", "lag", "error", "ivhac", "ols"),
#'         het = FALSE, verbose = FALSE, 
#'         na.action = na.fail,  HAC = FALSE, 
#'         distance = NULL, type =  c("Epanechnikov","Triangular",
#'                                    "Bisquare", "Parzen", "QS", "TH","Rectangular"), 
#'         bandwidth = "variable" , step1.c = FALSE,
#'         control = list(), Durbin = FALSE)
#' 
#' @param formula a description of the model to be fit 
#' @param data an object of class \link{data.frame}. An optional data frame containing the variables 
#' in the model
#' @param listw an object of class \code{listw}, \code{matrix}, or \code{Matrix}
#' @param listw2 an object of class \code{listw}, \code{matrix}, or \code{Matrix} specified only when \code{sarar} is true
#' @param endog additional endogenous variables. Default \code{NULL}. If not \code{NULL} should be specified
#' as a formula with no dependent variable (endog = ~ x1 + x2). Note the ~ before the expression. 
#' @param instruments external instruments. Default \code{NULL}. If not \code{NULL} should be specified 
#' as a formula with no dependent variable (instruments = ~ x1 + x2). Note the ~ before the expression.
#' @param lag.instr should the external instruments be spatially lagged?
#' @param initial.value The initial value for \eqn{\rho}. It can be either numeric (default is 0.2) or
#' set to \code{'SAR'}, in which case the optimization will start from the estimated coefficient of a regression of the 2SLS 
#' residuals over their spatial lag (i.e. a spatial AR model)
#' @param q default equal 2, if 1 the only WX is considered in matrix of instruments
#' @param model one of \code{lag}, \code{error}, \code{sarar}, \code{ivhac}, or \code{ols}. 
#' If HAC is TRUE, model should be one of \code{ivhac}, or \code{ols}.
#' @param het default FALSE: if TRUE uses the methods developed for heteroskedasticity
#' @param verbose print optimization details
#' @param na.action a function which indicates what should happen when the data contains missing values. See \link{lm} for details.
#' @param HAC perform the HAC estimator of Kelejian and Prucha, 2007.
#' @param distance an object of class \code{distance} created for example by \link{read.gwt2dist}
#' The object contains the specification of the distance measure 
#' to be employed in the estimation of the VC matrix. See Details. 
#' @param type One of \code{c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular")}.
#' The type of Kernel to be used. default = \code{"Epanechnikov"} See Details.
#' @param bandwidth "variable" (default) - or numeric when a fixed bandwidth is specified by the user.
#' @param step1.c Should step 1.c from Arraiz et al. 2012 be performed?
#' @param control A list of control arguments. See \link{nlminb}
#' @param Durbin Should (some of) the regressors be lagged? Default FALSE. If not \code{FALSE} should be specified 
#' as a formula with no dependent variable (Durbin = ~ x1 + x2) or set to \code{TRUE}. See details. 
#' @details
#' The procedure consists of two steps alternating GM and IV estimators. Each step consists of sub-steps.
#' In step one \eqn{\delta = [\beta',\lambda]'} is estimated by 2SLS. The 2SLS residuals are first employed
#' to obtain an consistent GM estimator of \eqn{\rho}.
#' 
#' In step two, the spatial Cochrane-Orcutt transformed model is estimated by 2SLS. This corresponds to a GS2SLS procedure. 
#' The GS2SLS residuals are used to obtain a consistent and efficient GM estimator for \eqn{\rho}. 
#' 
#' The initial value for the optimization in step 1b is taken to be \code{initial.value}. 
#' The initial value for the optimization of step 2b is the optimal parameter of step 1b.
#' 
#' Internally, the object of class \code{listw} is transformed into a \link{Matrix} 
#' using the function \link{listw2dgCMatrix}.
#' 
#' For the HAC estimator (Kelejian and Prucha, 2007),  there are four possibilities:
#' 
#'  \itemize{
#'        \item A model with only Wy
#'              \item A model with Wy and additional endogenous 
#'                    \item Additional endogenous variables but no Wy
#'                          \item No additional endogenous variables (A linear model with HAC estimation)
#'                          }
#' In the first two cases, the \code{model} should be \code{"ivhac"}, 
#' in the last two cases, the \code{model} should be \code{"ols"}. 
#'                          
#' Furthermore, the default sets the bandwidth for each observation to the maximum distance for that observation (i.e.
#' the max of each element of the list of distances). 
#' 
#' Six different kernel functions are implemented:
#'   \itemize{
#'         \item \code{'Epanechnikov'}: \eqn{K(z) = 1-z^2}
#'             \item \code{'Triangular'}: \eqn{K(z) = 1-z} 
#'                 \item \code{'Bisquare'}: \eqn{K(z) = (1-z^2)^2}
#'                     \item \code{'Parzen'}: \eqn{K(z) = 1-6z^2+6 |z|^3} if  \eqn{z \leq 0.5} and  
#'                           \eqn{ K(z) = 2(1-|z|)^3} if \eqn{0.5 < z \leq 1} 
#'                                 \item \code{'TH'} (Tukey - Hanning):  \eqn{ K(z) = \frac{1+ \cos(\pi z)}{2}}
#'                                       \item \code{'Rectangular'}:  \eqn{ K(z) = 1}
#'                                             \item \code{'QS'} (Quadratic Spectral): \eqn{K(z) = \frac{25}{12\pi^2z^2} 
#'                                                   (\frac{\sin(6\pi z)/5)}{6\pi z/5} - \cos(6\pi z)/5)}). 
#'                                                     }
#' If the kernel type is not one of the six implemented, the function will terminate with an error message.
#' @return  A list object of class \code{sphet}
#' \item{coefficients}{Generalized Spatial two stage least squares coefficient estimates of \eqn{\delta} and GM estimator for \eqn{\rho}. }
#' \item{var}{variance-covariance matrix of the estimated coefficients}
#'   \item{s2}{GS2SLS residuals variance}
#'     \item{residuals}{GS2SLS residuals}
#'       \item{yhat}{difference between GS2SLS residuals and response variable}
#'         \item{call}{the call used to create this object}
#'           \item{model}{the model matrix of data}
#'             \item{method}{\code{'gs2slshac'}}
#'             
#' @seealso  \code{\link{stslshac}}
#' @references
#' Arraiz, I. and Drukker, M.D. and Kelejian, H.H. and Prucha, I.R. (2010) 
#'   A spatial Cliff-Ord-type Model with Heteroskedastic Innovations: Small and Large Sample Results,
#'     \emph{Journal of Regional Sciences}, \bold{50}, pages 592--614.
#'     
#'         Drukker, D.M. and Egger, P. and Prucha, I.R. (2013)
#'           On Two-step Estimation of a Spatial Auto regressive Model with Autoregressive
#'             Disturbances and Endogenous Regressors, 
#'               \emph{Econometric Review}, \bold{32}, pages 686--733. 
#'                 
#'                 
#'                 Kelejian, H.H. and Prucha, I.R. (2010) 
#'                   Specification and Estimation of Spatial Autoregressive Models with Autoregressive and Heteroskedastic Disturbances,
#'                     \emph{Journal of Econometrics}, \bold{157}, pages 53--67. 
#'                     
#'                         Kelejian, H.H. and Prucha, I.R. (1999) 
#'                           A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model,
#'                             \emph{International Economic Review}, \bold{40}, pages 509--533.
#'                             
#'                                 Kelejian, H.H. and Prucha, I.R. (1998) 
#'                                   A Generalized Spatial Two Stage Least Square Procedure for Estimating a Spatial Autoregressive
#'                                     Model with Autoregressive Disturbances,
#'                                       \emph{Journal of Real Estate Finance and Economics}, \bold{17}, pages 99--121.
#'                                       
#'                                         
#'                                            Gianfranco Piras (2010). sphet: Spatial Models with Heteroskedastic Innovations in R. \emph{Journal of Statistical Software}, 35(1), 1-21. \doi{10.18637/jss.v035.i01}.
#'                                            
#'                                                Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of Estimation Methods for Spatial Econometrics. \emph{Journal of Statistical Software}, 63(18), 1-36. \doi{10.18637/jss.v063.i18}.
#'                                                
#'                                                  Gianfranco Piras, Paolo Postiglione (2022).  A deeper look at impacts in spatial Durbin model with sphet. \emph{Geographical Analysis}, 54(3), 664-684. 
#'@author  Gianfranco Piras \email{gpiras@mac.com}
#'@examples
#' data(columbus, package="spdep")
#' listw <- spdep::nb2listw(col.gal.nb)
#' res <- spreg(CRIME ~ HOVAL + INC, data = columbus , listw = listw,
#'              het = TRUE, verbose = FALSE, model = "sarar")
#' summary(res)
#' Effects <- impacts(res, listw = listw, R = 1000)
#' 
#' library(spdep)
#' data("baltimore", package = "spData")
#' mat <- nb2listw(knn2nb(knearneigh(cbind(baltimore$X,baltimore$Y), 3)))
#' 
#' knb10 <- knn2nb(knearneigh(cbind(baltimore$X,baltimore$Y), k=5))
#' dists <- nbdists(knb10, cbind(baltimore$X,baltimore$Y))
#' k10lw <- nb2listw(knb10, glist=dists, style="B")
#' class(k10lw) <- "distance"
#' 
#' # OLS MODEL
#' res <- spreg(PRICE ~ NROOM  +AGE, data = baltimore, listw = mat, 
#' verbose = FALSE, model = "ols", Durbin = TRUE, HAC = TRUE, 
#' distance = k10lw, type = "Triangular")
#' summary(res)
#' 
#' # note model = "ols" but with endogenous variables 
#' res <- spreg(PRICE ~ NROOM  +AGE, data = baltimore, listw = mat, 
#' verbose = FALSE, model = "ols", Durbin = TRUE, HAC = TRUE, 
#' distance = k10lw, type = "Triangular", endog = ~SQFT, 
#' instruments = ~GAR + PATIO)
#' summary(res)
#' 
#' # ERROR MODEL
#' res <- spreg(PRICE ~ NROOM  +AGE, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "error", Durbin = FALSE)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE + SQFT + NBATH, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "error", Durbin = TRUE)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE + SQFT + NBATH, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "error", Durbin = ~SQFT + NBATH)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE +  NBATH, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "error", Durbin = ~SQFT + NBATH)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "error", Durbin = ~SQFT + NBATH)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE -1, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "error", Durbin = ~SQFT + NBATH)
#' summary(res)
#' 
#' # LAG MODEL
#' res <- spreg(PRICE ~ NROOM  +AGE, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "lag", Durbin = FALSE)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE + SQFT + NBATH, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "lag", Durbin = TRUE)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE + SQFT + NBATH, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "lag", Durbin = ~SQFT + NBATH)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE +  NBATH, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "lag", Durbin = ~SQFT + NBATH)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "lag", Durbin = ~SQFT + NBATH)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE -1, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "lag", Durbin = ~SQFT + NBATH)
#' summary(res)
#' 
#' 
#' # IVHAC MODEL
#' res <- spreg(PRICE ~ NROOM  +AGE, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "ivhac", Durbin = FALSE,
#' HAC = TRUE, distance = k10lw, type = "Triangular", endog = ~SQFT, 
#' instruments = ~GAR + PATIO)
#' 
#' # SARAR MODEL
#' res <- spreg(PRICE ~ NROOM  +AGE, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "sarar", Durbin = FALSE, q = 1)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE + SQFT + NBATH, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "sarar", Durbin = TRUE, q = 1)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE + SQFT + NBATH, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "sarar", Durbin = ~SQFT + NBATH, q =2)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE +  NBATH, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "sarar", Durbin = ~SQFT + NBATH)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "sarar", Durbin = ~SQFT + NBATH)
#' summary(res)
#' 
#' res <- spreg(PRICE ~ NROOM + AGE -1, data = baltimore, listw = mat, 
#' het = TRUE, verbose = FALSE, model = "sarar", Durbin = ~SQFT + NBATH)
#' summary(res)
#' 
#' summary(res)
#' 
#' @keywords  models
#' @export
#' @importFrom methods as new
#' @importFrom stats model.matrix model.response aggregate  coef coefficients fitted lm sd model.extract na.fail nlminb pchisq
#' pnorm printCoefmat quantile residuals terms
#' @importFrom utils read.table write.table
#' @import nlme
#' @importFrom spdep  card  lag.listw  nb2listw  sym.attr.nb vi2mrc
#' @import spatialreg
#' @import Matrix
#' @import mvtnorm
#' @import sf 
#' @import spData
#' @importFrom sp spDists
#' @importFrom stringr str_remove
#' @importFrom coda as.mcmc










spreg<-function(formula, data=list(), listw, listw2=NULL, 
                endog = NULL, instruments= NULL, lag.instr = FALSE, 
                initial.value=0.2, q = 2,
                model = c("sarar", "lag", "error", "ivhac", "ols"), 
                het = FALSE, verbose=FALSE, na.action = na.fail,  
                HAC = FALSE, distance = NULL, 
                type =  c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"), bandwidth="variable", 
                step1.c = FALSE, control = list(), Durbin = FALSE){

  cl = match.call()
switch(match.arg(model),
       sarar = sarargmm(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                      instruments = instruments, lag.instr = lag.instr, initial.value = initial.value, 
                      q = q, het = het, verbose = verbose, na.action = na.action,
                      step1.c = step1.c, control = control, HAC = HAC, cl = cl, Durbin = Durbin),
       lag = laggmm(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                    instruments = instruments, lag.instr = lag.instr, q = q,
                     het = het, verbose = verbose, na.action = na.action, HAC = HAC, cl = cl, Durbin = Durbin),
        error = errorgmm(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                        instruments = instruments, lag.instr = lag.instr, initial.value = initial.value, q = q,
                        het = het, verbose = verbose, na.action = na.action,
                        step1.c = step1.c, control = control, HAC = HAC, cl = cl, Durbin = Durbin),
        ivhac = laghac(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog, 
                       instruments = instruments, lag.instr = lag.instr,  q = q, verbose = verbose, 
                       na.action = na.action, het = het, HAC = HAC, distance = distance, 
                       type = type, bandwidth = bandwidth, cl = cl, Durbin = Durbin),
        ols = olshac(formula = formula, data = data, listw = listw, 
                     endog = endog, instruments= instruments,  q = q,
                     na.action = na.action, het = het, HAC = HAC, distance = distance, 
                     type = type, bandwidth = bandwidth, cl = cl, Durbin = Durbin),
        stop("Argument model incorrectly specified")
  )

  }

# coefficients.sphet <- function(x) drop(as.numeric(x$coefficients))
# vcov.sphet <- function(x) x$var
