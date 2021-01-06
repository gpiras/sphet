#' @name stslshac
#' @aliases stslshac
#' @aliases tsls
#' @title Spatial two stages least square with HAC standard errors
#' @description Non-parametric heteroskedasticity and autocorrelation 
#'   consistent (HAC) estimator of the variance-covariance (VC) for a vector of sample moments within a 
#'   spatial context. The disturbance vector is generated as follows:
#'   \deqn{ u = R  \epsilon }
#'   where \eqn{R} is a non-stochastic matrix.
#' @usage stslshac(formula, data = list(), listw,
#'              na.action = na.fail, zero.policy = NULL, HAC = TRUE,
#'              distance = NULL, type = "Epanechnikov", 
#'              bandwidth = "variable", W2X = TRUE)
#' @param formula a description of the model to be fit 
#' @param data an object of class \link{data.frame}. An optional data frame containing the variables
#' in the model.
#' @param listw an object of class \code{listw} created for example by \code{nb2listw} 
#' @param distance an object of class \code{distance} created for example by \link{read.gwt2dist} 
#' The object contains the specification of the distance measure 
#' to be employed in the estimation of the VC matrix. See Details. 
#' @param type One of \code{c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH")}.
#' The type of Kernel to be used. See Details. 
#' @param na.action a function which indicates what should happen when the data contains missing values.
#' See \link{lm} for details.
#' @param zero.policy See \code{lagsarlm} for details
#' @param bandwidth "variable" (default) - or numeric when a fixed bandwidth is specified by the user.
#' @param HAC if FALSE traditional standard errors are provided.
#' @param W2X default TRUE. if FALSE only WX are used as instruments in the spatial two stage least squares.
#' @details 
#' The default sets the bandwith for each observation to the maximum distance for 
#' that observation (i.e.
#' the max of each element of the list of distances). 
#' 
#' Six different kernel functions are implemented:
#' \itemize{
#' \item \code{'Epanechnikov'}: \eqn{K(z) = 1-z^2}
#' \item \code{'Rectangular'}: \eqn{K(z) = 1}
#' \item \code{'Triangular'}: \eqn{K(z) = 1-z} 
#' \item \code{'Bisquare'}: \eqn{K(z) = (1-z^2)^2}
#' \item \code{'Parzen'}: \eqn{K(z) = 1-6z^2+6 |z|^3} if  \eqn{z \leq 0.5} and  
#' \eqn{ K(z) = 2(1-|z|)^3} if \eqn{0.5 < z \leq 1} 
#' \item \code{'TH'} (Tukey - Hanning):  \eqn{ K(z) = \frac{1+ \cos(\pi z)}{2}}
#' \item \code{'QS'} (Quadratic Spectral): \eqn{K(z) = \frac{25}{12\pi^2z^2} 
#' (\frac{\sin(6\pi z)/5)}{6\pi z/5} - \cos(6\pi z)/5)}). 
#' }
#' 
#' If the kernel type is not one of the six implemented, the function will terminate with an error message. 
#' The spatial two stage least square estimator is based on the matrix of instruments \eqn{H=[X,WX,W^2X^2]}.
#' 
#' @return   A list object of class \code{sphet}
#' \item{coefficients}{Spatial two stage least squares coefficient estimates}
#' \item{vcmat}{variance-covariance matrix of the estimated coefficients}
#' \item{s2}{S2sls residulas variance}
#' \item{residuals}{S2sls residuals}
#' \item{yhat}{difference between residuals and response variable}
#' \item{call}{the call used to create this object}
#' \item{model}{the model matrix of data}
#' \item{type}{the kernel employed in the estimation}
#' \item{bandwidth}{the type of bandwidth}
#' \item{method}{\code{'s2slshac'}}
#' @references
#' Kelejian, H.H. and Prucha, I.R. (2007) 
#'   HAC estimation in a spatial framework,
#'     \emph{Journal of Econometrics}, \bold{140}, pages 131--154.
#'     
#'     Kelejian, H.H. and Prucha, I.R. (1999) 
#'     A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model,
#'     \emph{International Economic Review}, \bold{40}, pages 509--533.
#'     
#'     Kelejian, H.H. and Prucha, I.R. (1998) 
#'     A Generalized Spatial Two Stage Least Square Procedure for Estimating a Spatial Autoregressive
#'     Model with Autoregressive Disturbances,
#'       \emph{Journal of Real Estate Finance and Economics}, \bold{17}, pages 99--121.
#'       
#' @author Gianfranco Piras \email{gpiras@mac.com}
#' @seealso \code{\link{gstslshet}}, \code{\link{distance}}, \code{\link{distance}}
#' @examples
#' library(spdep)
#' data(columbus)
#' listw <- nb2listw(col.gal.nb)
#' data(coldis)
#' res <- stslshac(CRIME ~ HOVAL + INC, data = columbus, listw = listw,  
#' distance = coldis, type = 'Triangular')
#' summary(res)
#' @keywords models
#' @export



stslshac<-function(formula, data=list(),listw,na.action=na.fail,zero.policy=NULL,HAC=TRUE, distance=NULL,type="Epanechnikov", bandwidth="variable",W2X=TRUE){

##functions that need to be sourced
	#source("twostagels.R")
	#source("kernelsfun.R")	
#extract model objects	
	
	if (is.null(zero.policy))
            zero.policy <- get.ZeroPolicyOption()
        stopifnot(is.logical(zero.policy))

	
	mt<-terms(formula,data=data)
	mf<-lm(formula, data, na.action=na.action,method="model.frame")
	na.act<-attr(mf,'na.action')


#preferences on missings values
if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
    }
    

#check that W ius an object of class listw
if(!inherits(listw,"listw")) 
	stop("The weights matrix is not a listw object")
	
##check that an exiting kernel is specified
if(HAC){
if(!(type %in% c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH","Rectangular"))) stop("Unknown kernel")
#check that the distance measure is specificed
if(is.null(distance) ) stop("No distance measure specified")

#check that dist is an object of class sphet distance
if(!inherits(distance,"distance")) 
	stop("The distance measure is not a distance object")

}
#generates x and y 
	y<-model.extract(mf,"response")
	x<-model.matrix(mt,mf)

cl<-match.call()
#checks on teh dimensions of x and y 	
if (length(y)!=nrow(x)) 
	stop("x and y have different length")

#check on the dimensions of x and W	
if (nrow(x) != length(listw$neighbours))
	stop("Input data and weights have different dimension")

#check that X and y does not have missing values	
if (any(is.na(y))) 
        stop("NAs in dependent variable")
if (any(is.na(x))) 
        stop("NAs in independent variable")

	
#fix the dimensions of the problem
	n<-nrow(x)
	k<-ncol(x)	
	xcolnames<-colnames(x)
	K<-ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)
	
	
	wy<-lag.listw(listw,y, zero.policy=zero.policy)
	wy<-array(wy,c(length(y),1))
	colnames(wy)<-("Wy")
	
	
if (any(is.na(wy)))
	stop("NAs in spatially lagged dependent variable")
	
if (k > 1) {
        WX <- matrix(nrow = n, ncol = (k  - (K - 1)))
        WWX <- matrix(nrow = n, ncol = (k  - (K - 1)))
        for (i in K:k) {
            wx <- lag.listw(listw, x[, i], zero.policy = zero.policy)
if(W2X)      wwx<- lag.listw(listw, wx, zero.policy = zero.policy)
            if (any(is.na(wx))) 
                stop("NAs in lagged independent variable")
            WX[, (i - (K - 1))] <- wx
				 WWX[, (i - (K - 1))] <- wwx
				         }
    }

#instr<-cbind(WX[,-c(1:8)],WWX[,-c(1:8)]) 
instr<-cbind(WX,WWX) 
#print(cbind(x,instr)[1,])

##spatial two stage least square of the initial model
#print(type)
results<-tsls(y=y,yend=wy, X=x, Zinst = instr, HAC=HAC, type=type, bandwidth=bandwidth, distance=distance, zero.policy=zero.policy)
model.data<-data.frame(cbind(y,x[,-1]))

results$call<-cl
results$model<-model.data
results$type<-type
results$bandwidth<-bandwidth
results$method<-"s2slshac"
results$HAC<-HAC
results$zero.policy<-zero.policy
class(results)<-c("sphet", "stsls")

return(results)
}
