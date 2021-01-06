#' GM estimation of a Cliff-Ord type model with Heteroskedastic Innovations 
#' @aliases gstslshet
#' @name gstslshet
#' 
#' @description Multi step GM/IV estimation of a linear Cliff and Ord -type of model
#' of the form:
#'    
#'    \deqn{y=\lambda W y + X \beta + u  }
#'  \deqn{u=\rho W u + e}  with \deqn{e ~ N(0,\sigma^2_i) }
#'  
#'  The model allows for spatial lag in the dependent variable
#'  and disturbances. The innovations in the disturbance process are assumed 
#'  heteroskedastic of an unknown form.
#'  
#' @usage  gstslshet(formula, data = list(), listw, na.action = na.fail, 
#'            zero.policy = NULL, initial.value = 0.2, abs.tol = 1e-20, 
#'            rel.tol = 1e-10, eps = 1e-5, inverse = T, sarar = T)
#'            
#' @param formula a description of the model to be fit              
#' @param data an object of class \link{data.frame}. An optional data frame containing the variables
#'  in the model.
#' @param listw an object of class \code{listw} created for example by \code{nb2listw}
#' @param na.action a function which indicates what should happen when the data contains missing values.
#'    See \link{lm} for details.
#' @param zero.policy See \code{lagsarlm} for details
#' @param initial.value The initial value for \eqn{\rho}. It can be either numeric (default is 0.2) or
#'    set to \code{'SAR'}, in which case the optimization will start from the estimated coefficient of a regression of the 2SLS 
#'    residuals over their spatial lag (i.e. a spatial AR model)
#' @param abs.tol Absolute tolerance. See \link{nlminb} for details.
#' @param rel.tol Relative tolerance. See \link{nlminb} for details.
#' @param eps Tolerance level for the approximation. See Details.
#' @param  inverse \code{TRUE}. If \code{FALSE}, an appoximated inverse is calculated. See Details.
#' @param sarar \code{TRUE}. If \code{FALSE}, a spatial error model is estimated. 
#' 
#' @details   The procedure consists of two steps alternating GM and IV estimators. Each step consists of sub-steps.
#' In step one \eqn{\delta = [\beta',\lambda]'} is estimated by 2SLS. The 2SLS residuals are first employed
#' to obtain an initial (consistent but not efficient) GM estimator of \eqn{\rho} and then a consistent and efficient 
#' estimator (involving the variance-covariance matrix of the limiting distribution of the normalized sample moments). 
#' In step two, the spatial Cochrane-Orcutt transformed model is estimated by 2SLS. This corresponds to a GS2SLS procedure. 
#' The GS2SLS residuals are used to obtain a consistent and efficient GM estimator for \eqn{\rho}. 
#'
#' The initial value for the optimization in step 1b is taken to be \code{initial.value}. The initial value in step 1c is the 
#' optimal parameter of step 1b. Finally, the initial value for the optimization of step 2b is the optimal parameter of step 1c.
#'
#' Internally, the object of class \code{listw} is transformed into a \link{Matrix} 
#' using the function \link{listw2dgCMatrix}.
#'
#'
#' The expression of the estimated variance covariance matrix of the limiting 
#' distribution of the normalized sample moments based on 2SLS residuals 
#' involves the inversion of \eqn{I-\rho W'}.
#' When \code{inverse} is \code{FALSE}, the inverse is calculated using the approximation 
#' \eqn{I +\rho W' + \rho^2 W'^2 + ...+ \rho^n W'^n}. 
#' The powers considered depend on a condition. 
#' The  function will keep adding terms until the absolute value of the \code{sum} of all elements 
#' of the matrix \eqn{\rho^i W^i} is greater than a fixed \eqn{\epsilon} (\code{eps}). By default \code{eps}
#' is set to 1e-5.
#' 
#' 
#' @return   A list object of class \code{sphet}
#' \item{coefficients}{Generalized Spatial two stage least squares coefficient estimates of \eqn{\delta} and GM estimator for \eqn{\rho}. }
#' \item{var}{variance-covariance matrix of the estimated coefficients}
#' \item{s2}{GS2SLS residuals variance}
#' \item{residuals}{GS2SLS residuals}
#' \item{yhat}{difference between GS2SLS residuals and response variable}
#' \item{call}{the call used to create this object}
#' \item{model}{the model matrix of data}
#' \item{method}{\code{'gs2slshac'}}
#' \item{W}{Wald test for both \eqn{\rho} and \eqn{\lambda} are zero}
#' 
#' @references  
#'   Arraiz, I. and Drukker, M.D. and Kelejian, H.H. and Prucha, I.R. (2007) 
#' A spatial Cliff-Ord-type Model with Heteroskedastic Innovations: Small and Large Sample Results,
#' \emph{Department of Economics, University of Maryland}'
#'
#'      Kelejian, H.H. and Prucha, I.R. (2007) 
#'Specification and Estimation of Spatial Autoregressive Models with Autoregressive and Heteroskedastic Disturbances,
#'    \emph{Journal of Econometrics}, forthcoming.
#'
#'  Kelejian, H.H. and Prucha, I.R. (1999) 
#' A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model,
#'    \emph{International Economic Review}, \bold{40}, pages 509--533.
#'    
#'      Kelejian, H.H. and Prucha, I.R. (1998) 
#'A Generalized Spatial Two Stage Least Square Procedure for Estimating a Spatial Autoregressive
#' Model with Autoregressive Disturbances,
#'    \emph{Journal of Real Estate Finance and Economics}, \bold{17}, pages 99--121.
#'    
#' @seealso \code{\link{stslshac}} 
#'    
#' @author  Gianfranco Piras \email{gpiras@mac.com}
#'    
#' @examples
#' data(columbus, package = "spdep")
#' listw <- spdep::nb2listw(col.gal.nb)
#' res <- gstslshet(CRIME ~ HOVAL + INC, data = columbus, listw = listw)
#' summary(res)
#' 
#' @keywords models
#' @export




gstslshet<-function(formula, data=list(),listw, na.action=na.fail, zero.policy=NULL, initial.value=0.2, abs.tol=1e-20, rel.tol=1e-10, eps=1e-5, inverse=T, sarar=T){

##functions that need to be sourced
	#source("twostagels.R")
	#source("utilities.R")
	#source("listw2dgCMatrix.R")
	#source("Omega.R")
	
	
	if (is.null(zero.policy))
            zero.policy <- get.ZeroPolicyOption()
        stopifnot(is.logical(zero.policy))
	
#extract model objects	
	
	mt<-terms(formula,data=data)
	mf<-lm(formula, data, na.action=na.action, method="model.frame")
	na.act<-attr(mf,'na.action')

# record call
cl<- match.call()


#if(!inverse) warning("Approximated inverse could be inaccurate")

#preferences on missings values
if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
    }

#check that W ius an object of class listw
if(!inherits(listw,"listw")) 
	stop("The weights matrix is not a listw object")

#generates x and y 
	y<-model.extract(mf,"response")
	x<-model.matrix(mt,mf)

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

if(sarar){	
	
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
            wwx<- lag.listw(listw, wx, zero.policy = zero.policy)
            if (any(is.na(wx))) 
                stop("NAs in lagged independent variable")
            WX[, (i - (K - 1))] <- wx
				 WWX[, (i - (K - 1))] <- wwx
				         }
    }

instr<-cbind(WX,WWX) 




##spatial two stage least square of the initial model
firststep<-tsls(y=y,yend=wy, X=x, Zinst = instr)

ubase<-residuals(firststep)


#GM step 1b
		int1<-Ggfastfast(listw=listw,ubase,n, zero.policy = zero.policy)

##One could start from this initial value but if uses optimize, the initial values are not an issue
if (initial.value=="SAR"){
		Wubase<-lag.listw(listw,ubase, zero.policy = zero.policy)
		pars<-coefficients(lm(ubase~Wubase-1))
		}
else pars<-initial.value
		
		
optres1<-nlminb(pars, arg, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1)

		param<-optres1$par
		# print(param)
		reg<-x
		u<-ubase
		toinst<-wy
		
		fi<-fifour(x,listw,instr,ubase,wy,param,n,inverse,eps, zero.policy= zero.policy)
		# print(fi$nlw)

GMMfeas3<-nlminb(param, arg1, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1, VC=fi$nlw)


rhotilde <- GMMfeas3$par
# print(rhotilde)

fi1<-fifour(x,listw,instr,ubase,wy,rhotilde,n,inverse=inverse,eps=eps, zero.policy = zero.policy)
Om<-Omega(n,gamma=fi1$res1,H=fi1$instr,param=rhotilde,G=int1$bigG,FIinv=fi1$nlw,a=fi1$a,FI=fi1$nl,P=fi1$P,gammas=fi1$gammas,Ws=fi1$Ws)

yt  <- y - rhotilde * wy
xt <- x - rhotilde * lag.listw(listw,x, zero.policy = zero.policy)
wyt<-wy - rhotilde * lag.listw(listw,wy, zero.policy = zero.policy)
colnames(xt)<-xcolnames
colnames(wyt)<-c('Wyt')

secstepb<-tsls(y=yt, yend=wyt, X=xt, Zinst = instr, reg=x, end=wy, yor=y,modified=TRUE)

utildeb<-secstepb$residuals
int1b<-Ggfastfast(listw, utildeb,n, zero.policy = zero.policy)
psippb<-fistslsfast(reg=xt, Ws=fi$Ws, instr=instr, resid=utildeb, toinst=wyt, param=rhotilde,solo=x,end=wy)
pars1<-rhotilde


GMMstsls1b<-nlminb(pars1, arg1, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1b, VC=psippb$nlw)

rhohat<-GMMstsls1b$par
#xt <- x - rhohat * lag.listw(listw,x)
psifinal<-fistslsfast(reg=xt, Ws=fi$Ws, instr=instr, resid=utildeb, toinst=wyt, param=rhohat,solo=x,end=wy)
Omfinal<- Omegabis(gammas=psifinal$res1, Hs=psifinal$instr, param=rhohat, G=int1b$bigG, FIinv=psifinal$nlw, as=psifinal$a, n, FI=psifinal$nl, P=psifinal$P) 
coef<-(secstepb$coefficients[2:(ncol(xt)+1)])
coeff<-as.matrix(c(coef, secstepb$coefficients[1], rhohat ))
rownames(coeff)<-c(xcolnames, 'lambda', 'rho')
s2<-crossprod(utildeb)/(n-k)


model.data<-data.frame(cbind(y,x[,-1]))

method<-"gs2slshac"

	k<-nrow(coeff)
	R<-matrix(0,1,k)
	R[,((k-1):k)]<-1
	Rbeta<-R%*%coeff
	Rvar<-R%*%Omfinal$VarCov%*%t(R)
	stat<-as.numeric(t(Rbeta)%*% Rbeta/Rvar)
	pval <- pchisq(stat,df=1,lower.tail=FALSE)
W<-list(stat=stat,pval=pval)



results<-list(coefficients=coeff,var=Omfinal$VarCov, s2=s2, call=cl, residuals=as.numeric(utildeb), model=model.data,method=method,W=W,yhat=secstepb$yhat )
}
else{
	##this is an error model estimated by ols
firststep<-lm(y~x-1)
ubase<-residuals(firststep)

int1<-Ggfastfast(listw=listw,ubase,n, zero.policy = zero.policy)

##One could start from this initial value but if uses optimize, the initial values are not an issue
if (initial.value=="SAR"){
		Wubase<-lag.listw(listw,ubase, zero.policy = zero.policy)
		pars<-coefficients(lm(ubase~Wubase-1))
		}
else pars<-initial.value
		
		
optres1<-nlminb(pars, arg, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1)

		param<-optres1$par
		reg<-x
		u<-ubase
				
		fi<-fierror(x,listw, ubase, param,n,inverse=inverse,eps=eps)

GMMfeas3<-nlminb(param, arg1, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1, VC=fi$nlw)

rhotilde<- GMMfeas3$par
# print(rhotilde)

fi1<-fierror(x,listw,ubase,rhotilde,n,inverse=inverse,eps=eps)
Om<-Omega(n,gamma=fi1$res1, H=x, param=rhotilde,G=int1$bigG,FIinv=fi1$nlw,a=fi1$a,FI=fi1$nl,P=fi1$P,gammas=fi1$gammas,Ws=fi1$Ws) ## in an error model H=X
#####################


yt  <- y - rhotilde * lag.listw(listw,y, zero.policy=zero.policy)
xt <- x - rhotilde * lag.listw(listw,x, zero.policy=zero.policy)
colnames(xt)<-xcolnames

secstepb<-lm(yt~xt-1)

utildeb<-y - x %*% as.matrix(coefficients(secstepb)) 

int1b<-Ggfastfast(listw, utildeb,n, zero.policy = zero.policy)

psippb<-fistslserror(reg=x, Ws=fi$Ws, resid=utildeb, param=rhotilde, solo = cbind(x,lag.listw(listw,x[,-1], zero.policy=zero.policy)))
pars1<-rhotilde


GMMstsls1b<-nlminb(pars1, arg1, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=int1b, VC=psippb$nlw)

rhohat<-GMMstsls1b$par
psifinal<-fistslserror(reg=x, Ws=fi$Ws, resid=utildeb, param=rhohat,solo=cbind(x,lag.listw(listw,x[,-1], zero.policy=zero.policy)))
Omfinal<- Omegabis(gammas=psifinal$res1, Hs = cbind(x,lag.listw(listw,x[,-1], zero.policy=zero.policy)) , param=rhohat, G=int1b$bigG, FIinv=psifinal$nlw, as=psifinal$a, n, FI=psifinal$nl, P=psifinal$P) 
coef<-coefficients(secstepb)
coeff<-as.matrix(c(coef, rhohat ))
rownames(coeff)<-c(xcolnames, 'rho')
s2<-crossprod(utildeb)/(n-k)


model.data<-data.frame(cbind(y,x[,-1]))

method<-"gs2slshac"


results<-list(coefficients=coeff,var=Omfinal$VarCov, s2=s2, call=cl, residuals=as.numeric(utildeb), model=model.data,method=method,W=NULL,yhat=fitted(secstepb) )

	
	}

class(results)<-c("sphet", "gstsls")

return(results)
}

impacts.gstsls <- function(obj, ..., tr=NULL, R=NULL, listw=NULL,
    evalues=NULL, tol=1e-6, empirical=FALSE, Q=NULL) {
    if (!is.null(obj$listw_style) && obj$listw_style != "W") 
        stop("Only row-standardised weights supported")
    coefs <- drop(obj$coefficients)
    p2 <- length(coefs)
    rho <- coefs[(p2-1)]
    beta <- coefs[1:(p2-2)]
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
    p <- length(beta)
    n <- length(obj$residuals)
    mu <- c(beta, rho, lambda)
    Sigma <- obj$var
    irho <- p2-1
    drop2beta <- c(p2-1, p2)
    if (!requireNamespace("spatialreg", quietly=TRUE))
        stop("install spatialreg")
    res <- spatialreg::intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu,
        Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
        interval=NULL, type="lag", tr=tr, R=R, listw=listw, evalues=evalues,
        tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "iClass") <- class(obj)
    res
}


