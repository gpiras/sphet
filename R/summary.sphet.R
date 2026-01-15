#' @rdname spreg
#' @method coef sphet
#' @export
coef.sphet <- function(object, ...) object$coefficients


#' @rdname spreg
#' @method vcov sphet
#' @export
vcov.sphet <- function(object, ...) object$var

#' @rdname spreg
#' @method summary sphet
#' @export
summary.sphet <- function(object, width=getOption("width"),digits=getOption("digits"),obsinfo=FALSE,...){

	if(inherits(object,"distance")){
		weights<-	attributes(object)$GeoDa$dist
		neighbors<- object
		n<- length(weights)
		maxdist<-unlist(lapply(weights,max))
		nneigh<-unlist(lapply(weights,length))
		object$weights<-weights
		object$n<-n
		object$maxdist <- maxdist
		object$nneigh <- nneigh
		class(object)<- c('summary.sphet','sphet','distance')
		}
	else{
	coeff<-object$coefficients
	vcmat<-object$var
	se<-sqrt(diag(vcmat))
   t <- as.numeric(coeff/se)
  # print(t)
   pval<-pnorm(abs(t),lower.tail=FALSE) * 2
   CoefTable <- cbind(coeff,se,t,pval)
if(object$method=="legacy")  { 
	 if (object$HAC) colnames(CoefTable) <- c("Estimate","SHAC St.Er.","t-value","Pr(>|t|)")
	 else colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
	}
else colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
   object$CoefTable <- CoefTable
	class(object)<- c('summary.sphet','sphet')
	}
	object
}

#' @rdname spreg
#' @method print summary.sphet
#' @export
print.summary.sphet <- function(x, digits= max(3, getOption("digits") - 2), width=getOption("width"), obsinfo=FALSE,...){
if(inherits(x,"distance")){
	cat("\n Number of observations:\n")
	cat(" n:", x$n, "\n")
if (obsinfo){
	cat("\n Maximum distance for each observation:\n") 
	cat(x$maxdist,  fill=TRUE)
	}
	cat("\n Distance summary:\n")
	print(summary(x$maxdist), width=width)
if (obsinfo){	
	cat("\n Number of non-zero elements for each observation:\n")
	cat( x$nneigh, fill=TRUE)
	}
	cat("\n Neighbors summary:\n")	
	print(summary(x$nneigh), width=width)
	}
			
			else{
			  
				if(x$method == "legacy"){
				  if(x$HAC){
				  cat("======================================================\n")
				  cat("======================================================\n")
				  cat("                   Spatial Lag Model \n")
				  cat("                  HAC standard errors\n")
				  cat("======================================================\n")
				  cat("======================================================\n")
				  }
				  else{
				    cat("======================================================\n")
				    cat("======================================================\n")
				    cat("                   Spatial Lag Model \n")
				    cat("======================================================\n")
				    cat("======================================================\n")
				  }
				}
if (x$method=="s2slshac" ){ 
  
  if((isTRUE(x$Durbin)  | inherits(x$Durbin, "formula")) && x$HAC){
    cat("======================================================\n")
    cat("======================================================\n")
    cat("                   SLX model \n")
    cat("            with endogenous variables\n")
    cat("              HAC standard errors\n")
    cat("======================================================\n")
    cat("======================================================\n")
  }
  if(!isTRUE(x$Durbin) && !inherits(x$Durbin, "formula") && x$HAC){
    cat("======================================================\n")
    cat("======================================================\n")
    cat("                  Linear model \n")
    cat("             with endogenous variables\n")
    cat("               HAC standard errors\n")
    cat("======================================================\n")
    cat("======================================================\n")
  }
  if((isTRUE(x$Durbin)  | inherits(x$Durbin, "formula")) && !x$HAC){
    if(x$het){
      cat("======================================================\n")
      cat("======================================================\n")
      cat("                    SLX model \n")
      cat("            with endogenous variables\n")
      cat("              White standard errors\n")
      cat("======================================================\n")
      cat("======================================================\n")
    }
    else{
      cat("======================================================\n")
      cat("======================================================\n")
      cat("                   SLX model  \n")
      cat("            with endogenous variables\n")
      cat("======================================================\n")
      cat("======================================================\n")
    }
  }
  if(!isTRUE(x$Durbin) && !inherits(x$Durbin, "formula") && !x$HAC){
    if(x$het){
      cat("======================================================\n")
      cat("======================================================\n")
      cat("                   Linear model \n")
      cat("            with endogenous variables\n")
      cat("              White standard errors\n")
      cat("======================================================\n")
      cat("======================================================\n")
    }
    else{
      cat("======================================================\n")
      cat("======================================================\n")
      cat("                 Linear model  \n")
      cat("            with endogenous variables\n")
      cat("======================================================\n")
      cat("======================================================\n")
    }
    
  }
}
			  
			  if ( x$method=="olshac" ){ 
			      if((isTRUE(x$Durbin)  | inherits(x$Durbin, "formula")) && x$HAC){
			    cat("======================================================\n")
			    cat("======================================================\n")
			    cat("                   SLX model \n")
			    cat("              HAC standard errors\n")
			    cat("======================================================\n")
			    cat("======================================================\n")
			      }
			       if(!isTRUE(x$Durbin) && !inherits(x$Durbin, "formula") && x$HAC){
			         cat("======================================================\n")
			         cat("======================================================\n")
			         cat("                  Linear model \n")
			         cat("              HAC standard errors\n")
			         cat("======================================================\n")
			         cat("======================================================\n")
			       }
			    if((isTRUE(x$Durbin)  | inherits(x$Durbin, "formula")) && !x$HAC){
			      if(x$het){
			      cat("======================================================\n")
			      cat("======================================================\n")
			      cat("                    SLX model \n")
			      cat("              White standard errors\n")
			      cat("======================================================\n")
			      cat("======================================================\n")
			      }
			      else{
			        cat("======================================================\n")
			        cat("======================================================\n")
			        cat("                     SLX model  \n")
			        cat("======================================================\n")
			        cat("======================================================\n")
			      }
			    }
			    if(!isTRUE(x$Durbin) && !inherits(x$Durbin, "formula") && !x$HAC){
			      if(x$het){
			        cat("======================================================\n")
			        cat("======================================================\n")
			        cat("                   Linear model  \n")
			        cat("              White standard errors\n")
			        cat("======================================================\n")
			        cat("======================================================\n")
			      }
			      else{
			        cat("======================================================\n")
			        cat("======================================================\n")
			        cat("                 Linear model  \n")
			        cat("======================================================\n")
			        cat("======================================================\n")
			      }
			      
			    }
			    }
			  
if (x$method %in% c("gmm error", "gmm lag", "gmm sarar")){
  if(x$method =="gmm error"){
    cat("======================================================\n")
    cat("======================================================\n")
    cat("               GMM Spatial Error Model\n")
    cat("======================================================\n")
    cat("======================================================\n")
  }
  if(x$method =="gmm lag"){
    if(is.null(x$endog)){
    
    if((isTRUE(x$Durbin)  | inherits(x$Durbin, "formula")) && x$HAC){
      cat("======================================================\n")
      cat("======================================================\n")
      cat("              Spatial Durbin Model \n")
      cat("              HAC standard errors\n")
      cat("======================================================\n")
      cat("======================================================\n")
    }
    if(!isTRUE(x$Durbin) && !inherits(x$Durbin, "formula") && x$HAC){
      cat("======================================================\n")
      cat("======================================================\n")
      cat("                Spatial Lag Model \n")
      cat("               HAC standard errors\n")
      cat("======================================================\n")
      cat("======================================================\n")
    }
    if((isTRUE(x$Durbin)  | inherits(x$Durbin, "formula")) && !x$HAC){
      if(x$het){
        cat("======================================================\n")
        cat("======================================================\n")
        cat("              Spatial Durbin Model\n")
        cat("              White standard errors\n")
        cat("======================================================\n")
        cat("======================================================\n")
      }
      else{
        cat("======================================================\n")
        cat("======================================================\n")
        cat("               Spatial Durbin Model\n")
        cat("======================================================\n")
        cat("======================================================\n")
      }
    }
    if(!isTRUE(x$Durbin) && !inherits(x$Durbin, "formula") && !x$HAC){
      if(x$het){
        cat("======================================================\n")
        cat("======================================================\n")
        cat("                Spatial Lag Model \n")
        cat("              White standard errors\n")
        cat("======================================================\n")
        cat("======================================================\n")
      }
      else{
        cat("======================================================\n")
        cat("======================================================\n")
        cat("                 Spatial Lag model  \n")
        cat("======================================================\n")
        cat("======================================================\n")
      }
      
    }
    }
    
    else{
      if((isTRUE(x$Durbin)  | inherits(x$Durbin, "formula")) && x$HAC){
        cat("======================================================\n")
        cat("======================================================\n")
        cat("              Spatial Durbin Model \n")
        cat("            with endogenous variables\n")
        cat("              HAC standard errors\n")
        cat("======================================================\n")
        cat("======================================================\n")
      }
      if(!isTRUE(x$Durbin) && !inherits(x$Durbin, "formula") && x$HAC){
        cat("======================================================\n")
        cat("======================================================\n")
        cat("                 Spatial Lag Model \n")
        cat("             with endogenous variables\n")
        cat("               HAC standard errors\n")
        cat("======================================================\n")
        cat("======================================================\n")
      }
      if((isTRUE(x$Durbin)  | inherits(x$Durbin, "formula")) && !x$HAC){
        if(x$het){
          cat("======================================================\n")
          cat("======================================================\n")
          cat("               Spatial Durbin Model \n")
          cat("            with endogenous variables\n")
          cat("              White standard errors\n")
          cat("======================================================\n")
          cat("======================================================\n")
        }
        else{
          cat("======================================================\n")
          cat("======================================================\n")
          cat("               Spatial Durbin Model  \n")
          cat("            with endogenous variables\n")
          cat("======================================================\n")
          cat("======================================================\n")
        }
      }
      if(!isTRUE(x$Durbin) && !inherits(x$Durbin, "formula") && !x$HAC){
        if(x$het){
          cat("======================================================\n")
          cat("======================================================\n")
          cat("                Spatial Lag Model \n")
          cat("            with endogenous variables\n")
          cat("              White standard errors\n")
          cat("======================================================\n")
          cat("======================================================\n")
        }
        else{
          cat("======================================================\n")
          cat("======================================================\n")
          cat("                Spatial Lag Model  \n")
          cat("            with endogenous variables\n")
          cat("======================================================\n")
          cat("======================================================\n")
        }
        
      }  
    }
  }
  if(x$method =="gmm sarar"){
    cat("====================================================\n")
    cat("====================================================\n")
    cat("               GS2SLS SARAR Model\n")
    cat("====================================================\n")
    cat("====================================================\n")
  }
} 
			  
	cat("\nCall:\n")
  print(x$call)

  cat("\nResiduals:\n")
  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))
  print(sumres(x))


		  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable,digits=digits)

if(x$method=="gmm sarar") {
	
	if(!is.null(x$W)){
		  cat("\nWald test that rho and lambda are both zero:\n")
cat(" Statistics:", x$W$stat,"p-val:", x$W$pval, "\n")
}
}
}
  cat("\n")  
  invisible(x)
}

sumres <- function(x){
  sr <- summary(residuals(x))
  srm <- mean(residuals(x))
  if (abs(srm) < 1e-10){
    sr <- sr[c(1:3,5:6)]
  }
  sr
}

#' @rdname spreg
#' @method print sphet
#' @export
print.sphet <- function(x, digits = max(3, getOption("digits") - 3),...) 
{
	if(inherits(x,"distance")){ 
		tmp<-summary.sphet(x)
		print.summary.sphet(tmp)		
		}

else{	
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(x$coefficients[,1], digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")
}
    invisible(x)
}
