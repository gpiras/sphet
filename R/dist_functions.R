dist.euclidean<-function(X,Y,...){
	X<-matrix(X,nrow(Y),2,byrow=TRUE)
	dif<-rowSums((X-Y)^2)
	return<-sqrt(dif)
	}
	

dist.chebyshev<-function(X,Y,...){
	X<-matrix(X,nrow(Y),2,byrow=TRUE)
	dif<-matrix(,nrow(Y),1)
for(i in 1:nrow(Y))dif[i,]<-max(abs(X-Y)[i,])
	return<- dif
	}

dist.braycur<-function(X,Y,...){
	X<-matrix(X,nrow(Y),2,byrow=TRUE)
	dif<-rowSums(abs((X-Y)))
	su<-rowSums(abs(X+Y))
	return<- dif /su
	}

dist.canberra<-function(X,Y,...){
	X<-matrix(X,nrow(Y),2,byrow=TRUE)
	dif<-abs((X-Y))
	su<-abs(X)+abs(Y) 
	return<- rowSums(dif /su)
	}
