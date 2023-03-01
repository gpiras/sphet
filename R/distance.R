#' Writes distance matrices
#' @aliases distance
#' @title Distance measures available in distance
#' @description Reads points coordinates and generates objects of class \code{distance.matrix}
#'
#' @usage distance(coord, region.id = NULL, output = TRUE,
#' type = c("NN", "distance", "inverse"),
#' measure = c("euclidean", "gcircle", "chebyshev", "braycur", "canberra"),
#' nn = 6, cutoff = FALSE, miles = TRUE, R = NULL, shape.name = NULL, region.id.name = NULL,
#' firstline = FALSE, file.name = NULL)
#' 
#' @details The object created is similar to the content of a \code{'GWT'} file. The output file can be of any format. 
#' In particular, it could be a \code{'GWT'} file. When \code{firstline} is \code{TRUE}, an header line is added to the \code{'GWT'} file. 
#' The first element is simply a place holder, the second is the number of observations. The name of the shape file 
#' and of the id variable can be specified by the options \code{shape.name} and \code{region.id.name} respectively. 
#' The function performs a series of test on the \code{region.id} variable. If a \code{region.id} variable is not specified and 
#' \code{coord} only has two columns, a sequence from 1 to the number of observations is generated and 
#' used as identification variable. If \code{region.id} is specified 
#' and the first column of \code{coord} contains an id variable they should be the same. 
#' 
#' The distance measures implemented in \code{sphet} are:
#' \itemize{
#'  \item \code{'euclidean'}: \eqn{\sqrt{\sum{(x_i - y_i)^2}}}
#'  \item \code{'chebyshev'}: \eqn{max(|x_i - y_i|)} 
#'  \item \code{'braycur'}: \eqn{ \frac{\sum{|x_i - y_i|}}{\sum{|x_i + y_i|}}}
#'  \item \code{'canberra'}:  \eqn{  \frac{\sum{|x_i - y_i|}}{\sum{|x_i| + |y_i|}} }
#'  \item \code{'gcircle'}: see \code{\link[sp]{spDists}}, which uses an approximation to the WGS84 spheroid. 
#' }
#' 
#' @param coord a matrix with the (X,Y)-coordinates of the points. The first column can be the 
#' region.id variable giving the ordering of the observations
#' @param region.id variable that defines the ordering of the observations
#' @param output when TRUE (default) writes the object to a file 
#' @param type one of \code{("NN","distance","inverse")}. Nearest neighbors, distance or inverse distance 
#' @param measure one of \code{("euclidean","gcircle","chebyshev","braycur","canberra")}.The distance measure to be employed in the calculations
#' (See Details)
#' @param nn the number of nearest neighbors 
#' @param cutoff If type is \code{distance} or \code{inverse}. Assumes values 1, 2 or 3. 
#' When 1, the cutoff is set to the first quantile of the distribution of distances. When 2 to the median, and when 3 to 
#' the third quantile. Only observations with distance less than cutoff distance are neighbors.  
#' @param miles If TRUE (default), distances are in miles, otherwise in Km. (See \code{\link[sp]{spDists}} wcich returns km, and are converted if required) 
#' @param R  deprecated, \code{\link[sp]{spDists}} uses an approximation to the WGS84 spheroid 
#' @param shape.name The name of the shape file. See Details 
#' @param region.id.name The name of the \code{region.id} variable. See Details 
#' @param firstline If \code{TRUE}, a first line is added to the output file. See Details
#' @param file.name If \code{output}, the name of the output file. See Details 
#'
#'  
#' @return A \code{matrix} of three columns: \code{from}, \code{to}, and \code{distance}
#' 
#' @examples  
#' set.seed("1234")
#' X <- runif(100, 0, 70)
#' Y <- runif(100, -30, 20)
#' coord1 <- cbind(seq(1,100), X, Y)
#' thm2 <- distance(coord1, region.id = NULL, 
#' output = FALSE, type = "NN", nn = 6)
#' thm2 <- distance(coord1, region.id = NULL, output = FALSE, type = "distance", cutoff = 1)
#'
#' @export
#' @author  Gianfranco Piras \email{gpiras@mac.com}
#' 
#' @keywords spatial




distance <- function(coord,region.id=NULL,
                     output=TRUE, type=c("NN","distance","inverse"),
                     measure=c("euclidean","gcircle","chebyshev","braycur","canberra"),
                     nn=6, cutoff=FALSE, miles=TRUE,R=NULL, 
                     shape.name=NULL,region.id.name=NULL,
                     firstline=FALSE,file.name=NULL) {
	
	if ( !type %in% c("NN","distance","inverse") ) stop("unknown type")
	if ( !measure[1] %in% c("euclidean","gcircle","chebyshev","braycur","canberra") ) stop("unknown measure")
	if(!inherits(coord,c("data.frame","matrix"))){
		coord <-  as.matrix(coord)
		}
	if (is.null(region.id) & ncol(coord)!=3) {
		id<-seq(1,nrow(coord))
		warning("region.id variable not specified")
		}
	if (is.null(region.id) & ncol(coord)==3) id<-coord[,1]
	if (!is.null(region.id) & ncol(coord) == 3){ 
		check<- all.equal(region.id,coord[,1])
		if ( check!=TRUE ) stop ("region.id and coord[,1] are different")
		id <- region.id
		}
	if (!is.null(region.id) && ncol(coord) != 3) id<-region.id 
	
	if (length(unique(id)) != length(id)) 
            stop("non-unique region.id given")

		
    k<- ifelse(ncol(coord)==3,2,1)
	 x<- coord[,k]
	 y<- coord[,k+1]
    h <- length(x)
    hh2 <- h*(h-1)/2
    hh <- h*(h-1)
    dd <- matrix(,nrow=hh2,ncol=3) 
	 vec<-id
	 vec1<-x
	 vec2<-y
		coord<-cbind(vec,vec1,vec2)
Weights<-cbind(0,0,0)
		Sq2<-c(1,seq(h-1,2))
		#print(Sq2)
		SSq2<-cumsum(Sq2)
		#print(SSq2)
						for (i in vec[-length(vec)]) {
									ref<-	vec[i]
									w<-which(vec!=i & vec > i)

									coordi<-cbind(vec1[i],vec2[i])
									coordj<-cbind(vec1[w],vec2[w])
####

		dist.fun<-switch(match.arg(measure), euclidean={
			dist.euclidean
			}, gcircle={
				dist.gcircle
				}, chebyshev= {
					dist.chebyshev
					}, braycur={
						dist.braycur
						}, canberra = {
							dist.canberra
								})

weights<-dist.fun(coordi,coordj,miles=TRUE,R=NULL)
if(measure[1]=="gcircle") weights<-t(weights)            
									
#####									
									
if(type=="inverse") tmp <- rbind(cbind(rep(vec[i],length(w)),vec[w],1/weights),cbind(vec[w],rep(vec[i],length(w)),1/weights))
else tmp <- rbind(cbind(rep(vec[i],length(w)),vec[w],weights),cbind(vec[w],rep(vec[i],length(w)),weights))

#print(tmp)
						if (!is.numeric(weights)) stop
									Weights<-rbind(Weights,tmp)
			}
			Weights<-Weights[-1,]

	if (type=="NN"){
		or<-order(Weights[,1],Weights[,3])
Weights<- Weights[or,]
final<-matrix(,h*nn,3)
for (i in 1:h){
	tmp<-Weights[Weights[,1]==i,][1:nn,]
	final[(((i*nn)-(nn-1)):(i*nn)),]<-as.matrix(tmp)
	}
Weights<-final
		}	
		
if (type=="distance" && cutoff){
			tmp<-Weights[,3]
		quant<-quantile(tmp)[cutoff+1]
tokeep<-which(tmp <= quant)
#print(tokeep)
Weights<-Weights[tokeep,]
}

if (type=="inverse" && cutoff){
			tmp<-Weights[,3]
		quant<-quantile(tmp)[cutoff+1]
tokeep<-which(tmp <= quant)
Weights<-Weights[tokeep,]
}

	    worder <- order(Weights[,1],Weights[,2])
   		 Weights<-Weights[worder,]
    colnames(Weights)<-c("from","to","distance")
    class(Weights)<- c("matrix","distance.matrix")
   
    if (output) {
   		 firstl<-c("0",h,shape.name,region.id.name)
   		 #print(firstline)
if(firstline) write(firstl,file.name,ncolumns=4)
    	write.table(Weights,file.name,row.names=FALSE,quote=FALSE,append=TRUE,col.names=FALSE)
    	}
Weights
    }
