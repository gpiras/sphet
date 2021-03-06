% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listw2dgCMatrix.R
\name{listw2dgCMatrix}
\alias{listw2dgCMatrix}
\title{Interface between Matrix class objects and weights list}
\usage{
listw2dgCMatrix(listw, zero.policy = NULL)
}
\arguments{
\item{listw}{a \code{listw} object created for example by \code{nb2listw}}

\item{zero.policy}{See \code{lagsarlm} for details}
}
\value{
Matrix class object: a sparse Matrix
}
\description{
Interface between Matrix class objects and weights list
}
\examples{
library(spdep)
data(columbus)
listw <- nb2listw(col.gal.nb)
spW <- listw2dgCMatrix(listw)
}
\author{
Gianfranco Piras \email{gpiras@mac.com}
}
\keyword{spatial}
