\name{calcPLSPCE}
\alias{calcPLSPCE}
\title{
Compute PLS-PCE Sensitivity Indexes
}
\description{
 Compute the optimal number of components, the PLS-PCE sensitivity indexes
  and related results.
}
\usage{
calcPLSPCE(pce, nc = 2)
}
\arguments{
  \item{pce}{
 an object of class \code{\linkS4class{PCEpoly}}.
    Design to analyze.
}
  \item{nc}{
integer. Required number of components.} 
}
\value{
  An object of class \code{\linkS4class{PLSPCE}}.}
\examples{
### Load the dataset
load(system.file("extdata",  "ishigami200.Rda", package="plspolychaos"))
X <- ishi200[, -ncol(ishi200)] # inputs
Y <- ishi200[,  ncol(ishi200)] # output
### Creation of the full polynomials
degree <- 6 # polynomial degree
pce <- polyLeg(X, Y, degree) 
### Compute the PLS-PCE sensitivity indexes for 25 components
ret <-  calcPLSPCE(pce, nc=25) 
}
\keyword{ regression }

