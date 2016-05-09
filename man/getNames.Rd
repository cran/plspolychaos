\name{getNames}
\alias{getNames}
\alias{getNames,PCEpoly-method}
\title{
Display Structure of a Class
}
\description{
 Display the names, class and 
 length of all the slots of a
 \code{\linkS4class{PCEpoly}} 
 object.
 
}
\usage{
getNames(object)
}
\arguments{
  \item{object}{
object from class \code{\linkS4class{PCEpoly}}.
}
}
\details{
  It is a generic function. Its  method is defined in classe
  \code{\linkS4class{PCEpoly}}.
}
\value{
Nothing. It is a display function.
}

\seealso{
 Classe
  \code{\linkS4class{PCEpoly}}.
}
\examples{
### Load the datase
load(system.file("extdata",  "ishigami200.Rda", package="plspolychaos"))
X <- ishi200[, -ncol(ishi200)] # inputs
Y <- ishi200[,  ncol(ishi200)] # output
pce <- polyLeg(X, Y, degree=6)
### Display recursively the names, class and 
### length of all the components.
getNames(pce)
}

