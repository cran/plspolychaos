\name{polyLeg}
\alias{polyLeg}
\title{
Calculate Legendre Polynomials on a Dataset
}
\description{
  This function calculates Legendre polynomials, optionally reducted
  to the most significant monomials, on
  a user dataset.
  
  Legendre polynomials are computed
  after calibration within the bounds [-1, +1].
}
\usage{
polyLeg(lhs, Y, degree, forward=NULL)
}
\arguments{
  \item{lhs}{ matrix with as many columns
    as inputs. Dataset of  inputs. Generally, a space filling design is used for
    forming
    this dataset. Typically, this is a simple LHS (see  McKay, 1979)
    or a modified  LHS.
}
  \item{Y}{
    vector of length equal to the number of rows in \code{lhs}.
    Model outputs.
}
\item{degree}{
  integer greater than 1 and less than 11. Degree of the polynomial. 
}
  \item{forward}{
    NULL or an integer equal to the required number of monomials.
    A null value (the default), or a value less than the number of
    inputs or greater than
    the total number of monomials, means that all the monomials are kept. See details.}
}
\value{
  An objet of class \code{\linkS4class{PCEpoly}}.
}
\details{
    When the value of the argument \code{forward} is non NULL,
    it should be an integer equal to the required
    number of the monomials (let say \code{q}). The \code{q} monomials are selected,
    among all the monomials of the full polynomial, by all the
    linear simple regressions of the output versus all the monomials.
    Those associated with the \code{q} largest \eqn{R^2}{R^2} values
    are kept.
}
\references{
 McKay, M.D. and Beckman, R.J. and Conover, W.J. 
    1979. \dQuote{A Comparison of Three Methods for Selecting Values of
      Input Variables in the Analysis of Output from a Computer Code}.In
    \emph{Technometrics},
 21 (2). 239-245p. 
  }    
\seealso{
 \itemize{ 
    \item  Function \code{\link{analyticsPolyLeg}}  builds
  Legendre polynomials from a simulated dataset.
  \item
  Function \code{\link{calcPLSPCE}} calculates PLS-PCE sensivity indexes
from the returned
  object.
  }
}
\examples{
### Load the dataset
load(system.file("extdata", "ishigami200.Rda", package="plspolychaos"))
X <- ishi200[, -ncol(ishi200)] # inputs
Y <- ishi200[,  ncol(ishi200)] # output
degree <- 6 # polynomial degree
### Creation of the full polynomials
pce <- polyLeg(X, Y, degree) 
print(pce)
### Selection of the 50 most significant monomials
pcef <- polyLeg(X, Y, degree, forward=50)
print(pcef)
}


\keyword{ regression }
