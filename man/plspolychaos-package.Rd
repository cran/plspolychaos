\name{plspolychaos-package}
\alias{plspolychaos-package}
\alias{plspolychaos}
\docType{package}
\title{
\packageTitle{plspolychaos}
}
\description{
\packageDescription{plspolychaos}
}
\details{

  The Legendre chaos polynomials are calculated, either on a
  user provided dataset  by function  \code{\link{polyLeg}},
  or on a simulated LHS by function \code{\link{analyticsPolyLeg}}.
  Then, function \code{\link{calcPLSPCE}} calculates
  PLS-regression coefficients,
  PLS-PCE sensitivity indexes and some other
  results.
}
\author{
\packageAuthor{plspolychaos}

Maintainer: \packageMaintainer{plspolychaos}
}
\references{
  \itemize{ 
    \item Metamodeling and global sensitivity analysis for computer models with
correlated inputs: a practical approach tested with a 3D light interception
computer model.
J.-P. Gauchi et al.
In \emph{Environmental Modelling & Software},
2016, submitted.
  }
  }
\seealso{
\code{\link[polychaosbasics:polychaosbasics-package]{polychaosbasics}}
  package.
  }
  \examples{
### First example: the dataset is simulated
nlhs <- 200 # number of rows
degree <- 6 # polynomial degree
set.seed(42)# fix the seed for reproductible results
# Generate data and calculate Legendre polynomials
# Independent inputs; response calculated by the Ishigami function
pce <- analyticsPolyLeg(nlhs, degree, 'ishigami')
# Compute the PLS-PCE sensitivity indexes for ten components
ret <- calcPLSPCE(pce, nc=10) 
# Plot the result
\dontrun{plot(ret, pce)}

### Second example: the dataset is provided and the
### most significant monomials are selected
# Load the dataset
load(system.file("extdata",  "ishigami200.Rda", package="plspolychaos"))
X <- ishi200[, -ncol(ishi200)] # inputs
Y <- ishi200[,  ncol(ishi200)] # output
# Build Legendre polynomial with the 50 most significant monomials
pce <- polyLeg(X, Y, degree=6, forward=50) 
# Compute the PLS-PCE sensitivity indexes 
ret <-  calcPLSPCE(pce, nc=10) 
print(ret, all=TRUE)
}
\keyword{ package }
\keyword{ regression}
