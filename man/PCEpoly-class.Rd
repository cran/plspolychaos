\name{PCEpoly-class}
\Rdversion{1.1}
\docType{class}
\alias{PCEpoly-class}
%\alias{show,PCEpoly-method}
\alias{print.PCEpoly}
%\alias{getNames,PCEpoly-method}
%\alias{getNames}

\title{Class \code{"PCEpoly"}}
\description{
  Container of the PCE design: the
  computed values of the monomials
  and the polynomial structure description.
}
\section{Objects from the Class}{
Objects  from this class are created by calls to  functions \code{\link{polyLeg}}
or  \code{\link{analyticsPolyLeg}}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{matrix.
      The computed values of the monomials of the Legendre polynomial.
      The number of rows is the number of
      rows of the dataset.
      The number of  columns  is the number of monomials plus one.
      The first column is equal to one:
      it is for the constant term.
    }
    \item{\code{STRUC}:}{object of class \code{\linkS4class{PCEdesign}}. Matrix
      coding the polynomial structure. }
    \item{\code{nvx}:}{integer equal to the number of inputs. }
    \item{\code{call}:}{expression of class \sQuote{call}. The command
      which creates the object.}
  }
}
\section{Methods}{
    \describe{
    \item{getNames}{\code{signature(object = "PCEpoly")}: display the names, class and 
      length of all the components. }
    \item{print}{\code{signature(object = "PCEpoly", all=FALSE, ...)}:
      method of function \code{\link[base:print.default]{print}}.
      The polynomial expression is printed 
       when option \code{all} is set to TRUE.
      The additional arguments are passed to the
       \code{\link[base:print.default]{print.default}} function. }
    \item{show}{\code{signature(object = "PCEpoly")}: same as
      function \code{print}, without any arguments. }
 }
   
}

\seealso{
\itemize{
\item
  Functions \code{\link{polyLeg}} and \code{\link{analyticsPolyLeg}},
  creators of objects from this class
\item the generic method \code{\link{getNames}}
}
}
\examples{
showClass("PCEpoly")
}
\keyword{classes}
