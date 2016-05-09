\name{PCEdesign-class}
\Rdversion{1.1}
\docType{class}
\alias{PCEdesign-class}
%\alias{show,PCEdesign-method}
\alias{print.PCEdesign}

\title{Class \code{"PCEdesign"}}
\description{
  Container of the polynomial description structure.
}
\section{Objects from the Class}{
  Objects from this class are created by calls to  functions
  \code{\link{polyLeg}} or
   \code{\link{analyticsPolyLeg}}.
  They are stored  in the slot \code{design}
  in the object of  class \code{\linkS4class{PCEpoly}}
  returned by these functions.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{matrix
      with as many columns as inputs and as many rows as
      monomials plus one.
      Element \code{(i,j)} is an integer
      equal to the degree of the input \code{j}
      in the monomial \code{i-1}.      
      The first row is equal to zero: it is for the constant term.
    }
    \item{\code{degree}:}{integer equal to the polynomial degree.
    }
    \item{\code{total.nmono}:}{integer equal to the number of monomials
      in the full polynomial.}
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "PCEdesign", all=FALSE, ...)}:
      method of function \code{\link[base:print.default]{print}}. If option \code{all} is set to TRUE,
      all the monomials are printed. The additional arguments are passed to the
      \code{\link[base:print.default]{print.default}} function. }
    \item{show}{\code{signature(object = "PCEdesign")}: same as
      function \code{print}, without any arguments. }
	 }
       }

       \seealso{
	 \itemize{
	   \item
 Functions \code{\link{polyLeg}}
 and   \code{\link{analyticsPolyLeg}}, creators of objects from this
 class.
 \item Class \code{\linkS4class{PCEpoly}} in which 
 objects from this class are stored.
 }
}

\keyword{classes}
