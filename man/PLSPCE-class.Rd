\name{PLSPCE-class}
\Rdversion{1.1}
\docType{class}
\alias{PLSPCE-class}
%\alias{show,PLSPCE-method}
\alias{print.PLSPCE}
\alias{plot.PLSPCE}

\title{Class \code{"PLSPCE"}}
\description{
  An S4 class  container for the result of the function
  \code{\link{calcPLSPCE}}. }
\section{Slots}{
  \describe{
    \item{\code{indexes}:}{ matrix with as many rows as inputs and three
      columns. Values of the PLS-PCE sensitivity indexes.
      The column labels are LE, PE, TPE.
      \itemize{
      \item
      \code{indexes[i, "LE"]} is the Linear Effect of the  input \code{i}.
      \item
      \code{indexes[i, "PE"]} is the Polynomial Effect. It is the  effect of the
      monomials in which only the input \code{i} appears.
      \item
      \code{indexes[i, "TPE"]} is the Total  Polynomial Effect. It is the  effect of all the
      monomials in which the input \code{i} appears.
    }
    }
    \item{\code{indexes.percent}:}{matrix.
      Percentages of the PLS-PCE sensitivity indexes, i.e
       values of  \code{indexes} expressed as percentages
      of the sums of their  columns.
    }
    \item{\code{ncopt}:}{number of the optimal component.}
    \item{\code{R2}:}{matrix
      with as many rows as components and three columns.
      The columns are labeled
      R2, \%R2 and \%R2cum. They store respectively the R2 values,
       percentages
       of their column sums, and cumulated values of their
       percentages.}
    \item{\code{Q2}:}{matrix
      with as many rows as components and two columns.
      The columns are labeled
      Q2 and Q2cum. They store respectively the Q2 and
      the Q2cum values. }      
    \item{\code{rmsep}:}{matrix
      with as many rows as components and one column.
      Root mean square error predictions.}
    \item{\code{COEF}:}{matrix
      with as many rows as monomials plus one and as many columns
      as components.
      PLS-regression coefficients. The first row is the constant term.
  }
  \item{\code{betaCR}:}{matrix
      with as many rows as monomials and as many columns
      as components. The centered-reducted PLS-regression coefficients,
      for all the components.
    }
    \item{\code{y.hat}:}{vector of length  equal to
    the number of rows of the dataset. Metamodel outputs.
  }
  \item{\code{STRUC}:}{object of class \code{\linkS4class{PCEdesign}}  coding the polynomial structure. }
}}
\note{
  The \code{COEF} and \code{STRUC} slots can be used for making
  predictions
  for Legendre-coded new input values.
}

\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "PLSPCE", all=FALSE, ...)}:
      method of function \code{\link[base:print.default]{print}}.
      When option \code{all} is FALSE, some
      slots are hidden :
      \code{RMSEP, coef, betaCR, y.hat, STRUC}.
      When option \code{all} is set to TRUE,
      the names of  the hidden slots are printed. The additional arguments are passed to the
      \code{\link[base:print.default]{print.default}} function. }
    \item{show}{\code{signature(object = "PLSPCE")}: same as
      function \code{print}, without any arguments. }
       \item{plot}{\code{signature(x = "PLSPCE", pce, options =c("fit", "bar", "compo"))}:
	 method of function \code{\link[graphics:plot]{plot}}.
	 Different plots are drawn, according to \code{options}.
	 When it includes \dQuote{fit}, 
	 computer model outputs against metamodel outputs.
	 When it includes \dQuote{bar},
	 barplots of the polynomial 
	 and total polynomial effects (PE and TPE).
	 When it includes \dQuote{compo},
	  TPE  against  components. This last plot requires the
  calculation of the TPE for each component, which may take some time.}
}}	   
       \seealso{
 Function \code{\link{calcPLSPCE}}, creator of objects from this
 class.
}

\keyword{classes}
