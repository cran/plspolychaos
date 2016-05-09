\name{descrdata}
\alias{descrdata}
\title{
Main Characteristics of the Dataset
}
\description{
Display the number of rows, the mean, standard deviation, range and
correlations of the inputs and output.
}
\usage{
descrdata(X, Y)
}
\arguments{
  \item{X}{
matrix with as many columns
    as inputs. Dataset of  inputs.
}
  \item{Y}{
    vector of length equal to the number of rows in \code{X}.
     Model outputs.
}}
\value{
  Nothing. It is a display function.
}
\examples{
### Load the dataset
load(system.file("extdata",  "ishigami200.Rda", package="plspolychaos"))
X <- ishi200[, -ncol(ishi200)] # inputs
Y <- ishi200[,  ncol(ishi200)] # output
### Data characteristics
descrdata(X, Y)
}
