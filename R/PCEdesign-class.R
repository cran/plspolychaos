###################################################################
# plspolychaos R package
# Copyright INRA 2017
# INRA, UR1404, Research Unit MaIAGE
# F78350 Jouy-en-Josas, France.
#
# URL: http://genome.jouy.inra.fr/logiciels/plspolychaos
#
# This file is part of plspolychaos R package.
# plspolychaos is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#
###################################################################
#########################################################
# Class to store the matrix which describes the polynomial
#########################################################
setClass("PCEdesign", slots=c(.Data="matrix", degree="numeric",
                        total.nmono="numeric"))
#########################################################
# print method
# INPUT
# x : return of Structure()
#########################################################
print.PCEdesign <- function(x, all = FALSE, ...) {
    # oter le terme constant
    planx <- x[-1, , drop = FALSE]
    # degree <- max(planx)
    degree <- x@degree
    cat("Total number of monomials:", x@total.nmono, "\n")
    nmono <- nrow(planx)
    if (nmono != x@total.nmono) {
      cat("Number of selected monomials:", nmono, "\n")
    }
    
    if (all) 
        {
            cat("Polynomial expression:\n")
            # inverser les colonnes planx <- planx[, seq(ncol(planx), 1, -1), drop=FALSE]
            descr <- paste(rownames(x), collapse = " + ")
            
            cat(descr, "\n")
            # Number of monomials per degree
            so <- apply(planx, 1, sum)
            for (i in 1:degree) {
                cat("Number of monomials of degree ", i, ": ", length(so[so == i]), 
                  "\n", sep = "")
            }  # fin i
        }  # fin all

    cat("Number of inputs: ", ncol(planx), "\n")
    cat("Polynomial degree: ", degree, "\n")
    return(invisible())
}  # end print.PCEdesign

#########################################################
# show method
# INPUT
#  return of Structure()
#########################################################
show.PCEdesign <- function(object) {
    print.PCEdesign(object)
    return(invisible())
}  # end show.PCEdesign


# --------------------------------------
setMethod("show", signature(object = "PCEdesign"), definition = show.PCEdesign) 
