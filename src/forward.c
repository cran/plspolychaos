/* ---------------------------------------------------------------
  plspolychaos R package
  Copyright INRA 2016
  INRA, UR1404, Research Unit MaIAGE
  F78352 Jouy-en-Josas, France.
 
  URL: http://cran.r-project.org/web/packages/plspolychaos
 
  This file is part of plspolychaos R package.
  plspolychaos is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  See the GNU General Public License at:
  http://www.gnu.org/licenses/
 
----------- --------------------------------------------------- */

#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h> // to allow user interrupt
#include "string.h" // for memcpy

/* ---------------------------------------------------------------
 Programmation in C of the computation of Legendre polynomials
 NOTE: in C, the values are stored by column
----------- --------------------------------------------------- */
  
void Cpolleg2(double *xin, int *struc,
	       int *nl, int *nvx, double *ret) {
  /* 
  # xin: matrix (nl x nvx) (number of rows x number of inputs)
  #   Calibrated inputs
  # struc : vector (nvx). Polynomial degree of each variable
  # of the monomial
  # RETURN
  # ret: vector (nl) product of the Legendre coded inputs of the monome
  */
/*
  The allocations are done before the call
 ret= (double *) S_alloc( *nl, sizeof(double));
  */

  int i,j, ind;
  double Leg=0;

  for (i=0; i<*nl; i++) {
    ret[i] = 1;
    for (j=0; j<*nvx; j++) {
      ind = (j* (*nl)) +i;
      switch( struc[j])  {
      case 1 :
	Leg = xin[ind];
	break;
      case 2:
	Leg = (1.0/2.0)*(3.0*xin[ind]*xin[ind]-1.0); 
	break;
      case 3:
	Leg = (1.0/2.0)*((5.0* R_pow_di(xin[ind],3))-3.0*xin[ind]); 
	break;
      case 4:
	Leg = (1.0/8.0)*((35.0 * R_pow_di(xin[ind],4))-(30.0 *xin[ind]*xin[ind])+3.0);
	break;
      case 5:
	Leg = (1.0/8.0)*((63.0 * R_pow_di(xin[ind],5))-(70.0 * R_pow_di(xin[ind],3))+15.0 *xin[ind]);
	break;
      case 6:
	Leg = (1.0/(16.0))*((231.0 * R_pow_di(xin[ind],6))-(315.0 * R_pow_di(xin[ind],4))
                                +(105.0 *xin[ind]*xin[ind])-5.0);
	break;
      case 7:
	Leg = (1.0/16.0)*(429.0 * R_pow_di(xin[ind],7)-693.0 * R_pow_di(xin[ind],5)
			+315.0 * R_pow_di(xin[ind],3)-35.0 *xin[ind]);
	break;
      case 8:
	Leg = (1.0/128.0)*(6435.0 * R_pow_di(xin[ind],8)-12012.0 * R_pow_di(xin[ind],6)
			 +6930.0 * R_pow_di(xin[ind],4)-1260.0 *xin[ind]*xin[ind]+35.0);
	break;
      case 9:
	Leg = (1.0/128.0)*(12155.0 * R_pow_di(xin[ind],9)-25740.0 * R_pow_di(xin[ind],7)
			 +18018.0 * R_pow_di(xin[ind],5)-4620.0 * R_pow_di(xin[ind],3)+315.0 *xin[ind]);
	break;
      case 10:
	Leg = (1.0/256.0)*(46189.0 * R_pow_di(xin[ind],10)-109395.0 * R_pow_di(xin[ind],8)
   +90090.0 * R_pow_di(xin[ind],6)-30030.0 * R_pow_di(xin[ind],4) +3465.0 *xin[ind]*xin[ind]-63.0);
	break;
      default:
	/* degree=0 */
	if (struc[j] !=0) {
	  Rprintf("degree=%d\n", struc[j]);
	  error("polleg2: bad argument struc\n");
	}
	Leg = 1;
      } /* fin switch */

        ret[i] = ret[i] * Leg;

    } // fin j


  } // fin i
  
} // fin Cpolleg2

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void TCpolleg2(double *xin, int *struc,
	      int *nmono, 
	       int *nl, int *nvx, int *trav, double *ret) {
  /* 
  # xin: matrix (nl x nvx) (number of rows x number of inputs)
  #   Calibrated inputs
  # struc : matrix (nmono, nvx). Polynomial degree of each variable
  # of the monomials
  # trav: vector (nvx) Working array
  # RETURN
  # ret: matrix (nl, nmono) product of the Legendre coded inputs of the monomials
  */
/*
  The allocations are done before the call
 ret= (double *) S_alloc( *nl * *nmono, sizeof(double));
 trav = (int *) S_alloc( *nvx, sizeof(int));
  */

/* Loop over Cpolleg2 */
  int j, imono;
  for (imono=0; imono<*nmono; imono++) {
    for (j=0; j< *nvx; j++) {
      trav[j]=struc[(j * (*nmono))+ imono];
    }
    Cpolleg2(xin, trav, nl, nvx, ret+(imono* (*nl)));
  }



} // fin TCpolleg2

