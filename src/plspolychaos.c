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
#include <R_ext/Rdynload.h>

#include "string.h" // for memcpy
#include "boucle.h"

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
/* --------------------------------------------------- */

/* calcHii
FUNCTION
calculate: B= t(T[1:hcur, , drop = FALSE]) %*% A
           Hii[i] = sum(B[i, ] * T[1:hcur, i]), i=1,n
INPUT
nc: nombre de lignes de T
hcur:  nombre de lignes a considerer dans T et
nombre de lignes et de colonnes de A
n: nombre de colonnes de T 
T: matrix nc X n
A: mtrix hcur X hcur
OUTPUT
Hii: vector n
*/
void calcHii(int *nc, int *hcur, int *n,
	     double *T, double *A, 
	     double *Hii) {
  int i,j,k;
  double B;

  for (i=0; i< *n; i++) {
    Hii[i] =0;
    for (j=0; j< *hcur; j++) {
      B = 0.0;
      for (k=0; k< *hcur; k++) {
	B += (T[i*(*nc)+k] * A[j* (*hcur)+k]);
      } /* fin k */
      Hii[i] += (T[i*(*nc)+j] * B);
    } /* fin j */
  } /* fin i */


} /* fin calcHii */

/* --------------------------------------------------- */
/*
boucleregpls
INPUT
Xold: matrix nlig X ncolX
YYold: matrix nlig X 1
INPUT/OUTPUT
somunew
WORKING
wold, wdif: matrix  ncolX (i.e=nbre de monomes) X 1
OUTPUT
wnew: matrix  ncolX (i.e=nbre de monomes) X 1
cnew:  matrix 1 X 1 (dim=1 car une seule reponse)
unew, tnew: matrix nlig X 1
pnew:   matrix ncolX X 1
 leRSS: scalar (colSums((YY.old - t.new %*% t(c.new))^2))
*/

void boucleregpls(int *nlig, int *ncolX, 
		  double *Xold, double * YYold, 
		  double *wnew, double *wold,
		  double *cnew, double *unew, double *tnew,
		  double *pnew,
		  double *wdif,
		  double *somunew,
		  double *leRSS) {
  /*  The allocations are done before the call */
  double somtnew2;
  int ncolYY=1; /* une seule reponse */


  somtnew2 = bouclewhile(*nlig, *ncolX,
		  Xold,  YYold, 
		  wnew, wold,
		  cnew, unew, tnew,
		  wdif,	 somunew);
  /* p.new <- t(X.old) %*% t.new/sum(t.new^2); P[hcur, ] <- p.new */
  multppsc(Xold, tnew, somtnew2, *nlig, *ncolX, pnew);
  /* c.new <- t(YY.old) %*% t.new/sum(t.new^2) */
  multppsc(YYold, tnew, somtnew2, *nlig, ncolYY,  cnew);
  /*  RSS[hcur + 1, ] <- colSums((YY.old - t.new %*% t(c.new))^2) */
  *leRSS = multsom( YYold, tnew, *cnew, *nlig);
  /* X.old <- X.old - t.new %*% t(p.new) */
  multXmoins(*nlig, *ncolX, Xold, tnew, pnew);
  /* YY.old <- YY.old - t.new %*% t(c.new) */
  multXmoins(*nlig, ncolYY, YYold, tnew, cnew);


} /* fin boucleregpls */



/* ++++++++++++++ INIT +++++++++++++++++++ */
static R_NativePrimitiveArgType Cpolleg2_t[] = {
  REALSXP, INTSXP, INTSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType TCpolleg2_t[] = {
  REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,REALSXP
};

static R_NativePrimitiveArgType boucleregpls_t[] = {
 INTSXP, INTSXP, 
 REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
 REALSXP, REALSXP, REALSXP,  REALSXP,  REALSXP
};

static R_NativePrimitiveArgType calcHii_t[] = {
  INTSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP
};


static const R_CMethodDef cMethods[] = {
  {"Cpolleg2", (DL_FUNC) &Cpolleg2, 5, Cpolleg2_t},
  {"TCpolleg2", (DL_FUNC) &TCpolleg2, 7, TCpolleg2_t},
  {"boucleregpls",  (DL_FUNC) &boucleregpls, 13, boucleregpls_t},
  {"calcHii",  (DL_FUNC) &calcHii, 6, calcHii_t},
  {NULL, NULL, 0}
};

void   R_init_plspolychaos(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL,  NULL, NULL);
  R_useDynamicSymbols(info, FALSE);

}
