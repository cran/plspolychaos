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
 Programmation in C of the loop of fastregpls2nomissing
 NOTE: in C, the values are stored by column
--------------------------------------------------------------- */
/* multsom
Internal function
FUNCTION
Compute:
RSS[hcur + 1, ] <- colSums((YY.old - t.new %*% t(c.new))^2)
YYold: matrix nlig X 1
tnew : vector of length nlig
cnew :  scalar
nlig: number of rows of YYold
OUTPUT
RSS[hcur + 1, ]
*/

    double multsom( double *YYold, double *tnew, double cnew,
		  int nlig) {

     int i; 
     double b, RSS;

      /*     RSS[h+1, j] =0.0; */
      RSS = 0.0;
      for (i=0; i < nlig; i++) {
	/* b = YYold[i,j] - ( tnew[i] * cnew[j]); */
	b= YYold[ i]  - ( tnew[i] * cnew);
	/* RSS[h+1, j] += (b*b); */
	 RSS += (b*b);
      } /* fin i */
      return(RSS);
    } /* fin multsom */

/* --------------------------------------------------- */
/* multXmoins
Internal function
FUNCTION
calculate
X.old <- X.old - t.new %*% t(p.new)
INPUT
nrow, ncol
Xold: matrix nrow X ncol
tnew: vector nrow
pnew: vector ncol
OUTPUT
Xold
*/

void multXmoins(int nrow, int ncol, double *Xold,
		double *tnew, double *pnew) {
  int i,j;
  for (i=0; i< nrow; i++) {
    for (j=0; j< ncol; j++) {
      Xold[j*nrow +i] -= (tnew[i] * pnew[j]);
    }
  }
} /* finmultXmoins */
/* --------------------------------------------------- */

/*  multmat
Internal function
FUNCTION
Multiplication of a  matrix (mat1) by the first column of t(mat2)
The number of columns of mat1 should be equal to the number 
of rows of mat2.
INPUT
mat1: nrow X ncol
mat2: nrow X ignored
OUTPUT
res:  ncol X 1
RETURN
somme(res*res)
*/

  double multmat(double *mat1, double *mat2,
	      int nrow, int ncol, double *res) {
  int i,j=0 ;
  double som2=0.0;

  for (i=0;  i<  nrow; i++) {
       res[i] =0.0;
       for (j=0; j< ncol; j++) {
	 /*	 Rprintf(" %g, %g ;", mat1[nrow * j+i] , mat2[j]); */
	 res[i]  += (mat1[nrow * j+i] * mat2[j]);
       }/* fin j */
       /*       Rprintf(" R[%d]=%g\n", l, res[i]); */
       som2 += (res[i]*res[i]);
   } /* fin i */
  return(som2);
} /* fin multmat */

/* --------------------------------------------------- */
/* multmat2
Internal function
Same as multmat, where mat2 is divided by 'somme'.
*/

double multmat2(double *mat1, double *mat2, double somme, 
	      int nrow, int ncol, double *res) {
  int i,j=0 ;
  double som2=0.0;
 
  for (i=0;  i<  nrow; i++) {
       res[i] =0.0;
       for (j=0; j< ncol; j++) {
	 /*	 Rprintf(" %g, %g ;", mat1[nrow * j+i] , mat2[j]); */
	 res[i]  += (mat1[nrow * j+i] * (mat2[j]/somme));
       }/* fin j */
       /*       Rprintf(" R[%d]=%g\n", l, res[i]); */
       som2 += ( res[i]* res[i]);

  } /* fin i */
  return(som2);

} /* fin multmat2 */



/* --------------------------------------------------- */
/* multpps
Internal function
FUNCTION
 Multiplication of a  matrix t(mat1) by the first column of mat2
Elements of mat2 are divided by 'somme'.
 mat1 should have the same number of rows as mat2.
Same as 'multpp', but with all the rows of mat1.
INPUT
mat1: nrow X ncol1
mat2: nrow X ignored
OUTPUT
res:  ncol1 X 1
RETURN
sum(res**2)
*/
  double multpps(double *mat1, double *mat2, double somme,
	      int nrow, int ncol1,
	      double *res) {
  int i,j,l=0 ;
  double som2=0.0;

  for (j=0; j< ncol1; j++) {
    res[j]=0.0;
    l=0;
    for (i=0; i< nrow; i++) {
	res[j] += (mat1[nrow * j+i] * (mat2[l++]/somme));
    }
    som2 += (res[j]*res[j]);
  } /* fin j */
  return(som2);
} /* fin multpps */
/* --------------------------------------------------- */
/* multppsc
Internal function
FUNCTION
 Multiplication of a  matrix t(mat1) by the first column of mat2
Elements of mat2 are divided by 'somme'.
 mat1 should have the same number of rows as mat2.
Same as 'multpps', but return nothing
INPUT
mat1: nrow X ncol1
mat2: nrow X ignored
OUTPUT
res:  ncol1 X 1
*/
  void multppsc(double *mat1, double *mat2, double somme,
	      int nrow, int ncol1,
	      double *res) {
  int i,j,l=0 ;

  for (j=0; j< ncol1; j++) {
    res[j]=0.0;
    l=0;
    for (i=0; i< nrow; i++) {
	res[j] += (mat1[nrow * j+i] * (mat2[l++]/somme));
    }
  } /* fin j */

} /* fin multppsc */

/* --------------------------------------------------- */
/*
bouclewhile
Internal function
INPUT
Xold: matrix nlig X ncolX
YYold: matrix nlig X 1
wnew, wold, wdif: matrix  ncolX (i.e=nbre de monomes) X 1
cnew:  matrix 1 X 1 (dim=1 car une seule reponse)
unew, tnew: matrix nlig X 1
RETURN
somtnew2: sum(tnew**2)
*/

double bouclewhile(int nlig, int ncolX,
		  double *Xold, double * YYold, 
		  double *wnew, double *wold,
		  double *cnew, double *unew, double *tnew,
		  double *wdif,
		  double *somunew) {
  /*  The allocations are done before the call */


  int i;
  double somwnew2, somcnew2,  somtnew2, somwdif2, somunew2;

  int iloop=0;
  int ncolYY=1; /* une seule reponse */
  somunew2=*somunew;

  while (1) {
    iloop++; /* pour s'assurer que le while se termine */

 /* w.new <- t(X.old) %*% u.new/sum(u.new^2) */
  somwnew2 = multpps(Xold, unew, somunew2, nlig, ncolX,   wnew);
/* w.new <- w.new/sqrt(sum(w.new^2)) */
  somwnew2 = R_pow(somwnew2, (double)0.5);

 for (i=0; i<ncolX; i++) {
	wnew[i] /= somwnew2;
      }
  /* t.new <- X.old %*% w.new */
  somtnew2 = multmat(Xold, wnew,  nlig, ncolX, tnew);

  /*  c.new <- t(YY.old) %*% t.new/sum(t.new^2) */
  somcnew2 = multpps(YYold, tnew, somtnew2, nlig, ncolYY,  cnew);

  /* u.new <- YY.old %*% c.new/sum(c.new^2) */
  somunew2 =  multmat2(YYold, cnew, somcnew2, nlig, ncolYY, unew);

  /* w.dif <- w.new - w.old ; w.old <- w.new */
  somwdif2=0.0;
  for (i=0; i< ncolX; i++) {
    wdif[i] = wnew[i] - wold[i];
    somwdif2 += ( wdif[i]* wdif[i]);
    wold[i] = wnew[i];
  }


  if ( (somwdif2 < 1.0e-12)|| (iloop >3000))	{ 
    break;
  }
  R_CheckUserInterrupt(); // check User interrupt
  } /* fin while */
  return(somtnew2);
		  
} /* fin bouclewhile */

