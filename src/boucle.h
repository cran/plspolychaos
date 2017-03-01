 double multsom( double *YYold, double *tnew, double cnew,
		 int nlig);
void multXmoins(int nrow, int ncol, double *Xold,
		double *tnew, double *pnew);
double multmat(double *mat1, double *mat2,
	       int nrow, int ncol, double *res);

double multmat2(double *mat1, double *mat2, double somme, 
		int nrow, int ncol, double *res);
double multpps(double *mat1, double *mat2, double somme,
	      int nrow, int ncol1,
	       double *res);
void multppsc(double *mat1, double *mat2, double somme,
	      int nrow, int ncol1,
	      double *res);
double bouclewhile(int nlig, int ncolX,
		  double *Xold, double * YYold, 
		  double *wnew, double *wold,
		  double *cnew, double *unew, double *tnew,
		  double *wdif,
		   double *somunew);
