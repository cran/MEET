#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rdefines.h>
#include "modelgetList.h"

SEXP readH(SEXP iicc,SEXP p, SEXP ncoltraining){
	int n;
	char *Pmychar[1];
	SEXP HX;
	ncoltraining=coerceVector(ncoltraining,INTSXP);
    int m=INTEGER(ncoltraining)[0];
    PROTECT(p=AS_CHARACTER(p));
    PROTECT(HX=allocVector(REALSXP,m));
	Pmychar[0]=R_alloc(strlen(CHAR(STRING_ELT(p, 0))), sizeof(char));
	
	for (n=0;n<m;n++){
		REAL(HX)[n]=0;
        strcpy(Pmychar[0], CHAR(STRING_ELT(p,n)));
		double *v = REAL(getListElement(iicc, Pmychar[0]));
		REAL(HX)[n]=v[n];	
	}
	UNPROTECT(2);
	return HX;
}
