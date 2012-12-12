#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rdefines.h>
#include "modelgetList.h"

SEXP MI(SEXP iicc, SEXP name,SEXP lengthname, SEXP ncoltraining, SEXP p,SEXP interA,SEXP interB){
    
    int i,j,k;
    char *Pmychar[1];
    SEXP Dval;
    lengthname=coerceVector(lengthname,INTSXP);
    ncoltraining=coerceVector(ncoltraining,INTSXP);
    interB=coerceVector(interB,INTSXP);
    interA=coerceVector(interA,INTSXP);
    p=coerceVector(p,INTSXP);
   
    int n=INTEGER(lengthname)[0];
    int m=INTEGER(ncoltraining)[0];
    PROTECT(name=AS_CHARACTER(name));
    PROTECT(Dval = allocVector(REALSXP, m*m));
    Pmychar[0] = R_alloc(strlen(CHAR(STRING_ELT(name, 0))), sizeof(char));
    
    for (k=0;k<(m*m);k++){
            REAL(Dval)[k]=0;
    }
    for (j=0;j<n;j++){
          strcpy(Pmychar[0], CHAR(STRING_ELT(name, j)));
          double *v = REAL(getListElement(iicc, Pmychar[0]));
          int A=INTEGER(interA)[j];
          int B=INTEGER(interB)[j];
          int AB=(B-1)+m*(A-1);
          int PB=INTEGER(p)[B-1];
          int PA=INTEGER(p)[A-1];
          int PAB= (PB-1)+4*(PA-1);
          REAL(Dval)[AB]=v[PAB];
    }
    UNPROTECT(2);
    return Dval;
}




