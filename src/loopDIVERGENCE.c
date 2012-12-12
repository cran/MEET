#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rdefines.h>
#include "modelgetList.h"


SEXP loopDIVERGENCE(SEXP validationsetx, SEXP lengthDNA,SEXP ncoltraining,SEXP lengthname, SEXP interA, SEXP interB, SEXP nameDivergence, SEXP MachineDouble,SEXP Divergence,SEXP D,SEXP Mperfil){
    
    int i,ii,j,jj,k;
    char *Pmychar[1];
    
    lengthDNA=coerceVector(lengthDNA,INTSXP);
    lengthname=coerceVector(lengthname,INTSXP);
    ncoltraining=coerceVector(ncoltraining,INTSXP);
    int n=INTEGER(ncoltraining)[0];
    int l=INTEGER(lengthDNA)[0];
    int nn=INTEGER(lengthname)[0];
    int m=l-n+1;
    
    MachineDouble=coerceVector(MachineDouble,REALSXP);
    interA=coerceVector(interA,INTSXP);
    interB=coerceVector(interB,INTSXP);
    
    SEXP Dval;
    SEXP Dout;
    PROTECT(Dout=allocVector(REALSXP,m));
    PROTECT(Dval = allocVector(REALSXP, n*n));
    PROTECT(nameDivergence=AS_CHARACTER(nameDivergence));
    Mperfil=coerceVector(Mperfil,REALSXP);
    D=coerceVector(D,REALSXP);
    Pmychar[0] = R_alloc(strlen(CHAR(STRING_ELT(nameDivergence, 0))), sizeof(char));
    double *mat;
    for (j=0;j<m;j++){
        char *seqrand[n];
        int p[n];
        double x=0;
        for (k=0;k<(n*n);k++){
            REAL(Dval)[k]=0;
        }
        for (i=0;i<n;i++){
            
            seqrand[i]=R_alloc(strlen(CHAR(STRING_ELT(validationsetx, 0))), sizeof(char));
            strcpy(seqrand[i], CHAR(STRING_ELT(validationsetx, i+j)));
           
            if (*seqrand[i]=='A') {
                p[i]=1;
            }else if (*seqrand[i]=='T') {
                p[i]=2;
            }else if (*seqrand[i]=='C') {
                p[i]=3;
            }else {
                p[i]=4;
            }
            
        }
        for (jj=0;jj<nn;jj++){
            strcpy(Pmychar[0], CHAR(STRING_ELT(nameDivergence, jj)));
            double *v = REAL(getListElement(Divergence, Pmychar[0]));
            int A=INTEGER(interA)[jj];
            int B=INTEGER(interB)[jj];
            int AB=(B-1)+n*(A-1);
            int PB=p[B-1];
            int PA=p[A-1];
            int PAB= (PB-1)+4*(PA-1);
            REAL(Dval)[AB]=v[PAB];
        }
        for (i=0;i<n;i++){
              for (ii=0;ii<n;ii++){
                  if (  REAL(Dval)[n*i+ii] > REAL(D)[n*i+ii]){
                      x= x+ (REAL(Dval)[n*i+ii] - REAL(D)[n*i+ii])*REAL(Mperfil)[n*i+ii];
                  }else{
                      x= x+(REAL(D)[n*i+ii] - REAL(Dval)[n*i+ii])*REAL(Mperfil)[n*i+ii];
                  }
              }
         }
        REAL(Dout)[j]=1/(x+REAL(MachineDouble)[0]);
    }
    UNPROTECT(3);
    return Dout;
}
