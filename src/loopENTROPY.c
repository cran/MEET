#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rdefines.h>
#include "modelgetList.h"


SEXP loopENTROPY(SEXP validationsetx, SEXP lengthDNA,SEXP ncoltraining,SEXP Redundanciacorregida,SEXP iicc, SEXP Entropy,SEXP nameEntropy, SEXP MachineDouble){
    
    int i,j;
    lengthDNA=coerceVector(lengthDNA,INTSXP);
    ncoltraining=coerceVector(ncoltraining,INTSXP);
    int n=INTEGER(ncoltraining)[0];
    int l=INTEGER(lengthDNA)[0];
    int m=l-n+1;
    
    double *HXmax=REAL(getListElement(iicc,"HXmax"));
    double *Herror=REAL(getListElement(iicc,"Herror"));
    double Rerror=(1-Herror[0]/HXmax[0]);
    MachineDouble=coerceVector(MachineDouble,REALSXP);
    Redundanciacorregida=coerceVector(Redundanciacorregida,REALSXP);
    
    SEXP R;
    SEXP Rout;
    PROTECT(R=allocVector(REALSXP,n));
    PROTECT(Rout=allocVector(REALSXP,m));
    PROTECT(nameEntropy=AS_CHARACTER(nameEntropy));

    for (j=0;j<m;j++){
        char *seqrand[n];
        char *Pmychar[1];
        double x=0;
        for (i=0;i<n;i++){
            
            int index;
            REAL(R)[i]=0;
            seqrand[i]=R_alloc(strlen(CHAR(STRING_ELT(validationsetx, 0))), sizeof(char));
            strcpy(seqrand[i], CHAR(STRING_ELT(validationsetx, i+j)));
        
            if (*seqrand[i]=='A') {
                index=0;
            }else if (*seqrand[i]=='T') {
                index=1;
            }else if (*seqrand[i]=='C') {
                index=2;
            }else {
                index=3;
            }
            
            Pmychar[0]=R_alloc(strlen(CHAR(STRING_ELT(nameEntropy, 0))), sizeof(char));
            strcpy(Pmychar[0], CHAR(STRING_ELT(nameEntropy,index)));
            double *v = REAL(getListElement(Entropy, Pmychar[0]));
            REAL(R)[i]=(1-v[i]/HXmax[0]);
           
            if ( Rerror > REAL(R)[i]){
                REAL(R)[i]=0;
            }
            
            if (REAL(R)[i] > REAL(Redundanciacorregida)[i]){
                x=x+(REAL(R)[i]-REAL(Redundanciacorregida)[i])*REAL(Redundanciacorregida)[i];
           
            }else{
                x=x+(REAL(Redundanciacorregida)[i]-REAL(R)[i])*REAL(Redundanciacorregida)[i];
            }
        }
        
        REAL(Rout)[j]=1/(x+REAL(MachineDouble)[0]);
    }
    UNPROTECT(3);
    return Rout;
}
