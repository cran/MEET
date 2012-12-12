#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include "modelgetList.h"


SEXP getListElement(SEXP list, char *str) {
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;
    for ( i = 0; i < length(list); i++ )
        if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) {
            elmt = VECTOR_ELT(list, i);
            break;
        }
    
    if ( elmt == R_NilValue )
        error("%s missing from list", str);
    
    return elmt;
}
