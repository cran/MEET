#include <R.h>
#include <Rmath.h>
#include <stdio.h>

/* Input Parameters */

/* HXmax  = Maximum Entropy */
/* Herror = Given N binding sites sequnces, Herror is entropy error */
/* lengthHerror = Length of the vector */
/* r = Redundancy vector */
 

void correctionredundancy( double *HXmax, double *Herror, double *r, double *lengthRedundancy)
{   
   double herror=Herror[0];
   double s=lengthRedundancy[0];
   int i=0;
   double Rerror;
   double Hmax=HXmax[0];
    
   Rerror=1-(herror/Hmax);

   for (i=0; i<s; i++){
        if ( Rerror>r[i]){
            r[i]=0;
        }
    }
}
