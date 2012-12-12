#include <R.h>
#include <Rmath.h>
#include <stdio.h>

/* Input parameters */
/* wind Probability matrix */

void entropyShannon(double *prob, double *lengthprob, double *H){
	double n=lengthprob[0];
    int i;
    double h=0;
    
	for (i=0;i<n;i++){
        if (prob[i]!=0){
            h=h-prob[i]*log2(prob[i]);
        }
    }
    H[0]=h;
}

