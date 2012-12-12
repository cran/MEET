#include <R.h>
#include <Rmath.h>
#include <stdio.h>

/* Input parameters */
/* prob = Probability matrix */
/* q = Renyi Order */
/* H = Entropy */

void entropyRenyi(double *prob, double *lengthprob, double *q, double *H){
		int i;
		double h=0;
		double m=q[0];
        double r=0;
        double n=lengthprob[0];
        
    for (i=0;i<n;i++){
			h=h+pow(prob[i],m);
		}
    if (h>0) {
        r=log2(h);
        H[0]=r/(1-m);
    } else {
        H[0]=0;
    }
}
