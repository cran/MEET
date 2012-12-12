#include <R.h>
#include <Rmath.h>
#include <stdio.h>

/* Input Parameters */
/* rowTFBS= rows TFBS */
/* symbols = nunmber  A, C, G and T/*
/* R = number of missing values */
/* Prob = probability */
/* pm = probability matrix */

void missingfun(int *symbols, double *Prob, int *R,int *rowTFBS, double *pm)
 	{
 	int i;	
 	int m=R[0];
 	double n;
 	n=1./(*rowTFBS);

	for(i=0;i<4;i++){
		pm[i]=n*(symbols[i]+(Prob[i]*m));
	}
}
