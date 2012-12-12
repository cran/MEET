#include <R.h>
#include <Rmath.h>
#include <stdio.h>

/* Input Parameters */

/* HXmax  = Maximum Entropy */
/* ErrorHX = Given N binding sites sequnces, Herror is entropy error */
/* lengthH = Length of the vector */
/* HXmax= Maximun entropy */


void entropycorrected(double *H, double *lengthH, double *ErrorHX,double *HXmax){
	
	double w[1];
	double n=lengthH[0];
	int i;
	w[0]=HXmax[0]-ErrorHX[0];
	
	for (i=0;i<n;i++){
		if (H[i]>w[0]){
				H[i]=HXmax[0];
		} 
		if (H[i]<ErrorHX[0]){
			H[i]=0;
			}	
		}
}
