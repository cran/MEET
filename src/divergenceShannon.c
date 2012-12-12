#include <R.h>
#include <Rmath.h>
#include <stdio.h>

/* Input Parameters */

/* HXmax  = Maximum Entropy */
/* Herror = Given N binding sites sequnces, Herror is entropy error */
/* lengthHerror = Length of the vector */
/* r = Redundancy vector */


void divergenceShannon(int *numpos,double *H,double Hxx[numpos[0]][numpos[0]], double HXY[numpos[0]][numpos[0]],
                        double mi[numpos[0]][numpos[0]],double *Herror,double *ErrorMI){
		
		int n=numpos[0];
		int i,j;
        double m=fabs(ErrorMI[0]);
		
        for (i=0; i<n;i++){
			for (j=0;j<n;j++){
				Hxx[i][j]=H[i]+H[j];
                mi[i][j]=Hxx[i][j]- HXY[i][j];
                if ((mi[i][j]+m)>Herror[0]){
                    mi[i][j]=Herror[0];
                    }
                if ((mi[i][j]-m)<0){
                    mi[i][j]=0;
                }
			}	
		}
}

	
