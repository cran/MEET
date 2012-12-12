#include <R.h>
#include <Rmath.h>
#include <stdio.h>

/* Input Parameters */
/* wind = TFBS vector */
/* Prob = probability */

void probability(char **wind, int *nrowTFBS, double *Prob, double *R, double *a, double *t, double *c, double *g){

    double m=nrowTFBS[0];
	int i=0,j=0;
	char seq;
    
    a[0]=0;
    t[0]=0;
    c[0]=0;
    g[0]=0;
    
    Prob[0]=0,Prob[1]=0,Prob[2]=0,Prob[3]=0;
	
		for (i=0; i<m ;i++)
			{
            seq=*wind[i];
            if (seq=='A'){
			       	a[0]=a[0]+1;
                }
            if (seq=='T'){
					t[0]=t[0]+1;
                }
            if (seq=='C'){
                    c[0]=c[0]+1;
                }
            if (seq=='G'){
                    g[0]=g[0]+1;
                }
		}
    
    *R=m-(*a+*c+*t+*g);
    
    if (R[0]==0){
        Prob[0]=a[0]/m;
        Prob[1]=t[0]/m;
        Prob[2]=c[0]/m;
        Prob[3]=g[0]/m;
    }
}
