#include <R.h>
#include <Rmath.h>
#include <stdio.h>

/* Input Parameters */

/* numfila  = number row */
/* numpos = number columns */
/* pmX1= probability matrix vector */
/* pmX2= probability matrix vector */
/* pmXY= joint probability */
/* q= Renyi Order */
/* rr= power Renyi */

void powRenyi( int *numfila, int *numpos, double *pmX1,double *pmX2, double pmXY[numfila[0]][numfila[0]],double *q,double rr[numfila[0]][numfila[0]]){
	
	int m=numfila[0];
	int k,p;
    double order=q[0];
    double L;
    L=1-order;
    double v,w,vw;
    
    Rprintf("L=%f\n",L) ;
    Rprintf("m=%i\n",m) ;
  
   
  	for (k=0;k<m;k++){
        for (p=0;p<m;p++){
            v= pmX1[k];
            w= pmX2[p];
            vw=pmXY[k][p];
            rr[k][p]=pow(vw,order)*pow(v*w+2.220446e-16,L);
          }
    }
}
