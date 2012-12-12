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

void powX2(int *numfila,int *numpos, double *pmX1,double *pmX2, double pmXY[numfila[0]][numfila[0]],double *q,double rr[numfila[0]][numfila[0]]){
	
	int m=numfila[0];
	int k,p;
    double order=q[0];
    double v,w,vw,ww;
    
  
   
  	for (k=0;k<m;k++){
        for (p=0;p<m;p++){
            v= pmX1[k];
            w= pmX2[p];
            vw=pmXY[k][p];
            ww=vw-(v*w);
            rr[k][p]=pow(ww,2)/(v*w+2.220446e-16);
 		   }
        }
}
