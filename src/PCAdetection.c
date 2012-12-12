/* Convert DNA to numeric DNA in C*/
#include <R.h>
#include <Rmath.h>
#include <stdio.h>


void PCAdetection(char **sequence, int *lengthDNA, int *lengthTFBS, double *Prob, double*center, double *loadings, int *nPCS, double *Scores){
  int i, j, k, l, ind;
  //char DNAm[lengthDNA][lengthTFBS]
  int n=lengthTFBS[0];
  int m=lengthDNA[0];
  char seq;
  double a[3], c[3], g[3], t[3], N[3];
  int dimensions=3*n;
  double seqNum[m][dimensions];
  int components=nPCS[0];
  double suma;
  double prediction[m][components];
  double xhat[m][dimensions];
  double residual[m][dimensions];
  double vector[dimensions];
  //double Scores[m];
 /*inicialitzem l'ADN*/
  a[0]=0;
  a[1]=0;
  a[2]=1;
  c[0]=-sqrt(2)/3;
  c[1]=sqrt(6)/3;
  c[2]=-1./3.;
  g[0]=-sqrt(2)/3;
  g[1]=-sqrt(6)/3;
  g[2]=-1./3.;
  t[0]=2.*sqrt(2)/3;
  t[1]=0.;	
  t[2]=-1./3.;
  for(i=0; i<3; i++){
    N[i]=Prob[0]*a[i]+Prob[1]*t[i]+Prob[2]*c[i]+Prob[3]*g[i];
    }
/*comencem la conversió a numèric*/
    for(i=0; i<(m); i++){
     // Rprintf("%c\n", *sequence[i]);
    for(j=0; j<n; j++){
       l=i+j;       
       seq=*sequence[l];
      //Rprintf("%c\n",seq);
	if(seq=='A'){
	  for(k=0; k<3;k++){
	  ind=3*j+k;
	  seqNum[i][ind]=a[k];
	  
	  }
	}else if(seq=='T'){
	  for(k=0; k<3;k++){
	  ind=3*j+k;
	  seqNum[i][ind]=t[k];
	  
	  }
	}else if(seq=='C'){
	  for(k=0; k<3;k++){
	  ind=3*j+k;
	  seqNum[i][ind]=c[k];
	
	  }
	}else if(seq=='G'){
	  for(k=0; k<3;k++){
	  ind=3*j+k;
	  seqNum[i][ind]=g[k];
	  ;
	  }
	}else if(seq=='-'){
	  for(k=0; k<3;k++){
	  ind=3*j+k;
	  seqNum[i][ind]=N[k];
	  
	  }
	} else{
	Rprintf("%c\n", seq);
	}
	
     
     }
    
  }
/*centrem la matriu*/
      for(i=0; i<m; i++){
	for(j=0; j<dimensions; j++){
	  seqNum[i][j]=seqNum[i][j]-center[j];
	}
      }
/*Projectem al subespai*/      
        for(i=0; i<m; i++){
 	for(j=0; j<components; j++){
 	    suma=0.;
 	  for(k=0; k<dimensions; k++){
 	    suma=suma+seqNum[i][k]*loadings[dimensions*j+k];
 	    //Rprintf("%f\n", suma);
 	 }
 	    
 	prediction[i][j]=suma;
 	}
       }
  
/*Calculem la nova matriu*/
	
     for(i=0; i<m; i++){
       for(j=0; j<dimensions; j++){
 	  suma=0.;
 	for(k=0; k<components; k++){
	    suma=suma+prediction[i][k]*loadings[dimensions*k+j];
 	    //Rprintf("m=%i",(dimensions*k+j)); 	    
 	}
// 	
  	  xhat[i][j]=suma;
  	  residual[i][j]=seqNum[i][j]-xhat[i][j];
	}
     }


      for(i=0; i<m; i++){
	suma=0;
        for(j=0; j<dimensions; j++){
	  vector[j]=residual[i][j];
	  suma=suma+vector[j]*vector[j];
        }
      Scores[i]=suma/dimensions;
      
      }
       


}

