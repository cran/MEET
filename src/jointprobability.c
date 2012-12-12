#include <R.h>
#include <Rmath.h>
#include <stdio.h>

/* Input Parameters */
/* trainingset = A set of Transcription Factor Binding Sites sequences */
/* nrowTFBS = Number of sequences */
/* Prob = Probability vector */
/* Probtrans = Probability couple */


void jointprobability(char **trainingset,int * indexA, int *indexB,int *nrowTFBS, double *Probtrans,double *background, double *Prob){
	
	int i=0,j=0;
    double m=nrowTFBS[0];
	char seq[2];
    double N;
    double AA=0,AT=0,AC=0,AG=0,AR=0;
    double TA=0,TT=0,TC=0,TG=0,TR=0;
    double CA=0,CT=0,CC=0,CG=0,CR=0;
    double GA=0,GT=0,GC=0,GG=0,GR=0;
    double RA=0,RT=0,RC=0,RG=0,RR=0;
    
    
    Prob[0]=0,Prob[1]=0,Prob[2]=0,Prob[3]=0;
    Prob[4]=0,Prob[5]=0,Prob[6]=0,Prob[7]=0;
    Prob[8]=0,Prob[9]=0,Prob[10]=0,Prob[11]=0;
    Prob[12]=0,Prob[13]=0,Prob[14]=0,Prob[15]=0;
    
    
    
	for (i=0; i<m; i++)
		{
			seq[0]=trainingset[i][0];
            seq[1]=trainingset[i][1];
            
            if (seq[0]=='A'){
                    if (seq[1]=='A'){
                        AA=AA+1;
                        }
                    if (seq[1]=='T'){
                        AT=AT+1;
                        }
                    if (seq[1]=='C'){
                        AC=AC+1;
                        }
                    if (seq[1]=='G'){
                        AG=AG+1;
                        }
                    if (seq[1]=='-'){
                        AR=AR+1;
                        }
                    }
            if (seq[0]=='T'){
                    if (seq[1]=='A'){
                        TA=TA+1;
                    }
                    if (seq[1]=='T'){
                        TT=TT+1;
                    }
                    if (seq[1]=='C'){
                        TC=TC+1;
                    }
                    if (seq[1]=='G'){
                        TG=TG+1;
                    }
                    if (seq[1]=='-'){
                        TR=TR+1;
                    }
                }
            if (seq[0]=='C'){
                    if (seq[1]=='A'){
                        CA=CA+1;
                    }
                    if (seq[1]=='T'){
                        CT=CT+1;
                    }
                    if (seq[1]=='C'){
                        CC=CC+1;
                    }
                    if (seq[1]=='G'){
                        CG=CG+1;
                    }
                    if (seq[1]=='-'){
                        CR=CR+1;
                    }
                }
                if (seq[0]=='G'){
                    if (seq[1]=='A'){
                        GA=GA+1;
                    }
                    if (seq[1]=='T'){
                        GT=GT+1;
                    }
                    if (seq[1]=='C'){
                        GC=GC+1;
                    }
                    if (seq[1]=='G'){
                        GG=GG+1;
                    }
                    if (seq[1]=='-'){
                    GR=GR+1;
                    }
                }
                if (seq[0]=='-'){
                    if (seq[1]=='A'){
                        RA=RA+1;
                    }
                    if (seq[1]=='T'){
                        RT=RT+1;
                    }
                    if (seq[1]=='C'){
                        RC=RC+1;
                    }
                    if (seq[1]=='G'){
                        RG=RG+1;
                    }
                    if (seq[1]=='-'){
                        RR=RR+1;
                }
            }
        }
    
    N=AA+AT+AC+AG+TA+TT+TC+TG+CA+CT+CC+CG+GA+GT+GC+GG;

        if (N==m) {
        
                Prob[0]=AA/m,Prob[1]=AT/m,Prob[2]=AC/m,Prob[3]=AG/m;
                Prob[4]=TA/m,Prob[5]=TT/m,Prob[6]=TC/m,Prob[7]=TG/m;
                Prob[8]=CA/m,Prob[9]=CT/m,Prob[10]=CC/m,Prob[11]=CG/m;
                Prob[12]=GA/m,Prob[13]=GT/m,Prob[14]=GC/m,Prob[15]=GG/m;
            
            } else {
                
                if (indexA[0]!=indexB[0]){

                    Prob[0]=(AA+RA*background[0]+AR*background[0]+RR*Probtrans[0])/m;
                    Prob[1]=(AT+RT*background[0]+AR*background[1]+RR*Probtrans[1])/m;
                    Prob[2]=(AC+RC*background[0]+AR*background[2]+RR*Probtrans[2])/m;
                    Prob[3]=(AG+RG*background[0]+AR*background[3]+RR*Probtrans[3])/m;
            
                    Prob[4]=(TA+RA*background[1]+TR*background[0]+RR*Probtrans[4])/m;
                    Prob[5]=(TT+RT*background[1]+TR*background[1]+RR*Probtrans[5])/m;
                    Prob[6]=(TC+RC*background[1]+TR*background[2]+RR*Probtrans[6])/m;
                    Prob[7]=(TG+RG*background[1]+TR*background[3]+RR*Probtrans[7])/m;
            
                    Prob[8]=(CA+RA*background[2]+CR*background[0]+RR*Probtrans[8])/m;
                    Prob[9]=(CT+RT*Prob[2]+CR*background[1]+RR*Probtrans[9])/m;
                    Prob[10]=(CC+RC*background[2]+CR*background[2]+RR*Probtrans[10])/m;
                    Prob[11]=(CG+RG*background[2]+CR*background[3]+RR*Probtrans[11])/m;
                
                    Prob[12]=(GA+RA*background[3]+GR*background[0]+RR*Probtrans[12])/m;
                    Prob[13]=(GT+RT*background[3]+GR*background[1]+RR*Probtrans[13])/m;
                    Prob[14]=(GC+RC*background[3]+GR*background[2]+RR*Probtrans[14])/m;
                    Prob[15]=(GG+RG*background[3]+GR*background[3]+RR*Probtrans[15])/m;
                   
                    
                } else{

                    Prob[0]=(AA+RR*Probtrans[0])/m;
                    Prob[5]=(TT+RR*Probtrans[5])/m;
                    Prob[10]=(CC+RR*Probtrans[10])/m;
                    Prob[15]=(GG+RR*Probtrans[15])/m;

                }
            }
            /*Rprintf("%.*s\n", 2, seq)*/ ;
}
