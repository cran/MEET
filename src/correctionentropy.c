#include <R.h>
#include <Rmath.h>
#include <stdio.h>

/* Input Parameters */
/* nrowTFBS = Number of sequences */
/* q = Renyi order*/
/* length */
/* Prob = Probability */
/* class= kind of Entropy */
/* Herror = Error entropy */
/* sderror = sd */

void correctionentropy(double *nrowTFBS, double *q, double *Prob, char **classentropy, double *Herror, double *VAR){
    
    int n=0;
    int ii;
    double a=0,t=0,c=0,g=0;
    double pp=0;
    double p=nrowTFBS[0];
    double probmat[4];
    double j=1,jj=1;
    double potA, potT, potC, potG;
    double Pn;
    double factA, factT, factC, factG, factPP;
    double HA=0,HT=0,HC=0,HG=0;
    double H;
    double E;
    double Etotal;
    double Var;
    double r=0;
    double m=q[0];
    double h=0;
    char class;
    
    
    class=*classentropy[0];
    
    probmat[0]=0, probmat[1]=0, probmat[2]=0, probmat[3]=0;
    
    for (pp=1; pp<=p; pp++){
        Etotal=0;
        Var=0;
        for (a=pp; a>=0;a--){
           
            for (t=(pp-a);t>=0; t--){
                 
                for (c=pp-(a+t);c>=0;c--){
                    
                    g=pp-(a+t+c);
                    
                    probmat[0]=a/pp, probmat[1]=t/pp;
                    probmat[2]=c/pp, probmat[3]=g/pp;
                    
                    factA=1,factT=1;
                    factC=1,factG=1;
                    factPP=1;
                    
                    if (a > 1){
                        for (j=2; j<=a; j++){
                            factA= factA*j;
                            }
                    }
                    if (t > 1){
                        for (j=2; j <= t; j++){
                            factT = factT*j;
                        }
                    }
                    if (c > 1){
                        for (j=2; j <= c; j++){
                            factC = factC*j;
                        }
                    }
                    if (g > 1){
                        for (j=2; j <= g; j++){
                            factG = factG*j;
                        }
                    }
                    if (pp > 1){
                        for (j=2; j <= pp; j++){
                            factPP = factPP*j;
                        }
                    }
                    
                    potA=pow(Prob[0],a),potT=pow(Prob[1],t);
                    potC=pow(Prob[2],c),potG=pow(Prob[3],g);
                    
                    Pn=(factPP*potA*potT*potC*potG)/(factA*factT*factC*factG);
                    
                    if (class=='S'){
                        
                        if (probmat[0]!=0){
                            HA=-probmat[0]*log2(probmat[0]);
                        }
                        if (probmat[1]!=0){
                            HT=-probmat[1]*log2(probmat[1]);
                        }
                        if (probmat[2]!=0){
                            HC=-probmat[2]*log2(probmat[2]);
                        }
                        if (probmat[3]!=0){
                            HG=-probmat[3]*log2(probmat[3]);
                        }
                       
                        H=HA+HT+HC+HG;
                        
                    } else {
                      
                       
                        h=pow(probmat[0],m)+pow(probmat[1],m)+pow(probmat[2],m)+pow(probmat[3],m);
                        
                        if (h>0) {
                            r=log2(h);
                            H=r/(1-m);
                        } else {
                            H=0;
                        }
                    }
                    E=H*Pn;
                    Etotal=E+Etotal;
                    Var=Var+E*H;
                    }
                }
            }
        n=pp-1;
        Herror[n]=Etotal;
        VAR[n]=Var-pow(Etotal,2);
        }
}
