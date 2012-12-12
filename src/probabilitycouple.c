#include <R.h>
#include <Rmath.h>
#include <stdio.h>

/* Input Parameters */
/* background= probability */


void probabilitycouple(double *background, double *Probtrans){
	
		Probtrans[0]=background[0]*background[0], Probtrans[1]=background[0]*background[1], Probtrans[2]=background[0]*background[2], Probtrans[3]=background[0]*background[3];
		Probtrans[4]=background[1]*background[0], Probtrans[5]=background[1]*background[1], Probtrans[6]=background[1]*background[2], Probtrans[7]=background[1]*background[3];
		Probtrans[8]=background[2]*background[0], Probtrans[9]=background[2]*background[1], Probtrans[10]=background[2]*background[2], Probtrans[11]=background[2]*background[3];
		Probtrans[12]=background[3]*background[0],Probtrans[13]=background[3]*background[1],Probtrans[14]=background[3]*background[2], Probtrans[15]=background[3]*background[3];
		
}
