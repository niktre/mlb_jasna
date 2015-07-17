#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defs.h"

/*********************************************************/
/**first derivative of f at x,y*/
/*double der1(DerivativeOf what, int x, int y, int direction){
	double der_f;
	double fm1,fm2,fp1,fp2; //fm1=function at position minus 1, fp2=function at position plus 2
	int xm1, xp1, xm2, xp2, ym1, yp1, ym2, yp2;
	int a, b;
	
	if(direction==0){
		xm1=x-1;
		xm2=x-2;
		xp1=x+1;
		xp2=x+2;
		yp1=yp2=ym1=ym2=y;		

		if (xm1 >= XDIM) xm1-=XDIM;
		if (xm1 < 0) xm1+=XDIM;
		if (xm2 >= XDIM) xm2-=XDIM;
		if (xm2 < 0) xm2+=XDIM;
		if (xp1 >= XDIM) xp1-=XDIM;
		if (xp1 < 0) xp1+=XDIM;
		if (xp2 >= XDIM) xp2-=XDIM;
		if (xp2 < 0) xp2+=XDIM;
	}
	else if(direction==1){
		ym1=y-1;
		ym2=y-2;
		yp1=y+1;
		yp2=y+2;
		xp1=xp2=xm1=xm2=x;		

		if (ym1 >= YDIM) ym1-=YDIM;
		if (ym1 < 0) ym1+=YDIM;
		if (ym2 >= YDIM) ym2-=YDIM;
		if (ym2 < 0) ym2+=YDIM;
		if (yp1 >= YDIM) yp1-=YDIM;
		if (yp1 < 0) yp1+=YDIM;
		if (yp2 >= YDIM) yp2-=YDIM;
		if (yp2 < 0) yp2+=YDIM;
	}

	switch (what){
		case DENSITY:
			fm1 = node[YDIM*xm1+ym1].rho;
			fm2 = node[YDIM*xm2+ym2].rho;
			fp1 = node[YDIM*xp1+yp1].rho;
			fp2 = node[YDIM*xp2+yp2].rho;
			break;
		case UX:
			fm1 = node[YDIM*xm1+ym1].u[0];
			fm2 = node[YDIM*xm2+ym2].u[0];
			fp1 = node[YDIM*xp1+yp1].u[0];
			fp2 = node[YDIM*xp2+yp2].u[0];
			break; 
		case UY:
			fm1 = node[YDIM*xm1+ym1].u[1];
			fm2 = node[YDIM*xm2+ym2].u[1];
			fp1 = node[YDIM*xp1+yp1].u[1];
			fp2 = node[YDIM*xp2+yp2].u[1];
			break; 
		case PUXUX:
			a=0; b=0;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].S*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].S*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].S*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].S*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			break;
		case PUYUY:
			a=1; b=1;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].S*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].S*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].S*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].S*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			break;
		case PUXUY:
			a=0; b=1;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].S*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].S*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].S*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].S*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			break;
		case PRESSURE:
			fm1=node[YDIM*xm1+ym1].S*node[YDIM*xm1+ym1].rho;
			fm2=node[YDIM*xm2+ym2].S*node[YDIM*xm2+ym2].rho;
			fp1=node[YDIM*xp1+yp1].S*node[YDIM*xp1+yp1].rho;
			fp2=node[YDIM*xp2+yp2].S*node[YDIM*xp2+yp2].rho;
			break;
		case JX:
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[0];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].u[0];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[0];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].u[0];
			break; 
		case JY:
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[1];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].u[1];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[1];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].u[1];
			break; 
		case RHOUXUX:
			a=0; b=0;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			break;
		case RHOUYUY:
			a=1; b=1;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			break;
		case RHOUXUY:
			a=0; b=1;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			break;
		case DERP:
			fm1=derP(node[YDIM*xm1+ym1].rho);
			fm2=derP(node[YDIM*xm2+ym2].rho);
			fp1=derP(node[YDIM*xp1+yp1].rho);
			fp2=derP(node[YDIM*xp2+yp2].rho);
			break;
		case SIGMA:
			fm1=node[YDIM*xm1+ym1].S;
			fm2=node[YDIM*xm2+ym2].S;
			fp1=node[YDIM*xp1+yp1].S;
			fp2=node[YDIM*xp2+yp2].S;
			break;

	}
	der_f = 2./3. * (fp1-fm1) - 1./12. * (fp2 - fm2);
//	der_f = (8.0*fp1 - 8.0*fm1 - fp2 + fm2)/12.0;
return der_f;
}
*/

/*********************************************************/
/**first derivative of f at x,y higher accuracy*/
double der1v2(DerivativeOf what, int x, int y, int direction){
	double der_f;
	double fm1,fm2,fp1,fp2,fm3,fp3; //fm1=function at position minus 1, fp2=function at position plus 2
	int xm1, xp1, xm2, xp2, ym1, yp1, ym2, yp2, xm3, xp3, ym3, yp3; 
	int a, b;
	
	if(direction==0){
		xm1=x-1;
		xm2=x-2;
		xp1=x+1;
		xp2=x+2;
		xp3=x+3;
		xm3=x-3;
		yp1=yp2=ym1=ym2=ym3=yp3=y;		

		if(xm1 >= XDIM) xm1 -= XDIM;
		if(xm1 < 0) xm1 += XDIM;
		if(xm2 >= XDIM) xm2 -= XDIM;
		if(xm2 < 0) xm2 += XDIM;
		if(xp1 >= XDIM) xp1 -= XDIM;
		if(xp1 < 0) xp1 += XDIM;
		if(xp2 >= XDIM) xp2 -= XDIM;
		if(xp2 < 0) xp2 += XDIM;
		if(xm3 >= XDIM) xm3 -= XDIM;
		if(xm3 < 0) xm3 += XDIM;
		if(xp3 >= XDIM) xp3 -= XDIM;
		if(xp3 < 0) xp3 += XDIM;
	}
	else if(direction==1){
		ym1=y-1;
		ym2=y-2;
		yp1=y+1;
		yp2=y+2;
		yp3=y+3;
		ym3=y-3;
		xp1=xp2=xm1=xm2=xm3=xp3=x;		

		if (ym1 >= YDIM) ym1-=YDIM;
		if (ym1 < 0) ym1+=YDIM;
		if (ym2 >= YDIM) ym2-=YDIM;
		if (ym2 < 0) ym2+=YDIM;
		if (yp1 >= YDIM) yp1-=YDIM;
		if (yp1 < 0) yp1+=YDIM;
		if (yp2 >= YDIM) yp2-=YDIM;
		if (yp2 < 0) yp2+=YDIM;
		if (ym3 >= YDIM) ym3-=YDIM;
		if (ym3 < 0) ym3+=YDIM;
		if (yp3 >= YDIM) yp3-=YDIM;
		if (yp3 < 0) yp3+=YDIM;
	}
		switch (what){
		case DENSITY:
			fm1=node[YDIM*xm1+ym1].rho;
			fm2=node[YDIM*xm2+ym2].rho;
			fp1=node[YDIM*xp1+yp1].rho;
			fp2=node[YDIM*xp2+yp2].rho;
			fm3=node[YDIM*xm3+ym3].rho;
			fp3=node[YDIM*xp3+yp3].rho;
			break;
		case UX:
			fm1=node[YDIM*xm1+ym1].u[0];
			fm2=node[YDIM*xm2+ym2].u[0];
			fp1=node[YDIM*xp1+yp1].u[0];
			fp2=node[YDIM*xp2+yp2].u[0];
			fm3=node[YDIM*xm3+ym3].u[0];
			fp3=node[YDIM*xp3+yp3].u[0];
			break; 
		case UY:
			fm1=node[YDIM*xm1+ym1].u[1];
			fm2=node[YDIM*xm2+ym2].u[1];
			fp1=node[YDIM*xp1+yp1].u[1];
			fp2=node[YDIM*xp2+yp2].u[1];
			fm3=node[YDIM*xm3+ym3].u[1];
			fp3=node[YDIM*xp3+yp3].u[1];
			break; 
		case PUXUX:
			a=0; b=0;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].S*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].S*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].S*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].S*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			fm3=node[YDIM*xm3+ym3].rho*node[YDIM*xm3+ym3].S*node[YDIM*xm3+ym3].u[a]*node[YDIM*xm3+ym3].u[b];
			fp3=node[YDIM*xp3+yp3].rho*node[YDIM*xp3+yp3].S*node[YDIM*xp3+yp3].u[a]*node[YDIM*xp3+yp3].u[b];
			break;
		case PUYUY:
			a=1; b=1;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].S*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].S*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].S*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].S*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			fm3=node[YDIM*xm3+ym3].rho*node[YDIM*xm3+ym3].S*node[YDIM*xm3+ym3].u[a]*node[YDIM*xm3+ym3].u[b];
			fp3=node[YDIM*xp3+yp3].rho*node[YDIM*xp3+yp3].S*node[YDIM*xp3+yp3].u[a]*node[YDIM*xp3+yp3].u[b];
			break;
		case PUXUY:
			a=0; b=1;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].S*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].S*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].S*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].S*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			fm3=node[YDIM*xm3+ym3].rho*node[YDIM*xm3+ym3].S*node[YDIM*xm3+ym3].u[a]*node[YDIM*xm3+ym3].u[b];
			fp3=node[YDIM*xp3+yp3].rho*node[YDIM*xp3+yp3].S*node[YDIM*xp3+yp3].u[a]*node[YDIM*xp3+yp3].u[b];
			break;
		case PRESSURE:
			fm1=node[YDIM*xm1+ym1].S*node[YDIM*xm1+ym1].rho;
			fm2=node[YDIM*xm2+ym2].S*node[YDIM*xm2+ym2].rho;
			fp1=node[YDIM*xp1+yp1].S*node[YDIM*xp1+yp1].rho;
			fp2=node[YDIM*xp2+yp2].S*node[YDIM*xp2+yp2].rho;
			fm3=node[YDIM*xm3+ym3].S*node[YDIM*xm3+ym3].rho;
			fp3=node[YDIM*xp3+yp3].S*node[YDIM*xp3+yp3].rho;
			break;
		case JX:
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[0];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].u[0];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[0];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].u[0];
			fm3=node[YDIM*xm3+ym3].rho*node[YDIM*xm3+ym3].u[0];
			fp3=node[YDIM*xp3+yp3].rho*node[YDIM*xp3+yp3].u[0];
			break; 
		case JY:
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[1];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].u[1];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[1];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].u[1];
			fm3=node[YDIM*xm3+ym3].rho*node[YDIM*xm3+ym3].u[1];
			fp3=node[YDIM*xp3+yp3].rho*node[YDIM*xp3+yp3].u[1];
			break; 
		case RHOUXUX:
			a=0; b=0;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			fm3=node[YDIM*xm3+ym3].rho*node[YDIM*xm3+ym3].u[a]*node[YDIM*xm3+ym3].u[b];
			fp3=node[YDIM*xp3+yp3].rho*node[YDIM*xp3+yp3].u[a]*node[YDIM*xp3+yp3].u[b];
			break;
		case RHOUYUY:
			a=1; b=1;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			fm3=node[YDIM*xm3+ym3].rho*node[YDIM*xm3+ym3].u[a]*node[YDIM*xm3+ym3].u[b];
			fp3=node[YDIM*xp3+yp3].rho*node[YDIM*xp3+yp3].u[a]*node[YDIM*xp3+yp3].u[b];
			break;
		case RHOUXUY:
			a=0; b=1;
			fm1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[a]*node[YDIM*xm1+ym1].u[b];
			fm2=node[YDIM*xm2+ym2].rho*node[YDIM*xm2+ym2].u[a]*node[YDIM*xm2+ym2].u[b];
			fp1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[a]*node[YDIM*xp1+yp1].u[b];
			fp2=node[YDIM*xp2+yp2].rho*node[YDIM*xp2+yp2].u[a]*node[YDIM*xp2+yp2].u[b];
			fm3=node[YDIM*xm3+ym3].rho*node[YDIM*xm3+ym3].u[a]*node[YDIM*xm3+ym3].u[b];
			fp3=node[YDIM*xp3+yp3].rho*node[YDIM*xp3+yp3].u[a]*node[YDIM*xp3+yp3].u[b];
			break;
	}
	der_f= fp3 - fm3 + 9.0* fm2 -9.0 * fm2 + 45.0 * fp1 -45.0* fm1;
	der_f/=60.0;
return der_f;
}


/*********************************************************/
/**Second - Laplacian derivative of f at x,y*/
double der22(DerivativeOf what, int x, int y){
	double der_f;
	double f00,fm10,fp10,f0p1,f0m1; //fm1=function at position minus 1, fp2=function at position plus 2
	int xm1, xp1, ym1, yp1;
	int a, b;
	
	xm1=x-1;
	xp1=x+1;
	ym1=y-1;
	yp1=y+1;

	if(xm1>(XDIM-1))xm1-=XDIM;
	else if(xm1<0)xm1+=XDIM;
	if(xp1>(XDIM-1))xp1-=XDIM;
	else if(xp1<0)xp1+=XDIM;

	if(ym1>(YDIM-1))ym1-=YDIM;
	else if(ym1<0)ym1+=YDIM;
	if(yp1>(YDIM-1))yp1-=YDIM;
	else if(yp1<0)yp1+=YDIM;

	switch (what){
		case DENSITY:
			f00=node[YDIM*x+y].rho;
			fm10=node[YDIM*xm1+y].rho;
			fp10=node[YDIM*xp1+y].rho;
			f0p1=node[YDIM*x+yp1].rho;
			f0m1=node[YDIM*x+ym1].rho;
			break;
		case UX:
			f00=node[YDIM*x+y].u[0];
			fm10=node[YDIM*xm1+y].u[0];
			fp10=node[YDIM*xp1+y].u[0];
			f0p1=node[YDIM*x+yp1].u[0];
			f0m1=node[YDIM*x+ym1].u[0];
			break; 
		case UY:
			f00=node[YDIM*x+y].u[1];
			fm10=node[YDIM*xm1+y].u[1];
			fp10=node[YDIM*xp1+y].u[1];
			f0p1=node[YDIM*x+yp1].u[1];
			f0m1=node[YDIM*x+ym1].u[1];
			break; 
		case PRESSURE:
			f00=node[YDIM*x+y].rho*node[YDIM*x+y].S;
			fm10=node[YDIM*xm1+y].rho*node[YDIM*xm1+y].S;
			fp10=node[YDIM*xp1+y].rho*node[YDIM*xp1+y].S;
			f0p1=node[YDIM*x+yp1].rho*node[YDIM*x+yp1].S;
			f0m1=node[YDIM*x+ym1].rho*node[YDIM*x+ym1].S;
			break;
		case RHOUXUX:
			f00=node[YDIM*x+y].rho*node[YDIM*x+y].u[0]*node[YDIM*x+y].u[0];
			fm10=node[YDIM*xm1+y].rho*node[YDIM*xm1+y].u[0]*node[YDIM*xm1+y].u[0];
			fp10=node[YDIM*xp1+y].rho*node[YDIM*xp1+y].u[0]*node[YDIM*xp1+y].u[0];
			f0p1=node[YDIM*x+yp1].rho*node[YDIM*x+yp1].u[0]*node[YDIM*x+yp1].u[0];
			f0m1=node[YDIM*x+ym1].rho*node[YDIM*x+ym1].u[0]*node[YDIM*x+ym1].u[0];
			break;
		case RHOUXUY:
			f00=node[YDIM*x+y].rho*node[YDIM*x+y].u[0]*node[YDIM*x+y].u[1];
			fm10=node[YDIM*xm1+y].rho*node[YDIM*xm1+y].u[0]*node[YDIM*xm1+y].u[1];
			fp10=node[YDIM*xp1+y].rho*node[YDIM*xp1+y].u[0]*node[YDIM*xp1+y].u[1];
			f0p1=node[YDIM*x+yp1].rho*node[YDIM*x+yp1].u[0]*node[YDIM*x+yp1].u[1];
			f0m1=node[YDIM*x+ym1].rho*node[YDIM*x+ym1].u[0]*node[YDIM*x+ym1].u[1];
			break;
		case RHOUYUY:
			f00=node[YDIM*x+y].rho*node[YDIM*x+y].u[1]*node[YDIM*x+y].u[1];
			fm10=node[YDIM*xm1+y].rho*node[YDIM*xm1+y].u[1]*node[YDIM*xm1+y].u[1];
			fp10=node[YDIM*xp1+y].rho*node[YDIM*xp1+y].u[1]*node[YDIM*xp1+y].u[1];
			f0p1=node[YDIM*x+yp1].rho*node[YDIM*x+yp1].u[1]*node[YDIM*x+yp1].u[1];
			f0m1=node[YDIM*x+ym1].rho*node[YDIM*x+ym1].u[1]*node[YDIM*x+ym1].u[1];
			break;
		case JX:
			f00=node[YDIM*x+y].rho*node[YDIM*x+y].u[0];
			fm10=node[YDIM*xm1+y].rho*node[YDIM*xm1+y].u[0];
			fp10=node[YDIM*xp1+y].rho*node[YDIM*xp1+y].u[0];
			f0p1=node[YDIM*x+yp1].rho*node[YDIM*x+yp1].u[0];
			f0m1=node[YDIM*x+ym1].rho*node[YDIM*x+ym1].u[0];
			break; 
		case JY:
			f00=node[YDIM*x+y].rho*node[YDIM*x+y].u[1];
			fm10=node[YDIM*xm1+y].rho*node[YDIM*xm1+y].u[1];
			fp10=node[YDIM*xp1+y].rho*node[YDIM*xp1+y].u[1];
			f0p1=node[YDIM*x+yp1].rho*node[YDIM*x+yp1].u[1];
			f0m1=node[YDIM*x+ym1].rho*node[YDIM*x+ym1].u[1];
			break; 

	}

	der_f=-4.0*f00 + fp10 + fm10 + f0m1 + f0p1;
return der_f;
}

/*********************************************************/
/**Second derivative of f at x,y*/
double der12(DerivativeOf what, int x, int y, int dir1, int dir2){
	double der_f=0.0;
	double fm1m1,fm1p1,fp1m1,fp1p1, fp10, fm10, f0p1, f0m1, fp20, fm20, f0p2, f0m2 ; //fm1=function at position minus 1, fp2=function at position plus 2
	int xm1, xp1, ym1, yp1, xm2, xp2, ym2, yp2;
	int a, b;
	
	fp10=0.0;
	fm10=0.0;
	f0p1=0.0;
	f0m1=0.0;
	fp20=0.0;
	fm20=0.0;
	f0p2=0.0;
	f0m2=0.0;
	
	//if(dir1 == dir2){
	//	der_f=0.5*der22(what, x, y);
	//}
	//else if(dir1 != dir2)	
	//{

	xm1=x-1;
	xm2=x-2;
	xp1=x+1;
	xp2=x+2;

	if(xm1>(XDIM-1))xm1-=XDIM;
	else if(xm1<0)xm1+=XDIM;
	if(xm2>(XDIM-1))xm2-=XDIM;
	else if(xm2<0)xm2+=XDIM;
	if(xp1>(XDIM-1))xp1-=XDIM;
	else if(xp1<0)xp1+=XDIM;
	if(xp2>(XDIM-1))xp2-=XDIM;
	else if(xp2<0)xp2+=XDIM;

	ym1=y-1;
	ym2=y-2;
	yp1=y+1;
	yp2=y+2;

	if(ym1>(YDIM-1))ym1-=YDIM;
	else if(ym1<0)ym1+=YDIM;
	if(ym2>(YDIM-1))ym2-=YDIM;
	else if(ym2<0)ym2+=YDIM;
	if(yp1>(YDIM-1))yp1-=YDIM;
	else if(yp1<0)yp1+=YDIM;
	if(yp2>(YDIM-1))yp2-=YDIM;
	else if(yp2<0)yp2+=YDIM;


	switch (what){
		case DENSITY:
			fm1m1=node[YDIM*xm1+ym1].rho;
			fm1p1=node[YDIM*xm1+yp1].rho;
			fp1p1=node[YDIM*xp1+yp1].rho;
			fp1m1=node[YDIM*xp1+ym1].rho;
			if(dir1==dir2){
				fp10=node[YDIM*xp1+y].rho;
				fm10=node[YDIM*xm1+y].rho;
				f0p1=node[YDIM*x+yp1].rho;
				f0m1=node[YDIM*x+ym1].rho;
			
				fp20=node[YDIM*xp2+y].rho;
				fm20=node[YDIM*xm2+y].rho;
				f0p2=node[YDIM*x+yp2].rho;
				f0m2=node[YDIM*x+ym2].rho;
			}
			break;
		case UX:
			fm1m1=node[YDIM*xm1+ym1].u[0];
			fm1p1=node[YDIM*xm1+yp1].u[0];
			fp1p1=node[YDIM*xp1+yp1].u[0];
			fp1m1=node[YDIM*xp1+ym1].u[0];
			if(dir1==dir2){
				fp10=node[YDIM*xp1+y].u[0];
				fm10=node[YDIM*xm1+y].u[0];
				f0p1=node[YDIM*x+yp1].u[0];
				f0m1=node[YDIM*x+ym1].u[0];
			
				fp20=node[YDIM*xp2+y].u[0];
				fm20=node[YDIM*xm2+y].u[0];
				f0p2=node[YDIM*x+yp2].u[0];
				f0m2=node[YDIM*x+ym2].u[0];
			}
			break; 
		case UY:
			fm1m1=node[YDIM*xm1+ym1].u[1];
			fm1p1=node[YDIM*xm1+yp1].u[1];
			fp1p1=node[YDIM*xp1+yp1].u[1];
			fp1m1=node[YDIM*xp1+ym1].u[1];
			if(dir1==dir2){
				fp10=node[YDIM*xp1+y].u[1];
				fm10=node[YDIM*xm1+y].u[1];
				f0p1=node[YDIM*x+yp1].u[1];
				f0m1=node[YDIM*x+ym1].u[1];
			
				fp20=node[YDIM*xp2+y].u[1];
				fm20=node[YDIM*xm2+y].u[1];
				f0p2=node[YDIM*x+yp2].u[1];
				f0m2=node[YDIM*x+ym2].u[1];
			}
			break; 
		case JX:
			fm1m1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[0];
			fm1p1=node[YDIM*xm1+yp1].rho*node[YDIM*xm1+yp1].u[0];
			fp1p1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[0];
			fp1m1=node[YDIM*xp1+ym1].rho*node[YDIM*xp1+ym1].u[0];
			if(dir1==dir2){
				fp10=node[YDIM*xp1+y].rho*node[YDIM*xp1+y].u[0];
				fm10=node[YDIM*xm1+y].rho*node[YDIM*xm1+y].u[0];
				f0p1=node[YDIM*x+yp1].rho*node[YDIM*x+yp1].u[0];
				f0m1=node[YDIM*x+ym1].rho*node[YDIM*x+ym1].u[0];
			
				fp20=node[YDIM*xp2+y].rho*node[YDIM*xp2+y].u[0];
				fm20=node[YDIM*xm2+y].rho*node[YDIM*xm2+y].u[0];
				f0p2=node[YDIM*x+yp2].rho*node[YDIM*x+yp2].u[0];
				f0m2=node[YDIM*x+ym2].rho*node[YDIM*x+ym2].u[0];
			}
			break; 
		case JY:
			fm1m1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].u[1];
			fm1p1=node[YDIM*xm1+yp1].rho*node[YDIM*xm1+yp1].u[1];
			fp1p1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].u[1];
			fp1m1=node[YDIM*xp1+ym1].rho*node[YDIM*xp1+ym1].u[1];
			if(dir1==dir2){
				fp10=node[YDIM*xp1+y].rho*node[YDIM*xp1+y].u[1];
				fm10=node[YDIM*xm1+y].rho*node[YDIM*xm1+y].u[1];
				f0p1=node[YDIM*x+yp1].rho*node[YDIM*x+yp1].u[1];
				f0m1=node[YDIM*x+ym1].rho*node[YDIM*x+ym1].u[1];
			
				fp20=node[YDIM*xp2+y].rho*node[YDIM*xp2+y].u[1];
				fm20=node[YDIM*xm2+y].rho*node[YDIM*xm2+y].u[1];
				f0p2=node[YDIM*x+yp2].rho*node[YDIM*x+yp2].u[1];
				f0m2=node[YDIM*x+ym2].rho*node[YDIM*x+ym2].u[1];
			}
			break; 

		case PRESSURE:
			fm1m1=node[YDIM*xm1+ym1].rho*node[YDIM*xm1+ym1].S;
			fm1p1=node[YDIM*xm1+yp1].rho*node[YDIM*xm1+yp1].S;
			fp1p1=node[YDIM*xp1+yp1].rho*node[YDIM*xp1+yp1].S;
			fp1m1=node[YDIM*xp1+ym1].rho*node[YDIM*xp1+ym1].S;
			if(dir1==dir2){
				fp10=node[YDIM*xp1+y].rho*node[YDIM*xp1+y].S;
				fm10=node[YDIM*xm1+y].rho*node[YDIM*xm1+y].S;
				f0p1=node[YDIM*x+yp1].rho*node[YDIM*x+yp1].S;
				f0m1=node[YDIM*x+ym1].rho*node[YDIM*x+ym1].S;
			
				fp20=node[YDIM*xp2+y].rho*node[YDIM*xp2+y].S;
				fm20=node[YDIM*xm2+y].rho*node[YDIM*xm2+y].S;
				f0p2=node[YDIM*x+yp2].rho*node[YDIM*x+yp2].S;
				f0m2=node[YDIM*x+ym2].rho*node[YDIM*x+ym2].S;
			}
			break;
		case RHOUXUX:
			fm1m1=node[YDIM*xm1+ym1].u[0]*node[YDIM*xm1+ym1].u[0]*node[YDIM*xm1+ym1].rho;
			fm1p1=node[YDIM*xm1+yp1].u[0]*node[YDIM*xm1+yp1].u[0]*node[YDIM*xm1+yp1].rho;
			fp1p1=node[YDIM*xp1+yp1].u[0]*node[YDIM*xp1+yp1].u[0]*node[YDIM*xp1+yp1].rho;
			fp1m1=node[YDIM*xp1+ym1].u[0]*node[YDIM*xp1+ym1].u[0]*node[YDIM*xp1+ym1].rho;
			if(dir1==dir2){
				fp10=node[YDIM*xp1+y].u[0]*node[YDIM*xp1+y].u[0]*node[YDIM*xp1+y].rho;
				fm10=node[YDIM*xm1+y].u[0]*node[YDIM*xm1+y].u[0]*node[YDIM*xm1+y].rho;
				f0p1=node[YDIM*x+yp1].u[0]*node[YDIM*x+yp1].u[0]*node[YDIM*x+yp1].rho;
				f0m1=node[YDIM*x+ym1].u[0]*node[YDIM*x+ym1].u[0]*node[YDIM*x+ym1].rho;
			
				fp20=node[YDIM*xp2+y].u[0]*node[YDIM*xp2+y].u[0]*node[YDIM*xp2+y].rho;
				fm20=node[YDIM*xm2+y].u[0]*node[YDIM*xm2+y].u[0]*node[YDIM*xm2+y].rho;
				f0p2=node[YDIM*x+yp2].u[0]*node[YDIM*x+yp2].u[0]*node[YDIM*x+yp2].rho;
				f0m2=node[YDIM*x+ym2].u[0]*node[YDIM*x+ym2].u[0]*node[YDIM*x+ym2].rho;
			}
			break;
		case RHOUXUY:
			fm1m1=node[YDIM*xm1+ym1].u[0]*node[YDIM*xm1+ym1].u[1]*node[YDIM*xm1+ym1].rho;
			fm1p1=node[YDIM*xm1+yp1].u[0]*node[YDIM*xm1+yp1].u[1]*node[YDIM*xm1+yp1].rho;
			fp1p1=node[YDIM*xp1+yp1].u[0]*node[YDIM*xp1+yp1].u[1]*node[YDIM*xp1+yp1].rho;
			fp1m1=node[YDIM*xp1+ym1].u[0]*node[YDIM*xp1+ym1].u[1]*node[YDIM*xp1+ym1].rho;
			if(dir1==dir2){
				fp10=node[YDIM*xp1+y].u[0]*node[YDIM*xp1+y].u[1]*node[YDIM*xp1+y].rho;
				fm10=node[YDIM*xm1+y].u[0]*node[YDIM*xm1+y].u[1]*node[YDIM*xm1+y].rho;
				f0p1=node[YDIM*x+yp1].u[0]*node[YDIM*x+yp1].u[1]*node[YDIM*x+yp1].rho;
				f0m1=node[YDIM*x+ym1].u[0]*node[YDIM*x+ym1].u[1]*node[YDIM*x+ym1].rho;
			
				fp20=node[YDIM*xp2+y].u[0]*node[YDIM*xp2+y].u[1]*node[YDIM*xp2+y].rho;
				fm20=node[YDIM*xm2+y].u[0]*node[YDIM*xm2+y].u[1]*node[YDIM*xm2+y].rho;
				f0p2=node[YDIM*x+yp2].u[0]*node[YDIM*x+yp2].u[1]*node[YDIM*x+yp2].rho;
				f0m2=node[YDIM*x+ym2].u[0]*node[YDIM*x+ym2].u[1]*node[YDIM*x+ym2].rho;
			}
			break; 
		case RHOUYUY:
			fm1m1=node[YDIM*xm1+ym1].u[1]*node[YDIM*xm1+ym1].u[1]*node[YDIM*xm1+ym1].rho;
			fm1p1=node[YDIM*xm1+yp1].u[1]*node[YDIM*xm1+yp1].u[1]*node[YDIM*xm1+yp1].rho;
			fp1p1=node[YDIM*xp1+yp1].u[1]*node[YDIM*xp1+yp1].u[1]*node[YDIM*xp1+yp1].rho;
			fp1m1=node[YDIM*xp1+ym1].u[1]*node[YDIM*xp1+ym1].u[1]*node[YDIM*xp1+ym1].rho;
			if(dir1==dir2){
				fp10=node[YDIM*xp1+y].u[1]*node[YDIM*xp1+y].u[1]*node[YDIM*xp1+y].rho;
				fm10=node[YDIM*xm1+y].u[1]*node[YDIM*xm1+y].u[1]*node[YDIM*xm1+y].rho;
				f0p1=node[YDIM*x+yp1].u[1]*node[YDIM*x+yp1].u[1]*node[YDIM*x+yp1].rho;
				f0m1=node[YDIM*x+ym1].u[1]*node[YDIM*x+ym1].u[1]*node[YDIM*x+ym1].rho;
			
				fp20=node[YDIM*xp2+y].u[1]*node[YDIM*xp2+y].u[1]*node[YDIM*xp2+y].rho;
				fm20=node[YDIM*xm2+y].u[1]*node[YDIM*xm2+y].u[1]*node[YDIM*xm2+y].rho;
				f0p2=node[YDIM*x+yp2].u[1]*node[YDIM*x+yp2].u[1]*node[YDIM*x+yp2].rho;
				f0m2=node[YDIM*x+ym2].u[1]*node[YDIM*x+ym2].u[1]*node[YDIM*x+ym2].rho;
			}
			break; 
	}

	if(dir1==0 && dir2==0)
		der_f= (fp1p1 + fm1m1 + fp1m1 + fm1p1)/4.0 - (fp10 + fm10) + 0.5*(fp20 + fm20) - 0.5 * der22(what, x, y);
	
	else if(dir1==1 && dir2 ==1)
		der_f= (fp1p1 + fm1m1 + fp1m1 + fm1p1)/4.0 - (f0p1 + f0m1) + 0.5*(f0p2 + f0m2) - 0.5 * der22(what, x, y);
	
	else{	
	der_f= (fp1p1 + fm1m1 - fp1m1 - fm1p1)/4.0;
	}
return der_f;
}


