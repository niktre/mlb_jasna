/*Streaming and random number generator functions*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defs.h"

/**********************************************************/
/**Streaming step - only 2D D2Q21 model*/
void stream(){
	int x,y,i;

	for (y = 0; y < YDIM; y++) {
		for (x = 0; x < XDIM; x++) {
			int xnew = 0;
			int ynew = 0;
			
			for (i = 0; i < VDIM; i++){
				xnew = x + C[2 * i];
				ynew = y + C[2 * i + 1];
				
				/* boundary conditions */
				if (xnew >= XDIM) xnew -= XDIM;
				if (xnew < 0) xnew += XDIM;
				if (ynew >= YDIM) ynew -= YDIM;
				if (ynew < 0) ynew += YDIM;
				
				new_node[YDIM*xnew+ynew].pop[i]=node[YDIM*x+y].pop[i];
			}
		}
	}
	
	/* swapping */
	for (y = 0; y < YDIM; y++) {
		for (x = 0; x < XDIM; x++) {
			for (i = 0; i < VDIM; i++){
				temp_node[YDIM*x+y].pop[i] = node[YDIM*x+y].pop[i];
				node[YDIM*x+y].pop[i] = new_node[YDIM*x+y].pop[i];
				new_node[YDIM*x+y].pop[i] = temp_node[YDIM*x+y].pop[i];
			}
		}
	}
	
	/* Nikita: I can believe that the following works and so on, but I really do not want to check it. Really
	 
	//    extern Velocity vel[VDIM];
	//	extern node[]
	//	double y_vel_1[YDIM], y_vel_2[YDIM], y_vel_5[YDIM], y_vel_6[YDIM], y_vel_7[YDIM], y_vel_8[YDIM],  y_vel_9[2*YDIM], y_vel_10[2*YDIM], y_vel_13[2*YDIM], y_vel_14[2*YDIM], y_vel_15[2*YDIM], y_vel_16[2*YDIM], y_vel_17[4*YDIM], y_vel_18[4*YDIM];
	//	double x_vel_3[XDIM], x_vel_4[XDIM], x_vel_5[XDIM], x_vel_6[XDIM], x_vel_7[XDIM], x_vel_8[XDIM], x_vel_11[2*XDIM], x_vel_12[2*XDIM], x_vel_13[2*XDIM], x_vel_14[2*XDIM], x_vel_15[2*XDIM], x_vel_16[2*XDIM], x_vel_19[4*XDIM], x_vel_20[4*XDIM];
	
	/*	printf("\n");
	 :	for(y=YDIM-1;y>-1;y--){printf("\n");
	 for(x=0;x<XDIM;x++){
	 vel[12].n[YDIM*x+y]=YDIM*x+y;
	 printf("%f    ",vel[12].n[YDIM*x+y]);}}
	 printf("\n\n");*/
	
	/*Using periodic boundary conditions. One lattice only, therefore the values on the boundaries have to be saved first*/
	/*  for(y=YDIM-1;y>-1;y--){printf("\n");
  for(x=0;x<XDIM;x++)
	 printf("%f    ",vel[12].n[YDIM*x+y]);}
	 */
	
//	&temp_node = &node;
//	&node = &new_node;
//	&new_node = &temp_node;
/*
	for(y=0;y<YDIM;y++){
		y_vel_1[y]=node[YDIM*(XDIM-1)+y].pop[1];// vel[1].n[YDIM*(XDIM-1)+y];
		y_vel_2[y]=node[y].pop[2];

		y_vel_9[y]=node[YDIM*(XDIM-1)+y].pop[9];
		y_vel_9[YDIM+y]=node[YDIM*(XDIM-2)+y].pop[9];

		y_vel_10[y]=node[y].pop[10];
		y_vel_10[YDIM+y]=node[YDIM+y].pop[10];

		y_vel_17[y]=node[YDIM*(XDIM-1)+y].pop[17];
		y_vel_17[YDIM+y]=node[YDIM*(XDIM-2)+y].pop[17];
		y_vel_17[YDIM*2+y]=node[YDIM*(XDIM-3)+y].pop[17];
		y_vel_17[YDIM*3+y]=node[YDIM*(XDIM-4)+y].pop[17];

		y_vel_18[y]=node[y].pop[18];
		y_vel_18[YDIM+y]=node[YDIM+y].pop[18];
		y_vel_18[YDIM*2+y]=node[YDIM*2+y].pop[18];
		y_vel_18[YDIM*3+y]=node[YDIM*3+y].pop[18];

		y_vel_8[y]=node[YDIM*(XDIM-1)+y].pop[8];
		y_vel_16[y]=node[YDIM*(XDIM-1)+y].pop[16];
		y_vel_16[YDIM+y]=node[YDIM*(XDIM-2)+y].pop[16];
	
		y_vel_7[y]=node[y].pop[7];
		y_vel_15[y]=node[y].pop[15];
		y_vel_15[YDIM+y]=node[YDIM+y].pop[15];

		y_vel_6[y]=node[y].pop[6];
		y_vel_14[y]=node[y].pop[14];
		y_vel_14[YDIM+y]=node[YDIM+y].pop[14];
	
		y_vel_5[y]=node[YDIM*(XDIM-1)+y].pop[5];
		y_vel_13[y]=node[YDIM*(XDIM-1)+y].pop[13];
		y_vel_13[YDIM+y]=node[YDIM*(XDIM-2)+y].pop[13];
	}
	for(x=0;x<XDIM;x++){
		x_vel_4[x]=vel[4].n[YDIM*x];
		x_vel_3[x]=vel[3].n[YDIM*x+YDIM-1];

		x_vel_12[x]=vel[12].n[YDIM*x];
		x_vel_12[XDIM+x]=vel[12].n[YDIM*x+1];
		x_vel_11[x]=vel[11].n[YDIM*x+YDIM-1];
		x_vel_11[XDIM+x]=vel[11].n[YDIM*x+YDIM-2];

		x_vel_20[x]=vel[20].n[YDIM*x];
		x_vel_20[XDIM+x]=vel[20].n[YDIM*x+1];
		x_vel_20[XDIM*2+x]=vel[20].n[YDIM*x+2];
		x_vel_20[XDIM*3+x]=vel[20].n[YDIM*x+3];

		x_vel_19[x]=vel[19].n[YDIM*x+YDIM-1];
		x_vel_19[XDIM+x]=vel[19].n[YDIM*x+YDIM-2];
		x_vel_19[XDIM*2+x]=vel[19].n[YDIM*x+YDIM-3];
		x_vel_19[XDIM*3+x]=vel[19].n[YDIM*x+YDIM-4];


		x_vel_8[x]=vel[8].n[YDIM*x];
		x_vel_16[x]=vel[16].n[YDIM*x];
		x_vel_16[XDIM+x]=vel[16].n[YDIM*x+1];

		x_vel_7[x]=vel[7].n[YDIM*x+YDIM-1];
		x_vel_15[x]=vel[15].n[YDIM*x+YDIM-1];
		x_vel_15[XDIM+x]=vel[15].n[YDIM*x+YDIM-2];

		x_vel_6[x]=vel[6].n[YDIM*x];
		x_vel_14[x]=vel[14].n[YDIM*x];
		x_vel_14[XDIM+x]=vel[14].n[YDIM*x+1];

		x_vel_5[x]=vel[5].n[YDIM*x+YDIM-1];
		x_vel_13[x]=vel[13].n[YDIM*x+YDIM-1];
		x_vel_13[XDIM+x]=vel[13].n[YDIM*x+YDIM-2];
	}
*/

	/*Moving the density distributions according to the velocity vectors*/
	/*direction (1,0)*/ 
//	memmove(&vel[1].n[YDIM],&vel[1].n[0],(XDIM-1)*YDIM*sizeof(double));
  /*direction (2,0)*/
//  memmove(&vel[9].n[YDIM*2],&vel[9].n[0],(XDIM-2)*YDIM*sizeof(double));
  /*direction (4,0)*/
//  memmove(&vel[17].n[YDIM*4],&vel[17].n[0],(XDIM-4)*YDIM*sizeof(double));


	/*direction (-1,0)*/
//  memmove(&vel[2].n[0],&vel[2].n[YDIM],(XDIM-1)*YDIM*sizeof(double));
  /*direction (-2,0)*/
//  memmove(&vel[10].n[0],&vel[10].n[YDIM*2],(XDIM-2)*YDIM*sizeof(double));
  /*direction (-4,0)*/
//  memmove(&vel[18].n[0],&vel[18].n[YDIM*4],(XDIM-4)*YDIM*sizeof(double));

  /*direction (0,-1)*/
//  memmove(&vel[4].n[0],&vel[4].n[1],(XDIM*YDIM-1)*sizeof(double));
  /*direction (0,-2)*/
 // memmove(&vel[12].n[0],&vel[12].n[2],(XDIM*YDIM-2)*sizeof(double));
  /*direction (0,-4)*/
//  memmove(&vel[20].n[0],&vel[20].n[4],(XDIM*YDIM-4)*sizeof(double));

  /*direction (0,1)*/
//  memmove(&vel[3].n[1],&vel[3].n[0],(XDIM*YDIM-1)*sizeof(double));
  /*direction (0,2)*/
//  memmove(&vel[11].n[2],&vel[11].n[0],(XDIM*YDIM-2)*sizeof(double));
  /*direction (0,4)*/
//  memmove(&vel[19].n[4],&vel[19].n[0],(XDIM*YDIM-4)*sizeof(double));

  /*direction (1,-1)*/
//  memmove(&vel[8].n[YDIM],&vel[8].n[1],((XDIM-1)*YDIM-1)*sizeof(double));
  /*direction (2,-2)*/
//  memmove(&vel[16].n[YDIM*2],&vel[16].n[2],((XDIM-2)*YDIM-2)*sizeof(double));

  /*direction (-1,1)*/
//  memmove(&vel[7].n[1],&vel[7].n[YDIM],((XDIM-1)*YDIM-1)*sizeof(double));
  /*direction (-2,2)*/
 // memmove(&vel[15].n[2],&vel[15].n[YDIM*2],((XDIM-2)*YDIM-2)*sizeof(double));

  /*direction (-1,-1)*/
//  memmove(&vel[6].n[0],&vel[6].n[YDIM+1],((XDIM-1)*YDIM-1)*sizeof(double));
  /*direction (-2,-2)*/
//  memmove(&vel[14].n[0],&vel[14].n[2*YDIM+2],((XDIM-2)*YDIM-2)*sizeof(double));

  /*direction (1,1)*/
//  memmove(&vel[5].n[YDIM+1],&vel[5].n[0],((XDIM-1)*YDIM-1)*sizeof(double));
  /*direction (2,2)*/
//  memmove(&vel[13].n[2*YDIM+2],&vel[13].n[0],((XDIM-2)*YDIM-2)*sizeof(double));

	
	/*Fixing the boundary values*/
/*	for(y=0;y<YDIM;y++){
		vel[1].n[y]=y_vel_1[y];
		vel[2].n[YDIM*(XDIM-1)+y]=y_vel_2[y];
		
		vel[9].n[y]=y_vel_9[YDIM+y];
		vel[9].n[YDIM+y]=y_vel_9[y];
		vel[10].n[YDIM*(XDIM-1)+y]=y_vel_10[YDIM+y];
		vel[10].n[YDIM*(XDIM-2)+y]=y_vel_10[y];
		
		vel[17].n[y]=y_vel_17[YDIM*3+y];
		vel[17].n[YDIM+y]=y_vel_17[YDIM*2+y];
		vel[17].n[YDIM*2+y]=y_vel_17[YDIM+y];
		vel[17].n[YDIM*3+y]=y_vel_17[y];

		vel[18].n[YDIM*(XDIM-1)+y]=y_vel_18[YDIM*3+y];
		vel[18].n[YDIM*(XDIM-2)+y]=y_vel_18[YDIM*2+y];
		vel[18].n[YDIM*(XDIM-3)+y]=y_vel_18[YDIM+y];
		vel[18].n[YDIM*(XDIM-4)+y]=y_vel_18[y];

		vel[8].n[y]=y_vel_8[(y+1)%YDIM];
		vel[16].n[y]=y_vel_16[YDIM+(y+2)%YDIM];
		vel[16].n[YDIM+y]=y_vel_16[(y+2)%YDIM];

		vel[7].n[YDIM*(XDIM-1)+(y+1)%YDIM]=y_vel_7[y];
		vel[15].n[YDIM*(XDIM-1)+(y+2)%YDIM]=y_vel_15[YDIM+y];
		vel[15].n[YDIM*(XDIM-2)+(y+2)%YDIM]=y_vel_15[y];

		vel[6].n[YDIM*(XDIM-1)+y]=y_vel_6[(y+1)%YDIM];
		vel[14].n[YDIM*(XDIM-1)+y]=y_vel_14[YDIM+(y+2)%YDIM];
		vel[14].n[YDIM*(XDIM-2)+y]=y_vel_14[(y+2)%YDIM];

		vel[5].n[(y+1)%(YDIM)]=y_vel_5[y];
		vel[13].n[(y+2)%(YDIM)]=y_vel_13[YDIM+y];
		vel[13].n[YDIM+(y+2)%(YDIM)]=y_vel_13[y];

	}
	for(x=0;x<XDIM;x++){
		vel[4].n[YDIM*x+YDIM-1]=x_vel_4[x];
		vel[3].n[YDIM*x]=x_vel_3[x];

		vel[12].n[YDIM*x+YDIM-1]=x_vel_12[XDIM+x];
		vel[12].n[YDIM*x+YDIM-2]=x_vel_12[x];
		vel[11].n[YDIM*x]=x_vel_11[XDIM+x];
		vel[11].n[YDIM*x+1]=x_vel_11[x];

		vel[20].n[YDIM*x+YDIM-1]=x_vel_20[XDIM*3+x];
		vel[20].n[YDIM*x+YDIM-2]=x_vel_20[XDIM*2+x];
		vel[20].n[YDIM*x+YDIM-3]=x_vel_20[XDIM+x];
		vel[20].n[YDIM*x+YDIM-4]=x_vel_20[x];

		vel[19].n[YDIM*x]=x_vel_19[XDIM*3+x];
		vel[19].n[YDIM*x+1]=x_vel_19[XDIM*2+x];
		vel[19].n[YDIM*x+2]=x_vel_19[XDIM+x];
		vel[19].n[YDIM*x+3]=x_vel_19[x];

		vel[8].n[YDIM*((x+1)%XDIM)+YDIM-1]=x_vel_8[x];
		vel[16].n[YDIM*((x+2)%XDIM)+YDIM-2]=x_vel_16[x];
		vel[16].n[YDIM*((x+2)%XDIM)+YDIM-1]=x_vel_16[XDIM+x];

		vel[7].n[YDIM*x]=x_vel_7[(x+1)%XDIM];
		vel[15].n[YDIM*x]=x_vel_15[XDIM+(x+2)%XDIM];
		vel[15].n[YDIM*x+1]=x_vel_15[(x+2)%XDIM];

		vel[6].n[YDIM*x+YDIM-1]=x_vel_6[(x+1)%XDIM];
		vel[14].n[YDIM*x+YDIM-1]=x_vel_14[XDIM+(x+2)%XDIM];
		vel[14].n[YDIM*x+YDIM-2]=x_vel_14[(x+2)%XDIM];
	
		vel[5].n[YDIM*((x+1)%XDIM)]=x_vel_5[x];
		vel[13].n[YDIM*((x+2)%XDIM)]=x_vel_13[XDIM+x];
		vel[13].n[YDIM*((x+2)%XDIM)+1]=x_vel_13[x];
	}
*/	/*for(y=YDIM-1;y>-1;y--){printf("\n");
	for(x=0;x<XDIM;x++)
printf("%f    ",vel[12].n[YDIM*x+y]);}*/
return;
}


/*********************************************************/
/*Random number [0,1] generator*/

double ran2(long *idum)
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

{
  int j;
  long k;
  static long idum2=243456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;
  //for(i=1;i<=3;i++)
    //{

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
  //}
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/***********************************************************************/
/* Added by uschille */
/***********************************************************************/

/** Generate Gaussian random numbers via Box-Muller transformation */
double gaussian_random() {
  static double x1, x2, r2, fac;
  static int calc_new = 1;

  /* On every second call two gaussian random numbers are calculated
     via the Box-Muller transformation. One is returned as the result
     and the second one is stored for use on the next call. */

  /* Note: Use variable seeds for production runs! */

  if (calc_new) {

    calc_new = 0;

    /* draw two uniform random numbers in the unit circle */
    do {
      x1 = 2.0*rand()/(double)RAND_MAX-1.0;
      x2 = 2.0*rand()/(double)RAND_MAX-1.0;
      r2 = x1*x1 + x2*x2;
    } while (r2 >= 1.0);

    /* perform Box-Muller transformation */
    fac = sqrt(-2.0*log(r2)/r2);

    /* return the first gaussian random number */
    return x1*fac;

  } else {

    calc_new = 1;

    /* return the stored gaussian random number */
    return x2*fac;

  }

}

/*********************************************************/
