#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defs.h"

double cs2_eq;

double tau[TDIM]={TAU1, TAU1, TAU1, TAU1,TAU2, TAU2, TAU2, TAU2, TAU3, TAU3, TAU3, TAU3};

//Velocity vel[VDIM];

int D[DIM*TDIM]={
		 1,  0,
		-1,  0,
		 0,  1,
		 0, -1,
		 1,  1,
		-1, -1,
		-1,  1,
		 1, -1,
		 2,  0,
		-2,  0,
		 0,  2,
		 0, -2,
		};

double C[DIM*VDIM]={
		 0.,  0.,
		 1.,  0.,
		-1.,  0.,
		 0.,  1.,
		 0., -1.,
		 1.,  1.,
		-1., -1.,
		-1.,  1.,
		 1., -1.,
		 2.,  0.,
		-2.,  0.,
		 0.,  2.,
		 0., -2.,
		 2.,  2.,
		-2., -2.,
		-2.,  2.,
		 2., -2.,
		 4.,  0.,
		-4.,  0.,
		 0.,  4.,
		 0., -4.,
		};

/***********************************************************/
/***********************************************************/
/**Calculate sigma2 as a function of density from equation of state */
double eq_state(double rho){
	double sigma=0.0;
	double integral=0.0;

	if(rho <= R1 || rho > R5 ){
		integral = 0.0;
	}
	else if(R1 < rho && rho <= R2){
		integral = (rho-R1) * (rho-R1) / (R2-R1);
	}
	else if(R2 < rho && rho<= R3){
		integral = (rho - R2) * (2. - (rho - R2) / (R3 - R2)) + (R2 - R1);
	}
	else if(R3 < rho && rho<= R4){
		integral = R3 - R1 - (rho - R3) * (rho - R3) / (R4 - R3);
	}
	else if(R4 < rho && rho <= R5){
		integral = 2. * R3 - R1 - R4 - (rho - R4) * (2. - (rho - R4) / (R5 - R4));
	} else {
		printf ("I've lost myself in the EOS!\n");
		exit (0);
	}
	sigma= S1 * exp(0.5 * R0 * integral);
return sigma;
}

/******************************************************/
/**Derivative of sigma_2 with respect to density*/
double derS(double rho){
	double der;

	if((rho <= R1) || (rho > R5) ){
		der = 0.0;
	}
	else if(R1 < rho && rho <= R2){
		der = eq_state(rho) * R0 * (rho-R1) / (R2-R1);
	}
	else if(R2 < rho && rho<= R3){
		der = eq_state(rho) * R0 * (R3 - rho) / (R3 - R2);
	}
	else if(R3 < rho && rho<= R4){
		der = eq_state(rho) * R0 * (R3 - rho) / (R4 - R3);
	}
	else if(R4 < rho && rho <= R5){
		der = eq_state(rho) * R0 * (rho - R5) / (R5 - R4);
	} else {
		printf ("Something is wrong with the derS\n");
		exit (0);
	}
return der;
}


/******************************************************/
/**Second derivative of sigma_2 with respect to density*/
double der2S(double rho){
	double der;
	
	if((rho <= R1) || (rho > R5) ){
		der = 0.0;
	}
	else if(R1 < rho && rho <= R2){
		der = derS(rho) * R0 * (rho - R1) / (R2 - R1) + eq_state(rho) * R0 / (R2 - R1);
	}
	else if(R2 < rho && rho <= R3){
		der = derS(rho) * R0 * (R3 - rho) / (R3 - R2) - eq_state(rho) * R0 / (R3 - R2);
	}
	else if(R3 < rho && rho <= R4){
		der = derS(rho) * R0 * (R3 - rho) / (R4 - R3) - eq_state(rho) * R0 / (R4 - R3);
	}
	else if(R4 < rho && rho <= R5){
		der = derS(rho) * R0 * (rho - R5) / (R5 - R4) + eq_state(rho) * R0 / (R5 - R4);
	} else {
		printf ("Something is wrong with the der2S\n");
		exit (0);
	}
return der;
}

/******************************************************/
/**Derivative of p with respect to density*/
double derP(double rho){
	double der;

	if((rho <= R1) || (rho > R5) ){
		der = S1;
	} else {
		der = eq_state(rho) + rho * derS(rho); 
	}

return der;
}


/******************************************************/
/**Second derivative of p with respect to density*/
double der2P(double rho){
	double der;

	if((rho <= R1) || (rho > R5) ){
		der = 0.0;
	} else {
		der = 2.0*derS(rho) + rho * der2S(rho); 
	}

return der;
}

/***********************************************************/
/**Calculate c_s^2 and weights as a function of density*/
void weights(int pos){
	int i;
	double rho, sigma_2;
	/**set c_s^2 at each node as a function of density at that node*/
	rho = node[pos].rho;
	sigma_2 = eq_state(rho);
	node[pos].S= sigma_2;
	//if(rho<RHO_HIGH && rho>RHO_LOW)node[pos].S=(5.824*RHO_LOW-8.6272*RHO_LOW*RHO_LOW+3.408*RHO_LOW*RHO_LOW*RHO_LOW)*RHO_LOW/rho;
		
	/**set weights at each node*/
	node[pos].w[0]= 1. - sigma_2 * 45. * 0.5 * (7./60. - 7./48. * sigma_2 + sigma_2 * sigma_2 / 16.);
	node[pos].w[1]= sigma_2 / 3. * (32./15. - 4. * sigma_2 + 2. * sigma_2 * sigma_2);
	node[pos].w[5]= sigma_2 * sigma_2 * (1./3. - 0.25 * sigma_2);
	node[pos].w[9]= sigma_2 * (-1./18. + 3./16. * sigma_2 - sigma_2 * sigma_2 / 12.);
	node[pos].w[13]= sigma_2 * sigma_2 / 96. * (-0.5 + 1.5 * sigma_2);
	node[pos].w[17]= sigma_2 / 384. * (4./15. - sigma_2 + sigma_2 * sigma_2);

	for(i=1;i<4;i++){
		node[pos].w[1+i]=node[pos].w[1];
		node[pos].w[5+i]=node[pos].w[5];
		node[pos].w[9+i]=node[pos].w[9];		
		node[pos].w[13+i]=node[pos].w[13];
		node[pos].w[17+i]=node[pos].w[17];
	}

return;
}

/*******************************************************/
/**calculate density*/
void density(int pos){
	int v;
	double density=0.0;
	for(v=0;v<VDIM;v++){
		density += node[pos].pop[v];
	} 
	node[pos].rho=density;
//if(density<0.5 || density>1.4 ){printf("pos=%d rho=%e\n",pos,density); for(v=0;v<VDIM;v++){printf("node[%d].n=%e\n",v,vel[v].n[pos]);} exit(0);}
	
	return;
}
//

/*********************************************************/
/**Initialize macroscopic properties of the initial condition and density distributions at each node to equilibrium distributions - only 2D*/
void init(){
	int limit=8;
	int chx=-1, chy=-1;
	int x,y,v,i;
	long int seed=-15789138;
	double ux, uy, uxx, uyy, uxy, usq,r;
	FILE *velFile;
	for(x=0;x<XDIM;x++){
		r=ran2(&seed);
		for(y=0;y<YDIM;y++){
			
			/*2D small square initial condition*/
			int sqSide;
			sqSide = (int) (0.6*XDIM);
			int sqCentX, sqCentY;
			sqCentX = (int) (0.45*XDIM);
			sqCentY = (int) (0.45*YDIM);
			if (x >= sqCentX - (int)(sqSide * .5) &&
					x < sqCentX + (int)(sqSide * .5) &&
					y >= sqCentY - (int)(sqSide * .5) &&
					y < sqCentY + (int)(sqSide * .5)){
				node[YDIM*x+y].rho=RHO_MEAN + ran2(&seed) * (RHO_HIGH-RHO_MEAN);
			} else {
				node[YDIM*x+y].rho=RHO_LOW + ran2(&seed) * (RHO_HIGH-RHO_MEAN);
			}
			
			node[YDIM*x+y].u[0]=ux=0.0;//UAMP*sin(y*2.0*M_PI*N/(YDIM));
			node[YDIM*x+y].u[1]=uy=0.0;
			
			weights(YDIM*x+y);
			
			/**calculate equilibrium distributions at each node*/
			uxx=ux*ux;
			uyy=uy*uy;
			uxy=2*ux*uy;
			usq=uxx+uyy;
			
			for(v=0;v<VDIM;v++){
				node[YDIM*x+y].pop[v] = node[YDIM*x+y].w[v]*node[YDIM*x+y].rho;
			}
			
		}      
	}

	velFile = fopen("velocity_init","w");
	for(x=0;x<XDIM;x++) {
		for(y=0;y<YDIM;y++){
			fprintf(velFile, "%-*d %-*d %-*f  %-*f \n", 5,x,5,y,10,node[YDIM*x+y].u[0],10,node[YDIM*x+y].u[1]);
		}
	}
	fclose(velFile);

	return;
}

/**********************************************************/
/**Calculate total momentum: \vec j=\sum_r \sum n_i \vec c_i - only 2D*/
void total_momentum(){

	int x,y,v,pos;
	double Mx=0.0, My=0.0, M=0.0;
	double Fx, Fy;
	for(x=0;x<XDIM;x++){
		for(y=0;y<YDIM;y++){
			pos=YDIM*x+y;
			for(v=0;v<VDIM;v++){
				Mx += (node[pos].pop[v]*C[2*v]);
				My += (node[pos].pop[v]*C[2*v+1]);
				M += node[pos].pop[v];
			}
			Mx += (H*node[pos].Fx/2.0);
			Mx -= (der_g_pi1minus(0,x,y)/12.0);  //!!! ce dodamo ta clen je momentum pri 2D random init ful velik
			My += (H*node[pos].Fy/2.0); 
			My -= (der_g_pi1minus(1,x,y)/12.0);
		}
	}
printf("total Mom=(%e, %e) total Mass=%e\n\n",Mx,My,M); 

return;
}

/*************************************************************/
/**Interface force calculation*/
void interface_force(){
	int t, pos;
	int x, y, xnew, ynew;
	double Fx,Fy;

	for(x=0;x<XDIM;x++){
		for(y=0;y<YDIM;y++){
			node[pos].Fx = 0.;
			node[pos].Fy = 0.;
		}
	}
	
	for(x=0;x<XDIM;x++){
		for(y=0;y<YDIM;y++){
			pos=YDIM*x+y;
			Fx=0.0; Fy=0.0;
			for(t=0;t<TDIM;t++){
			//	check+=tau[t]*D[2*t]*D[2*t]*D[2*t]*D[2*t];
				
				xnew=x+D[2*t];
				ynew=y+D[2*t+1];
				if (xnew >= XDIM) xnew -= XDIM;
				if (xnew < 0)	xnew += XDIM;
				if (ynew >= YDIM) ynew -= YDIM;
				if (ynew < 0)	ynew += YDIM;
				
				/* formula 70, page 9 - Nikita: does not correspond to the paper anymore.
				 Probably, it should be page 9, Eq. 67 */
				Fx += tau[t] * D[2*t]   * node[YDIM * xnew + ynew].rho;
				Fy += tau[t] * D[2*t+1] * node[YDIM * xnew + ynew].rho;
			}
			Fx *= KAPPA * node[YDIM * xnew + ynew].rho / (A*A*A);
			Fy *= KAPPA * node[YDIM * xnew + ynew].rho / (A*A*A);
			node[pos].Fx += Fx;
			node[pos].Fy += Fy;
		}
	}
return;
}

/********************************************************/
/**Delta function*/
double delta(int a, int b){
	if (a == b) return 1.0;
	else return 0.0;
}

/******************************************************/
/** $\partial_\g (\pi_{\a\g}^{*(1)} - \pi_{\a\g}^{(1)})$. Nikita: probably, Appendix C.1. */
double der_g_pi1minus(int a, int x, int y){
	double der=0.0,derv1,derv2;
	double rho;
	rho = node[YDIM*x+y].rho;
	int g;
	DerivativeOf whata, whatg;

	if(a==0) whata = UX; else if (a==1) whata = UY;
	for(g=0;g<DIM;g++){
		if(g==0) whatg = UX; else if (g==1) whatg = UY;
		der += ((der1v2(PRESSURE,x,y,g) * (der1v2(whatg,x,y,a) + der1v2(whata,x,y,g))+
			rho*node[YDIM*x+y].S *der12(whatg,x,y,a,g)) );
		der -= (derS(rho) * rho * rho *  der12(whatg,x,y,a,g));

	}
	der += (rho * node[YDIM*x+y].S * der22(whata,x,y)); 
//	derv1 = (rho * der2P(rho) * der1(DENSITY,x,y,a) * (der1(UX,x,y,0)+der1(UY,x,y,1))); //derv1 in derv2 bi mogla dat isti rezultat pa ga ne
	derv2 = (rho *  der1v2(DERP,x,y,a) * (der1v2(UX,x,y,0)+der1v2(UY,x,y,1)));

	der-=derv2;
	
return der;
}

/*******************************************************/
/**$\partial_\g (\pi_{\a\g}^{(0)})$. Nikita: probably, Appendix C.4.*/
double der_g_pi0(int a, int x, int y){
	double der=0.0;
	double rho;
	rho = node[YDIM*x+y].rho;
	DerivativeOf what;
	double ua,ux,uy;
	ux=node[YDIM*x+y].u[0];
	uy=node[YDIM*x+y].u[1];

	if(a==0){
		der = der1v2(PRESSURE,x,y,a) + der1v2(RHOUXUX,x,y,0) + der1v2(RHOUXUY,x,y,1);
	}
	else if (a==1){
		der = der1v2(PRESSURE,x,y,a) + der1v2(RHOUXUY,x,y,0) + der1v2(RHOUYUY,x,y,1);
	}

return der;
}


/******************************************************/
/**$\partial_{\b\g} \pi_{\a\g}^{(0)}$. Nikita: no idea so far */
double der_b_g_pi0(int a, int x, int y, int b){
	double der=0.0;
	double rho;
	rho = node[YDIM*x+y].rho;
	DerivativeOf what;

		if(a==0)
			der =  der12(PRESSURE,x,y,a,b) + der12(RHOUXUX,x,y,0,b) + der12(RHOUXUY,x,y,1,b);
		else if(a==1)
			der = der12(PRESSURE,x,y,a,b) + der12(RHOUXUY,x,y,0,b) + der12(RHOUYUY,x,y,1,b);

return der;
}

/********************************************************/
/**$\partial_{\a}(frac{u_{\b}}{\rho})$ Nikita: no idea so far */
double der_u_rho(int b, int x, int y, int a){
	double der=0.0;
	double invRho;
	invRho = 1./node[YDIM*x+y].rho;
	if (b == 0){
		der = der1v2(UX,x,y,a) * invRho - node[YDIM*x+y].u[b] * der1v2(DENSITY,x,y,a) * invRho * invRho;
	} else if (b == 1){
		der = der1v2(UY,x,y,a) * invRho - node[YDIM*x+y].u[b] * der1v2(DENSITY,x,y,a) * invRho * invRho;
	} else {
	}

return der;
}


/*********************************************************/
/**First order correction term $\Ksi_{\alpha\beta\gamma}$*/
double correction_ksi(int a, int b, int g, int x, int y){
	double Xi=0.0;                         
	int a1,b1,g1,i;
	double rho;
	rho=node[YDIM*x+y].rho;
	DerivativeOf what;

	for(i=0;i<3;i++){
		if(i==0){a1=a; b1=b; g1=g;}
		else if(i==1){a1=b; b1=g; g1=a;}
		else if(i==2){a1=g; b1=a; g1=b;}
	
		if(b1==0 && g1==0) what=PUXUX;
		else if (b1==1 && g1==1) what=PUYUY;
		else if ((b1==1 && g1==0)||(b1==0 && g1==1)) what=PUXUY;

		Xi+= der1v2(what,x,y,a1);
		if(b1==g1){
			Xi+= (rho*node[YDIM*x+y].S * der1v2(SIGMA,x,y,a1) + rho*rho*derS(rho)*node[YDIM*x+y].u[a1] *
											 (der1v2(UX,x,y,0)+der1v2(UY,x,y,1)));
		}
	}
	Xi*=((GAMMA*GAMMA+4.0*GAMMA+1.0)/(3.0*(GAMMA+1.0)));

return Xi;
}

/**********************************************************/
/**Second order correction term $Sigma_{alpha\beta}$*/
double correction_sigma(int a, int b, int x, int y){
	double Sigma, Sigma1, Sigma2, Sigma3, Sigma4, Sigma5, S, rho, Sigma0, Sigma00;
	int a1=0;
	DerivativeOf whata,whatb,whata0,whata1,whatb0,whatb1;
	Sigma0=0.0; Sigma00=0.0;
	Sigma=0.0;
	Sigma1=0.0;
	Sigma2=0.0;
	Sigma3=0.0;
	Sigma4=0.0;
	Sigma5=0.0;

	if (a==0) whata=UX;  else if (a==1)whata=UY;
	if (b==0) whatb=UX;  else if (b==1)whatb=UY;

	S=node[YDIM*x+y].S;
	rho=node[YDIM*x+y].rho;

	/*first line*/
	Sigma1= ((1.0+GAMMA)*(1.0+GAMMA)/(4.0*(1.0-GAMMA)))* (node[YDIM*x+y].u[a]*der_g_pi1minus(b,x,y) + node[YDIM*x+y].u[b]*der_g_pi1minus(a,x,y));   //first line  $\part_{t_2} \pi_{\a\b}^{(0)}$

	//second line
	Sigma2= (der1v2(JX,x,y,0) + der1v2(JY,x,y,1)) * (- derP(rho) * (der1v2(whatb,x,y,a) + der1v2(whata,x,y,b)) + delta(a,b) * der2P(rho) * rho * (der1v2(UX, x,y,0) + der1v2(UY, x,y,1))); //second line

	//third line
	Sigma3= der1v2(DENSITY, x, y, a) *der_g_pi0(b,x,y) / (rho) - der_b_g_pi0(b,x,y,a)
		+der1v2(DENSITY, x, y, b) *der_g_pi0(a,x,y) / (rho) - der_b_g_pi0(a,x,y,b);
	Sigma3*= S;

	//fourth line
	Sigma4= (der1v2(JX,x,y,0) + der1v2(JY,x,y,1)) * (der_u_rho(b,x,y,a) + der_u_rho(a,x,y,b))
		+ node[YDIM*x+y].u[b] * (der12(JX,x,y,a,0) +der12(JY,x,y,a,1))/(rho) +  node[YDIM*x+y].u[a] * (der12(JX,x,y,b,0) +der12(JY,x,y,b,1))/(rho);

	Sigma4*= (S*rho);

	//fifth line
	if(a==b){
		Sigma5= (der1v2(DENSITY,x,y,0) * der_g_pi0(0,x,y) + der1v2(DENSITY,x,y,1) * der_g_pi0(1,x,y))/(rho*rho)
			- (der_b_g_pi0(0,x,y,0) +der_b_g_pi0(1,x,y,1))/(rho)
			+ (der_u_rho(0,x,y,0) + der_u_rho(1,x,y,1)) * (der1v2(JX,x,y,0) + der1v2(JY,x,y,1))
			+ (node[YDIM*x+y].u[0] * (der12(JX,x,y,0,0) + der12(JY,x,y,0,1)) + node[YDIM*x+y].u[1] *  (der12(JX,x,y,0,1) + der12(JY,x,y,1,1)))/(rho);

		Sigma5*= (S*rho - derP(rho) * rho);
	}
	Sigma=Sigma1 + (Sigma2 +Sigma3+Sigma4+Sigma5)*(GAMMA*GAMMA + 4.0*GAMMA +1.0)/(6.0*(GAMMA - 1.0));

return Sigma;
}

/**********************************************************/
/**Implicit calculation of velocity*/
void velocity(int rep){

	double ux0[XDIM*YDIM], uy0[XDIM*YDIM];
	double dif=0.0, difference=0.0;
	int x,y,t,v,xnew, ynew;
	double ux=0.0,uy=0.0;
	int count,pos;

	for(x=0;x<XDIM;x++){
		for(y=0;y<YDIM;y++){
			pos=YDIM*x+y;

			/*determine velocity*/
			ux0[pos] = 0.0;
			uy0[pos] = 0.0;

			for(v=0;v<VDIM;v++){
				ux0[pos] += (node[pos].pop[v]*C[2*v]);
				uy0[pos] += (node[pos].pop[v]*C[2*v+1]);
			}

			ux0[pos] = (ux0[pos] + node[pos].Fx*0.5)/(node[pos].rho);
			uy0[pos] = (uy0[pos] + node[pos].Fy*0.5)/(node[pos].rho);
			node[pos].u[0] = ux0[pos];
			node[pos].u[1] = uy0[pos];

//if((fabs(ux0[pos])>0.5) || (fabs(uy0[pos])>0.5)){printf("ERR:velocity too big r=(%d,%d)    v=(%f,%e)\n",x,y,ux0[pos],uy0[pos]);}
		}
	}	


	/*Implicit algorithm for the calculation of velocity*/  //!!!ce to dodam momentum ni vec ohranjen!!!!
	dif=1.0; count=0;
	while(dif > 1e-15){	
		count++;
		dif=0.0;
		for(x=0;x<XDIM;x++){
			for(y=0;y<YDIM;y++){
				ux=0.0;
				uy=0.0;
				ux = ux0[YDIM*x+y] - (der_g_pi1minus(0,x,y)/(12.0 * node[YDIM*x+y].rho));
				uy = uy0[YDIM*x+y] - (der_g_pi1minus(1,x,y)/(12.0 * node[YDIM*x+y].rho));
				
				difference = fabs(node[YDIM*x+y].u[0]-ux);
				if(difference > dif) dif=difference;
				difference = fabs(node[YDIM*x+y].u[1]-uy);
				if(difference > dif) dif=difference;
			
				node[YDIM*x+y].u[0]=ux;
				node[YDIM*x+y].u[1]=uy;

			}
		}
	}

return;
}


/**********************************************************/
/**Collision step - only 2D*/
void collide(int rep, double *W){

	int x,y,v,a,b,g,pos,t; 
	double delta_i[VDIM],delta_j2[VDIM], delta_c2[VDIM];  //collision operator due to external force and corrections terms
	double corr=0.0,delta_c1=0.0, n_v_eq=0.0;
	double ux=0.0,uy=0.0,S;
	double ux0,r,fx0,check,correction;
	for(v=0;v<VDIM;v++){
		delta_i[v]=0.0;
		delta_j2[v]=0.0;
		delta_c2[v]=0.0;
	}

	/*calculate total density and velocity at every node*/
	for(x=0;x<XDIM;x++){
		for(y=0;y<YDIM;y++){
			pos=YDIM*x+y;
			density(pos);  
			weights(pos);
		}
	}
	interface_force();
	velocity(rep);
	

	for(x=0;x<XDIM;x++){
		for(y=0;y<YDIM;y++){
			pos=YDIM*x+y;
			//if(y==11 && (x==100 || x==101 ))printf("ux0=%e, rho=%e %e\n\n",ux0/r,r, node[YDIM*x+y].rho);
			ux=node[pos].u[0];
			uy=node[pos].u[1];
			
			
			/*calculate interface collision operator delta_i[v]*/
			for(v=1;v<VDIM;v++){
				check=0.0;
				/* formula 71, page 9 */
				delta_i[v]=((GAMMA+1.0)/(2.0*cs2_eq))*(C[2*v]*node[pos].Fx + C[2*v+1]*node[pos].Fy) * W[v];
			}
			
			/*new density distributions - with interface force*/
			ux0=0.0; fx0=0.0; r=0.0;
			S=node[YDIM*x+y].S;
			for(v=0;v<VDIM;v++){
				/*calculate equilibrium distribution*/
				n_v_eq = node[pos].w[v] * node[pos].rho * (1.0 + (1. / S) * (ux*C[2*v]+uy*C[2*v+1]) +
																								 0.5 * (ux*C[2*v] + uy*C[2*v+1]) * (ux*C[2*v] + uy*C[2*v+1]) / S -
																								 0.5 * (ux*ux+uy*uy));

				node[pos].pop[v] = GAMMA*node[pos].pop[v] + (1.0-GAMMA)*n_v_eq + delta_i[v];
			}
		}
	}
	/*calculate correction collision operator of first order*/

	S=cs2_eq;//node[YDIM*x+y].S;
	for(x=0;x<XDIM;x++){
		for(y=0;y<YDIM;y++){
			for(v=0;v<VDIM;v++){
				delta_c1=0.0;
				for(a=0;a<DIM;a++){
					for(b=0;b<DIM;b++){
						for(g=0;g<DIM;g++){
							delta_c1 += (correction_ksi(a,b,g,x,y)*(C[2*v+a]*C[2*v+b]*C[2*v+g] - S*(delta(a,b)*C[2*v+g] + delta(a,g)*C[2*v+b] + delta(b,g)*C[2*v+a]) ));
						}
					}
				}
				delta_c1 *= (W[v]/(6.0*S*S*S));
				node[YDIM*x+y].pop[v] += delta_c1;
			}
			//printf("%f %f\n",node[YDIM*x+y].rho, derS(node[YDIM*x+y].rho));
		}
		check=0;
		for(v=0;v<VDIM;v++){
			if(node[YDIM*x+1].pop[v] < 0.0) {
				check=1;
			} break;
		}

		if(check==1) {
			for(v=0;v<VDIM;v++){
				printf("(%d,%d)    v=%d,    n=%e      rho=%f\n",x,1,v,node[YDIM*x+1].pop[v],node[YDIM*x+1].rho  );
			}
			exit(0);
		}
	}

	/*calculate correction collision operator of second order*/
	for(x=0;x<XDIM;x++){
		for(y=0;y<YDIM;y++){
			S=cs2_eq;//node[YDIM*x+y].S;
			for(v=0;v<VDIM;v++){
				delta_c2[v]=0.0;
				delta_j2[v]=0.0;
			}
			for(a=0;a<DIM;a++){
				for(b=0;b<DIM;b++){
					correction=correction_sigma(a,b,x,y);
					for(v=0;v<VDIM;v++){
						delta_c2[v]+=(correction*(C[2*v+a]*C[2*v+b] - S*(delta(a,b))));
						
					}
				}
			}
			
			for(a=0;a<DIM;a++){
				correction=der_g_pi1minus(a,x,y);
				for(v=0;v<VDIM;v++){
					delta_j2[v]+=(correction * C[2*v+a]);
				}
			}
			
			for(v=0;v<VDIM;v++){
				delta_c2[v]*=(W[v]/(2.0*S*S));//printf("delta_c2[%d]=%e\n",v,delta_c2);
				delta_j2[v]*=((1.0 - GAMMA)*W[v]/(S * 12.0)); //if(fabs(delta_j2)>0.5)printf("delta_j2[%d]=%e\n",v,delta_j2);
				//	if(delta_c2[v]>1e-1)printf("v=%d	delta_c2=%e	x=%d y=%d\n",v, delta_c2[v],x,y );
				
				node[YDIM*x+y].pop[v] += (delta_c2[v] + delta_j2[v]); //!!! ce dodamo j2 je pri 2D random init ful velik momentum
			}
		}
	}
	return;
}

/**********************************************************/
/**********************************************************/
#if 0
int main(){
	int x,y,rep;
	char filename[20], filename2[20];
	FILE *density, *velocityX, *velocityY, *decay, *density2, *state, *pressure, *speed, *forceX, *forceY;
	double di;

	int xnew, ynew,t,i;
	double Fx,Fy;
	double W[VDIM];
	double test;

	cs2_eq = eq_state(RHO_MEAN);
	
	W[0]=(1.0-cs2_eq*45.0*0.5*(7.0/60.0 - cs2_eq * 7.0/48.0 +  cs2_eq *  cs2_eq/16.0));
	W[1]=cs2_eq/3.0 *(32.0/15.0 - 4.0*cs2_eq +2.0*cs2_eq * cs2_eq);
	W[5]=cs2_eq * cs2_eq * (1.0/3.0 - cs2_eq/4.0);
	W[9]=cs2_eq * (-1.0/18.0 + (3.0/16.0) * cs2_eq - cs2_eq * cs2_eq/12.0 );
	W[13]=(cs2_eq * cs2_eq/96.0) * (-0.5 + 1.5 * cs2_eq);
	W[17]=(cs2_eq /384.0) * (4.0/15.0 - cs2_eq + cs2_eq * cs2_eq);
	for(i=1;i<4;i++){
		W[1+i]=W[1];
		W[5+i]=W[5];
		W[9+i]=W[9];		
		W[13+i]=W[13];
		W[17+i]=W[17];
	}
	printf("cs2_eq=%f\n",eq_state(RHO_MEAN));

	printf("parameters are:\n KAPPA_LB=%f\n GAMMA=%f\n TAU2=%e\n MU_SIM=%f\n cs2_eq=%f\n",KAPPA,GAMMA,TAU2, (1+GAMMA)*RHO_MEAN*cs2_eq/(2*(1-GAMMA)),cs2_eq);

	state = fopen("output/eq_state","w");
	for(di=0.01;di<2.0;di+=0.001)
	fprintf(state, "%-*f %-*f  %-*f  %-*f  %-*f  %-*f  %-*f  \n", 5,di,5, eq_state(di),5, derS(di), 5, der2S(di), 5, di*eq_state(di), 5, derP(di), 5, der2P(di));
	fclose(state);

	init();
	printf ("FINISHED init!\n");
	for(rep=0;rep<REP;rep++){
		collide(rep,W);
		stream();

		if(rep%10==0){
			printf("%d	",rep);
			total_momentum();

			sprintf(filename, "output/density%d",rep);
			density = fopen(filename,"w");
			for(y=0;y<XDIM;y++){
				fprintf(density,"\n");
				for(x=0;x<YDIM;x++){ 
					fprintf(density, "%-*f", 10,node[YDIM*x+y].rho);
				}
			}
			fclose(density);

			sprintf(filename, "output/pressure%d",rep);
			pressure = fopen(filename,"w");
			for(y=0;y<XDIM;y++){
				fprintf(pressure,"\n");
				for(x=0;x<YDIM;x++){ 
					fprintf(pressure, "%-*f", 10,node[YDIM*x+y].rho*node[YDIM*x+y].S);
				}
			}
			fclose(pressure);
			
			sprintf(filename, "output/S%d",rep);
			speed = fopen(filename,"w");
			for(y=0;y<XDIM;y++){
				fprintf(speed,"\n");
				for(x=0;x<YDIM;x++){ 
					fprintf(speed, "%-*f", 10,node[YDIM*x+y].S);
				}
			}
			fclose(speed);

			sprintf(filename, "output/Ux%d",rep);
			velocityX = fopen(filename,"w");
			for(y=0;y<XDIM;y++){
				fprintf(velocityX,"\n");
				for(x=0;x<YDIM;x++){ 
					fprintf(velocityX, "%-*f", 10,node[YDIM*x+y].u[0]);
				}	
			}
			fclose(velocityX);
			
			sprintf(filename, "output/Uy%d",rep);
			velocityY = fopen(filename,"w");
			for(y=0;y<XDIM;y++){
				fprintf(velocityY,"\n");
				for(x=0;x<YDIM;x++){ 
					fprintf(velocityY, "%-*f", 10,node[YDIM*x+y].u[1]);
				}	
			}
			fclose(velocityY);
			
			sprintf(filename, "output/Fx%d",rep);
			sprintf(filename2, "output/Fy%d",rep);
			forceX = fopen(filename,"w");
			forceY = fopen(filename2,"w");
			for(y=0;y<XDIM;y++){
				fprintf(forceX,"\n");
				fprintf(forceY,"\n");
				for(x=0;x<YDIM;x++){ 
					Fx=node[YDIM*x+y].Fx;
					Fy=node[YDIM*x+y].Fy;
				
					fprintf(forceY, "%-*f", 10,Fy);
					fprintf(forceX, "%-*f", 10,Fx);
				}	
			}
			fclose(forceX);	
			fclose(forceY);
		}
	}

return 0;
}
#endif
