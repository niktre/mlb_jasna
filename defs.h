#ifndef _DEFS_H
#define _DEFS_H

/*****************************************/
/*PARAMETERS TO BE DETERMINED*/
/* Reynolds and Weber numbers */
#define RE 1.
#define WE 2.0 // must be < XDIM and it doesnt work if we<2??--> too big force, the lower limit is determined by tau2, which should be <1e-1! If size of lattice is increased, lower limit for We increases!

/* lattice size */
#define XDIM 50
#define YDIM XDIM

/* number of iterations */
#define REP 100000

/*
 RE=1.000000
 MU_LB=0.500000
 WE=2.000000
 KAPPA_LB=0.0625001
 AA=2.000000e-02
 HH=2.000000e-04
 GAMMA=0.0641433
 TAU2=1.662723e-02
 MU_SIM=0.500000
 SIGMA2=0.676498
 */

/*****************************************/
/*PARAMETERS OF THE SYSTEM AND SIMULATION*/
/* dimension of the problem: 2D, 3D */
#define DIM 2

/* number of velocity vectors */
#define VDIM 21

/* number of velocity vectors for interface force calculation */
#define TDIM 12

/* constants that define equation of state: p/rho */
#define R0 1.2 
#define R1 0.5
#define R2 0.7
#define R3 1.0
#define R4 1.3
#define R5 1.5
#define S1 0.6

/* phase separation densities */
#define RHO_LOW 1.13 //0.577348
#define RHO_HIGH 1.47 //1.325965
#define RHO_MEAN 1.3//(RHO_LOW+RHO_HIGH)/2.0

/* time step: 1 because we do simulation in lattice units- viscosity and surface tension parameter kappa have to be rescaled accordingly! */
#define H 1.0

/* lattice size */
#define A 1.0

/* time step and lattice size used to convert constants from dimensionless to lattice units condition for small mach number: HH<<AA */
#define KK 0.01
#define AA (1.0/(XDIM)) //lattice
#define HH (KK*AA)    //time

/* 2D velocity vectors (21) */
extern double C[DIM*VDIM];

/* 2D velocity vectors (12) for interface force calculation */
extern int D[DIM*TDIM];

/* surface tension constant */
#define KAPPA (KK*KK)/(AA*AA*WE*WE)		// 4e-6

/* Initialization constants */
#define AMP 1.0
#define WIDTH 10

/* SIGMA2=c_s^2 at mean density (rho=1.3) */
#define SIGMA2 0.676498//(0.3*sin(RHO_MEAN*10.0)+0.8)//(A1*RHO_MEAN + A2*RHO_MEAN*RHO_MEAN + A3*RHO_MEAN*RHO_MEAN*RHO_MEAN)

/* BGK collision operator constant OMEGA=1/tau, GAMMA has to be between -1 and 1 */
#define CC (2.0*KK)/(RE*AA*RHO_MEAN*SIGMA2)
#define GAMMA (CC-1.0)/(CC+1.0)
#define OMEGA (1.0-GAMMA)

/* weight tau for velocity set D */
//#define TAU2 1.662723e-02
#define TAU2 ((GAMMA+1.0)*KAPPA/4.0) //(GAMMA+1)*H*H*KAPPA/(4*A*A*A*A*SIGMA2)
#define TAU1 (-4.0*TAU2)
#define TAU3 (0.5*TAU2)

/* weights tau */
extern double tau[TDIM];

/* Force and initial velocity amplitude */
#define FAMP 0.0
#define UAMP 0.0
#define N 1

typedef struct
{
double rho;
double u[DIM];
double S;        /*\sigma_2 = S is any function of density rho*/
double w[VDIM]; /*weights corresponding to the velocity vectors*/
double Fx;
double Fy;
}Lattice;
Lattice node[XDIM*YDIM];

typedef struct
{
double n[XDIM*YDIM];
}Velocity;

typedef enum{
	UX,
	UY,
	PUXUX,
	PUYUY,
	PUXUY,
	PRESSURE,
	DENSITY,
	JX,
	JY,
	RHOUXUX,
	RHOUXUY,
	RHOUYUY,
	DERP,
	SIGMA
}DerivativeOf;

/********************************************/
/*MAIN FUNCTIONS OF LATTICE BOLTZMANN*/

double eq_state(double rho);/**Calculate sigma2 as a function of density from equation of state */

void init();/**Initialize macroscopic properties of the initial condition and density distributions at each node to equilibrium distributions - only 2D*/

void weights(int pos);/**Calculate c_s^2 and weights as a function of density*/

void density(int pos);/**calculate density*/

void total_momentum();/**Calculate total momentum: \vec j=\sum_r \sum n_i \vec c_i - only 2D*/

void velocity(int rep); /**Implicit calculation of velocity*/

void stream();/**Streaming step - only 2D*/

void collide(int rep, double *W);/**Collision step - only 2D*/

/******************************************/
/*FUNCTIONS FOR DERIVATIVES*/

double der1(DerivativeOf what, int x, int y, int direction);/**first derivative of f at x,y*/

double der1v2(DerivativeOf what, int x, int y, int direction);/**first derivative of f at x,y with higher precission*/

double der22(DerivativeOf what, int x, int y);/**Second - Laplacian derivative of f at x,y*/

double der12(DerivativeOf what, int x, int y, int dir1, int dir2);/**Second derivative of f at x,y*/

double derS(double rho);/**Derivative of sigma_2 with respect to density*/

double der2S(double rho);/**Second derivative of sigma_2 with respect to density*/

double derP(double rho);/**Derivative of p with respect to density*/

double der2P(double rho);/**Second derivative of p with respect to density*/

/******************************************/
/*FUNCTIONS FOR COLLISION OPERATOR*/

double der_g_pi1minus(int a, int x, int y);/** $\partial_\g (\pi_{\a\g}^{*(1)} - \pi_{\a\g}^{(1)})$ */

double der_g_pi0(int a, int x, int y); /**$\partial_\g (\pi_{\a\g}^{(0)})$*/

double der_b_g_pi0(int a, int x, int y, int b);/**$\partial_{\b\g} \pi_{\a\g}^{(0)}$*/

double correction_ksi(int a, int b, int g, int x, int y);/**First order correction term $\Ksi_{\alpha\beta\gamma}$*/

double correction_sigma(int a, int b, int x, int y);/**Second order correction term $Sigma_{alpha\beta}$*/

void interface_force();/**Interface collision operator calculation*/

/*****************************************/

double ran2(long *idum);/*Random number [0,1] generator*/

double delta(int a, int b); /**Delta function*/


#endif
