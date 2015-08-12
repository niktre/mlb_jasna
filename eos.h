/***********************************************************************
 *
 * eos.h
 *
 * Copyright (c) 2015 Ulf D. Schiller <ulf@lattice-boltzmann.de>
 * All rights reserved.
 *
 ***********************************************************************/

#ifndef EOS_H
#define EOS_H

#include <stdlib.h>

#define SCALE (1.0)

/* phase separation densities */
#define RHO_LOW  (1.1*pow(SCALE,3)) //1.13 //0.577348
#define RHO_HIGH (1.2*pow(SCALE,3)) //1.47 //1.325965
#define RHO_MEAN ((RHO_LOW+RHO_HIGH)/2.0)

/* constants that define equation of state: p/rho */
#define A0 (1.2/pow(SCALE,3)) /* specific volume */
#define S1 0.5                /* reference speed of sound squared */
#define S2 0.7
#define R1 (0.5*pow(SCALE,3)) /* densities */
#define R2 (1.0*pow(SCALE,3))
#define R3 (2.0*R2-R1)
#define W 0.05

/***********************************************************************/

#if defined EOS_PSI

static double psi(double rho) {
  return A0*sin(M_PI*(rho-R1)/(R2-R1));
}

static double dpsi(double rho) {
  return A0*M_PI/(R2-R1)*cos(M_PI*(rho-R1)/(R2-R1));
}

static double intpsi(double rho) {
  return A0*(R2-R1)/M_PI*(1.0 - cos(M_PI*(rho-R1)/(R2-R1)));
}

/***********************************************************************/

static double eq_state(double rho) {
  double cs2;
  
  if (rho <= R1 || rho >= R3) {
    fprintf(stderr, "Don't like to be here with rho = %f\n", rho);
    cs2 = S1;
  } else if (rho < R3) {
    cs2 = S1 * exp(intpsi(rho));
  } else {
    printf ("I've lost myself in the EOS (rho=%f)!\n", rho);
    exit (0);
  }

  return cs2;
    
}

/***********************************************************************/

static double derP(double rho) {
  double dp;

  if (rho <= R1 || rho >= R3) {
    dp = S1;
  } else {
    /* dp = eq_state(rho) + rho*derS(rho); */
    dp = eq_state(rho)*(1. + rho*psi(rho));
  }

  return dp;

}
  
/***********************************************************************/

static double der2P(double rho) {
  double d2p;
  
  if (rho <= R1 || rho >= R3) {
    d2p = 0.0;
  } else {
    /* d2p = 2.0*derS(rho) + rho*der2S(rho); */
    d2p = eq_state(rho)*(psi(rho)*(2. + rho*psi(rho)) + rho*dpsi(rho));
  }

  return d2p;

}

#else

/***********************************************************************/

static double eq_state(double rho) {
  double cs2;

  cs2 = (S2 - S1) * exp( -(rho-R2)*(rho-R2)/W ) + S1;

  return cs2;
}

static double derS(double rho) {
  double ds;

  ds = -2./W*(rho-R2)*(eq_state(rho) - S1);

  return ds;
}

static double der2S(double rho) {
  double d2s;

  d2s = -2./W*(eq_state(rho) - S1 + (rho-R2)*derS(rho));

  return d2s;
}

static double derP(double rho) {
  double dp;

  dp = eq_state(rho) + rho*derS(rho);

  return dp;
}

static double der2P(double rho) {
  double d2p;

  d2p = 2.*derS(rho) + rho*der2S(rho);

  return d2p;
}

#endif

/***********************************************************************/

inline static void write_eos() {
  FILE* file;
  double rho, cs2, p, dp, d2p;

  file = fopen("eos.dat", "w");

  for (rho=0.4; rho<=1.6; rho+=0.01) {

    cs2 = eq_state(rho);
    p = rho*cs2;
    dp = derP(rho);
    d2p = der2P(rho);

    fprintf(file, "%f %f %f %f %f\n", rho, cs2, p, dp, d2p);

  }

  fclose(file);

}

/***********************************************************************/

#endif
