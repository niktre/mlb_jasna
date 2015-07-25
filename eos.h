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

/* phase separation densities */
#define RHO_LOW 1.13 //0.577348
#define RHO_HIGH 1.47 //1.325965
#define RHO_MEAN ((RHO_LOW+RHO_HIGH)/2.0)

/* constants that define equation of state: p/rho */
#define A  1.2
#define S1 0.6
#define R1 0.5
#define R2 1.0
#define R3 (2.0*R2-R1)

/***********************************************************************/

static double eq_state(double rho) {
  double cs2, integral;
  
  if (rho <= R1 || rho >= R3) {
    integral = 0.0;
  } else if (rho < R3) {
    integral = A*(R2-R1)/M_PI*(1.0 - cos(M_PI*(rho-R1)/(R2-R1)));
  } else {
    printf ("I've lost myself in the EOS!\n");
    exit (0);
  }

  cs2 = S1 * exp(integral);
  
  return cs2;
    
}

/***********************************************************************/

static double psi(double rho) {
  return A*sin(M_PI*(rho-R1)/(R2-R1)); 
}

static double dpsi(double  rho) {
  return A*M_PI/(R2-R1)*cos(M_PI*(rho-R1)/(R2-R1));
}

/***********************************************************************/

static double derS(double rho) {
  double ds;

  if (rho <= R1 || rho >= R3) {
    ds = 0.0;
  } else if (rho < R3) {
    ds = eq_state(rho)*psi(rho);
  } else {
    printf ("I've lost myself in derS!\n");
    exit (0);
  }

  return ds;
  
}

/***********************************************************************/

static double der2S(double rho) {
  double d2s;
  
  if (rho <= R1 || rho >= R3) {
    d2s = 0.0;
  } else if (rho < R3) {
    d2s = derS(rho) * psi(rho) + eq_state(rho) * dpsi(rho);
  } else {
    printf ("I've lost myself in der2S!\n");
    exit (0);
  }

  return d2s;

}

/***********************************************************************/

static double derP(double rho) {
  double dp;

  if (rho <= R1 || rho >= R3) {
    dp = S1;
  } else {
    dp = eq_state(rho) + rho*derS(rho);
  }

  return dp;

}
  
/***********************************************************************/

static double der2P(double rho) {
  double d2p;
  
  if (rho <= R1 || rho >= R3) {
    d2p = 0.0;
  } else {
    d2p = 2.0*derS(rho) + rho*der2S(rho);
  }

  return d2p;

}

/***********************************************************************/

#endif
