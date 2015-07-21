/***********************************************************************
 *
 * mlb.c
 *
 * Copyright (c) 2015 Ulf D. Schiller <ulf@lattice-boltzmann.de>
 * All rights reserved.
 *
 ***********************************************************************/

#include "defs.h"

/***********************************************************************/

void mlb_calc_force(double *force, double *f, int x, int y) {
}

void mlb_interface_collisions(double *f, double *force) {
  int i;
  double rho, cs2, fc;
  double w[lbmodel.n_vel];

  rho = RHO_MEAN;
  cs2 = eq_state(rho);
  lb_weights(w, cs2);

  for (i=0; i<lbmodel.n_vel; ++i) {
    fc = lbmodel.c[i][0]*force[0] + lbmodel.c[i][1]*force[1];
    f[i] += 0.5*(1. + lbpar.gamma)*w[i]/cs2*fc;
  }

}

/***********************************************************************/
