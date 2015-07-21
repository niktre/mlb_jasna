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
  const double *w = lbmodel.fd_weights[3];
  const double (*c)[lbmodel.n_dim] = lbmodel.c;
  double *m = f + lblattice.halo_grid_volume*lbmodel.n_vel;
  int i;
  double rho, nb_rho;

  rho = m[0];

  for (i=0; i<lbmodel.n_fd; ++i) {
    nb_rho = m[lblattice.nb_offset[i]*lbmodel.n_vel];
    force[0] += w[i]*c[i][0]*nb_rho;
    force[1] += w[i]*c[i][1]*nb_rho;
  }

  force[0] *= lbpar.kappa*rho;
  force[1] *= lbpar.kappa*rho;

}

/***********************************************************************/

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
