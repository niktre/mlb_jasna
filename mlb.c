/***********************************************************************
 *
 * mlb.c
 *
 * Copyright (c) 2015 Ulf D. Schiller <ulf@lattice-boltzmann.de>
 * All rights reserved.
 *
 ***********************************************************************/

//#include "defs.h"
#include <math.h>
#include "d2q21.h"
/***********************************************************************/

void mlb_calc_force(double *force, double *f, int x, int y) {
  const double *w = lbmodel.fd_weights[3];
  const double (*c)[lbmodel.n_dim] = lbmodel.c;
  double *m = f + lblattice.halo_grid_volume*lbmodel.n_vel;
  int i;
  double rho, nb_rho;

  rho = m[0];

  for (i=0; i<lbmodel.n_fd; ++i) {
    nb_rho = *(m + lblattice.nb_offset[i]*lbmodel.n_vel);
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

  rho = lbpar.rho;
  cs2 = eq_state(rho);
  lb_weights(w, cs2);

  for (i=0; i<lbmodel.n_vel; ++i) {
    fc = lbmodel.c[i][0]*force[0] + lbmodel.c[i][1]*force[1];
    f[i] += 0.5*(1. + lbpar.gamma)*w[i]/cs2*fc;
  }

}

/***********************************************************************/

double eq_state(double rho){
	double sigma=0.0;
	double integral=0.0;
	
	if(rho <= R1 || rho >= R3 ){
		integral = 0.0;
	}
	else if(R1 < rho && rho < R3){
		integral = 2*A*(R2-R1)/M_PI*sin(0.5*M_PI*(rho-2*R1+R0)/(R2-R1))*sin(0.5*M_PI*(rho-R0)/(R2-R1));
	} else {
		printf ("I've lost myself in the EOS!\n");
		exit (0);
	}
	sigma= S1 * exp(integral);
	return sigma;
}

/***********************************************************************/
