/***********************************************************************
 *
 * mlb.c
 *
 * Copyright (c) 2015 Ulf D. Schiller <ulf@lattice-boltzmann.de>
 * All rights reserved.
 *
 ***********************************************************************/

#include <math.h>
#include <stdio.h>
#include "defs.h"

#define TOLERANCE 1.e-9

/***********************************************************************/

void mlb_calc_force(double *force, double *m, int x, int y) {
  const double *w = lbmodel.fd_weights[3];
  const double (*c)[lbmodel.n_dim] = lbmodel.c;
  int i;
  double rho, nb_rho;

  force[0] = force[1] = 0.0;

  rho = m[0];

  for (i=0; i<lbmodel.n_fd; ++i) {
    nb_rho = (m + lblattice.nb_offset[i]*lbmodel.n_vel)[0];
    force[0] += w[i]*c[i][0]*nb_rho;
    force[1] += w[i]*c[i][1]*nb_rho;
  }

  force[0] *= lbpar.kappa*rho;
  force[1] *= lbpar.kappa*rho;

}

/***********************************************************************/

static void mlb_calc_current(double *jc, double *m, int x, int y) {
  jc[0] = 0.0;
  jc[1] = 0.0;
}

/***********************************************************************/

static void mlb_init_current(double *m) {
  int x, y, xl, xh, yl, yh, xoff;
  double rho, p, dp, d2p, force[lbmodel.n_dim];

  xl = lblattice.halo_size[0];
  xh = lblattice.halo_size[0] + lblattice.grid[0];
  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1];

  xoff = lblattice.halo_grid[0] - (xh - xl);

  m += lbmodel.n_vel*(xl*lblattice.stride[0]+yl);

  for (x=xl; x<xh; ++x, m+=lbmodel.n_vel*xoff) {
    for (y=yl; y<yh; ++y, m+=lbmodel.n_vel) {

      rho = m[0];
      p   = rho*eq_state(rho);
      dp  = derP(rho);
      d2p = der2P(rho);

      /* Step 1 */
      m[6] = p;
      m[7] = dp;
      m[8] = p - rho*dp;
      m[9] = rho*d2p;

      /* Step 2 */
      mlb_calc_force(force, m, x, y);
      m[10] = force[0];
      m[11] = force[1];

      /* Step 3 */
      m[12] = (m[1] + 0.5*force[0])/m[0];
      m[13] = (m[2] + 0.5*force[1])/m[0];

      //m[14] = 0.0; /* correction current */
      //m[15] = 0.0;

    }
  }

}

/***********************************************************************/

static void ic_read(double *dmax, double *m, int x, int y) {
  double jnew[lbmodel.n_dim], *jc = m + 14, d;

  mlb_calc_current(jnew, m, x, y);

  d = fabs(jnew[0] - jc[0]);
  if (d > *dmax) *dmax = d;
  d = fabs(jnew[1] - jc[1]);
  if (d > *dmax) *dmax = d;

  jc[0] = jnew[0];
  jc[1] = jnew[1];

}

/***********************************************************************/

static void ic_write(double *m, int x, int y) {

  double rho = m[0];
  double *j  = m + 1;
  double *g  = m + 10;
  double *u  = m + 12;
  double *jc = m + 14;

  u[0] = (j[0] + 0.5*g[0] + jc[0])/rho;
  u[1] = (j[1] + 0.5*g[1] + jc[1])/rho;

}

/***********************************************************************/

static void ic_read_column(double *dmax, double *m, int x) {
  int y, yl, yh;

  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1] - 1;

  for (y=yl, m+=yl*lbmodel.n_vel; y<=yh; ++y, m+=lbmodel.n_vel) {
    ic_read(dmax, m, x, y);
  }

}

/***********************************************************************/

static void ic_write_column(double *m, int x) {
  int y, yl, yh;

  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1] - 1;

  for (y=yl, m+=yl*lbmodel.n_vel; y<=yh; ++y, m+=lbmodel.n_vel) {
    ic_write(m, x, y);
  }

}

/***********************************************************************/

void mlb_correction_current(double *m0) {
  double *m = m0;
  int x, xl, xh;
  int xstride = lblattice.stride[0]*lbmodel.n_vel;

  int niter = 0;
  double dmax = 0.0;

  lb_halo_copy(); /* need up to date densities in halo */

  mlb_init_current(m);

  do {

    m = m0;

    ++niter;

    lb_halo_copy(); /* need up to date currents in halo */

    /* no information is `streaming' so we loop the internal region */
    xl = lblattice.halo_size[0];
    xh = lblattice.halo_size[0]+lblattice.grid[0]+VMAX;

    /* Columns in the lower range will read only */
    for (x=xl, m+=xl*xstride; x<xl+VMAX; ++x, m+=xstride) {
      ic_read_column(&dmax, m, x);
    }

    /* Calculate the correction current but do not overwrite the column yet
     * x-MAXV can be overwritten with the updated current */
    for (x=xl+VMAX; x<xh-VMAX; ++x, m+=xstride) {
      ic_read_column(&dmax, m, x);
      ic_write_column(m-VMAX*xstride, x-VMAX);
    }

    /* Columns in the higher range will write only */
    for (x=xh-VMAX; x<xh; ++x, m+=xstride) {
      ic_write_column(m-VMAX*xstride, x-VMAX);
    }

  } while (dmax > TOLERANCE);

  //fprintf(stderr, "Implicit algorithm converged after %d iteration(s).\n", niter);

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
