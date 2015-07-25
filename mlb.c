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
#include "derivFD.h"

#define TOLERANCE 1.e-3

/***********************************************************************/

void mlb_calc_force(double *force, double *m, int x, int y) {
  const double *w = lbmodel.fd_weights[3];
  const double (*c)[lbmodel.n_dim] = lbmodel.c;
  int i;
  double rho, nb_rho;

  force[0] = force[1] = 0.0;

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

static void mlb_calc_current(double *jc, double *m, int x, int y) {
  double Dp[lbmodel.n_dim],
    Dpmrdp[lbmodel.n_dim],
    Du[lbmodel.n_dim][lbmodel.n_dim],
    D2u[lbmodel.n_dim][lbmodel.n_dim][lbmodel.n_dim],
    divu;

  double *p     = m + 6;
  double *pmrdp = m + 8;
  double *u     = m + 12;

  firstDer(Dp, p);
  firstDer(Dpmrdp, pmrdp);
  firstDer(Du[0], &u[0]);
  firstDer(Du[1], &u[1]);
  secDerAB(D2u[0], &u[0]);
  secDerAB(D2u[1], &u[1]);

  divu = Du[0][0] + Du[1][1];

  jc[0] = Dpmrdp[0] * divu;
  jc[1] = Dpmrdp[1] * divu;

  jc[0] += Dp[0] * (Du[0][0] + Du[0][0]) + Dp[1] * (Du[0][1] + Du[1][0]);
  jc[1] += Dp[0] * (Du[1][0] + Du[0][1]) + Dp[1] * (Du[1][1] + Du[1][1]);

  jc[0] += (*pmrdp + *p) * (D2u[0][0][0] + D2u[1][0][1]);
  jc[1] += (*pmrdp + *p) * (D2u[0][1][0] + D2u[1][1][1]);

  jc[0] += *p * (D2u[0][0][0] + D2u[0][1][1]);
  jc[1] += *p * (D2u[1][0][0] + D2u[1][1][1]);

  jc[0] *= - 1./12.;
  jc[1] *= - 1./12.;

  //jc[0] = 0.0;
  //jc[1] = 0.0;

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

    m = m0; dmax = 0.0; ++niter;

    lb_halo_copy(); /* need up to date currents in halo */

    //fprintf(stderr, "Starting iteration #%d of implicit algorithm...\n", niter);

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

    //fprintf(stderr, "Iteration #%d: dmax = %f\n", niter, dmax);

  } while (dmax > TOLERANCE);

  fprintf(stderr, "Implicit algorithm converged after %d iteration(s).\n", niter);

}

/***********************************************************************/

void mlb_interface_collisions(double *f) {
  int i;
  double *m = f + lblattice.halo_grid_volume*lbmodel.n_vel;
  double rho, cs2, fc, *force = m + 10;
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

inline static void mlb_calc_sigma(double Sigma[][lbmodel.n_dim], double *m) {
  const double gamma = lbpar.gamma;
  int i, j, k;

  double Dr[lbmodel.n_dim],
    Dp[lbmodel.n_dim],
    Dpdu[lbmodel.n_dim],
    Dpmrdp[lbmodel.n_dim],
    Dpmrdpdivu[lbmodel.n_dim],
    Drdpudu[lbmodel.n_dim][lbmodel.n_dim],
    Du[lbmodel.n_dim][lbmodel.n_dim],
    Dj[lbmodel.n_dim][lbmodel.n_dim],
    D2p[lbmodel.n_dim][lbmodel.n_dim],
    D2u[lbmodel.n_dim][lbmodel.n_dim][lbmodel.n_dim],
    divu, divj;

  double *rho   = m;
  double *p     = m + 6;
  double *dp    = m + 7;
  double *pmrdp = m + 8;
  double *rd2p  = m + 9;
  double *u     = m + 12;

  firstDer(Dr, rho);
  firstDer(Dp, p);
  firstDer(Dpmrdp, pmrdp);
  firstDer(Du[0], &u[0]);
  firstDer(Du[1], &u[1]);
  firstDer(Dj[0], &m[1]);
  firstDer(Dj[1], &m[2]);
  secDerAB(D2p, p);
  secDerAB(D2u[0], &u[0]);
  secDerAB(D2u[1], &u[1]);

  divu = Du[0][0] + Du[1][1];
  divj = Dj[0][0] + Dj[1][1];

  /* Step 10 and 11 */
  for (i=0; i<lbmodel.n_dim; ++i) {
    Dpdu[i] = ( Dp[0] * ( Du[0][1] + Du[1][0] )
		+ *p * ( D2u[0][0][i] + D2u[1][1][i]
			 + D2u[i][0][0] + D2u[i][1][1] ) );
    Dpmrdpdivu[i] = Dpmrdp[i]*divu - *pmrdp*(D2u[0][i][0] + D2u[1][i][1]);
    for (j=0; j<lbmodel.n_dim; ++j) {
      Drdpudu[i][j] = ( D2p[i][j]/(*rho) - Dr[i]*Dp[j]/(*rho**rho)
		       + Du[0][i]*Du[j][0] + Du[1][i]*Du[j][1]
		       + u[0]*D2u[j][i][0] + u[1]*D2u[j][i][1] );
      for (k=0; k<lbmodel.n_dim; ++k) {
      }
    }
  }

  /* calculate Sigma */
  for (i=0; i<lbmodel.n_dim; ++i) {
    for (j=0; j<lbmodel.n_dim; ++j) {
      Sigma[i][j] = ( -0.25*(gamma+1.)*(gamma+1.)/(gamma-1.)
		      * ( u[i]*Dpmrdpdivu[j] + u[j]*Dpmrdpdivu[i]
			  + u[i]*Dpdu[j] + u[j]*Dpdu[i] )
		      + (gamma*gamma + 4.*gamma + 1.)/(gamma-1.)/6.
		      * ( *dp * divj * ( Du[i][j] + Du[j][i] )
			  - *p * ( Drdpudu[i][j] + Drdpudu[j][i] ) ) );
    }
    Sigma[i][i] += ( (gamma*gamma + 4.*gamma + 1.)/(gamma-1.)/6.
		     * ( - *pmrdp * ( Drdpudu[0][0] + Drdpudu[1][1] )
			 + *rd2p * divu * divj ) );
  }

}

/***********************************************************************/

inline static void mlb_calc_xi(double Xi[][lbmodel.n_dim][lbmodel.n_dim],
			       double *m) {
  const double gamma = lbpar.gamma;
  int i, j, k;

  double Dr[lbmodel.n_dim],
    Dp[lbmodel.n_dim],
    Dpr[lbmodel.n_dim],
    Du[lbmodel.n_dim][lbmodel.n_dim],
    Dpuu[lbmodel.n_dim][lbmodel.n_dim][lbmodel.n_dim],
    divu;

  double *rho   = m;
  double *p     = m + 6;
  double *pmrdp = m + 8;
  double *u     = m + 12;
  double *jcorr = m + 14;

  firstDer(Dr, rho);
  firstDer(Dp, p);
  firstDer(Du[0], &u[0]);
  firstDer(Du[1], &u[1]);

  divu = Du[0][0] + Du[1][1];

  /* Step 10 and 11 */
  for (i=0; i<lbmodel.n_dim; ++i) {
    Dpr[i] = Dp[i]/(*rho) - *p/(*rho**rho)*Dr[i];
    for (j=0; j<lbmodel.n_dim; ++j) {
      for (k=0; k<lbmodel.n_dim; ++k) {
	Dpuu[i][j][k] = Dp[i]*u[j]*u[k] + *p*Du[j][i]*u[k] + *p*u[j]*Du[k][i];
      }
    }
  }

  /* calculate Xi */
  for (i=0; i<lbmodel.n_dim; ++i) {
    for (j=0; j<lbmodel.n_dim; ++j) {
      for (k=0; k<lbmodel.n_dim; ++k) {
	Xi[i][j][k] = ( (gamma*gamma + 4.*gamma + 1.)/(gamma + 1.)/3.
			* ( Dpuu[k][i][j] + Dpuu[i][j][k] + Dpuu[j][k][i]
			    + delta(i,j) * ( *p*Dpr[k] - *pmrdp*u[k]*divu )
			    + delta(j,k) * ( *p*Dpr[i] - *pmrdp*u[i]*divu )
			    + delta(k,i) * ( *p*Dpr[j] - *pmrdp*u[j]*divu ) )
			+ (1. - gamma) * *p / *rho
			* ( delta(i,j)*jcorr[k]
			    + delta(j,k)*jcorr[i]
			    + delta(k,i)*jcorr[j] ) );
      }
    }
  }

}

/***********************************************************************/

void mlb_correction_collisions(double *f) {
  const double (*c)[lbmodel.n_dim] = lbmodel.c;
  int i, j, k, l;
  double *m = f + lblattice.halo_grid_volume*lbmodel.n_vel;
  double rho, cs2, jc;
  double w[lbmodel.n_vel];
  double Sigma[lbmodel.n_dim][lbmodel.n_dim];
  double Xi[lbmodel.n_dim][lbmodel.n_dim][lbmodel.n_dim];

  double *force = m + 10;
  double *jcorr = m + 14;

  rho = m[0];
  cs2 = eq_state(rho);
  lb_weights(w, cs2);

  m[1] += 0.5*force[0] + jcorr[0];
  m[2] += 0.5*force[1] + jcorr[1];

  mlb_calc_sigma(Sigma, m);
  mlb_calc_xi(Xi, m);

   for (i=0; i<lbmodel.n_vel; ++i) {
     jc = c[i][0]*jcorr[0] + c[i][1]*jcorr[1];
     f[i] += (lbpar.gamma - 1.)*w[i]/cs2*jc;
     for (j=0; j<lbmodel.n_dim; ++j) {
       f[i] -= w[i]/(2.*cs2)*Sigma[j][j];
       for (k=0; k<lbmodel.n_dim; ++k) {
	 f[i] += w[i]/(2.*cs2*cs2)*Sigma[j][k]*(c[i][j]*c[i][k]);
	 for (l=0; l<lbmodel.n_dim; ++l) {
	   f[i] += ( w[i]/(6.*cs2*cs2*cs2) * Xi[j][k][l] * ( c[i][j]*c[i][k]*c[i][l] - cs2 * delta(j,k) * c[i][l]- cs2 * delta(k,l) * c[i][j] - cs2 * delta(l,j) * c[i][k] ) );
	 }
       }
     }
   }

}

/***********************************************************************/
