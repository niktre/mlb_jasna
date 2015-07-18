/***********************************************************************
 *
 * d2q21.h
 *
 * Copyright (c) 2009-2015 Ulf D. Schiller <ulf@lattice-boltzmann.de>
 * All rights reserved.
 *
 * Please cite the following publications when using this code:
 *
 * U. D. Schiller 
 * "A unified operator splitting approach for multi-scale
 * fluid–particle coupling in the lattice Boltzmann method"
 * Comp. Phys. Comm. 185, 2586-2597 (2014)
 * http://dx.doi.org/10.1016/j.cpc.2014.06.005
 * 
 * B. Duenweg, U. D. Schiller, A. J. C. Ladd
 * "Statistical mechanics of the fluctuating lattice Boltzmann equation"
 * Phys. Rev. E 76, 036704 (2007)
 * http://dx.doi.org/10.1103/PhysRevE.76.036704
 * 
 ***********************************************************************/

#define NDIM 2
#define NVEL 21
#define HALO 4
#define WGRID (2*HALO+1)

/***********************************************************************/

#define VELS(D,Q)    d##D##q##Q##_velocities
#define WEIGHTS(D,Q) d##D##q##Q##_weights
#define NORMS(D,Q)   d##D##q##Q##_norms
/* these macros are necessary to force prescan in the macro below */
/* because prescan does not occur for stringify and concat */
#define DnQm(D,Q)  { D, Q, Q, VELS(D,Q) }

/***********************************************************************/

typedef struct _LBpar {
  double rho;
  double gamma;
} LB_Parameters;

typedef const struct _LBmodel {
  const int n_dim;
  const int n_vel;
  const int n_mom;
  const double (*c)[NDIM];
} LB_Model;

typedef struct _Lattice {
  int grid[NDIM];
  int halo_grid[NDIM];
  int halo_size[NDIM];
  int stride[NDIM];
  int halo_grid_volume;
} LB_Lattice;

/***********************************************************************/

static const double d2q21_velocities[21][2] = { {  0.,  0. },
						{  1.,  0. },
						{ -1.,  0. },
						{  0.,  1. },
						{  0., -1. },
						{  1.,  1. },
						{ -1., -1. },
						{ -1.,  1. },
						{  1., -1. },
						{  2.,  0. },
						{ -2.,  0. },
						{  0.,  2. },
						{  0., -2. },
						{  2.,  2. },
						{ -2., -2. },
						{ -2.,  2. },
						{  2., -2. },
						{  4.,  0. },
						{ -4.,  0. },
						{  0.,  4. },
						{  0., -4. } };

/***********************************************************************/

void lb_weights(double *w, double sigma2) {
  int i;

  w[ 0] = 1.0 - 45./2.*sigma2*(7./60. - 7./48.*sigma2 + sigma2*sigma2/16.);
  w[ 1] = sigma2/3.*(32./15. - 4.*sigma2 + 2.*sigma2*sigma2);
  w[ 5] = sigma2*(sigma2/3. - sigma2*sigma2/4.);
  w[ 9] = sigma2*(-1./18. + 3./16.*sigma2 - sigma2*sigma2/12.);
  w[13] = sigma2/96.*(-sigma2/2. + 3./2.*sigma2*sigma2);
  w[17] = sigma2/384.*(4./15. - sigma2 + sigma2*sigma2);

  for(i=1; i<4; ++i) {
    w[ 1+i] = w[1];
    w[ 5+i] = w[5];
    w[ 9+i] = w[9];
    w[13+i] = w[13];
    w[17+i] = w[17];
  }

}
