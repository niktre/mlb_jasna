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

/***********************************************************************/

#define VELS(D,Q)    d##D##q##Q##_velocities
#define WEIGHTS(D,Q) d##D##q##Q##_weights
#define NORMS(D,Q)   d##D##q##Q##_norms
/* these macros are necessary to force prescan in the macro below */
/* because prescan does not occur for stringify and concat */
#define DnQm(D,Q)  { D, Q, Q, VELS(D,Q) } // WEIGHTS(D,Q), NORMS(D,Q) }

/***********************************************************************/

typedef struct _LBpar {
  double rho;
  double viscosity;
  double mu;
  double grad_p;
  double ext_force[NDIM];
} LB_Parameters;

typedef const struct _LBmodel {
  const int n_dim;
  const int n_vel;
  const int n_mom;
  const double (*c)[NDIM];
  const double *w;
  const double *b;
} LB_Model;

typedef struct _Lattice {
  int grid[NDIM];
  int halo_grid[NDIM];
  int halo_size[NDIM];
  int stride[NDIM];
  int halo_grid_volume;
} LB_Lattice;

/***********************************************************************/

static const double d3q19_velocities[19][3] = { {  0.,  0.,  0. },
						{  1.,  0.,  0. },
						{ -1.,  0.,  0. },
						{  0.,  1.,  0. },
						{  0., -1.,  0. },
						{  0.,  0.,  1. },
						{  0.,  0., -1. },
						{  1.,  1.,  0. },
						{ -1., -1.,  0. },
						{  1., -1.,  0. },
						{ -1.,  1.,  0. },
						{  1.,  0.,  1. },
						{ -1.,  0., -1. },
						{  1.,  0., -1. },
						{ -1.,  0.,  1. },
						{  0.,  1.,  1. },
						{  0., -1., -1. },
						{  0.,  1., -1. },
						{  0., -1.,  1. } };

static const double d3q19_weights[19] = { 1./3.,
					  1./18., 1./18., 1./18.,
					  1./18., 1./18., 1./18.,
					  1./36., 1./36., 1./36.,
					  1./36., 1./36., 1./36.,
					  1./36., 1./36., 1./36.,
					  1./36., 1./36., 1./36. };

static const double d3q19_norms[19] = { 1.0, 
					1./3., 1./3., 1./3., 
					2./3., 4./9., 4./3., 
					1./9., 1./9., 1./9., 
					2./3., 2./3., 2./3., 
					2./9., 2./9., 2./9., 
					2.0, 4./9., 4./3. };

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
