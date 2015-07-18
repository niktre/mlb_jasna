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

typedef double vector[3];

typedef struct _LBpar {
  double rho;
  double viscosity;
  double mu;
  double grad_p;
  double ext_force[3];
} LB_Parameters;

typedef const struct _LBmodel {
  const int n_vel;
  const int n_mom;
  const vector *c;
  const double *b;
} LB_Model;

typedef struct _Lattice {
  int grid[3];
  int halo_grid[3];
  int halo_size[3];
  int stride[3];
  int halo_grid_volume;
} LB_Lattice;

/***********************************************************************/

static const double d3q19_velocities[19][3] = { {  1,  0,  0 },
						{ -1,  0,  0 },
						{  0,  1,  0 },
						{  0, -1,  0 },
						{  0,  0,  1 },
						{  0,  0, -1 },
						{  1,  1,  0 },
						{ -1, -1,  0 },
						{  1, -1,  0 },
						{ -1,  1,  0 },
						{  1,  0,  1 },
						{ -1,  0, -1 },
						{  1,  0, -1 },
						{ -1,  0,  1 },
						{  0,  1,  1 },
						{  0, -1, -1 },
						{  0,  1, -1 },
						{  0, -1,  1 } };

static const double d3q19_norms[19] = { 1.0, 
					1./3., 1./3., 1./3., 
					2./3., 4./9., 4./3., 
					1./9., 1./9., 1./9., 
					2./3., 2./3., 2./3., 
					2./9., 2./9., 2./9., 
					2.0, 4./9., 4./3. };

/***********************************************************************/
