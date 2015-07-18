/***********************************************************************
 *
 * d2q21.c
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

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "d2q21.h"

/***********************************************************************/

#define NVEL  19

#define D3Q19 { NVEL, NVEL, d3q19_velocities, d3q19_norms }

/* temporary/secondary grid */
#define MAX_Y 102
#define PFI *FI(NVEL,3,MAX_Y+2)
#define FI(i,x,y) fi[y][i][x]

/***********************************************************************/

static const LB_Model lbmodel = D3Q19;

static LB_Lattice lblattice = { {0,0,0}, {0,0,0}, {1,1,1}, {0,0,0}, 0 };

static LB_Parameters lbpar = { 1.0, 1./12., 0.0, 0.0, { 0.0, 0.0, 0.0} };

static double *lbmom = NULL;
static double PFI;

static double gamma_even;
static double gamma_odd;

static double lb_phi[NVEL];

static int fluct = 0;

/***********************************************************************/

void lb_halo_copy() {

  int x, y;
  int lsrc, rsrc, ldst, rdst, size;

  int xstride = lbmodel.n_mom*lblattice.stride[0];
  int ystride = lbmodel.n_mom*lblattice.stride[1];
  int zstride = lbmodel.n_mom*lblattice.stride[2];

  /***************
   * X direction *
   ***************/

  size = lblattice.halo_size[0]*xstride;

  ldst = 0;
  lsrc = size;
  rsrc = lblattice.grid[0]*xstride;
  rdst = rsrc + size;

  memcpy(lbmom+ldst, lbmom+rsrc, size*sizeof(*lbmom));
  memcpy(lbmom+rdst, lbmom+lsrc, size*sizeof(*lbmom));

  /***************
   * Y direction *
   ***************/

  size = lblattice.halo_size[1]*ystride;

  ldst = 0;
  lsrc = size;
  rsrc = lblattice.grid[1]*ystride;
  rdst = rsrc + size;

  for (x=0; x<lblattice.halo_grid[0]; x++) {
    memcpy(lbmom+ldst, lbmom+rsrc, size*sizeof(*lbmom));
    memcpy(lbmom+rdst, lbmom+lsrc, size*sizeof(*lbmom));

    ldst += xstride;
    lsrc += xstride;
    rsrc += xstride;
    rdst += xstride;
  }

#if 0
  /***************
   * Z direction *
   ***************/

  size = lblattice.halo_size[2]*zstride;

  ldst = 0;
  lsrc = size;
  rsrc = lblattice.grid[2]*zstride;
  rdst = rsrc + size;

  for (x=0; x<lblattice.halo_grid[0]; x++) {
    for (y=0; y<lblattice.halo_grid[1]; y++) {
      memcpy(lbmom+ldst, lbmom+rsrc, size);
      memcpy(lbmom+rdst, lbmom+lsrc, size);

      ldst += ystride;
      lsrc += ystride;
      rsrc += ystride;
      rdst += ystride;
    }
  }
#endif

}

/***********************************************************************/

static void lb_inlet_x(double PFI, int x, int y, int z) {  
  double rho, f, m0, my, mz, jx, jy, jz, nxy, nxz;

  int xc, yc;
  xc = x%3;
  yc = y+1;

  rho = 3.*lbpar.grad_p*lblattice.grid[0];

  f = FI( 0, xc, yc)[z]; m0  = f;  
  f = FI( 2, xc, yc)[z]; m0 += 2*f; 
  f = FI( 3, xc, yc)[z]; m0 += f;   my  = f;
  f = FI( 4, xc, yc)[z]; m0 += f;   my -= f;
  f = FI( 5, xc, yc)[z]; m0 += f;            mz  = f;
  f = FI( 6, xc, yc)[z]; m0 += f;            mz -= f;
  f = FI( 8, xc, yc)[z]; m0 += 2*f; 
  f = FI(10, xc, yc)[z]; m0 += 2*f;
  f = FI(12, xc, yc)[z]; m0 += 2*f;
  f = FI(14, xc, yc)[z]; m0 += 2*f;
  f = FI(15, xc, yc)[z]; m0 += f;   my += f; mz += f;
  f = FI(16, xc, yc)[z]; m0 += f;   my -= f; mz -= f;
  f = FI(17, xc, yc)[z]; m0 += f;   my += f; mz -= f;
  f = FI(18, xc, yc)[z]; m0 += f;   my -= f; mz += f;

  jx = (rho - m0)/6.;
  jy = 0.0;
  jz = 0.0;

  nxy = 2.*jy - 0.5*my;
  nxz = 2.*jz - 0.5*mz;
  
  FI( 1, xc, yc)[z] = FI( 2, xc, yc)[z] + 2.*jx;
  FI( 7, xc, yc)[z] = FI( 8, xc, yc)[z] + jx + jy + nxy;
  FI( 9, xc, yc)[z] = FI(10, xc, yc)[z] + jx - jy - nxy;
  FI(11, xc, yc)[z] = FI(12, xc, yc)[z] + jx + jz + nxz;
  FI(13, xc, yc)[z] = FI(14, xc, yc)[z] + jx - jz - nxz;

}

/***********************************************************************/

static void lb_outlet_x(double PFI, int x, int y, int z) {  
  double rho, f, m0, my, mz, jx, jy, jz, nxy, nxz;

  int xc, yc;
  xc = x%3;
  yc = y+1;

  rho = 0.0;

  f = FI( 0, xc, yc)[z]; m0  = f;  
  f = FI( 1, xc, yc)[z]; m0 += 2*f; 
  f = FI( 3, xc, yc)[z]; m0 += f;   my  = f;
  f = FI( 4, xc, yc)[z]; m0 += f;   my -= f;
  f = FI( 5, xc, yc)[z]; m0 += f;            mz  = f;
  f = FI( 6, xc, yc)[z]; m0 += f;            mz -= f;
  f = FI( 7, xc, yc)[z]; m0 += 2*f; 
  f = FI( 9, xc, yc)[z]; m0 += 2*f;
  f = FI(11, xc, yc)[z]; m0 += 2*f;
  f = FI(13, xc, yc)[z]; m0 += 2*f;
  f = FI(15, xc, yc)[z]; m0 += f;   my += f; mz += f;
  f = FI(16, xc, yc)[z]; m0 += f;   my -= f; mz -= f;
  f = FI(17, xc, yc)[z]; m0 += f;   my += f; mz -= f;
  f = FI(18, xc, yc)[z]; m0 += f;   my -= f; mz += f;

  jx = (m0 - rho)/6.;
  jy = 0.0;
  jz = 0.0;

  nxy = 2.*jy - 0.5*my;
  nxz = 2.*jz - 0.5*mz;
  
  FI( 2, xc, yc)[z] = FI( 1, xc, yc)[z] - 2.*jx;
  FI( 8, xc, yc)[z] = FI( 7, xc, yc)[z] - jx - jy - nxy;
  FI(10, xc, yc)[z] = FI( 9, xc, yc)[z] - jx + jy + nxy;
  FI(12, xc, yc)[z] = FI(11, xc, yc)[z] - jx - jz - nxz;
  FI(14, xc, yc)[z] = FI(13, xc, yc)[z] - jx + jz + nxz;

}

/***********************************************************************/

static void lb_bounce_back_read(double PFI, int x, int y, int z) {
  int xc, xm, yc, yp, ym, zp, zm;

  xc = x%3; xm = (xc+2)%3;
  yc = y+1; yp = yc+1; ym = yc-1;
  zp = z+1; zm = z-1;

  FI( 2, xm, yc)[z]  = FI( 1, xc, yc)[z];
  FI( 4, xc, ym)[z]  = FI( 3, xc, yc)[z];
  FI( 6, xc, yc)[zm] = FI( 5, xc, yc)[z];
  FI( 8, xm, ym)[z]  = FI( 7, xc, yc)[z];
  FI(10, xm, yp)[z]  = FI( 9, xc, yc)[z];
  FI(12, xm, yc)[zm] = FI(11, xc, yc)[z];
  FI(14, xm, yc)[zp] = FI(13, xc, yc)[z];
  FI(16, xc, ym)[zm] = FI(15, xc, yc)[z];
  FI(18, xc, ym)[zp] = FI(17, xc, yc)[z];

}

/***********************************************************************/

static void lb_bounce_back_write(double PFI, int x, int y, int z) {
  int xc, xp, yc, yp, ym, zp, zm;

  xc = x%3; xp = (xc+1)%3;
  yc = y+1; yp = yc+1; ym = yc-1;
  zp = z+1; zm = z-1;

  FI( 1, xp, yc)[z]  = FI( 2, xc, yc)[z];
  FI( 3, xc, yp)[z]  = FI( 4, xc, yc)[z];
  FI( 5, xc, yc)[zp] = FI( 6, xc, yc)[z];
  FI( 7, xp, yp)[z]  = FI( 8, xc, yc)[z];
  FI( 9, xp, ym)[z]  = FI(10, xc, yc)[z];
  FI(11, xp, yc)[zp] = FI(12, xc, yc)[z];
  FI(13, xp, yc)[zm] = FI(14, xc, yc)[z];
  FI(15, xc, yp)[zp] = FI(16, xc, yc)[z];
  FI(17, xc, yp)[zm] = FI(18, xc, yc)[z];

}

/***********************************************************************/

void lb_add_forces(double *m) {
  double rho, u[3], f[3], C[6];

  /* halfstep force */
  f[0] = 0.5*lbpar.ext_force[0];
  f[1] = 0.5*lbpar.ext_force[1];
  f[2] = 0.5*lbpar.ext_force[2];

  /* mass density */
  rho = m[0] + lbpar.rho;

  /* hydrodynamic velocity */
  u[0] = m[1]/rho;
  u[1] = m[2]/rho;
  u[2] = m[3]/rho;

  C[0] = 2.0*u[0]*f[0];
  C[2] = 2.0*u[1]*f[1];
  C[5] = 2.0*u[2]*f[2];
  C[1] = u[0]*f[1] + f[0]*u[1];
  C[3] = u[0]*f[2] + f[0]*u[2];
  C[4] = u[1]*f[2] + f[1]*u[2];

  /* update momentum density */
  m[1] += f[0];
  m[2] += f[1];
  m[3] += f[2];

  /* update stress modes */
  m[4] += C[0] + C[2] + C[5];
  m[5] += C[0] - C[2];
  m[6] += C[0] + C[2] - 2.*C[5];
  m[7] += C[1];
  m[8] += C[3];
  m[9] += C[4];

  /* kinetic modes have no contribution in the basis used here */

}

/***********************************************************************/

static void lb_relax_modes(double *m) {

  double rho, u[3], m4eq, m5eq, m6eq, m7eq, m8eq, m9eq;

  /* mass density */
  rho = m[0] + lbpar.rho;

  /* momentum density */
  u[0] = m[1]/rho;
  u[1] = m[2]/rho;
  u[2] = m[3]/rho;

  /* equilibrium part of the stress modes */
  m4eq = (u[0]*u[0] + u[1]*u[1] + u[2]*u[2])*rho;
  m5eq = (u[0]*u[0] - u[1]*u[1])*rho;
  m6eq = (u[0]*u[0] + u[1]*u[1] - 2.*u[2]*u[2])*rho;
  m7eq = u[0]*u[1]*rho;
  m8eq = u[0]*u[2]*rho;
  m9eq = u[1]*u[2]*rho;
  
  /* relax stress modes */  
  m[4] = m4eq + gamma_even * (m[4] - m4eq);
  m[5] = m5eq + gamma_even * (m[5] - m5eq);
  m[6] = m6eq + gamma_even * (m[6] - m6eq);
  m[7] = m7eq + gamma_even * (m[7] - m7eq);
  m[8] = m8eq + gamma_even * (m[8] - m8eq);
  m[9] = m9eq + gamma_even * (m[9] - m9eq);
  
  /* ghost modes have no equilibrium part due to orthogonality */
  m[10] = gamma_odd  * m[10];
  m[11] = gamma_odd  * m[11];
  m[12] = gamma_odd  * m[12];
  m[13] = gamma_odd  * m[13];
  m[14] = gamma_odd  * m[14];
  m[15] = gamma_odd  * m[15];
  m[16] = gamma_even * m[16];
  m[17] = gamma_even * m[17];
  m[18] = gamma_even * m[18];
  
}

/***********************************************************************/

void lb_thermalize_modes(double *m) {
  
  double rootrho = sqrt(m[0]+lbpar.rho);

  /* stress modes */
  m[4] += rootrho*lb_phi[4]*gaussian_random();
  m[5] += rootrho*lb_phi[5]*gaussian_random();
  m[6] += rootrho*lb_phi[6]*gaussian_random();
  m[7] += rootrho*lb_phi[7]*gaussian_random();
  m[8] += rootrho*lb_phi[8]*gaussian_random();
  m[9] += rootrho*lb_phi[9]*gaussian_random();
    
  /* ghost modes */
  m[10] += rootrho*lb_phi[10]*gaussian_random();
  m[11] += rootrho*lb_phi[11]*gaussian_random();
  m[12] += rootrho*lb_phi[12]*gaussian_random();
  m[13] += rootrho*lb_phi[13]*gaussian_random();
  m[14] += rootrho*lb_phi[14]*gaussian_random();
  m[15] += rootrho*lb_phi[15]*gaussian_random();
  m[16] += rootrho*lb_phi[16]*gaussian_random();
  m[17] += rootrho*lb_phi[17]*gaussian_random();
  m[18] += rootrho*lb_phi[18]*gaussian_random();

}

/***********************************************************************/

static void lb_collisions(double *m) {

  lb_relax_modes(m);
  if (fluct) lb_thermalize_modes(m);

}

/***********************************************************************/

static void lb_calc_modes(double *m, double PFI, int x, int y, int z) {
  int xc, yc;
  double f, mc0, mc1, mc2;
  double mx1, my1, mz1, mx2, my2, mz2, mx3, my3, mz3;
  double mxx1, myy1, mzz1, mxx2, myy2, mzz2;
  double mxy, mxz, myz;

  xc = x%3;
  yc = y+1;

  f = FI( 0, xc, yc)[z]; mc0  = f;
  f = FI( 1, xc, yc)[z]; mx1  = f; mxx1  = f;
  f = FI( 2, xc, yc)[z]; mx1 -= f; mxx1 += f;
  f = FI( 3, xc, yc)[z]; my1  = f; myy1  = f;
  f = FI( 4, xc, yc)[z]; my1 -= f; myy1 += f;
  f = FI( 5, xc, yc)[z]; mz1  = f; mzz1  = f;
  f = FI( 6, xc, yc)[z]; mz1 -= f; mzz1 += f;
  f = FI( 7, xc, yc)[z]; mx2  = f; my3  = f; mxy  = f; mzz2  = f;
  f = FI( 8, xc, yc)[z]; mx2 -= f; my3 -= f; mxy += f; mzz2 += f;
  f = FI( 9, xc, yc)[z]; mx2 += f; my3 -= f; mxy -= f; mzz2 += f;
  f = FI(10, xc, yc)[z]; mx2 -= f; my3 += f; mxy -= f; mzz2 += f;
  f = FI(11, xc, yc)[z]; mz2  = f; mx3  = f; mxz  = f; mxx2  = f;
  f = FI(12, xc, yc)[z]; mz2 -= f; mx3 -= f; mxz += f; mxx2 += f;
  f = FI(13, xc, yc)[z]; mz2 -= f; mx3 += f; mxz -= f; mxx2 += f;
  f = FI(14, xc, yc)[z]; mz2 += f; mx3 -= f; mxz -= f; mxx2 += f;
  f = FI(15, xc, yc)[z]; my2  = f; mz3  = f; myz  = f; myy2  = f;
  f = FI(16, xc, yc)[z]; my2 -= f; mz3 -= f; myz += f; myy2 += f;
  f = FI(17, xc, yc)[z]; my2 += f; mz3 -= f; myz -= f; myy2 += f;
  f = FI(18, xc, yc)[z]; my2 -= f; mz3 += f; myz -= f; myy2 += f;

  mc1 = mxx1 + myy1 + mzz1;
  mc2 = mxx2 + myy2 + mzz2;

  m[ 0] = mc0 + mc1 + mc2;
  m[ 1] = mx1 + mx2 + mx3;
  m[ 2] = my1 + my2 + my3;
  m[ 3] = mz1 + mz2 + mz3;
  m[ 4] = mc2 - mc0;
  m[ 5] = mxx1 - myy1 + mxx2 - myy2;
  m[ 6] = mc1 - 3.*mzz1 - mc2 + 3.*mzz2;
  m[ 7] = mxy;
  m[ 8] = mxz;
  m[ 9] = myz;
  m[10] = m[1] - 3.*mx1;
  m[11] = m[2] - 3.*my1;
  m[12] = m[3] - 3.*mz1;
  m[13] = mx2 - mx3;
  m[14] = my3 - my2;
  m[15] = mz2 - mz3;
  m[16] = m[0] - 3.*mc1;
  m[17] = myy1 - mxx1 + mxx2 - myy2;
  m[18] = 3.*mzz1 - mc1 + 3.*mzz2 - mc2;

}

/***********************************************************************/

static void lb_calc_fi_stream(double *m, double PFI, int x, int y, int z) {
  int xc, xp, xm, yc, yp, ym, zp, zm;
  double mc0, mc1, mc2;
  double mx1, my1, mz1, mx2, my2, mz2, mx3, my3, mz3;
  double mxx1, myy1, mzz1, mxy, mxz, myz, mxy2, mxz2, myz2;

  xc = x%3; xp = (xc+1)%3; xm = (xc+2)%3;
  yc = y+1; yp = yc+1; ym = yc-1;
  zp = z+1; zm = z-1;

  m[ 0] /= 36.;
  m[ 1] /= 12.;
  m[ 2] /= 12.;
  m[ 3] /= 12.;
  m[ 4] /= 24.;
  m[ 5] /= 16.;
  m[ 6] /= 48.;
  m[ 7] /= 4.;
  m[ 8] /= 4.;
  m[ 9] /= 4.;
  m[10] /= 24.;
  m[11] /= 24.;
  m[12] /= 24.;
  m[13] /= 8.;
  m[14] /= 8.;
  m[15] /= 8.;
  m[16] /= 72.;
  m[17] /= 16.;
  m[18] /= 48.;

  mc0 = 12.*(m[0] - m[4] + m[16]);
  mc1 =  2.*(m[0] - 2.*m[16]);
  mc2 = m[0] + m[4] + m[16];

  mx1 = 2.*(m[1] - 2.*m[10]);
  my1 = 2.*(m[2] - 2.*m[11]);
  mz1 = 2.*(m[3] - 2.*m[12]);

  mx2 = m[1] + m[10] + m[13];
  my2 = m[2] + m[11] - m[14];
  mz2 = m[3] + m[12] + m[15];

  mx3 = m[1] + m[10] - m[13];
  my3 = m[2] + m[11] + m[14];
  mz3 = m[3] + m[12] - m[15];

  mxx1 = mc1 + 2.*(m[5] + m[6]) - 2.*(m[17] + m[18]);
  myy1 = mc1 - 2.*(m[5] - m[6]) + 2.*(m[17] - m[18]);
  mzz1 = mc1 - 4.*(m[6] - m[18]);
  
  mxy2 = mc2 + 2.*(m[6] + m[18]);
  mxz2 = mc2 + (m[5] + m[17]) - (m[6] + m[18]);
  myz2 = mc2 - (m[5] + m[17]) - (m[6] + m[18]);

  mxy = m[7];
  mxz = m[8];
  myz = m[9];

  FI( 0, xc, yc)[z]  = mc0;
  FI( 1, xp, yc)[z]  = mxx1 + mx1;
  FI( 2, xm, yc)[z]  = mxx1 - mx1;
  FI( 3, xc, yp)[z]  = myy1 + my1;
  FI( 4, xc, ym)[z]  = myy1 - my1;
  FI( 5, xc, yc)[zp] = mzz1 + mz1;
  FI( 6, xc, yc)[zm] = mzz1 - mz1;
  FI( 7, xp, yp)[z]  = mxy2 + mx2 + my3 + mxy;
  FI( 8, xm, ym)[z]  = mxy2 - mx2 - my3 + mxy;
  FI( 9, xp, ym)[z]  = mxy2 + mx2 - my3 - mxy;
  FI(10, xm, yp)[z]  = mxy2 - mx2 + my3 - mxy;
  FI(11, xp, yc)[zp] = mxz2 + mz2 + mx3 + mxz;
  FI(12, xm, yc)[zm] = mxz2 - mz2 - mx3 + mxz;
  FI(13, xp, yc)[zm] = mxz2 - mz2 + mx3 - mxz;
  FI(14, xm, yc)[zp] = mxz2 + mz2 - mx3 - mxz;
  FI(15, xc, yp)[zp] = myz2 + my2 + mz3 + myz;
  FI(16, xc, ym)[zm] = myz2 - my2 - mz3 + myz;
  FI(17, xc, yp)[zm] = myz2 + my2 - mz3 - myz;
  FI(18, xc, ym)[zp] = myz2 - my2 + mz3 - myz;

}

/***********************************************************************/

static void lb_read(double *m, double PFI, int x, int y, int z) {

  lb_add_forces(m);
  lb_collisions(m);
  lb_calc_fi_stream(m, fi, x, y, z);    

}

/***********************************************************************/

static void lb_write(double *m, double PFI, int x, int y, int z) {

  lb_calc_modes(m, fi, x, y, z);
  lb_collisions(m);
  lb_add_forces(m);

}

/***********************************************************************/

static void lb_read_column(double *m, double PFI, int x, int y) {
  int z, zl, zh;

  zl = lblattice.halo_size[2];
  zh = lblattice.halo_size[2] + lblattice.grid[2] - 1;

  for (z=zl, m+=zl*lbmodel.n_mom; z<=zh; ++z, m+=lbmodel.n_mom) {
    lb_read(m, fi, x, y, z);
  }

  /* boundary conditions in z */
  int xc, xp, xm, yc, yp, ym;
  xc = x%3; xp = (xc+1)%3; xm = (xc+2)%3;
  yc = y+1; yp =  yc+1;    ym =  yc-1;
#ifdef ZPERIODIC
  FI( 5, xc, yc)[zl] = FI( 5, xc, yc)[zh+1];
  FI(11, xp, yc)[zl] = FI(11, xp, yc)[zh+1];
  FI(14, xm, yc)[zl] = FI(14, xm, yc)[zh+1];
  FI(15, xc, yp)[zl] = FI(15, xc, yp)[zh+1];
  FI(18, xc, ym)[zl] = FI(18, xc, ym)[zh+1];

  FI( 6, xc, yc)[zh] = FI( 6, xc, yc)[zl-1];
  FI(12, xm, yc)[zh] = FI(12, xm, yc)[zl-1];
  FI(13, xp, yc)[zh] = FI(13, xp, yc)[zl-1];
  FI(16, xc, ym)[zh] = FI(16, xc, ym)[zl-1];
  FI(17, xc, yp)[zh] = FI(17, xc, yp)[zl-1];
#else /* bounce back */
  FI( 5, xc, yc)[zl] = FI( 6, xc, yc)[zl-1];
  FI(11, xc, yc)[zl] = FI(12, xm, yc)[zl-1];
  FI(14, xc, yc)[zl] = FI(13, xp, yc)[zl-1];
  FI(15, xc, yc)[zl] = FI(16, xc, ym)[zl-1];
  FI(18, xc, yc)[zl] = FI(17, xc, yp)[zl-1];

  FI( 6, xc, yc)[zh] = FI( 5, xc, yc)[zh+1];
  FI(12, xc, yc)[zh] = FI(11, xp, yc)[zh+1];
  FI(13, xc, yc)[zh] = FI(14, xm, yc)[zh+1];
  FI(16, xc, yc)[zh] = FI(15, xc, yp)[zh+1];
  FI(17, xc, yc)[zh] = FI(18, xc, ym)[zh+1];
#endif

}

/***********************************************************************/

static void lb_write_column(double *m, double PFI, int x, int y) {
  int z, zl, zh;

  zl = lblattice.halo_size[2];
  zh = lblattice.halo_size[2] + lblattice.grid[2] - 1;

  for (z=zl, m+=zl*lbmodel.n_mom; z<=zh; ++z, m+=lbmodel.n_mom) {
    lb_write(m, fi, x, y, z);
  }
  
}

/***********************************************************************/

void lb_update() {
  int x, y;
  int ystride = lblattice.stride[1]*lbmodel.n_mom;
  int ioff    = lblattice.stride[0] + lblattice.stride[1];
  int moff    = ioff*lbmodel.n_mom;
  double *m;

  lb_halo_copy(); /* need up to date moments in halo */

  m = lbmom;

  for (y=0; y<lblattice.halo_grid[1]; ++y) {
    lb_read_column(m, fi, 0, y); 
    m += ystride;
  }

  for (x=1; x<lblattice.halo_grid[0]; ++x) {
    
    lb_read_column(m, fi, x, 0);
    m += ystride;

    for (y=1; y<lblattice.halo_grid[1]; ++y) {
      lb_read_column(m, fi, x, y);
      lb_write_column(m-moff, fi, x-1, y-1);
      m += ystride;
    }
  
  }

}

/***********************************************************************/

static void lb_init_fluid() {
  int x, y, z, zl, zh;
  double *m = lbmom;

  zl = lblattice.halo_size[2] - 1;
  zh = lblattice.halo_size[2] + lblattice.grid[2];

  for (x=0; x<lblattice.halo_grid[0]; ++x) {
    for (y=0; y<lblattice.halo_grid[1]; ++y) {
      for (z=0; z<lblattice.halo_grid[2]; ++z, m+=lbmodel.n_mom) {
	if (z==zl || z==zh) {
	  m[0] = - lbpar.rho;
	} else {
	  m[0] = 0.0;
	}
      }
    }
  }

}

/***********************************************************************/

static void lb_init_lattice(int *grid) {
  int i, x, y, hgrid[3], hsize[3];

  lblattice.grid[0] = grid[0];
  lblattice.grid[1] = grid[1];
  lblattice.grid[2] = grid[2];

  lblattice.halo_size[0] = hsize[0] = 1;
  lblattice.halo_size[1] = hsize[1] = 1;
  lblattice.halo_size[2] = hsize[2] = 1;

  lblattice.halo_grid[0] = hgrid[0] = lblattice.grid[0] + 2*hsize[0];
  lblattice.halo_grid[1] = hgrid[1] = lblattice.grid[1] + 2*hsize[1];
  lblattice.halo_grid[2] = hgrid[2] = lblattice.grid[2] + 2*hsize[2];

  lblattice.stride[2] = 1;
  lblattice.stride[1] = hgrid[2];
  lblattice.stride[0] = hgrid[2]*hgrid[1];

  lblattice.halo_grid_volume = hgrid[0]*hgrid[1]*hgrid[2];

  if (lblattice.halo_grid[1] > MAX_Y) {
    fprintf(stderr,"y grid size too large (%d > %d)!\n",lblattice.halo_grid[1],MAX_Y);
    exit(EXIT_FAILURE);
  }

  lbmom   = calloc(lblattice.halo_grid_volume*lbmodel.n_mom, sizeof(*lbmom));

  for (i=0; i<lbmodel.n_vel; ++i) {
    for (x=0; x<3; ++x) {
      for (y=0; y<lblattice.halo_grid[1]+2; ++y) {
	FI(i,x,y) = calloc(lblattice.halo_grid[2],sizeof(*FI(0,0,0)));
      }
    }
  }

}

/***********************************************************************/

static void lb_init_parameters(double rho, double viscosity, double temperature,
			       double grad_p, double *force) {

  const double *b = lbmodel.b;

  lbpar.rho       = rho;
  lbpar.viscosity = viscosity; /* kinematic viscosity */

  if (force) {
    lbpar.ext_force[0]  = force[0];
    lbpar.ext_force[1]  = force[1];
    lbpar.ext_force[2]  = force[2];
  } else {
    lbpar.ext_force[0] = 0.0;
    lbpar.ext_force[1] = 0.0;
    lbpar.ext_force[2] = 0.0;
  }

  if (grad_p > 0.0) {
    lbpar.grad_p = grad_p;
  } else {
    lbpar.grad_p = 0.0;
  }
  
  /* Eq. (33) Schiller, CPC 185, 2586-2597 (2014) */
  /* gamma_even = 1. - 1./(6.*lbpar.viscosity+0.5); */

  /* exact square root allows viscosities down to 1/6 */
  gamma_even = sqrt(1. - 1./(3.*lbpar.viscosity+0.5));

  /* Ginzburg and d'Humieres, PRE 68, 066614 (2003) */
  /* gamma_odd  = 1. - 8.*(1. + gamma_even)/(7. + gamma_even); */

  /* Note: this needs yet to be verified for collide-stream-collide */
  gamma_odd = 1. - 4.*(1. + gamma_even*gamma_even)/(7. + gamma_even*gamma_even);

  if (temperature > 0.0) {

    fluct = 1;

    /* Eq. (51) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007).
     * Note that the modes are not normalized here! */
    lbpar.mu = 3.0*temperature;
    
    lb_phi[0] = 0.0;
    lb_phi[1] = 0.0;
    lb_phi[2] = 0.0;
    lb_phi[3] = 0.0;
    lb_phi[4] = sqrt(lbpar.mu*b[4]*(1.-gamma_even*gamma_even));
    lb_phi[5] = sqrt(lbpar.mu*b[5]*(1.-gamma_even*gamma_even));
    lb_phi[6] = sqrt(lbpar.mu*b[6]*(1.-gamma_even*gamma_even));
    lb_phi[7] = sqrt(lbpar.mu*b[7]*(1.-gamma_even*gamma_even));
    lb_phi[8] = sqrt(lbpar.mu*b[8]*(1.-gamma_even*gamma_even));
    lb_phi[9] = sqrt(lbpar.mu*b[9]*(1.-gamma_even*gamma_even));

    lb_phi[10] = sqrt(lbpar.mu*b[10]*(1.-gamma_odd*gamma_odd));
    lb_phi[11] = sqrt(lbpar.mu*b[11]*(1.-gamma_odd*gamma_odd));
    lb_phi[12] = sqrt(lbpar.mu*b[12]*(1.-gamma_odd*gamma_odd));
    lb_phi[13] = sqrt(lbpar.mu*b[13]*(1.-gamma_odd*gamma_odd));
    lb_phi[14] = sqrt(lbpar.mu*b[14]*(1.-gamma_odd*gamma_odd));
    lb_phi[15] = sqrt(lbpar.mu*b[15]*(1.-gamma_odd*gamma_odd));
    
    lb_phi[16] = sqrt(lbpar.mu*b[16]*(1.-gamma_even*gamma_even));
    lb_phi[17] = sqrt(lbpar.mu*b[17]*(1.-gamma_even*gamma_even));
    lb_phi[18] = sqrt(lbpar.mu*b[18]*(1.-gamma_even*gamma_even));

  } else {

    fluct = 0;

  }
 
}

/***********************************************************************/

void lb_init(int *grid, double rho, double viscosity, double temperature,
	     double grad_p, double *force) {

  lb_init_parameters(rho, viscosity, temperature, grad_p, force);

  lb_init_lattice(grid);

  lb_init_fluid();

}

/***********************************************************************/

void lb_finalize() {
  int i, x, y;

  for (i=0; i<lbmodel.n_vel; ++i) {
    for (x=0; x<3; ++x) {
      for (y=0; y<lblattice.halo_grid[1]+2; ++y) {
	free(FI(i,x,y));
      }
    }
  }    

  free(lbmom);

}

/***********************************************************************/

void write_profile() {
  int x, y, z, zl, zh;  
  double rho, j[3], *m;
  FILE *file;

  x = lblattice.halo_grid[0]>>1;
  y = lblattice.halo_grid[1]>>1;

  zl = lblattice.halo_size[2];
  zh = lblattice.grid[2]+lblattice.halo_size[2];

  m = lbmom + lbmodel.n_mom*(lblattice.stride[0]*x+lblattice.stride[1]*y+zl);

  file = fopen("profile.dat","w");
  for (z=zl; z<zh; ++z, m+=lbmodel.n_mom) {
    rho  = m[0] + lbpar.rho;
    if (rho > 0.0) {
      j[0] = m[1]/rho;
    } else {
      j[0] = 0.0;
    }
    fprintf(file,"%f %f\n",(double)z,j[0]/rho);
  }
  fclose(file);

}

/***********************************************************************/

#if 0
int main(int argc, char *argv[]) {
  int i, n_steps, grid[3], vol;
  double rho, viscosity, temperature, pgrad, force[3];
  double start, finish, elapsed, mups;

  if (argc!=2) {
    fprintf(stderr, "Usage: lb-walking <nsteps>\n");
    return -1;
  }	   

  n_steps = atoi(argv[1]);

  grid[0] = 10;
  grid[1] = 10;
  grid[2] = 99;

  vol = grid[0]*grid[1]*grid[2];  

  rho         = 1.0;
  viscosity   = 1.0;
  temperature = 0.0;
  pgrad       = 0.0;
  force[0]    = 0.008*rho*viscosity/(grid[2]*grid[2]);
  force[1]    = 0.0;
  force[2]    = 0.0;

  lb_init(grid,rho,viscosity,temperature,pgrad,force);

  fprintf(stdout, "Running  %d iterations\n", n_steps); fflush(stdout);

  start = (double) clock();
  for (i=0; i<n_steps; ++i) {
    lb_update();
  }
  finish = (double) clock();

  elapsed = (finish-start)/CLOCKS_PER_SEC;
  mups = vol*n_steps/elapsed/1e6;

  fprintf(stdout, "Elapsed time: %.3f s (%.3e MUPS)\n", elapsed, mups); fflush(stdout); 

  write_profile();

  lb_finalize();

  return EXIT_SUCCESS;

}
#endif

/***********************************************************************/
