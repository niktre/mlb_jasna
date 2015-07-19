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
#include "defs.h"

/***********************************************************************/

/* temporary/secondary grid */
#define PFI *FI(NVEL,WGRID)
#define FI(i,x) fi[i][x]

/***********************************************************************/

static const LB_Model lbmodel = DnQm(NDIM,NVEL);

static LB_Lattice lblattice;

static LB_Parameters lbpar = { 1.0, 0.0 };

static double *lbmom = NULL;
static double PFI;

/***********************************************************************/

void lb_halo_copy() {

  int x, y;
  int lsrc, rsrc, ldst, rdst, size;

  int xstride = lbmodel.n_mom*lblattice.stride[0];
  int ystride = lbmodel.n_mom*lblattice.stride[1];

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

  /* Y direction is handled directly in the update loop */

}

/***********************************************************************/

static void lb_calc_equilibrium(double *f, double *f_eq) {

  int i;
  double w[lbmodel.n_vel];
  double rho, cs2, uc, u2;
  double j[lbmodel.n_dim], u[lbmodel.n_dim];

  rho = j[0] = j[1] = 0.0;
  for (i=0; i<lbmodel.n_vel; ++i) {
    rho  += f[i];
    j[0] += f[i]*lbmodel.c[i][0];
    j[1] += f[i]*lbmodel.c[i][1];
  }

  u[0] = j[0]/rho;
  u[1] = j[1]/rho;

  cs2 = eq_state(rho);

  lb_weights(w, cs2);

  u2   = (u[0]*u[0] + u[1]*u[1])/cs2;

  for (i=0; i<lbmodel.n_vel; ++i) {
    uc = (u[0]*lbmodel.c[i][0] + u[1]*lbmodel.c[i][1])/cs2;
    f_eq[i] = w[i]*rho*(1.0 + uc + 0.5*uc*uc - 0.5*u2);
  }

}

/***********************************************************************/

static void lb_collisions(double *f) {
  
  int i;
  double f_eq[lbmodel.n_vel];
  
  lb_calc_equilibrium(f, f_eq);

  for (i=0; i<lbmodel.n_vel; ++i) {
    f[i] = f_eq[i] + lbpar.gamma * ( f[i] - f_eq[i] );
  }
  
}

/***********************************************************************/

static void lb_stream(double *f, double PFI, int x, int y) {
  int xc, xp, xm, xp2, xm2, xp4, xm4;
  xc = x%WGRID;
  xp  = (x+1)%WGRID; xm  = (x-1+WGRID)%WGRID;
  xp2 = (x+2)%WGRID; xm2 = (x-2+WGRID)%WGRID;
  xp4 = (x+4)%WGRID; xm4 = (x-4+WGRID)%WGRID;

  FI( 0, xc )[y]   = f[0];
  FI( 1, xp )[y]   = f[1];
  FI( 2, xm )[y]   = f[2];
  FI( 3, xc )[y+1] = f[3];
  FI( 4, xc )[y-1] = f[4];
  FI( 5, xp )[y+1] = f[5];
  FI( 6, xm )[y-1] = f[6];
  FI( 7, xm )[y+1] = f[7];
  FI( 8, xp )[y-1] = f[8];
  FI( 9, xp2)[y]   = f[9];
  FI(10, xm2)[y]   = f[10];
  FI(11, xc )[y+2] = f[11];
  FI(12, xc )[y-2] = f[12];
  FI(13, xp2)[y+2] = f[13];
  FI(14, xm2)[y-2] = f[14];
  FI(15, xm2)[y+2] = f[15];
  FI(16, xp2)[y-2] = f[16];
  FI(17, xp4)[y]   = f[17];
  FI(18, xm4)[y]   = f[18];
  FI(19, xc )[y+4] = f[19];
  FI(20, xc )[y-4] = f[20];
  
}

/***********************************************************************/

static void lb_read(double *f, double PFI, int x, int y) {

  lb_collisions(f);
  lb_stream(f, fi, x, y);

}

/***********************************************************************/

static void lb_write(double *f, double PFI, int x, int y) {

  int xc = x%WGRID;

  f[ 0] = FI( 0, xc)[y];
  f[ 1] = FI( 1, xc)[y];
  f[ 2] = FI( 2, xc)[y];
  f[ 3] = FI( 3, xc)[y];
  f[ 4] = FI( 4, xc)[y];
  f[ 5] = FI( 5, xc)[y];
  f[ 6] = FI( 6, xc)[y];
  f[ 7] = FI( 7, xc)[y];
  f[ 8] = FI( 8, xc)[y];
  f[ 9] = FI( 9, xc)[y];
  f[10] = FI(10, xc)[y];
  f[11] = FI(11, xc)[y];
  f[12] = FI(12, xc)[y];
  f[13] = FI(13, xc)[y];
  f[14] = FI(14, xc)[y];
  f[15] = FI(15, xc)[y];
  f[16] = FI(16, xc)[y];
  f[17] = FI(17, xc)[y];
  f[18] = FI(18, xc)[y];
  f[19] = FI(19, xc)[y];
  f[20] = FI(20, xc)[y];

}

/***********************************************************************/

static void lb_read_column(double *m, double PFI, int x) {
  int y, yl, yh;

  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1] - 1;

  for (y=yl, m+=yl*lbmodel.n_mom; y<=yh; ++y, m+=lbmodel.n_mom) {
    lb_read(m, fi, x, y);
  }

  /* boundary conditions in y */
  int xc, xp, xm, xp2, xm2;
  xc  = x%WGRID;
  xp  = (x+1)%WGRID; xm  = (x-1+WGRID)%WGRID;
  xp2 = (x+2)%WGRID; xm2 = (x-2+WGRID)%WGRID;

  FI( 3, xc )[yl]   = FI( 3, xc )[yh+1];
  FI( 5, xp )[yl]   = FI( 5, xp )[yh+1];
  FI( 7, xm )[yl]   = FI( 7, xm )[yh+1];
  FI( 4, xc )[yh]   = FI( 4, xc )[yl-1];
  FI( 6, xm )[yh]   = FI( 6, xm )[yl-1];
  FI( 8, xp )[yh]   = FI( 8, xp )[yl-1];

  FI(11, xc )[yl]   = FI(11, xc )[yh+1];
  FI(13, xp2)[yl]   = FI(13, xp2)[yh+1];
  FI(15, xm2)[yl]   = FI(15, xm2)[yh+1];
  FI(12, xc )[yh]   = FI(12, xc )[yl-1];
  FI(14, xm2)[yh]   = FI(14, xm2)[yl-1];
  FI(16, xp2)[yh]   = FI(16, xp2)[yl-1];
  FI(11, xc )[yl+1] = FI(11, xc )[yh+2];
  FI(13, xp2)[yl+1] = FI(13, xp2)[yh+2];
  FI(15, xm2)[yl+1] = FI(15, xm2)[yh+2];
  FI(12, xc )[yh-1] = FI(12, xc )[yl-2];
  FI(14, xm2)[yh-1] = FI(14, xm2)[yl-2];
  FI(16, xp2)[yh-1] = FI(16, xp2)[yl-2];

  FI(19, xc )[yl]   = FI(19, xc )[yh+1];
  FI(19, xc )[yl+1] = FI(19, xc )[yh+2];
  FI(19, xc )[yl+2] = FI(19, xc )[yh+3];
  FI(19, xc )[yl+3] = FI(19, xc )[yh+4];
  FI(20, xc )[yh]   = FI(20, xc )[yl-1];
  FI(20, xc )[yh-1] = FI(20, xc )[yl-2];
  FI(20, xc )[yh-2] = FI(20, xc )[yl-3];
  FI(20, xc )[yh-3] = FI(20, xc )[yl-4];

}

/***********************************************************************/

static void lb_write_column(double *m, double PFI, int x) {
  int y, yl, yh;

  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1] - 1;

  for (y=yl, m+=yl*lbmodel.n_mom; y<=yh; ++y, m+=lbmodel.n_mom) {
    lb_write(m, fi, x, y);
  }
  
}

/***********************************************************************/

void lb_update() {
  int x;
  int xstride = lblattice.stride[0]*lbmodel.n_mom;
  double *m   = lbmom;

  lb_halo_copy(); /* need up to date moments in halo */

  /* Columns in the lower halo will be read only */
  for (x=0; x<lblattice.halo_size[0]; ++x, m+=xstride) {
    lb_read_column(m, fi, x);
  }

  /* Collide and stream column x, read back column x-HALO
   * x-HALO can be overwritten and all info is available now */
  for (x=lblattice.halo_size[0]; x<lblattice.halo_grid[0]; ++x, m+=xstride) {
    lb_read_column(m, fi, x);
    lb_write_column(m-HALO*xstride, fi, x-HALO);
  }

}

/***********************************************************************/

static void lb_init_fluid() {
  int x, y, z, i;
  double  cs2, w[lbmodel.n_vel];
  double *m = lbmom;

  cs2 = eq_state(lbpar.rho);

  lb_weights(w, cs2);

  for (x=0; x<lblattice.halo_grid[0]; ++x) {
    for (y=0; y<lblattice.halo_grid[1]; ++y, m+=lbmodel.n_mom) {
      for (i=0; i<lbmodel.n_vel; ++i) {
        m[i] = w[i]*lbpar.rho;
      }
    }
  }

}

/***********************************************************************/

static void lb_init_lattice(int *grid) {
  int i, x, y, hgrid[lbmodel.n_dim], hsize[lbmodel.n_dim];

  lblattice.grid[0] = grid[0];
  lblattice.grid[1] = grid[1];

  lblattice.halo_size[0] = hsize[0] = HALO;
  lblattice.halo_size[1] = hsize[1] = HALO;

  lblattice.halo_grid[0] = hgrid[0] = lblattice.grid[0] + 2*hsize[0];
  lblattice.halo_grid[1] = hgrid[1] = lblattice.grid[1] + 2*hsize[1];

  lblattice.stride[1] = 1;
  lblattice.stride[0] = hgrid[1];

  lblattice.halo_grid_volume = hgrid[0]*hgrid[1];

  lbmom = calloc(lblattice.halo_grid_volume*lbmodel.n_mom, sizeof(*lbmom));

  for (i=0; i<lbmodel.n_vel; ++i) {
    for (x=0; x<WGRID; ++x) {
      FI(i,x) = calloc(lblattice.halo_grid[1],sizeof(*FI(0,0)));
    }
  }

}

/***********************************************************************/

static void lb_init_parameters(double rho, double gamma) {

  lbpar.rho   = rho;
  lbpar.gamma = gamma;

}

/***********************************************************************/

void lb_init(int *grid, double rho, double gamma) {

  lb_init_parameters(rho, gamma);

  lb_init_lattice(grid);

  lb_init_fluid();

}

/***********************************************************************/

void lb_finalize() {
  int i, x, y;

  for (i=0; i<lbmodel.n_vel; ++i) {
    for (x=0; x<WGRID; ++x) {
      free(FI(i,x));
    }
  }    

  free(lbmom);

}

/***********************************************************************/

void write_profile(int write_halo) {
  int i, x, y, xl, xh, yl, yh, xoff;
  double rho, j[lbmodel.n_dim], *m;
  FILE *file;

  if (write_halo) {
    lb_halo_copy();
    xl = 0;
    xh = lblattice.halo_grid[0];
    yl = 0;
    yh = lblattice.halo_grid[1];
  } else {
    xl = lblattice.halo_size[0];
    xh = lblattice.halo_size[0] + lblattice.grid[0];
    yl = lblattice.halo_size[1];
    yh = lblattice.halo_size[1] + lblattice.grid[1];
  }

  xoff = lblattice.halo_grid[0] - (xh - xl);

  m = lbmom + lbmodel.n_mom*(lblattice.stride[0]*xl+yl);

  file = fopen("profile.dat","w");
  for (x=xl; x<xh; ++x, m+=lbmodel.n_mom*xoff) {
    for (y=yl; y<yh; ++y, m+=lbmodel.n_mom) {
      rho = j[0] = j[1] = 0.0;
      for (i=0; i<lbmodel.n_vel; ++i) {
	rho  += m[i];
	j[0] += m[i]*lbmodel.c[i][0];
	j[1] += m[i]*lbmodel.c[i][1];
      }
      fprintf(file,"%f %f %f %f\n",(double)x-xl,(double)y-yl,rho,j[0]/rho);
    }
  }
  fclose(file);

}

/***********************************************************************/

#if 1
int main(int argc, char *argv[]) {
  int i, n_steps, grid[lbmodel.n_dim], vol;
  double rho, gamma;
  double start, finish, elapsed, mups;

  if (argc!=2) {
    fprintf(stderr, "Usage: ./run <nsteps>\n");
    return -1;
  }	   

  n_steps = atoi(argv[1]);

  grid[0] = 10;
  grid[1] = 10;

  vol = grid[0]*grid[1];

  rho   = 1.0;
  gamma = 1.0;

  lb_init(grid,rho,gamma);

  fprintf(stdout, "Running  %d iterations\n", n_steps); fflush(stdout);

  start = (double) clock();
  for (i=0; i<n_steps; ++i) {
    lb_update();
  }
  finish = (double) clock();

  elapsed = (finish-start)/CLOCKS_PER_SEC;
  mups = vol*n_steps/elapsed/1e6;

  fprintf(stdout, "Elapsed time: %.3f s (%.3e MUPS)\n", elapsed, mups); fflush(stdout); 

  write_profile(0);

  lb_finalize();

  return EXIT_SUCCESS;

}
#endif

/***********************************************************************/
