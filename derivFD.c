#include <stdlib.h>
#include "derivFD.h"
#include "d2q21.h"
//#include "defs.h"

void firstDer (double res, double *f, int field_offset, int x, int y, int dir) {
  const double *tau = lbmodel.fd_weights[0];
	const double (*c)[lbmodel.n_dim] = lbmodel.c;
	double *m = f + lblattice.halo_grid_volume*lbmodel.n_vel;
	int i;
	double field;

	res = 0.;
	
	for (i=0; i<lbmodel.n_fd; ++i) {
//		field = (m + lblattice.nb_offset[i]*lbmodel.n_vel)[field_offset];
		field = *(m + lblattice.nb_offset[i]*lbmodel.n_vel + field_offset);
		res += tau[i]*c[i][dir]*field;
	}
}

/***********************************************************************/

void secDerAA (double res, double *f, int field_offset, int x, int y) {
  const double *tau = lbmodel.fd_weights[1];
	double *m = f + lblattice.halo_grid_volume*lbmodel.n_vel;
	int i;
	double field;
	res = 0.;
	
	for (i=0; i<lbmodel.n_fd; ++i) {
		field = *(m + lblattice.nb_offset[i]*lbmodel.n_vel + field_offset);
		res += tau[i]*field;
	}
}

/***********************************************************************/

void secDerAB (double res, double *f, int field_offset, int x, int y, int dir1, int dir2) {
	const double *tau = lbmodel.fd_weights[2];
	const double (*c)[lbmodel.n_dim] = lbmodel.c;
	double *m = f + lblattice.halo_grid_volume*lbmodel.n_vel;
	int i;
	double field;
	
	res = 0.;
	
	for (i=0; i<lbmodel.n_fd; ++i) {
		field = *(m + lblattice.nb_offset[i]*lbmodel.n_vel + field_offset);
		res += tau[i]*c[i][dir1]*c[i][dir2]*field;
	}
}

/***********************************************************************/

void thirdDer (double res, double *f, int field_offset, int x, int y, int dir) {
	const double *tau = lbmodel.fd_weights[3];
	const double (*c)[lbmodel.n_dim] = lbmodel.c;
	double *m = f + lblattice.halo_grid_volume*lbmodel.n_vel;
	int i;
	double field;
	
	res = 0.;
	
	for (i=0; i<lbmodel.n_fd; ++i) {
		field = *(m + lblattice.nb_offset[i]*lbmodel.n_vel + field_offset);
		res += tau[i]*c[i][dir]*field;
	}
}

/***********************************************************************/