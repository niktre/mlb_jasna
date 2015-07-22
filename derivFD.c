#include <stdlib.h>
#include "derivFD.h"
#include "defs.h"

void firstDer (double res, double *f, int field_offset, int x, int y, int dir) {
	const double (*c)[lbmodel.n_dim] = lbmodel.c;
	double *m = f + lblattice.halo_grid_volume*lbmodel.n_vel;
	int i;

	double field;
	res = 0.;
	
	for (i=0; i<lbmodel.n_fd; ++i) {
		// do not know what field will be needed. Therefore, field_offset is added.
		// if rho is needed then field_offset is 0, ux 1, uy 2, etc.
		// NOT sure that the next pointer is correct!!!
		field = *(m + (lblattice.nb_offset[i]+field_offset)*lbmodel.n_vel);
		res += tau_derFirst[i]*c[i][dir]*field;
	}
}

/***********************************************************************/

void secDerAA (double res, double *field, int x, int y, int dir) {
	int i;
	res = 0.;
	
	for (i=0; i<lbmodel.n_fd; ++i) {
		res += tau_derFirst[i]*lbmodel.c[i][dir]*field[i];
	}
}

/***********************************************************************/

void secDerAB (double res, double *field, int x, int y, int dir) {
	int i;
	res = 0.;
	
	for (i=0; i<lbmodel.n_fd; ++i) {
		res += tau_derFirst[i]*lbmodel.c[i][dir]*field[i];
	}
}

/***********************************************************************/

void thirdDer (double res, double *field, int x, int y, int dir) {
	int i;
	res = 0.;
	
	for (i=0; i<lbmodel.n_fd; ++i) {
		res += tau_derFirst[i]*lbmodel.c[i][dir]*field[i];
	}
}

/***********************************************************************/