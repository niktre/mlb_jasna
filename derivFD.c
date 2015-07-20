#include "derivFD.h"
#include "defs.h"

void firstDer (double res, double *field, int x, int y, int dir) {
	int i;
	res = 0.;

	for (i=0; i<NFORCE; ++i) {
	  res += tau_derFirst[i]*lbmodel.c[i][dir]*field[i];
	}
}
