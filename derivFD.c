#include "derivFD.h"
//#include "d2q21.h"

void firstDer (double res, double *field, int x, int y, int dir) {
	int i;
	res = 0.;
//	LB_Model lbmodel = DnQm(NDIM,NVEL);
	
	for (i=0; i<NFORCE; ++i) {
//		res += tau_derFirst[i]*lbmodel.c[i][dir]*field[i];
	}
}