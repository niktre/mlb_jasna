#include <stdlib.h>
#include "derivFD.h"
#include "defs.h"

void firstDer (double *res, double *m) {
  const double *tau = lbmodel.fd_weights[0];
  const double (*c)[lbmodel.n_dim] = lbmodel.c;
  int i, j;
  double field;

  for (j=0; j<lbmodel.n_dim; ++j) res[j] = 0.;

  for (i=0; i<lbmodel.n_fd; ++i) {
    field = m[lblattice.nb_offset[i]*lbmodel.n_vel];
    for (j=0; j<lbmodel.n_dim; ++j) {
      res[j] += tau[i]*c[i][j]*field;
    }
  }

}

/***********************************************************************/

void secDerAA (double *res, double *m) {
  const double *tau = lbmodel.fd_weights[1];
  int i;
  double field;

  *res = 0.;
	
  for (i=0; i<lbmodel.n_fd; ++i) {
    field = m[lblattice.nb_offset[i]*lbmodel.n_vel];
    *res += tau[i]*field;
  }

}

/***********************************************************************/

void secDerAB (double *res, double *m) {
  const double *tau = lbmodel.fd_weights[2];
  const double (*c)[lbmodel.n_dim] = lbmodel.c;
  int i, j, k;
  double field, trace;

  for (j=0; j<lbmodel.n_dim*(lbmodel.n_dim+1)/2; ++j) res[j] = 0.;
	
  for (i=0; i<lbmodel.n_fd; ++i) {
    field = m[lblattice.nb_offset[i]*lbmodel.n_vel];
    for (j=0; j<lbmodel.n_dim; ++j) {
      for (k=0; k<=j; ++k) {
	res[j*lbmodel.n_dim+k] += tau[i]*c[i][j]*c[i][k]*field;
      }
    }
  }

  secDerAA(&trace, m);

  for (j=0; j<lbmodel.n_dim; ++j) {
    res[j*lbmodel.n_dim+j] += trace/lbmodel.n_dim;
  }

}

/***********************************************************************/

void thirdDer (double *res, double *m) {
  const double *tau = lbmodel.fd_weights[3];
  const double (*c)[lbmodel.n_dim] = lbmodel.c;
  int i, j;
  double field;
	
  for (j=0; j<lbmodel.n_dim; ++j) res[j] = 0.;

  for (i=0; i<lbmodel.n_fd; ++i) {
    field = m[lblattice.nb_offset[i]*lbmodel.n_vel];
    for (j=0; j<lbmodel.n_dim; ++j) {
      res[j] += tau[i]*c[i][j]*field;
    }
  }

}

/***********************************************************************/
