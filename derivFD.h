#ifndef _DERIVFD_H
#define _DERIVFD_H

#include "d2q21.h"

void firstDer (double *, double *);
void secDerAA (double *, double *);
void secDerAB (double [][NDIM], double *);
void thirdDer (double *, double *);

void firstDer_coeff(double res[NDIM][NDIM]);
void secDer_coeff(double res[NFD][NDIM][NDIM]);


/***********************************************************************/

static const double tau_derFirst[13] = { 0.,
	2./3.,			2./3.,			2./3.,			2./3.,
	0.,				0.,				0.,				0.,
	-1./24.,		-1./24.,		-1./24.,		-1./24. };

/***********************************************************************/

static const double tau_derSecAA[13] = { -4.,
	1.,		1.,		1.,		1.,
	0.,		0.,		0.,		0.,
	0.,		0.,		0.,		0. };

/***********************************************************************/

static const double tau_derSecAB[13] = { 0.,
	1.,			1.,			1.,			1.,
	1./4.,		1./4.,		1./4.,		1./4.,
	0.,			0.,			0.,			0. };

/***********************************************************************/

static const double tau_derThird[13] = { 0.,
	-2.,			-2.,			-2.,			-2.,
	1./2.,		1./2.,		1./2.,		1./2.,
	1./4.,		1./4.,		1./4.,		1./4. };

/***********************************************************************/

#endif
