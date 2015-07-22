#ifndef _DERIVFD_H
#define _DERIVFD_H

void firstDer (double, double *, int, int, int, int);
void secDerAA (double, double *, int, int, int);
void secDerAB (double, double *, int, int, int, int, int);
void thirdDer (double, double *, int, int, int, int);

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