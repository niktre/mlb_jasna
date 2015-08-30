/***********************************************************************
 *
 * mlb.h
 *
 * Copyright (c) 2015 Ulf D. Schiller <ulf@lattice-boltzmann.de>
 * All rights reserved.
 *
 ***********************************************************************/

#define TOLERANCE 1.e-6

void mlb_calc_force(double *force, double *m, int x, int y);

void mlb_correction_current(double *m);

void mlb_interface_collisions(double *f);

void mlb_correction_collisions(double *f);

/***********************************************************************/
