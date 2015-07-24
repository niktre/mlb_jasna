/***********************************************************************
 *
 * mlb.h
 *
 * Copyright (c) 2015 Ulf D. Schiller <ulf@lattice-boltzmann.de>
 * All rights reserved.
 *
 ***********************************************************************/

void mlb_calc_force(double *force, double *m, int x, int y);

void mlb_correction_current(double *m);

void mlb_interface_collisions(double *f);

void mlc_correction_collisions(double *f);

/***********************************************************************/
