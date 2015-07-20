/***********************************************************************
 *
 * mlb.h
 *
 * Copyright (c) 2015 Ulf D. Schiller <ulf@lattice-boltzmann.de>
 * All rights reserved.
 *
 ***********************************************************************/

void mlb_calc_force(double *force, double *f, int x, int y);

void mlb_interface_collisions(double *f, double *force);

/***********************************************************************/
