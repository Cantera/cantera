// declarations for Fortran simplex routines from Numerical Recipes.

#ifndef NUM_RECIPES_H
#define NUM_RECIPES_H

extern "C" {

integer simplx_(doublereal *a, integer *m, integer *n, integer *
	mp, integer *np, integer *m1, integer *m2, integer *m3, integer *
	icase, integer *izrov, integer *iposv);

integer splin2_(doublereal *x1a, doublereal *x2a, doublereal *ya, 
	    doublereal *y2a, integer *m, integer *n,
	    doublereal *x1, doublereal *x2, doublereal *y);

integer splie2_(doublereal *x1a, doublereal *x2a, doublereal *ya, 
	    integer *m, integer *n, doublereal *y2a);
}


#endif
