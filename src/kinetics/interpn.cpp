#include <stdlib.h>

namespace Cantera
{

/*! Given an array xx[0..n-1] , and given a value x ,
 * returns a value j such that x is between xx[j] and xx[j+1].
 * xx must be monotonic, either increasing or decreasing.
 * j=0 or j=n is returned to indicate that x is out of range.
 * adapted from Numerical Recipes in C, 2nd Ed., p. 117.
 */
int locate(double const *xx, double x, int n) {
  xx -= 1;  // zero-offset to unit-offset

  int j;
  int ju, jm, jl;
  int ascnd;

  jl = 0;
  ju = n + 1;
  ascnd = (xx[n] >= xx[1]);
  while (ju - jl > 1) {
    jm = (ju + jl) >> 1;
    if (x >= xx[jm] == ascnd)
      jl = jm;
    else
      ju = jm;
  }
  if (x == xx[1])
    j = 1;
  else if (x == xx[n])
    j = n - 1;
  else
    j = jl;

  j -= 1;  // unit-offset to zero-offset
  return j;
}

/*! Multidimensional linear interpolation
 * val[0..nval-1]   : output values
 * coor[0..ndim-1]  : coordinate of the interpolation point
 * data[...]        : points to the start position of a multidimensional data
 * table. len[0..ndim-1]   : length of each dimension axis[...]        :
 * coordinates of each dimesnion is placed sequentially in axis
 */
void interpn(double *val, double const *coor, double const *data,
             double const *axis, size_t const *len, int ndim, int nval) {
  int i1, i2;
  i1 = locate(axis, *coor, *len);

  // if the interpolation value is out of bound
  // use the closest value
  if (i1 == -1) {
    i1 = 0;
    i2 = 0;
  } else if (i1 == *len - 1) {
    i1 = *len - 1;
    i2 = *len - 1;
  } else
    i2 = i1 + 1;

  double *v1 = (double *)malloc(nval * sizeof(double));
  double *v2 = (double *)malloc(nval * sizeof(double));

  double x1 = axis[i1];
  double x2 = axis[i2];

  if (ndim == 1) {
    for (int j = 0; j < nval; ++j) {
      v1[j] = data[i1 * nval + j];
      v2[j] = data[i2 * nval + j];
    }
  } else {
    int s = nval;
    for (int j = 1; j < ndim; ++j) s *= len[j];
    interpn(v1, coor + 1, data + i1 * s, axis + *len, len + 1, ndim - 1, nval);
    interpn(v2, coor + 1, data + i2 * s, axis + *len, len + 1, ndim - 1, nval);
  }

  if (x2 != x1)
    for (int j = 0; j < nval; ++j)
      val[j] = ((*coor - x1) * v2[j] + (x2 - *coor) * v1[j]) / (x2 - x1);
  else
    for (int j = 0; j < nval; ++j) val[j] = (v1[j] + v2[j]) / 2.;

  free(v1);
  free(v2);
}

}  // namespace Cantera
