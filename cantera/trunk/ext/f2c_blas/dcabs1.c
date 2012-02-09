#include "blaswrap.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

doublereal dcabs1_(doublecomplex *z__)
{
    /* System generated locals */
    doublereal ret_val;
    static doublecomplex equiv_0[1];
    /* Local variables */
#define t ((doublereal *)equiv_0)
#define zz (equiv_0)
    zz->r = z__->r, zz->i = z__->i;
    ret_val = abs(t[0]) + abs(t[1]);
    return ret_val;
} /* dcabs1_ */
#undef zz
#undef t
#ifdef __cplusplus
}
#endif

