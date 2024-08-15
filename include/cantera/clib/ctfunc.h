/**
 * @file ctfunc.h
 *
 * @warning  This module is an experimental part of the %Cantera API and
 *      may be changed or removed without notice.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_FUNC1_H
#define CTC_FUNC1_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    //! @implements checkFunc1
    //! @since New in %Cantera 3.1
    CANTERA_CAPI int func_check(const char* type, size_t len, char* buf);
    //! @implements newFunc1(const string&, double)
    CANTERA_CAPI int func_new_basic(const char* type, double c);
    //! @implements newFunc1(const string&, const vector<double>&)
    CANTERA_CAPI int func_new_advanced(const char* type, size_t lenp, const double* p);
    //! @implements newFunc1(const string&, const shared_ptr<Func1>, const shared_ptr<Func1>)
    CANTERA_CAPI int func_new_compound(const char* type, int a, int b);
    //! @implements newFunc1(const string&, const shared_ptr<Func1>, double)
    CANTERA_CAPI int func_new_modified(const char* type, int a, double c);

    CANTERA_CAPI int func_new_sum(int a, int b);  //!< @implements newSumFunction
    CANTERA_CAPI int func_new_diff(int a, int b);  //!< @implements newDiffFunction
    CANTERA_CAPI int func_new_prod(int a, int b);  //!< @implements newProdFunction
    CANTERA_CAPI int func_new_ratio(int a, int b);  //!< @implements newRatioFunction
    CANTERA_CAPI int func_del(int i);

    //! @implements Func1::type
    CANTERA_CAPI int func_type(int i, size_t lennm, char* nm);
    CANTERA_CAPI double func_value(int i, double t);  //!< @implements Func1::eval
    CANTERA_CAPI int func_derivative(int i);  //!< @implements Func1::derivative
    CANTERA_CAPI int func_duplicate(int i);
    //! @implements Func1::write
    //! @since Changed signature in %Cantera 3.1
    CANTERA_CAPI int func_write(int i, const char* arg, size_t lennm, char* nm);

    CANTERA_CAPI int ct_clearFunc();

#ifdef __cplusplus
}
#endif

#endif
