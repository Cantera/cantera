/**
 * @file ctrpath.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_RXNPATH_H
#define CTC_RXNPATH_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    CANTERA_CAPI int rdiag_new();
    CANTERA_CAPI int rdiag_del(int i);
    CANTERA_CAPI int rdiag_detailed(int i);
    CANTERA_CAPI int rdiag_brief(int i);
    CANTERA_CAPI int rdiag_setThreshold(int i, double v);
    CANTERA_CAPI int rdiag_setBoldColor(int i, const char* color);
    CANTERA_CAPI int rdiag_setNormalColor(int i, const char* color);
    CANTERA_CAPI int rdiag_setDashedColor(int i, const char* color);
    CANTERA_CAPI int rdiag_setDotOptions(int i, const char* opt);
    CANTERA_CAPI int rdiag_setBoldThreshold(int i, double v);
    CANTERA_CAPI int rdiag_setNormalThreshold(int i, double v);
    CANTERA_CAPI int rdiag_setLabelThreshold(int i, double v);
    CANTERA_CAPI int rdiag_setScale(int i, double v);
    CANTERA_CAPI int rdiag_setFlowType(int i, int iflow);
    CANTERA_CAPI int rdiag_setArrowWidth(int i, double v);
    CANTERA_CAPI int rdiag_setTitle(int i, const char* title);
    CANTERA_CAPI int rdiag_write(int i, int fmt, const char* fname);
    CANTERA_CAPI int rdiag_add(int i, int n);
    CANTERA_CAPI int rdiag_findMajor(int i, double threshold, size_t lda, double* a);
    CANTERA_CAPI int rdiag_setFont(int i, const char* font);
    CANTERA_CAPI int rdiag_displayOnly(int i, int k);

    CANTERA_CAPI int rbuild_new();
    CANTERA_CAPI int rbuild_del(int i);
    CANTERA_CAPI int rbuild_init(int i, const char* logfile, int k);
    CANTERA_CAPI int rbuild_build(int i, int k, const char* el, const char* dotfile,
                                  int idiag, int iquiet);

    CANTERA_CAPI int ct_clearReactionPath();

#ifdef __cplusplus
}
#endif

#endif
