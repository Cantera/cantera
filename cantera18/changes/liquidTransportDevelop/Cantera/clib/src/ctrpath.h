/**
 * @file ctrpath.h
 */
/*
 *      $Id: ctrpath.h,v 1.2 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifndef CTC_RXNPATH_H
#define CTC_RXNPATH_H

#include "clib_defs.h"

extern "C" {  

    EEXXTT int DLL_CPREFIX rdiag_new();
    EEXXTT int DLL_CPREFIX rdiag_del(int i);
    EEXXTT int DLL_CPREFIX rdiag_detailed(int i);
    EEXXTT int DLL_CPREFIX rdiag_brief(int i);
    EEXXTT int DLL_CPREFIX rdiag_setThreshold(int i, double v);
    EEXXTT int DLL_CPREFIX rdiag_setBoldColor(int i, char* color);
    EEXXTT int DLL_CPREFIX rdiag_setNormalColor(int i, char* color);
    EEXXTT int DLL_CPREFIX rdiag_setDashedColor(int i, char* color);
    EEXXTT int DLL_CPREFIX rdiag_setDotOptions(int i, char* opt);
    EEXXTT int DLL_CPREFIX rdiag_setBoldThreshold(int i, double v);
    EEXXTT int DLL_CPREFIX rdiag_setNormalThreshold(int i, double v);
    EEXXTT int DLL_CPREFIX rdiag_setLabelThreshold(int i, double v);
    EEXXTT int DLL_CPREFIX rdiag_setScale(int i, double v);
    EEXXTT int DLL_CPREFIX rdiag_setFlowType(int i, int iflow);
    EEXXTT int DLL_CPREFIX rdiag_setArrowWidth(int i, double v);
    EEXXTT int DLL_CPREFIX rdiag_setTitle(int i, char* title);
    EEXXTT int DLL_CPREFIX rdiag_write(int i, int fmt, char* fname);
    EEXXTT int DLL_CPREFIX rdiag_add(int i, int n);
    EEXXTT int DLL_CPREFIX rdiag_findMajor(int i, double threshold, int lda, double* a);
    EEXXTT int DLL_CPREFIX rdiag_setFont(int i, char* font);
    EEXXTT int DLL_CPREFIX rdiag_displayOnly(int i, int k);

    EEXXTT int DLL_CPREFIX rbuild_new();
    EEXXTT int DLL_CPREFIX rbuild_del(int i);
    EEXXTT int DLL_CPREFIX rbuild_init(int i, char* logfile, int k);
    EEXXTT int DLL_CPREFIX rbuild_build(int i, int k, char* el, char* dotfile, 
        int idiag, int iquiet);
}

#endif
