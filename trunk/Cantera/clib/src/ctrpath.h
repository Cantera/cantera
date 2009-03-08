#ifndef CTC_RXNPATH_H
#define CTC_RXNPATH_H

#include "clib_defs.h"

extern "C" {  

    int DLL_IMPORT rdiag_new();
    int DLL_IMPORT rdiag_del(int i);
    int DLL_IMPORT rdiag_detailed(int i);
    int DLL_IMPORT rdiag_brief(int i);
    int DLL_IMPORT rdiag_setThreshold(int i, double v);
    int DLL_IMPORT rdiag_setBoldColor(int i, char* color);
    int DLL_IMPORT rdiag_setNormalColor(int i, char* color);
    int DLL_IMPORT rdiag_setDashedColor(int i, char* color);
    int DLL_IMPORT rdiag_setDotOptions(int i, char* opt);
    int DLL_IMPORT rdiag_setBoldThreshold(int i, double v);
    int DLL_IMPORT rdiag_setNormalThreshold(int i, double v);
    int DLL_IMPORT rdiag_setLabelThreshold(int i, double v);
    int DLL_IMPORT rdiag_setScale(int i, double v);
    int DLL_IMPORT rdiag_setFlowType(int i, int iflow);
    int DLL_IMPORT rdiag_setArrowWidth(int i, double v);
    int DLL_IMPORT rdiag_setTitle(int i, char* title);
    int DLL_IMPORT rdiag_write(int i, int fmt, char* fname);
    int DLL_IMPORT rdiag_add(int i, int n);
    int DLL_IMPORT rdiag_findMajor(int i, double threshold, int lda, double* a);
    int DLL_IMPORT rdiag_setFont(int i, char* font);
    int DLL_IMPORT rdiag_displayOnly(int i, int k);

    int DLL_IMPORT rbuild_new();
    int DLL_IMPORT rbuild_del(int i);
    int DLL_IMPORT rbuild_init(int i, char* logfile, int k);
    int DLL_IMPORT rbuild_build(int i, int k, char* el, char* dotfile, 
        int idiag, int iquiet);
}

#endif
