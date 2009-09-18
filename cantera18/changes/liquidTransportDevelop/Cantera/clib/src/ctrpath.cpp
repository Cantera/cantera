/**
 * @file ctrpath.cpp
 */
/*
 *      $Id: ctrpath.cpp,v 1.4 2009/07/11 17:16:09 hkmoffa Exp $
 */

#define CANTERA_USE_INTERNAL
#include "ctrpath.h"


// Cantera includes

#include "ReactionPath.h"
#include "Cabinet.h"
#include "Storage.h"

using namespace Cantera;
using namespace std;

typedef ReactionPathDiagram diag_t;
typedef ReactionPathBuilder builder_t;


template<> Cabinet<ReactionPathDiagram>*   Cabinet<ReactionPathDiagram>::__storage = 0;

template<> Cabinet<builder_t>*             Cabinet<builder_t>::__storage = 0;

inline ReactionPathDiagram* _diag(int i) {
    return Cabinet<ReactionPathDiagram>::cabinet()->item(i);
}

inline builder_t* _builder(int i) {
    return Cabinet<builder_t>::cabinet()->item(i);
}

inline Kinetics* _kin(int n) {
    return Storage::__storage->__ktable[n];
}
 
extern "C" {  

    int DLL_EXPORT rdiag_new() {
        diag_t* d = new ReactionPathDiagram();
        return Cabinet<diag_t>::cabinet()->add(d);
    }

    int DLL_EXPORT rdiag_del(int i) {
        Cabinet<diag_t>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT rdiag_copy(int i) {
        return Cabinet<diag_t>::cabinet()->newCopy(i);
    }

    int DLL_EXPORT rdiag_assign(int i, int j) {
        return Cabinet<diag_t>::cabinet()->assign(i,j);
    }

    int DLL_EXPORT rdiag_detailed(int i) {
        _diag(i)->show_details = true;
        return 0;
    }

    int DLL_EXPORT rdiag_brief(int i) {
        _diag(i)->show_details = false;
        return 0;
    }

    int DLL_EXPORT rdiag_setThreshold(int i, double v) {
        _diag(i)->threshold = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setBoldColor(int i, char* color) {
        _diag(i)->bold_color = string(color);
        return 0;
    }

    int DLL_EXPORT rdiag_setNormalColor(int i, char* color) {
        _diag(i)->normal_color = string(color);
        return 0;
    }

    int DLL_EXPORT rdiag_setDashedColor(int i, char* color) {
        _diag(i)->dashed_color = string(color);
        return 0;
    }

    int DLL_EXPORT rdiag_setDotOptions(int i, char* opt) {
        _diag(i)->dot_options = string(opt);
        return 0;
    }

    int DLL_EXPORT rdiag_setFont(int i, char* font) {
        _diag(i)->setFont(string(font));
        return 0;
    }

    int DLL_EXPORT rdiag_setBoldThreshold(int i, double v) {
        _diag(i)->bold_min = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setNormalThreshold(int i, double v) {
        _diag(i)->dashed_max = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setLabelThreshold(int i, double v) {
        _diag(i)->label_min = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setScale(int i, double v) {
        _diag(i)->scale = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setFlowType(int i, int iflow) {
        if (iflow == 0) 
            _diag(i)->flow_type = OneWayFlow;
        else
            _diag(i)->flow_type = NetFlow;
        return 0;
    }

    int DLL_EXPORT rdiag_setArrowWidth(int i, double v) {
        _diag(i)->arrow_width = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setTitle(int i, char* title) {
        _diag(i)->title = string(title);
        return 0;
    }

    int DLL_EXPORT rdiag_add(int i, int n) {
        _diag(i)->add(*_diag(n));
        return 0;
    }

    int DLL_EXPORT rdiag_findMajor(int i, double threshold, 
        int lda, double* a) {
        _diag(i)->findMajorPaths(threshold, lda, a);
        return 0;
    }

    int DLL_EXPORT rdiag_write(int i, int fmt, char* fname) {
        ofstream f(fname);
        if (fmt == 0) 
            _diag(i)->exportToDot(f);
        else
            _diag(i)->writeData(f);
        f.close();
        return 0;
    }

    int DLL_EXPORT rdiag_displayOnly(int i, int k) {
        _diag(i)->displayOnly(k);
        return 0;
    }

    int DLL_EXPORT rbuild_new() {
        builder_t* d = new ReactionPathBuilder();
        return Cabinet<builder_t>::cabinet()->add(d);
    }

    int DLL_EXPORT rbuild_del(int i) {
        Cabinet<builder_t>::cabinet()->del(i);
        return 0;
    }

    int DLL_EXPORT rbuild_init(int i, char* logfile, int k) {
        ofstream flog(logfile);
        _builder(i)->init(flog, *_kin(k));
        return 0;
    }

    int DLL_EXPORT rbuild_build(int i, int k, char* el, char* dotfile, 
        int idiag, int iquiet) {
        ofstream fdot(dotfile);
        bool quiet = false;
        if (iquiet > 0) quiet = true;
        _builder(i)->build(*_kin(k), string(el), fdot, *_diag(idiag), quiet);
        return 0;
    }

}
