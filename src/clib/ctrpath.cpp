/**
 * @file ctrpath.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ctrpath.h"

// Cantera includes
#include "cantera/kinetics/ReactionPath.h"
#include "clib_utils.h"

#include <fstream>

using namespace Cantera;
using namespace std;

typedef SharedCabinet<ReactionPathBuilder> BuilderCabinet;
typedef SharedCabinet<ReactionPathDiagram> DiagramCabinet;
template<> DiagramCabinet* DiagramCabinet::s_storage = 0;
template<> BuilderCabinet* BuilderCabinet::s_storage = 0;

typedef SharedCabinet<Kinetics> KineticsCabinet;
template<> KineticsCabinet* KineticsCabinet::s_storage; // defined in ct.cpp

extern "C" {

    int rdiag_new()
    {
        try {
            return DiagramCabinet::add(make_shared<ReactionPathDiagram>());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_del(int i)
    {
        try {
            DiagramCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_detailed(int i)
    {
        try {
            DiagramCabinet::at(i)->show_details = true;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_brief(int i)
    {
        try {
            DiagramCabinet::at(i)->show_details = false;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setThreshold(int i, double v)
    {
        try {
            DiagramCabinet::at(i)->threshold = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setBoldColor(int i, const char* color)
    {
        try {
            DiagramCabinet::at(i)->bold_color = color;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setNormalColor(int i, const char* color)
    {
        try {
            DiagramCabinet::at(i)->normal_color = color;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setDashedColor(int i, const char* color)
    {
        try {
            DiagramCabinet::at(i)->dashed_color = color;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setDotOptions(int i, const char* opt)
    {
        try {
            DiagramCabinet::at(i)->dot_options = opt;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setFont(int i, const char* font)
    {
        try {
            DiagramCabinet::at(i)->setFont(font);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setBoldThreshold(int i, double v)
    {
        try {
            DiagramCabinet::at(i)->bold_min = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setNormalThreshold(int i, double v)
    {
        try {
            DiagramCabinet::at(i)->dashed_max = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setLabelThreshold(int i, double v)
    {
        try {
            DiagramCabinet::at(i)->label_min = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setScale(int i, double v)
    {
        try {
            DiagramCabinet::at(i)->scale = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setFlowType(int i, int iflow)
    {
        try {
            if (iflow == 0) {
                DiagramCabinet::at(i)->flow_type = OneWayFlow;
            } else {
                DiagramCabinet::at(i)->flow_type = NetFlow;
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setArrowWidth(int i, double v)
    {
        try {
            DiagramCabinet::at(i)->arrow_width = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setTitle(int i, const char* title)
    {
        try {
            DiagramCabinet::at(i)->title = title;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_add(int i, int n)
    {
        try {
            DiagramCabinet::at(i)->add(*DiagramCabinet::at(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_findMajor(int i, double threshold,
                        size_t lda, double* a)
    {
        try {
            DiagramCabinet::at(i)->findMajorPaths(threshold, lda, a);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_write(int i, int fmt, const char* fname)
    {
        try {
            ofstream f(fname);
            if (fmt == 0) {
                DiagramCabinet::at(i)->exportToDot(f);
            } else {
                DiagramCabinet::at(i)->writeData(f);
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_displayOnly(int i, int k)
    {
        try {
            DiagramCabinet::at(i)->displayOnly(k);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rbuild_new()
    {
        try {
            return BuilderCabinet::add(make_shared<ReactionPathBuilder>());
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rbuild_del(int i)
    {
        try {
            BuilderCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rbuild_init(int i, const char* logfile, int k)
    {
        try {
            ofstream flog(logfile);
            BuilderCabinet::at(i)->init(flog, *KineticsCabinet::at(k));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rbuild_build(int i, int k, const char* el, const char* dotfile,
                     int idiag, int iquiet)
    {
        try {
            ofstream fdot(dotfile);
            bool quiet = false;
            if (iquiet > 0) {
                quiet = true;
            }
            BuilderCabinet::at(i)->build(*KineticsCabinet::at(k), el, fdot,
                                         *DiagramCabinet::at(idiag), quiet);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_clearReactionPath()
    {
        try {
            DiagramCabinet::clear();
            BuilderCabinet::clear();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }
}
