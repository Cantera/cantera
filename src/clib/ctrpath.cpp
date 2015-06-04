/**
 * @file ctrpath.cpp
 */
#define CANTERA_USE_INTERNAL
#include "ctrpath.h"

// Cantera includes
#include "cantera/kinetics/ReactionPath.h"
#include "Cabinet.h"

#include <fstream>

using namespace Cantera;
using namespace std;

typedef Cabinet<ReactionPathBuilder> BuilderCabinet;
typedef Cabinet<ReactionPathDiagram> DiagramCabinet;
template<> DiagramCabinet* DiagramCabinet::s_storage = 0;
template<> BuilderCabinet* BuilderCabinet::s_storage = 0;

typedef Cabinet<Kinetics> KineticsCabinet;

extern "C" {

    int rdiag_new()
    {
        try {
            ReactionPathDiagram* d = new ReactionPathDiagram();
            return DiagramCabinet::add(d);
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

    int rdiag_copy(int i)
    {
        try {
            return DiagramCabinet::newCopy(i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_assign(int i, int j)
    {
        try {
            return DiagramCabinet::assign(i,j);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_detailed(int i)
    {
        try {
            DiagramCabinet::item(i).show_details = true;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_brief(int i)
    {
        try {
            DiagramCabinet::item(i).show_details = false;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setThreshold(int i, double v)
    {
        try {
            DiagramCabinet::item(i).threshold = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setBoldColor(int i, char* color)
    {
        try {
            DiagramCabinet::item(i).bold_color = color;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setNormalColor(int i, char* color)
    {
        try {
            DiagramCabinet::item(i).normal_color = color;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setDashedColor(int i, char* color)
    {
        try {
            DiagramCabinet::item(i).dashed_color = color;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setDotOptions(int i, char* opt)
    {
        try {
            DiagramCabinet::item(i).dot_options = opt;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setFont(int i, char* font)
    {
        try {
            DiagramCabinet::item(i).setFont(font);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setBoldThreshold(int i, double v)
    {
        try {
            DiagramCabinet::item(i).bold_min = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setNormalThreshold(int i, double v)
    {
        try {
            DiagramCabinet::item(i).dashed_max = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setLabelThreshold(int i, double v)
    {
        try {
            DiagramCabinet::item(i).label_min = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setScale(int i, double v)
    {
        try {
            DiagramCabinet::item(i).scale = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setFlowType(int i, int iflow)
    {
        try {
            if (iflow == 0) {
                DiagramCabinet::item(i).flow_type = OneWayFlow;
            } else {
                DiagramCabinet::item(i).flow_type = NetFlow;
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setArrowWidth(int i, double v)
    {
        try {
            DiagramCabinet::item(i).arrow_width = v;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_setTitle(int i, char* title)
    {
        try {
            DiagramCabinet::item(i).title = title;
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_add(int i, int n)
    {
        try {
            DiagramCabinet::item(i).add(DiagramCabinet::item(n));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_findMajor(int i, double threshold,
                        size_t lda, double* a)
    {
        try {
            DiagramCabinet::item(i).findMajorPaths(threshold, lda, a);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_write(int i, int fmt, char* fname)
    {
        try {
            ofstream f(fname);
            if (fmt == 0) {
                DiagramCabinet::item(i).exportToDot(f);
            } else {
                DiagramCabinet::item(i).writeData(f);
            }
            f.close();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rdiag_displayOnly(int i, int k)
    {
        try {
            DiagramCabinet::item(i).displayOnly(k);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rbuild_new()
    {
        try {
            ReactionPathBuilder* d = new ReactionPathBuilder();
            return BuilderCabinet::add(d);
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

    int rbuild_init(int i, char* logfile, int k)
    {
        try {
            ofstream flog(logfile);
            BuilderCabinet::item(i).init(flog, KineticsCabinet::item(k));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int rbuild_build(int i, int k, char* el, char* dotfile,
                     int idiag, int iquiet)
    {
        try {
            ofstream fdot(dotfile);
            bool quiet = false;
            if (iquiet > 0) {
                quiet = true;
            }
            BuilderCabinet::item(i).build(KineticsCabinet::item(k), el, fdot,
                                          DiagramCabinet::item(idiag), quiet);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int clear_rxnpath()
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
