/**
 * @file ctrpath.cpp
 */
#define CANTERA_USE_INTERNAL
#include "ctrpath.h"

// Cantera includes
#include "cantera/kinetics/ReactionPath.h"
#include "Cabinet.h"

using namespace Cantera;
using namespace std;

typedef Cabinet<ReactionPathBuilder> BuilderCabinet;
typedef Cabinet<ReactionPathDiagram> DiagramCabinet;
template<> DiagramCabinet* DiagramCabinet::__storage = 0;
template<> BuilderCabinet* BuilderCabinet::__storage = 0;

typedef Cabinet<Kinetics> KineticsCabinet;

extern "C" {

    int DLL_EXPORT rdiag_new()
    {
        ReactionPathDiagram* d = new ReactionPathDiagram();
        return DiagramCabinet::add(d);
    }

    int DLL_EXPORT rdiag_del(int i)
    {
        DiagramCabinet::del(i);
        return 0;
    }

    int DLL_EXPORT rdiag_copy(int i)
    {
        return DiagramCabinet::newCopy(i);
    }

    int DLL_EXPORT rdiag_assign(int i, int j)
    {
        return DiagramCabinet::assign(i,j);
    }

    int DLL_EXPORT rdiag_detailed(int i)
    {
        DiagramCabinet::item(i).show_details = true;
        return 0;
    }

    int DLL_EXPORT rdiag_brief(int i)
    {
        DiagramCabinet::item(i).show_details = false;
        return 0;
    }

    int DLL_EXPORT rdiag_setThreshold(int i, double v)
    {
        DiagramCabinet::item(i).threshold = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setBoldColor(int i, char* color)
    {
        DiagramCabinet::item(i).bold_color = string(color);
        return 0;
    }

    int DLL_EXPORT rdiag_setNormalColor(int i, char* color)
    {
        DiagramCabinet::item(i).normal_color = string(color);
        return 0;
    }

    int DLL_EXPORT rdiag_setDashedColor(int i, char* color)
    {
        DiagramCabinet::item(i).dashed_color = string(color);
        return 0;
    }

    int DLL_EXPORT rdiag_setDotOptions(int i, char* opt)
    {
        DiagramCabinet::item(i).dot_options = string(opt);
        return 0;
    }

    int DLL_EXPORT rdiag_setFont(int i, char* font)
    {
        DiagramCabinet::item(i).setFont(string(font));
        return 0;
    }

    int DLL_EXPORT rdiag_setBoldThreshold(int i, double v)
    {
        DiagramCabinet::item(i).bold_min = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setNormalThreshold(int i, double v)
    {
        DiagramCabinet::item(i).dashed_max = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setLabelThreshold(int i, double v)
    {
        DiagramCabinet::item(i).label_min = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setScale(int i, double v)
    {
        DiagramCabinet::item(i).scale = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setFlowType(int i, int iflow)
    {
        if (iflow == 0) {
            DiagramCabinet::item(i).flow_type = OneWayFlow;
        } else {
            DiagramCabinet::item(i).flow_type = NetFlow;
        }
        return 0;
    }

    int DLL_EXPORT rdiag_setArrowWidth(int i, double v)
    {
        DiagramCabinet::item(i).arrow_width = v;
        return 0;
    }

    int DLL_EXPORT rdiag_setTitle(int i, char* title)
    {
        DiagramCabinet::item(i).title = string(title);
        return 0;
    }

    int DLL_EXPORT rdiag_add(int i, int n)
    {
        DiagramCabinet::item(i).add(DiagramCabinet::item(n));
        return 0;
    }

    int DLL_EXPORT rdiag_findMajor(int i, double threshold,
                                   size_t lda, double* a)
    {
        DiagramCabinet::item(i).findMajorPaths(threshold, lda, a);
        return 0;
    }

    int DLL_EXPORT rdiag_write(int i, int fmt, char* fname)
    {
        ofstream f(fname);
        if (fmt == 0) {
            DiagramCabinet::item(i).exportToDot(f);
        } else {
            DiagramCabinet::item(i).writeData(f);
        }
        f.close();
        return 0;
    }

    int DLL_EXPORT rdiag_displayOnly(int i, int k)
    {
        DiagramCabinet::item(i).displayOnly(k);
        return 0;
    }

    int DLL_EXPORT rbuild_new()
    {
        ReactionPathBuilder* d = new ReactionPathBuilder();
        return BuilderCabinet::add(d);
    }

    int DLL_EXPORT rbuild_del(int i)
    {
        BuilderCabinet::del(i);
        return 0;
    }

    int DLL_EXPORT rbuild_init(int i, char* logfile, int k)
    {
        ofstream flog(logfile);
        BuilderCabinet::item(i).init(flog, KineticsCabinet::item(k));
        return 0;
    }

    int DLL_EXPORT rbuild_build(int i, int k, char* el, char* dotfile,
                                int idiag, int iquiet)
    {
        ofstream fdot(dotfile);
        bool quiet = false;
        if (iquiet > 0) {
            quiet = true;
        }
        BuilderCabinet::item(i).build(KineticsCabinet::item(k), string(el), fdot,
                                      DiagramCabinet::item(idiag), quiet);
        return 0;
    }

}
