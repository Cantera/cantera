/**
 * @file ctmethods.cpp
 *
 * The interface between the MATLAB environment and the C++ Cantera
 * kernel is through a single MEX file. This is top-level driver for
 * the MEX file.
 *
 * This file handles the methods of all Cantera MATLAB classes. The
 * class is indicated by the first parameter in the call from MATLAB.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include <string>

#include "cantera/clib/ct.h"
#include "ctmatutils.h"
#include "mllogger.h"

const int NO_CLASS = 0;
const int XML_CLASS = 10;
const int THERMO_CLASS = 20;
const int PHASE_CLASS = 30;
const int KINETICS_CLASS = 40;
const int TRANSPORT_CLASS = 50;
const int REACTOR_CLASS = 60;
const int REACTORNET_CLASS = 65;
const int WALL_CLASS = 70;
const int REACTORSURFACE_CLASS = 75;
const int FLOWDEVICE_CLASS = 80;
const int ONEDIM_CLASS = 90;
const int SURF_CLASS = 100;
const int FUNC_CLASS = 110;
const int MIXTURE_CLASS = 120;

void ctfunctions(int nlhs, mxArray* plhs[], int nrhs,
                 const mxArray* prhs[]);

void xmlmethods(int nlhs, mxArray* plhs[], int nrhs,
                const mxArray* prhs[]);

void thermomethods(int nlhs, mxArray* plhs[], int nrhs,
                   const mxArray* prhs[]);

void phasemethods(int nlhs, mxArray* plhs[], int nrhs,
                  const mxArray* prhs[]);

void mixturemethods(int nlhs, mxArray* plhs[], int nrhs,
                    const mxArray* prhs[]);

void surfmethods(int nlhs, mxArray* plhs[], int nrhs,
                 const mxArray* prhs[]);

void kineticsmethods(int nlhs, mxArray* plhs[], int nrhs,
                     const mxArray* prhs[]);

void transportmethods(int nlhs, mxArray* plhs[], int nrhs,
                      const mxArray* prhs[]);

void reactormethods(int nlhs, mxArray* plhs[], int nrhs,
                    const mxArray* prhs[]);

void reactornetmethods(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]);

void wallmethods(int nlhs, mxArray* plhs[], int nrhs,
                 const mxArray* prhs[]);

void reactorsurfacemethods(int nlhs, mxArray* plhs[], int nrhs,
                           const mxArray* prhs[]);

void flowdevicemethods(int nlhs, mxArray* plhs[], int nrhs,
                       const mxArray* prhs[]);

void onedimmethods(int nlhs, mxArray* plhs[], int nrhs,
                   const mxArray* prhs[]);

void funcmethods(int nlhs, mxArray* plhs[], int nrhs,
                 const mxArray* prhs[]);

static Cantera::ML_Logger* _logger = 0;

void initLogger()
{
    if (!_logger) {
        _logger = new Cantera::ML_Logger;
        // Call the DLL program to set the logger
        ct_setLogWriter(_logger);
    }
}

extern "C" {

    void mexFunction(int nlhs, mxArray* plhs[],
                     int nrhs, const mxArray* prhs[])
    {
        // create a log writer for error messages if this is the
        // first MATLAB function call
        initLogger();

        // flag specifying the class
        int iclass = getInt(prhs[0]);

        // Hand off to the appropriate routine, based on the
        // value of the first parameter
        switch (iclass) {
        case NO_CLASS:
            ctfunctions(nlhs, plhs, nrhs, prhs);
            break;
        case XML_CLASS:
            xmlmethods(nlhs, plhs, nrhs, prhs);
            break;
        case THERMO_CLASS:
            thermomethods(nlhs, plhs, nrhs, prhs);
            break;
        case PHASE_CLASS:
            phasemethods(nlhs, plhs, nrhs, prhs);
            break;
        case MIXTURE_CLASS:
            mixturemethods(nlhs, plhs, nrhs, prhs);
            break;
        case KINETICS_CLASS:
            kineticsmethods(nlhs, plhs, nrhs, prhs);
            break;
        case TRANSPORT_CLASS:
            transportmethods(nlhs, plhs, nrhs, prhs);
            break;
        case REACTOR_CLASS:
            reactormethods(nlhs, plhs, nrhs, prhs);
            break;
        case REACTORNET_CLASS:
            reactornetmethods(nlhs, plhs, nrhs, prhs);
            break;
        case WALL_CLASS:
            wallmethods(nlhs, plhs, nrhs, prhs);
            break;
        case REACTORSURFACE_CLASS:
            reactorsurfacemethods(nlhs, plhs, nrhs, prhs);
            break;
        case FLOWDEVICE_CLASS:
            flowdevicemethods(nlhs, plhs, nrhs, prhs);
            break;
        case ONEDIM_CLASS:
            onedimmethods(nlhs, plhs, nrhs, prhs);
            break;
        case SURF_CLASS:
            surfmethods(nlhs, plhs, nrhs, prhs);
            break;
        case FUNC_CLASS:
            funcmethods(nlhs, plhs, nrhs, prhs);
            break;
        default:
            mexPrintf("iclass = %d",iclass);
        }
    }
}
