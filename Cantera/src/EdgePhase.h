/**
 *
 *  @file EdgePhase.h
 *
 */

/*  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2002 California Institute of Technology
 *
 */


#ifndef CT_EDGEPHASE_H
#define CT_EDGEPHASE_H

#include "mix_defs.h"
#include "ThermoPhase.h"

namespace Cantera {

    class EdgePhase : public SurfPhase  {

    public:

        EdgePhase(doublereal n0 = 0.0);
        virtual ~EdgePhase() {}
        virtual int eosType() const { return cEdge; }
        virtual void setParametersFromXML(const XML_Node& eosdata);
    };
}
        
#endif





