
DEPRECATED

/**
 *  @file GRI30.h
 *
 *  GRI-Mech 3.0
 *
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_GRI30_H
#define CT_GRI30_H

#include <string>

#include "Phase.h"
#include "IdealGasThermo.h"
#include "GRI_30_Kinetics.h"

#include "global.h"
#include "import.h"
#include "SpeciesThermoFactory.h"

#include "speciesThermoTypes.h"


namespace Cantera {

    /**
     * Implements reaction mechanism GRI-Mech 3.0
     */
    class GRI30 : 
        public Phase, public IdealGasThermo, public GRI_30_Kinetics
    {
    public:

//         GRI30(map<string, string>& params) {
//             setSpThermo(NASA);
//             initThermo(*this);
//             GRI_30_Kinetics::setThermo(*this);
//             m_ok = importFromFile(this, this, params);
//         }

        GRI30() {
            setSpThermo(NASA);
            initThermo(*this);
            GRI_30_Kinetics::setThermo(*this);
            map<string, string> params;
            params["input"] = "gri30.xml";
            params["ID"] = "gri30";
            //if (validate) params["validate"] = "yes";
            m_ok = importFromFile(this, this, params);
            //m_ok = import("gri30.xml", "", false);
        }

        /**
         *  Destructor. Does nothing.
         */
        virtual ~GRI30() {}

        bool valid() const { return m_ok; }

        bool operator!() const {return !m_ok;}

    protected:

//         bool import(string infile, string dbfile="", bool validate=false) {
//             map<string, string> params;
//             params["input"] = infile;
//             params["database"] = dbfile;
//             if (validate) params["validate"] = "yes";
//             m_ok = importFromFile(this, this, params);
//             return m_ok;
//         }

        void setSpThermo(int paramType, SpeciesThermoFactory* fsp = 0) {
            if (fsp == 0) fsp = SpeciesThermoFactory::factory();
            setSpeciesThermo(fsp->newSpeciesThermo(paramType));
        }

        bool m_ok;

    private:
    };
}

#endif
