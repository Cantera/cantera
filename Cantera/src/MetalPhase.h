/**
 *
 *  @file MetalPhase.h
 *
 */

/*  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2003 California Institute of Technology
 *
 */


#ifndef CT_METALPHASE_H
#define CT_METALPHASE_H

//#include "ct_defs.h"
#include "mix_defs.h"
#include "ThermoPhase.h"
#include "SpeciesThermo.h"

namespace Cantera {

    /**
     * @ingroup thermoprops
     *
     * Class MetalPhase represents electrons in a metal.
     *
     */
    class MetalPhase : public ThermoPhase  {

    public:

        MetalPhase() {}

        virtual ~MetalPhase() {}

        /**
         * Equation of state flag.
         */
        virtual int eosType() const { return cMetal; }

        virtual doublereal enthalpy_mole() const { return 0.0; }
        virtual doublereal intEnergy_mole() const { return 0.0; }
        virtual doublereal entropy_mole() const { return 0.0; }
        virtual doublereal gibbs_mole() const { return 0.0; }
        virtual doublereal cp_mole() const { return 0.0; }
        virtual doublereal cv_mole() const { return 0.0; }

        virtual void getChemPotentials(doublereal* mu) const {
            mu[0] = 0.0;
        }

        virtual void getStandardChemPotentials(doublereal* mu0) const {
            mu0[0] = 0.0;
        }

        virtual void getActivityConcentrations(doublereal* c) const {
            c[0] = 1.0;
        }

         virtual doublereal standardConcentration(int k=0) const {
             return 1.0;
        }

         virtual doublereal logStandardConc(int k=0) const {
             return 0.0;
        }

    protected:

    private:

    };
}
        
#endif





