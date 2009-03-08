/**
 *
 *  @file MetalPhase.h
 *
 */

/*  $Author: dggoodwin $
 *  $Date: 2005/08/18 14:44:35 $
 *  $Revision: 1.3 $
 *
 *  Copyright 2003 California Institute of Technology
 *
 */


#ifndef CT_METALPHASE_H
#define CT_METALPHASE_H


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

        // Overloaded methoods of class ThermoPhase

        virtual int eosType() const { return cMetal; }

        virtual doublereal enthalpy_mole() const { return 0.0; }
        virtual doublereal intEnergy_mole() const { return 0.0; }
        virtual doublereal entropy_mole() const { return 0.0; }
        virtual doublereal gibbs_mole() const { return 0.0; }
        virtual doublereal cp_mole() const { return 0.0; }
        virtual doublereal cv_mole() const { return 0.0; }

        virtual void setPressure(doublereal pres) { m_press = pres; }
        virtual doublereal  pressure() const { return m_press; }
 
        virtual void getChemPotentials(doublereal* mu) const {
            int n, nsp = nSpecies();
            for (n = 0; n < nsp; n++) mu[n] = 0.0;
        }

        virtual void getStandardChemPotentials(doublereal* mu0) const {
            int n, nsp = nSpecies();
            for (n = 0; n < nsp; n++) mu0[n] = 0.0;
        }

        virtual void getActivityConcentrations(doublereal* c) const {
            int n, nsp = nSpecies();
            for (n = 0; n < nsp; n++) c[n] = 1.0;
        }

         virtual doublereal standardConcentration(int k=0) const {
             return 1.0;
        }

         virtual doublereal logStandardConc(int k=0) const {
             return 0.0;
        }

        virtual void setParametersFromXML(const XML_Node& eosdata) {
            eosdata._require("model","Metal");
            doublereal rho = getFloat(eosdata, "density", "-");
            setDensity(rho);
        }

    protected:

    private:
        doublereal m_press;
    };
}
        
#endif





