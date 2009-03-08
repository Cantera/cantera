/**
 *
 *  @file LatticeSolidPhase.h
 */

/*  $Author: hkmoffa $
 *  $Date: 2006/06/13 00:24:26 $
 *  $Revision: 1.5 $
 *
 *  Copyright 2005 California Institute of Technology
 *
 */

#ifndef CT_LATTICESOLID_H
#define CT_LATTICESOLID_H

#include "ct_defs.h"
#include "mix_defs.h"
#include "ThermoPhase.h"
#include "SpeciesThermo.h"
#include "utilities.h"

namespace Cantera {

    class LatticePhase;

    class LatticeSolidPhase : public ThermoPhase  {

    public:

        LatticeSolidPhase() : m_tlast(0.0) {}

        virtual ~LatticeSolidPhase() {}

        virtual int eosType() const { return cLatticeSolid; }

        virtual doublereal enthalpy_mole() const;

        virtual doublereal intEnergy_mole() const;

        virtual doublereal entropy_mole() const;

        virtual doublereal gibbs_mole() const;

        virtual doublereal cp_mole() const;

        virtual doublereal cv_mole() const {
            return cp_mole();
        }

        virtual doublereal pressure() const {
            return m_press;
        }

        virtual void setPressure(doublereal p) {
            m_press = p;
            setMolarDensity(m_molar_density);
        }

        virtual void getActivityConcentrations(doublereal* c) const;

	virtual void getActivityCoefficients(doublereal* ac) const;

        virtual void getChemPotentials(doublereal* mu) const;
        virtual void getStandardChemPotentials(doublereal* mu0) const;
        virtual doublereal standardConcentration(int k=0) const;
        virtual doublereal logStandardConc(int k=0) const;

        virtual void initThermo();

        virtual void setParametersFromXML(const XML_Node& eosdata);

        void setLatticeMoleFractions(int n, string x);

    protected:
        
        int m_mm;
        int m_kk;
        mutable doublereal        m_tlast;
        doublereal                m_press;
        doublereal                m_molar_density;
        int                       m_nlattice;
        vector<LatticePhase*>     m_lattice;
        mutable vector_fp                 m_x;

    private:

        void _updateThermo() const;
    };
}
        
#endif
