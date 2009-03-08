/**
 *
 *  @file SurfPhase.h
 *
 */

/*  $Author: dggoodwin $
 *  $Date: 2005/06/18 17:01:10 $
 *  $Revision: 1.11 $
 *
 *  Copyright 2002 California Institute of Technology
 *
 */


#ifndef CT_SURFPHASE_H
#define CT_SURFPHASE_H

#include "mix_defs.h"
#include "ThermoPhase.h"


namespace Cantera {


    /**
     * A simple model for a surface phase. The surface consists of a
     * grid of equivalent sites. Surface species may be defined that
     * occupy one or more sites. The surface species are assumed to be
     * independent, and thus the species form an ideal solution.
     * The definitions of the member functions are located in 
     * InterfaceKinetics.cpp.
     */
    class SurfPhase : public ThermoPhase  {

    public:

        /// Constructor.
        SurfPhase(doublereal n0 = 0.0);

        /// Destructor.
        virtual ~SurfPhase();

        //----- reimplimented methods of class ThermoPhase ------

        virtual int eosType() const { return cSurf; }
        virtual doublereal enthalpy_mole() const;
        virtual doublereal intEnergy_mole() const;
        virtual void getStandardChemPotentials(doublereal* mu0) const;
        virtual void getChemPotentials(doublereal* mu) const;
        virtual void getActivityConcentrations(doublereal* c) const;
        virtual doublereal standardConcentration(int k = 0) const;
        virtual doublereal logStandardConc(int k=0) const;
        virtual void setParameters(int n, doublereal* c);
        virtual void setParametersFromXML(const XML_Node& eosdata);
        virtual void initThermo();
        virtual void setStateFromXML(const XML_Node& state);
        doublereal siteDensity(){ return m_n0; }
        void setPotentialEnergy(int k, doublereal pe);
        doublereal potentialEnergy(int k) {return m_pe[k];}
        void setSiteDensity(doublereal n0);

        void getEnthalpy_RT(doublereal* hrt) const;
        void getEntropy_R(doublereal* sr) const;

        virtual doublereal pressure() const {
            return m_press;
        }

        virtual void setPressure(doublereal p) {
            m_press = p;
        }


        //------- new methods defined in this class ----------

        /**
         * Set the surface site fractions to a specified 
         * state. This routine converts to concentrations
         * in kmol/m2, using m_n0, the surface site density,
         * and size(k), which is defined to be the number of
         * surface sites occupied by the kth molecule.
         * It then calls State::setConcentrations to set the
         * internal concentration in the object.
	 *
	 * @param theta[k] This is the surface site fraction
	 *                 for the kth species in the surface phase.
	 *                 This is a dimensionless quantity.
         */
        void setCoverages(const doublereal* theta);

        /**
         * Set the coverages without normalizing them to sum to 1.0.
         * This may be used when the normalization condition is part
         * of the system of equations being solved.
         */
        void setCoveragesNoNorm(const doublereal* theta);

        /**
         * Set the coverages from a string of colon-separated
         * name:value pairs.
         */
        void setCoveragesByName(string cov);

        /**
         * Get the coverages. Array theta must be at least as long as
         * the number of species.
         */
        void getCoverages(doublereal* theta) const;

    protected:

        doublereal m_n0;
        doublereal m_logn0;
        doublereal m_tmin, m_tmax;
	doublereal m_press;

        mutable doublereal    m_tlast;
        mutable array_fp      m_h0;
        mutable array_fp      m_s0;
        mutable array_fp      m_cp0;
        mutable array_fp      m_mu0;
        mutable array_fp      m_work;
        mutable array_fp      m_pe;
        mutable array_fp      m_logsize;

    private:

        void _updateThermo(bool force=false) const;

    };
}
        
#endif





