/**
 *  @file ChemEquil.h
 *
 *  Chemical equilibrium.
 *
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2001 California Institute of Technology
 *
 */


#ifndef CT_CHEM_EQUIL_H
#define CT_CHEM_EQUIL_H


// Cantera includes
#include "ct_defs.h"
#include "vec_functions.h"
#include "ctexceptions.h"
#include "ThermoPhase.h"
#include "DenseMatrix.h"

namespace Cantera {

    int _equilflag(const char* xy);

    /**
     *  Chemical equilibrium options. Used internally by class ChemEquil.
     */
    class EquilOpt {
    public:
        EquilOpt() : relTolerance(1.e-10), maxIterations(1000), iterations(0), 
            maxStepSize(10.0), propertyPair(TP), contin(false) {}
        
        doublereal relTolerance;      ///< Relative tolerance 
        int maxIterations;            ///< Maximum number of iterations
        int iterations;               ///< Iteration counter

        /**
         * Maximum step size. Largest change in any element potential or
         * in log(T) allowed in one Newton step. Default: 10.0
         */
        doublereal maxStepSize;       

        /** 
         * Property pair flag. Determines which two thermodynamic properties
         * are fixed.
         */
        int propertyPair;

        /** 
         * Continuation flag. Set true if the calculation should be
         * initialized from the last calculation. Otherwise, the
         * calculation will be started from scratch and the initial
         * composition and element potentials estimated. (Not Implemented.)
         */
        bool contin;
    };

    template<class M>
    class PropertyCalculator;

    /**
     *  Chemical equilibrium processor. Sets a mixture to a state of
     *  chemical equilibrium.
     */
    class ChemEquil {

    public:
        ChemEquil();
        virtual ~ChemEquil();

        int equilibrate(thermo_t& s, int XY = 0);
        int equilibrate(thermo_t& s, int XY, vector_fp& elMoles);
        const vector_fp& elementPotentials() const { return m_lambda; }

        /**
         * Options controlling how the calculation is carried out. 
         * @see EquilOptions
         */
        EquilOpt options;


    protected:

        thermo_t*  m_phase;
        thermo_t* m_thermo;

        /// number of atoms of element m in species k.
        doublereal nAtoms(int k, int m) const { return m_comp[k*m_mm + m]; }

        void initialize(thermo_t& s);

        void setToEquilState(thermo_t& s, 
            const vector_fp& x, doublereal t);

        int setInitialMoles(thermo_t& s, vector_fp& elementMoles);

        int estimateElementPotentials(thermo_t& s, 
            vector_fp& lambda);
        
        int dampStep(thermo_t& s, vector_fp& oldx, 
            double oldf, vector_fp& grad, vector_fp& step, vector_fp& x, 
            double& f, vector_fp& elmols, int XY, double xval, double yval );

        void equilResidual(thermo_t& s, const vector_fp& x, 
            const vector_fp& elmtotal, vector_fp& resid, 
            int XY, double xval, double yval);

        void equilJacobian(thermo_t& s, vector_fp& x,  
            const vector_fp& elmols, DenseMatrix& jac, 
            int XY, double xval, double yval);

        void update(const thermo_t& s);

        int m_mm;
        int m_kk;
        int m_skip;

        PropertyCalculator<thermo_t> *m_p1, *m_p2;

        vector_fp m_molefractions;
        vector_fp m_lambda;
        vector_fp m_elementmolefracs;
        vector_fp m_reswork;
        vector_fp m_jwork1;
        vector_fp m_jwork2;
        vector_fp m_comp;
        doublereal m_temp, m_dens;
        doublereal m_p0;
        int m_eloc;
        doublereal m_abscharge;

        doublereal m_startTemp, m_startDens;
        vector_fp m_startSoln;

        vector_fp m_grt;
        vector_fp m_mu_RT;
    };


    //-----------------------------------------------------------
    //              convenience functions
    //-----------------------------------------------------------

    /**
     * Set a mixture to a state of chemical equilibrium. The flag 'XY'
     * determines the two properties that will be held fixed in the
     * calculation.
     */
    inline void equilibrate(thermo_t& s, int XY) {
        ChemEquil e;
        e.equilibrate(s,XY);
    }

    /**
     * Set a mixture to a state of chemical equilibrium. The flag 'XY'
     * determines the two properties that will be held fixed in the
     * calculation.
     */
    inline void equilibrate(thermo_t& s, const char* XY) {
        ChemEquil e;
        e.equilibrate(s,_equilflag(XY));
        s.setElementPotentials(e.elementPotentials());
    }

}


#endif

