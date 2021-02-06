/**
 *  @file ConstDensityThermo.h
 * Header for a Thermo manager for incompressible ThermoPhases
 * (see \ref thermoprops and \link Cantera::ConstDensityThermo ConstDensityThermo\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CONSTRHOTHERMO_H
#define CT_CONSTRHOTHERMO_H

#include "ThermoPhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

//! Overloads the virtual methods of class ThermoPhase to implement the
//! incompressible equation of state.
/**
 * ## Specification of Solution Thermodynamic Properties
 *
 * The density is assumed to be constant, no matter what the concentration of
 * the solution.
 *
 * @deprecated To be removed after Cantera 2.5.0. Replaceable with LatticePhase
 * or IdealSolidSolnPhase
 *
 * @ingroup thermoprops
 */
class ConstDensityThermo : public ThermoPhase
{
public:
    //! Constructor.
    ConstDensityThermo() {
      warn_deprecated("class ConstDensityThermo", "To be removed after Cantera "
        "2.5.0. Consider replacing with LatticePhase or IdealSolidSolnPhase\n");
    }

    virtual std::string type() const {
        return "ConstDensity";
    }

    virtual bool isCompressible() const {
        return false;
    }

    virtual doublereal enthalpy_mole() const;
    virtual doublereal entropy_mole() const;
    virtual doublereal cp_mole() const;
    virtual doublereal cv_mole() const;

    //! Return the thermodynamic pressure (Pa).
    virtual doublereal pressure() const;

    //! Set the internally stored pressure (Pa) at constant temperature and
    //! composition
    /*!
     *  @param p input Pressure (Pa)
     */
    virtual void setPressure(doublereal p);

    virtual void getActivityConcentrations(doublereal* c) const;
    virtual void getActivityCoefficients(doublereal* ac) const;

    virtual void getChemPotentials(doublereal* mu) const;
    virtual void getStandardChemPotentials(doublereal* mu0) const;

    //! Returns the standard Concentration in units of m3 kmol-1.
    //! @copydoc ThermoPhase::standardConcentration
    virtual doublereal standardConcentration(size_t k=0) const;

    virtual void getPureGibbs(doublereal* gpure) const {
        const vector_fp& gibbsrt = gibbs_RT();
        scale(gibbsrt.begin(), gibbsrt.end(), gpure, RT());
    }

    void getEnthalpy_RT(doublereal* hrt) const {
        const vector_fp& _h = enthalpy_RT();
        std::copy(_h.begin(), _h.end(), hrt);
    }

    void getEntropy_R(doublereal* sr) const {
        const vector_fp& _s = entropy_R();
        std::copy(_s.begin(), _s.end(), sr);
    }

    virtual void getGibbs_RT(doublereal* grt) const {
        const vector_fp& gibbsrt = gibbs_RT();
        std::copy(gibbsrt.begin(), gibbsrt.end(), grt);
    }

    void getCp_R(doublereal* cpr) const {
        const vector_fp& _cpr = cp_R();
        std::copy(_cpr.begin(), _cpr.end(), cpr);
    }

    //! Returns a reference to the vector of nondimensional enthalpies of the
    //! reference state at the current temperature of the solution and the
    //! reference pressure for the species.
    const vector_fp& enthalpy_RT() const {
        _updateThermo();
        return m_h0_RT;
    }

    //! Returns a reference to the vector of nondimensional Gibbs Free Energies
    //! of the reference state at the current temperature of the solution and
    //! the reference pressure for the species.
    const vector_fp& gibbs_RT() const {
        _updateThermo();
        return m_g0_RT;
    }

    //! Returns a reference to the vector of nondimensional entropies of the
    //! reference state at the current temperature of the solution and the
    //! reference pressure for each species.
    const vector_fp& entropy_R() const {
        _updateThermo();
        return m_s0_R;
    }

    //! Returns a reference to the vector of nondimensional constant pressure
    //! heat capacities of the reference state at the current temperature of the
    //! solution and reference pressure for each species.
    const vector_fp& cp_R() const {
        _updateThermo();
        return m_cp0_R;
    }

    virtual bool addSpecies(shared_ptr<Species> spec);

    virtual void setParameters(int n, doublereal* const c) {
        assignDensity(c[0]);
    }

    virtual void getParameters(int& n, doublereal* const c) const {
        double d = density();
        c[0] = d;
        n = 1;
    }

    virtual void initThermo();
    virtual void setParametersFromXML(const XML_Node& eosdata);

protected:
    //! Temporary storage for dimensionless reference state enthalpies
    mutable vector_fp m_h0_RT;

    //! Temporary storage for dimensionless reference state heat capacities
    mutable vector_fp m_cp0_R;

    //! Temporary storage for dimensionless reference state Gibbs energies
    mutable vector_fp m_g0_RT;

    //! Temporary storage for dimensionless reference state entropies
    mutable vector_fp m_s0_R;

    //! Current pressure (Pa)
    doublereal m_press;

private:
    //! Function to update the reference state thermo functions
    void _updateThermo() const;
};
}

#endif
