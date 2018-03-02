/**
 *  @file ConstDensityThermo.h
 * Header for a Thermo manager for incompressible ThermoPhases
 * (see \ref thermoprops and \link Cantera::ConstDensityThermo ConstDensityThermo\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

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
 * @ingroup thermoprops
 */
class ConstDensityThermo : public ThermoPhase
{
public:
    //! Constructor.
    ConstDensityThermo() {}

    virtual std::string type() const {
        return "ConstDensity";
    }

    virtual double enthalpy_mole() const;
    virtual double entropy_mole() const;
    virtual double cp_mole() const;
    virtual double cv_mole() const;

    //! Return the thermodynamic pressure (Pa).
    virtual double pressure() const;

    //! Set the internally stored pressure (Pa) at constant temperature and
    //! composition
    /*!
     *  @param p input Pressure (Pa)
     */
    virtual void setPressure(double p);

    virtual void getActivityConcentrations(double* c) const;
    virtual void getActivityCoefficients(double* ac) const;

    virtual void getChemPotentials(double* mu) const;
    virtual void getStandardChemPotentials(double* mu0) const;

    //! Returns the standard Concentration in units of m3 kmol-1.
    //! @copydoc ThermoPhase::standardConcentration
    virtual double standardConcentration(size_t k=0) const;

    virtual void getPureGibbs(double* gpure) const {
        const vector_fp& gibbsrt = gibbs_RT();
        scale(gibbsrt.begin(), gibbsrt.end(), gpure, RT());
    }

    void getEnthalpy_RT(double* hrt) const {
        const vector_fp& _h = enthalpy_RT();
        std::copy(_h.begin(), _h.end(), hrt);
    }

    void getEntropy_R(double* sr) const {
        const vector_fp& _s = entropy_R();
        std::copy(_s.begin(), _s.end(), sr);
    }

    virtual void getGibbs_RT(double* grt) const {
        const vector_fp& gibbsrt = gibbs_RT();
        std::copy(gibbsrt.begin(), gibbsrt.end(), grt);
    }

    void getCp_R(double* cpr) const {
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

    virtual void setParameters(int n, double* const c) {
        setDensity(c[0]);
    }

    virtual void getParameters(int& n, double* const c) const {
        double d = density();
        c[0] = d;
        n = 1;
    }

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
    double m_press;

private:
    //! Function to update the reference state thermo functions
    void _updateThermo() const;
};
}

#endif
