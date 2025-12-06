//! @file Reactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTOR_H
#define CT_REACTOR_H

#include "ReactorBase.h"
#include "cantera/numerics/eigen_sparse.h"


namespace Cantera
{

class Solution;
class AnyMap;

/**
 * Class Reactor is a general-purpose class for stirred reactors. The reactor
 * may have an arbitrary number of inlets and outlets, each of which may be
 * connected to a "flow device" such as a mass flow controller, a pressure
 * regulator, etc. Additional reactors may be connected to the other end of
 * the flow device, allowing construction of arbitrary reactor networks.
 *
 * The reactor class integrates the same governing equations no matter what
 * type of reactor is simulated. The differences among reactor types are
 * completely specified by the attached flow devices and the time-dependent
 * user-specified boundary conditions.
 *
 * If an instance of class Reactor is used directly, it will simulate an
 * adiabatic, constant volume reactor with gas-phase chemistry but no surface
 * chemistry. Other reactor types may be simulated by deriving a class from
 * Reactor.  This method allows specifying the following in terms of the
 * instantaneous reactor state:
 *
 *  - rate of change of the total volume (m^3/s)
 *  - surface heat loss rate (W)
 *  - species surface production rates (kmol/s)
 *
 * See the [Science Reference](../reference/reactors/controlreactor.html) for
 * the governing equations of class Reactor.
 *
 * @ingroup reactorGroup
 */
class Reactor : public ReactorBase
{
public:
    Reactor(shared_ptr<Solution> sol, const string& name="(none)");
    Reactor(shared_ptr<Solution> sol, bool clone, const string& name="(none)");

    string type() const override {
        return "Reactor";
    }

    //! Indicate whether the governing equations for this reactor type are a system of
    //! ODEs or DAEs. In the first case, this class implements the eval() method. In the
    //! second case, this class implements the evalDae() method.
    virtual bool isOde() const {
        return true;
    }

    void setInitialVolume(double vol) override {
        m_vol = vol;
    }

    void setChemistryEnabled(bool cflag=true) override {
        m_chem = cflag;
    }

    bool chemistryEnabled() const override {
        return m_chem;
    }

    void setEnergyEnabled(bool eflag=true) override {
        m_energy = eflag;
    }

    bool energyEnabled() const override {
        return m_energy;
    }

    void getState(double* y) override;
    void initialize(double t0=0.0) override;
    void eval(double t, double* LHS, double* RHS) override;
    vector<size_t> steadyConstraints() const override;
    void updateState(double* y) override;
    void addSensitivityReaction(size_t rxn) override;
    void addSensitivitySpeciesEnthalpy(size_t k) override;

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass", "volume",
    //! "int_energy", the name of a homogeneous phase species, or the name of a
    //! surface species.
    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;
    double upperBound(size_t k) const override;
    double lowerBound(size_t k) const override;
    void resetBadValues(double* y) override;

    //! Set absolute step size limits during advance
    //! @param limits array of step size limits with length neq
    void setAdvanceLimits(const double* limits);

    //! Check whether Reactor object uses advance limits
    //! @returns           True if at least one limit is set, False otherwise
    bool hasAdvanceLimits() const {
        return !m_advancelimits.empty();
    }

    //! Retrieve absolute step size limits during advance
    //! @param[out] limits array of step size limits with length neq
    //! @returns           True if at least one limit is set, False otherwise
    bool getAdvanceLimits(double* limits) const;

    //! Set individual step size limit for component name *nm*
    //! @param nm component name
    //! @param limit value for step size limit
    void setAdvanceLimit(const string& nm, const double limit);

    //! Calculate the reactor-specific Jacobian using a finite difference method.
    //!
    //! This method is used only for informational purposes. Jacobian calculations
    //! for the full reactor system are handled internally by CVODES.
    //!
    //! @warning  This method is an experimental part of the %Cantera
    //! API and may be changed or removed without notice.
    Eigen::SparseMatrix<double> finiteDifferenceJacobian();

    void setDerivativeSettings(AnyMap& settings) override;

    void applySensitivity(double* params) override;
    void resetSensitivity(double* params) override;

protected:
    //! Update #m_sdot to reflect current production rates of bulk phase species due to
    //! reactions on adjacent surfaces.
    void updateSurfaceProductionRates();

    //! Evaluate terms related to Walls. Calculates #m_vdot and #m_Qdot based on
    //! wall movement and heat transfer.
    //! @param t     the current time
    virtual void evalWalls(double t);

    //! Pointer to the homogeneous Kinetics object that handles the reactions
    Kinetics* m_kin = nullptr;

    double m_vdot = 0.0; //!< net rate of volume change from moving walls [m^3/s]
    double m_Qdot = 0.0; //!< net heat transfer into the reactor, through walls [W]
    vector<double> m_wdot; //!< Species net molar production rates
    vector<double> m_uk; //!< Species molar internal energies

    //! Total production rate of bulk phase species on surfaces [kmol/s]
    vector<double> m_sdot;

    bool m_chem = false;
    bool m_energy = true;

    vector<double> m_advancelimits; //!< Advance step limit

    //! Vector of triplets representing the jacobian
    vector<Eigen::Triplet<double>> m_jac_trips;
};
}

#endif
