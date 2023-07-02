//! @file StFlow.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_STFLOW_H
#define CT_STFLOW_H

#include "Domain1D.h"
#include "cantera/base/Array.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"

namespace Cantera
{

//------------------------------------------
//   constants
//------------------------------------------

//! Offsets of solution components in the 1D solution array.
enum offset
{
    c_offset_U   //! axial velocity
    , c_offset_V //! strain rate
    , c_offset_T //! temperature
    , c_offset_L //! (1/r)dP/dr
    , c_offset_E //! electric poisson's equation
    , c_offset_Y //! mass fractions
};

class Transport;

/**
 *  This class represents 1D flow domains that satisfy the one-dimensional
 *  similarity solution for chemically-reacting, axisymmetric flows.
 *  @ingroup onedim
 */
class StFlow : public Domain1D
{
public:
    //--------------------------------
    // construction and destruction
    //--------------------------------

    //! Create a new flow domain.
    //! @param ph Object representing the gas phase. This object will be used
    //!     to evaluate all thermodynamic, kinetic, and transport properties.
    //! @param nsp Number of species.
    //! @param points Initial number of grid points
    StFlow(ThermoPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    //! Delegating constructor
    StFlow(shared_ptr<ThermoPhase> th, size_t nsp = 1, size_t points = 1);

    //! Create a new flow domain.
    //! @param sol  Solution object used to evaluate all thermodynamic, kinetic, and
    //!     transport properties
    //! @param id  name of flow domain
    //! @param points  initial number of grid points
    StFlow(shared_ptr<Solution> sol, const std::string& id="", size_t points=1);

    ~StFlow();

    virtual string type() const;

    //! @name Problem Specification
    //! @{

    virtual void setupGrid(size_t n, const doublereal* z);

    virtual void resetBadValues(double* xg);

    ThermoPhase& phase() {
        return *m_thermo;
    }

    Kinetics& kinetics() {
        return *m_kin;
    }

    /**
     * Set the thermo manager.
     *
     * @deprecated  To be removed after Cantera 3.0 (unused)
     */
    void setThermo(ThermoPhase& th);

    virtual void setKinetics(shared_ptr<Kinetics> kin);

    //! Set the kinetics manager.
    //! @deprecated  To be removed after Cantera 3.0; replaced by Domain1D::setKinetics
    void setKinetics(Kinetics& kin);

    virtual void setTransport(shared_ptr<Transport> trans);

    //! Set transport model to existing instance
    //! @deprecated  To be removed after Cantera 3.0; replaced by Domain1D::setKinetics
    void setTransport(Transport& trans);

    //! Set the transport model
    //! @since  New in Cantera 3.0.
    void setTransportModel(const std::string& trans);

    //! Retrieve transport model
    //! @since  New in Cantera 3.0.
    std::string transportModel() const;

    //! Enable thermal diffusion, also known as Soret diffusion.
    //! Requires that multicomponent transport properties be
    //! enabled to carry out calculations.
    void enableSoret(bool withSoret) {
        m_do_soret = withSoret;
    }
    bool withSoret() const {
        return m_do_soret;
    }

    //! Set the pressure. Since the flow equations are for the limit of small
    //! Mach number, the pressure is very nearly constant throughout the flow.
    void setPressure(doublereal p) {
        m_press = p;
    }

    //! The current pressure [Pa].
    doublereal pressure() const {
        return m_press;
    }

    //! Write the initial solution estimate into array x.
    virtual void _getInitialSoln(double* x);

    virtual void _finalize(const doublereal* x);

    //! Sometimes it is desired to carry out the simulation using a specified
    //! temperature profile, rather than computing it by solving the energy
    //! equation. This method specifies this profile.
    void setFixedTempProfile(vector_fp& zfixed, vector_fp& tfixed) {
        m_zfix = zfixed;
        m_tfix = tfixed;
    }

    /*!
     * Set the temperature fixed point at grid point j, and disable the energy
     * equation so that the solution will be held to this value.
     */
    void setTemperature(size_t j, doublereal t) {
        m_fixedtemp[j] = t;
        m_do_energy[j] = false;
    }

    //! The fixed temperature value at point j.
    doublereal T_fixed(size_t j) const {
        return m_fixedtemp[j];
    }

    //! @}

    virtual std::string componentName(size_t n) const;

    virtual size_t componentIndex(const std::string& name) const;

    //! Returns true if the specified component is an active part of the solver state
    virtual bool componentActive(size_t n) const;

    //! Print the solution.
    virtual void show(const double* x);

    virtual shared_ptr<SolutionArray> asArray(const double* soln) const;
    virtual void fromArray(SolutionArray& arr, double* soln);

    //! Set flow configuration for freely-propagating flames, using an internal point
    //! with a fixed temperature as the condition to determine the inlet mass flux.
    void setFreeFlow() {
        m_type = cFreeFlow;
        m_dovisc = false;
        m_isFree = true;
        m_usesLambda = false;
    }

    //! Set flow configuration for axisymmetric counterflow flames, using specified
    //! inlet mass fluxes.
    void setAxisymmetricFlow() {
        m_type = cAxisymmetricStagnationFlow;
        m_dovisc = true;
        m_isFree = false;
        m_usesLambda = true;
    }

    //! Set flow configuration for burner-stabilized flames, using specified inlet mass
    //! fluxes.
    void setUnstrainedFlow() {
        m_type = cAxisymmetricStagnationFlow;
        m_dovisc = false;
        m_isFree = false;
        m_usesLambda = false;
    }

    //! Return the type of flow domain being represented, either "Free Flame" or
    //! "Axisymmetric Stagnation".
    //! @see setFreeFlow setAxisymmetricFlow
    //! @deprecated  To be removed after Cantera 3.0; replaced by type.
    virtual string flowType() const;

    void solveEnergyEqn(size_t j=npos);

    //! Get the solving stage (used by IonFlow specialization)
    //! @since New in Cantera 3.0
    virtual size_t getSolvingStage() const;

    //! Solving stage mode for handling ionized species (used by IonFlow specialization)
    //! - \c stage=1: the fluxes of charged species are set to zero
    //! - \c stage=2: the electric field equation is solved, and the drift flux for
    //!     ionized species is evaluated
    virtual void setSolvingStage(const size_t stage);

    //! Set to solve electric field in a point (used by IonFlow specialization)
    virtual void solveElectricField(size_t j=npos);

    //! Set to fix voltage in a point (used by IonFlow specialization)
    virtual void fixElectricField(size_t j=npos);

    //! Retrieve flag indicating whether electric field is solved or not (used by
    //! IonFlow specialization)
    virtual bool doElectricField(size_t j) const;

    //! Turn radiation on / off.
    /*!
     *  The simple radiation model used was established by Y. Liu and B. Rogg
     *  [Y. Liu and B. Rogg, Modelling of thermally radiating diffusion flames
     *  with detailed chemistry and transport, EUROTHERM Seminars, 17:114-127,
     *  1991]. This model considers the radiation of CO2 and H2O.
     */
    void enableRadiation(bool doRadiation) {
        m_do_radiation = doRadiation;
    }

    //! Returns `true` if the radiation term in the energy equation is enabled
    bool radiationEnabled() const {
        return m_do_radiation;
    }

    //! Return radiative heat loss at grid point j
    double radiativeHeatLoss(size_t j) const {
        return m_qdotRadiation[j];
    }

    //! Set the emissivities for the boundary values
    /*!
     * Reads the emissivities for the left and right boundary values in the
     * radiative term and writes them into the variables, which are used for the
     * calculation.
     */
    void setBoundaryEmissivities(double e_left, double e_right);

    //! Return emissivity at left boundary
    double leftEmissivity() const { return m_epsilon_left; }

    //! Return emissivity at right boundary
    double rightEmissivity() const { return m_epsilon_right; }

    void fixTemperature(size_t j=npos);

    bool doEnergy(size_t j) {
        return m_do_energy[j];
    }

    //! Change the grid size. Called after grid refinement.
    virtual void resize(size_t components, size_t points);

    //! Set the gas object state to be consistent with the solution at point j.
    void setGas(const doublereal* x, size_t j);

    //! Set the gas state to be consistent with the solution at the midpoint
    //! between j and j + 1.
    void setGasAtMidpoint(const doublereal* x, size_t j);

    doublereal density(size_t j) const {
        return m_rho[j];
    }

    virtual bool fixed_mdot();

    void setViscosityFlag(bool dovisc) {
        m_dovisc = dovisc;
    }

    /*!
     *  Evaluate the residual function for axisymmetric stagnation flow. If
     *  j == npos, the residual function is evaluated at all grid points.
     *  Otherwise, the residual function is only evaluated at grid points
     *  j-1, j, and j+1. This option is used to efficiently evaluate the
     *  Jacobian numerically.
     */
    virtual void eval(size_t j, doublereal* x, doublereal* r,
                      integer* mask, doublereal rdt);

    //! Evaluate all residual components at the right boundary.
    virtual void evalRightBoundary(double* x, double* res, int* diag,
                                   double rdt);

    //! Evaluate the residual corresponding to the continuity equation at all
    //! interior grid points.
    virtual void evalContinuity(size_t j, double* x, double* r,
                                int* diag, double rdt);

    //! Index of the species on the left boundary with the largest mass fraction
    size_t leftExcessSpecies() const {
        return m_kExcessLeft;
    }

    //! Index of the species on the right boundary with the largest mass fraction
    size_t rightExcessSpecies() const {
        return m_kExcessRight;
    }

protected:
    virtual AnyMap getMeta() const;
    virtual void setMeta(const AnyMap& state);

    doublereal wdot(size_t k, size_t j) const {
        return m_wdot(k,j);
    }

    //! Write the net production rates at point `j` into array `m_wdot`
    void getWdot(doublereal* x, size_t j) {
        setGas(x,j);
        m_kin->getNetProductionRates(&m_wdot(0,j));
    }

    //! Update the properties (thermo, transport, and diffusion flux).
    //! This function is called in eval after the points which need
    //! to be updated are defined.
    virtual void updateProperties(size_t jg, double* x, size_t jmin, size_t jmax);

    //! Evaluate the residual function. This function is called in eval
    //! after updateProperties is called.
    virtual void evalResidual(double* x, double* rsd, int* diag,
                              double rdt, size_t jmin, size_t jmax);

    /**
     * Update the thermodynamic properties from point j0 to point j1
     * (inclusive), based on solution x.
     */
    void updateThermo(const doublereal* x, size_t j0, size_t j1) {
        for (size_t j = j0; j <= j1; j++) {
            setGas(x,j);
            m_rho[j] = m_thermo->density();
            m_wtm[j] = m_thermo->meanMolecularWeight();
            m_cp[j] = m_thermo->cp_mass();
            m_thermo->getPartialMolarEnthalpies(&m_hk(0, j));
        }
    }

    //! @name Solution components
    //! @{

    doublereal T(const doublereal* x, size_t j) const {
        return x[index(c_offset_T, j)];
    }
    doublereal& T(doublereal* x, size_t j) {
        return x[index(c_offset_T, j)];
    }
    doublereal T_prev(size_t j) const {
        return prevSoln(c_offset_T, j);
    }

    doublereal rho_u(const doublereal* x, size_t j) const {
        return m_rho[j]*x[index(c_offset_U, j)];
    }

    doublereal u(const doublereal* x, size_t j) const {
        return x[index(c_offset_U, j)];
    }

    doublereal V(const doublereal* x, size_t j) const {
        return x[index(c_offset_V, j)];
    }
    doublereal V_prev(size_t j) const {
        return prevSoln(c_offset_V, j);
    }

    doublereal lambda(const doublereal* x, size_t j) const {
        return x[index(c_offset_L, j)];
    }

    doublereal Y(const doublereal* x, size_t k, size_t j) const {
        return x[index(c_offset_Y + k, j)];
    }

    doublereal& Y(doublereal* x, size_t k, size_t j) {
        return x[index(c_offset_Y + k, j)];
    }

    doublereal Y_prev(size_t k, size_t j) const {
        return prevSoln(c_offset_Y + k, j);
    }

    doublereal X(const doublereal* x, size_t k, size_t j) const {
        return m_wtm[j]*Y(x,k,j)/m_wt[k];
    }

    doublereal flux(size_t k, size_t j) const {
        return m_flux(k, j);
    }
    //! @}

    //! @name convective spatial derivatives.
    //!
    //! These use upwind differencing, assuming u(z) is negative
    //! @{
    doublereal dVdz(const doublereal* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (V(x,jloc) - V(x,jloc-1))/m_dz[jloc-1];
    }

    doublereal dYdz(const doublereal* x, size_t k, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (Y(x,k,jloc) - Y(x,k,jloc-1))/m_dz[jloc-1];
    }

    doublereal dTdz(const doublereal* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (T(x,jloc) - T(x,jloc-1))/m_dz[jloc-1];
    }
    //! @}

    doublereal shear(const doublereal* x, size_t j) const {
        doublereal c1 = m_visc[j-1]*(V(x,j) - V(x,j-1));
        doublereal c2 = m_visc[j]*(V(x,j+1) - V(x,j));
        return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    doublereal divHeatFlux(const doublereal* x, size_t j) const {
        doublereal c1 = m_tcon[j-1]*(T(x,j) - T(x,j-1));
        doublereal c2 = m_tcon[j]*(T(x,j+1) - T(x,j));
        return -2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    size_t mindex(size_t k, size_t j, size_t m) {
        return m*m_nsp*m_nsp + m_nsp*j + k;
    }

    //! Update the diffusive mass fluxes.
    virtual void updateDiffFluxes(const doublereal* x, size_t j0, size_t j1);

    //! Get the gradient of species specific molar enthalpies
    virtual void grad_hk(const double* x, size_t j);

    //---------------------------------------------------------
    //             member data
    //---------------------------------------------------------

    double m_press = -1.0; // pressure

    // grid parameters
    vector_fp m_dz;

    // mixture thermo properties
    vector_fp m_rho;
    vector_fp m_wtm;

    // species thermo properties
    vector_fp m_wt;
    vector_fp m_cp;

    // transport properties
    vector_fp m_visc;
    vector_fp m_tcon;
    vector_fp m_diff;
    vector_fp m_multidiff;
    Array2D m_dthermal;
    Array2D m_flux;

    //! Array of size #m_nsp by #m_points for saving molar enthalpies
    Array2D m_hk;

    //! Array of size #m_nsp by #m_points -1 for saving enthalpy fluxes
    Array2D m_dhk_dz;

    // production rates
    Array2D m_wdot;

    size_t m_nsp;

    ThermoPhase* m_thermo = nullptr;
    Kinetics* m_kin = nullptr;
    Transport* m_trans = nullptr;

    // boundary emissivities for the radiation calculations
    double m_epsilon_left = 0.0;
    double m_epsilon_right = 0.0;

    //! Indices within the ThermoPhase of the radiating species. First index is
    //! for CO2, second is for H2O.
    std::vector<size_t> m_kRadiating;

    // flags
    std::vector<bool> m_do_energy;
    bool m_do_soret = false;
    std::vector<bool> m_do_species;
    bool m_do_multicomponent = false;

    //! flag for the radiative heat loss
    bool m_do_radiation = false;

    //! radiative heat loss vector
    vector_fp m_qdotRadiation;

    // fixed T and Y values
    vector_fp m_fixedtemp;
    vector_fp m_zfix;
    vector_fp m_tfix;

    //! Index of species with a large mass fraction at each boundary, for which
    //! the mass fraction may be calculated as 1 minus the sum of the other mass
    //! fractions
    size_t m_kExcessLeft = 0;
    size_t m_kExcessRight = 0;

    bool m_dovisc;
    bool m_isFree;
    bool m_usesLambda;

    //! Update the transport properties at grid points in the range from `j0`
    //! to `j1`, based on solution `x`.
    virtual void updateTransport(doublereal* x, size_t j0, size_t j1);

public:
    //! Location of the point where temperature is fixed
    double m_zfixed = Undef;

    //! Temperature at the point used to fix the flame location
    double m_tfixed = -1.0;

private:
    vector_fp m_ybar;
};

}

#endif
