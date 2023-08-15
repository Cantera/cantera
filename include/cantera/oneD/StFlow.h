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

//! @defgroup flowGroup Flow Domains
//! One-dimensional flow domains.
//! @ingroup onedGroup

/**
 *  This class represents 1D flow domains that satisfy the one-dimensional
 *  similarity solution for chemically-reacting, axisymmetric flows.
 *  @ingroup flowGroup
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
    StFlow(shared_ptr<Solution> sol, const string& id="", size_t points=1);

    ~StFlow();

    string type() const override;

    //! @name Problem Specification
    //! @{

    void setupGrid(size_t n, const double* z) override;

    void resetBadValues(double* xg) override;

    ThermoPhase& phase() {
        return *m_thermo;
    }

    Kinetics& kinetics() {
        return *m_kin;
    }

    /**
     * Set the thermo manager.
     *
     * @deprecated To be removed after %Cantera 3.0 (unused)
     */
    void setThermo(ThermoPhase& th);

    void setKinetics(shared_ptr<Kinetics> kin) override;

    //! Set the kinetics manager.
    //! @deprecated To be removed after %Cantera 3.0;
    //!     replaced by Domain1D::setKinetics()
    void setKinetics(Kinetics& kin);

    void setTransport(shared_ptr<Transport> trans) override;

    //! Set transport model to existing instance
    //! @deprecated To be removed after %Cantera 3.0;
    //!     replaced by Domain1D::setKinetics()
    void setTransport(Transport& trans);

    //! Set the transport model
    //! @since New in %Cantera 3.0.
    void setTransportModel(const string& trans);

    //! Retrieve transport model
    //! @since New in %Cantera 3.0.
    string transportModel() const;

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
    void setPressure(double p) {
        m_press = p;
    }

    //! The current pressure [Pa].
    double pressure() const {
        return m_press;
    }

    //! Write the initial solution estimate into array x.
    void _getInitialSoln(double* x) override;

    void _finalize(const double* x) override;

    //! Sometimes it is desired to carry out the simulation using a specified
    //! temperature profile, rather than computing it by solving the energy
    //! equation. This method specifies this profile.
    void setFixedTempProfile(vector<double>& zfixed, vector<double>& tfixed) {
        m_zfix = zfixed;
        m_tfix = tfixed;
    }

    /**
     * Set the temperature fixed point at grid point j, and disable the energy
     * equation so that the solution will be held to this value.
     */
    void setTemperature(size_t j, double t) {
        m_fixedtemp[j] = t;
        m_do_energy[j] = false;
    }

    //! The fixed temperature value at point j.
    double T_fixed(size_t j) const {
        return m_fixedtemp[j];
    }

    //! @}

    string componentName(size_t n) const override;

    size_t componentIndex(const string& name) const override;

    //! Returns true if the specified component is an active part of the solver state
    virtual bool componentActive(size_t n) const;

    //! Print the solution.
    void show(const double* x) override;

    shared_ptr<SolutionArray> asArray(const double* soln) const override;
    void fromArray(SolutionArray& arr, double* soln) override;

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
    //! @deprecated To be removed after %Cantera 3.0; replaced by type().
    virtual string flowType() const;

    void solveEnergyEqn(size_t j=npos);

    //! Get the solving stage (used by IonFlow specialization)
    //! @since New in %Cantera 3.0
    virtual size_t getSolvingStage() const;

    //! Solving stage mode for handling ionized species (used by IonFlow specialization)
    //! - @c stage=1: the fluxes of charged species are set to zero
    //! - @c stage=2: the electric field equation is solved, and the drift flux for
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
     * The simple radiation model used was established by Liu and Rogg
     * @cite liu1991. This model considers the radiation of CO2 and H2O.
     *
     * This model uses the optically thin limit and the gray-gas approximation to
     * simply calculate a volume specified heat flux out of the Planck absorption
     * coefficients, the boundary emissivities and the temperature. Polynomial lines
     * calculate the species Planck coefficients for H2O and CO2. The data for the
     * lines are taken from the RADCAL program @cite RADCAL.
     * The coefficients for the polynomials are taken from
     * [TNF Workshop](https://tnfworkshop.org/radiation/) material.
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
    void resize(size_t components, size_t points) override;

    //! Set the gas object state to be consistent with the solution at point j.
    void setGas(const double* x, size_t j);

    //! Set the gas state to be consistent with the solution at the midpoint
    //! between j and j + 1.
    void setGasAtMidpoint(const double* x, size_t j);

    double density(size_t j) const {
        return m_rho[j];
    }

    //! @deprecated To be removed after %Cantera 3.0. Superseded by isFree()
    virtual bool fixed_mdot();

    /**
     * Retrieve flag indicating whether flow is freely propagating.
     * The flow is unstrained and the axial mass flow rate is not specified.
     * For free flame propagation, the axial velocity is determined by the solver.
     * @since New in %Cantera 3.0
     */
    bool isFree() const {
        return m_isFree;
    }

    /**
     * Retrieve flag indicating whether flow uses radial momentum.
     * If `true`, radial momentum equation for @f$ V @f$ as well as
     * @f$ d\Lambda/dz = 0 @f$ are solved; if `false`, @f$ \Lambda(z) = 0 @f$ and
     * @f$ V(z) = 0 @f$ by definition.
     * @since New in %Cantera 3.0
     */
    bool isStrained() const {
        return m_usesLambda;
    }

    void setViscosityFlag(bool dovisc) {
        m_dovisc = dovisc;
    }

    /**
     *  Evaluate the residual function for axisymmetric stagnation flow. If
     *  j == npos, the residual function is evaluated at all grid points.
     *  Otherwise, the residual function is only evaluated at grid points
     *  j-1, j, and j+1. This option is used to efficiently evaluate the
     *  Jacobian numerically.
     */
    void eval(size_t j, double* x, double* r, integer* mask, double rdt) override;

    //! Evaluate all residual components at the right boundary.
    virtual void evalRightBoundary(double* x, double* res, int* diag, double rdt);

    //! Evaluate the residual corresponding to the continuity equation at all
    //! interior grid points.
    virtual void evalContinuity(size_t j, double* x, double* r, int* diag, double rdt);

    //! Index of the species on the left boundary with the largest mass fraction
    size_t leftExcessSpecies() const {
        return m_kExcessLeft;
    }

    //! Index of the species on the right boundary with the largest mass fraction
    size_t rightExcessSpecies() const {
        return m_kExcessRight;
    }

protected:
    AnyMap getMeta() const override;
    void setMeta(const AnyMap& state) override;

    double wdot(size_t k, size_t j) const {
        return m_wdot(k,j);
    }

    //! Write the net production rates at point `j` into array `m_wdot`
    void getWdot(double* x, size_t j) {
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
    void updateThermo(const double* x, size_t j0, size_t j1) {
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

    double T(const double* x, size_t j) const {
        return x[index(c_offset_T, j)];
    }
    double& T(double* x, size_t j) {
        return x[index(c_offset_T, j)];
    }
    double T_prev(size_t j) const {
        return prevSoln(c_offset_T, j);
    }

    double rho_u(const double* x, size_t j) const {
        return m_rho[j]*x[index(c_offset_U, j)];
    }

    double u(const double* x, size_t j) const {
        return x[index(c_offset_U, j)];
    }

    double V(const double* x, size_t j) const {
        return x[index(c_offset_V, j)];
    }
    double V_prev(size_t j) const {
        return prevSoln(c_offset_V, j);
    }

    double lambda(const double* x, size_t j) const {
        return x[index(c_offset_L, j)];
    }

    double Y(const double* x, size_t k, size_t j) const {
        return x[index(c_offset_Y + k, j)];
    }

    double& Y(double* x, size_t k, size_t j) {
        return x[index(c_offset_Y + k, j)];
    }

    double Y_prev(size_t k, size_t j) const {
        return prevSoln(c_offset_Y + k, j);
    }

    double X(const double* x, size_t k, size_t j) const {
        return m_wtm[j]*Y(x,k,j)/m_wt[k];
    }

    double flux(size_t k, size_t j) const {
        return m_flux(k, j);
    }
    //! @}

    //! @name convective spatial derivatives.
    //!
    //! These use upwind differencing, assuming u(z) is negative
    //! @{
    double dVdz(const double* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (V(x,jloc) - V(x,jloc-1))/m_dz[jloc-1];
    }

    double dYdz(const double* x, size_t k, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (Y(x,k,jloc) - Y(x,k,jloc-1))/m_dz[jloc-1];
    }

    double dTdz(const double* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (T(x,jloc) - T(x,jloc-1))/m_dz[jloc-1];
    }
    //! @}

    double shear(const double* x, size_t j) const {
        double c1 = m_visc[j-1]*(V(x,j) - V(x,j-1));
        double c2 = m_visc[j]*(V(x,j+1) - V(x,j));
        return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    double divHeatFlux(const double* x, size_t j) const {
        double c1 = m_tcon[j-1]*(T(x,j) - T(x,j-1));
        double c2 = m_tcon[j]*(T(x,j+1) - T(x,j));
        return -2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    size_t mindex(size_t k, size_t j, size_t m) {
        return m*m_nsp*m_nsp + m_nsp*j + k;
    }

    //! Update the diffusive mass fluxes.
    virtual void updateDiffFluxes(const double* x, size_t j0, size_t j1);

    //! Get the gradient of species specific molar enthalpies
    virtual void grad_hk(const double* x, size_t j);

    //---------------------------------------------------------
    //             member data
    //---------------------------------------------------------

    double m_press = -1.0; // pressure

    // grid parameters
    vector<double> m_dz;

    // mixture thermo properties
    vector<double> m_rho;
    vector<double> m_wtm;

    // species thermo properties
    vector<double> m_wt;
    vector<double> m_cp;

    // transport properties
    vector<double> m_visc;
    vector<double> m_tcon;
    vector<double> m_diff;
    vector<double> m_multidiff;
    Array2D m_dthermal;
    Array2D m_flux;

    //! Array of size #m_nsp by #m_points for saving molar enthalpies
    Array2D m_hk;

    //! Array of size #m_nsp by #m_points-1 for saving enthalpy fluxes
    Array2D m_dhk_dz;

    // production rates
    Array2D m_wdot;

    size_t m_nsp; //!< Number of species in the mechanism

    ThermoPhase* m_thermo = nullptr;
    Kinetics* m_kin = nullptr;
    Transport* m_trans = nullptr;

    // boundary emissivities for the radiation calculations
    double m_epsilon_left = 0.0;
    double m_epsilon_right = 0.0;

    //! Indices within the ThermoPhase of the radiating species. First index is
    //! for CO2, second is for H2O.
    vector<size_t> m_kRadiating;

    // flags
    vector<bool> m_do_energy;
    bool m_do_soret = false;
    vector<bool> m_do_species;
    bool m_do_multicomponent = false;

    //! flag for the radiative heat loss
    bool m_do_radiation = false;

    //! radiative heat loss vector
    vector<double> m_qdotRadiation;

    // fixed T and Y values
    vector<double> m_fixedtemp;
    vector<double> m_zfix;
    vector<double> m_tfix;

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
    virtual void updateTransport(double* x, size_t j0, size_t j1);

public:
    //! Location of the point where temperature is fixed
    double m_zfixed = Undef;

    //! Temperature at the point used to fix the flame location
    double m_tfixed = -1.0;

private:
    vector<double> m_ybar;
};

}

#endif
