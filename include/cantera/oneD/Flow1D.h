//! @file Flow1D.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOW1D_H
#define CT_FLOW1D_H

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
    c_offset_U   //! axial velocity [m/s]
    , c_offset_V //! strain rate
    , c_offset_T //! temperature [kelvin]
    , c_offset_L //! (1/r)dP/dr
    , c_offset_E //! electric field
    , c_offset_Uo //! oxidizer axial velocity [m/s]
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
class Flow1D : public Domain1D
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
    Flow1D(ThermoPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    //! Delegating constructor
    Flow1D(shared_ptr<ThermoPhase> th, size_t nsp = 1, size_t points = 1);

    //! Create a new flow domain.
    //! @param sol  Solution object used to evaluate all thermodynamic, kinetic, and
    //!     transport properties
    //! @param id  name of flow domain
    //! @param points  initial number of grid points
    Flow1D(shared_ptr<Solution> sol, const string& id="", size_t points=1);

    ~Flow1D();

    string domainType() const override;

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

    void setKinetics(shared_ptr<Kinetics> kin) override;

    void setTransport(shared_ptr<Transport> trans) override;

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

    //! Compute species diffusive fluxes with respect to
    //! their mass fraction gradients (fluxGradientBasis = ThermoBasis::mass)
    //! or mole fraction gradients (fluxGradientBasis = ThermoBasis::molar, default)
    //! when using the mixture-averaged transport model.
    //! @param fluxGradientBasis set flux computation to mass or mole basis
    //! @since New in %Cantera 3.1.
    void setFluxGradientBasis(ThermoBasis fluxGradientBasis);

    //! Compute species diffusive fluxes with respect to
    //! their mass fraction gradients (fluxGradientBasis = ThermoBasis::mass)
    //! or mole fraction gradients (fluxGradientBasis = ThermoBasis::molar, default)
    //! when using the mixture-averaged transport model.
    //! @return the basis used for flux computation (mass or mole fraction gradients)
    //! @since New in %Cantera 3.1.
    ThermoBasis fluxGradientBasis() const {
        return m_fluxGradientBasis;
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

    void show(const double* x) override;

    shared_ptr<SolutionArray> asArray(const double* soln) const override;
    void fromArray(SolutionArray& arr, double* soln) override;

    //! Set flow configuration for freely-propagating flames, using an internal point
    //! with a fixed temperature as the condition to determine the inlet mass flux.
    void setFreeFlow() {
        m_dovisc = false;
        m_isFree = true;
        m_usesLambda = false;
    }

    //! Set flow configuration for axisymmetric counterflow flames, using specified
    //! inlet mass fluxes.
    void setAxisymmetricFlow() {
        m_dovisc = true;
        m_isFree = false;
        m_usesLambda = true;
    }

    //! Set flow configuration for burner-stabilized flames, using specified inlet mass
    //! fluxes.
    void setUnstrainedFlow() {
        m_dovisc = false;
        m_isFree = false;
        m_usesLambda = false;
    }

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
    double leftEmissivity() const {
        return m_epsilon_left;
    }

    //! Return emissivity at right boundary
    double rightEmissivity() const {
        return m_epsilon_right;
    }

    void fixTemperature(size_t j=npos);

    /**
     * @name Two-Point control method
     *
     * In this method two control points are designated in the 1D domain, and the value
     * of the temperature at these points is fixed. The values of the control points are
     * imposed and thus serve as a boundary condition that affects the solution of the
     * governing equations in the 1D domain. The imposition of fixed points in the
     * domain means that the original set of governing equations' boundary conditions
     * would over-determine the problem. Thus, the boundary conditions are changed to
     * reflect the fact that the control points are serving as internal boundary
     * conditions.
     *
     * The imposition of the two internal boundary conditions requires that two other
     * boundary conditions be changed. The first is the boundary condition for the
     * continuity equation at the left boundary, which is changed to be a value that is
     * derived from the solution at the left boundary. The second is the continuity
     * boundary condition at the right boundary, which is also determined from the flow
     * solution by using the oxidizer axial velocity equation variable to compute the
     * mass flux at the right boundary.
     *
     * This method is based on the work of Nishioka et al. @cite nishioka1996 .
     */
    //! @{

    //! Returns the temperature at the left control point
    double leftControlPointTemperature() const;

    //! Returns the z-coordinate of the left control point
    double leftControlPointCoordinate() const;

    //! Sets the temperature of the left control point
    void setLeftControlPointTemperature(double temperature);

    //! Sets the coordinate of the left control point
    void setLeftControlPointCoordinate(double z_left);

    //! Returns the temperature at the right control point
    double rightControlPointTemperature() const;

    //! Returns the z-coordinate of the right control point
    double rightControlPointCoordinate() const;

    //! Sets the temperature of the right control point
    void setRightControlPointTemperature(double temperature);

    //! Sets the coordinate of the right control point
    void setRightControlPointCoordinate(double z_right);

    //! Sets the status of the two-point control
    void enableTwoPointControl(bool twoPointControl);

    //! Returns the status of the two-point control
    bool twoPointControlEnabled() const {
        return m_twoPointControl;
    }
    //! @}

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
     * Evaluate the residual functions for axisymmetric stagnation flow.
     * If jGlobal == npos, the residual function is evaluated at all grid points.
     * Otherwise, the residual function is only evaluated at grid points j-1, j,
     * and j+1. This option is used to efficiently evaluate the Jacobian numerically.
     *
     * These residuals at all the boundary grid points are evaluated using a default
     * boundary condition that may be modified by a boundary object that is attached
     * to the domain. The boundary object connected will modify these equations by
     * subtracting the boundary object's values for V, T, mdot, etc. As a result,
     * these residual equations will force the solution variables to the values of
     * the connected boundary object.
     *
     *  @param jGlobal  Global grid point at which to update the residual
     *  @param[in] xGlobal  Global state vector
     *  @param[out] rsdGlobal  Global residual vector
     *  @param[out] diagGlobal  Global boolean mask indicating whether each solution
     *      component has a time derivative (1) or not (0).
     *  @param[in] rdt Reciprocal of the timestep (`rdt=0` implies steady-state.)
     */
    void eval(size_t jGlobal, double* xGlobal, double* rsdGlobal,
              integer* diagGlobal, double rdt) override;

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

    //! @name Updates of cached properties
    //! These methods are called by eval() to update cached properties and data that are
    //! used for the evaluation of the governing equations.
    //! @{

    /**
     * Update the thermodynamic properties from point j0 to point j1
     * (inclusive), based on solution x.
     *
     * The gas state is set to be consistent with the solution at the
     * points from j0 to j1.
     *
     * Properties that are computed and cached are:
     * * #m_rho (density)
     * * #m_wtm (mean molecular weight)
     * * #m_cp (specific heat capacity)
     * * #m_hk (species specific enthalpies)
     * * #m_wdot (species production rates)
     */
    void updateThermo(const double* x, size_t j0, size_t j1) {
        for (size_t j = j0; j <= j1; j++) {
            setGas(x,j);
            m_rho[j] = m_thermo->density();
            m_wtm[j] = m_thermo->meanMolecularWeight();
            m_cp[j] = m_thermo->cp_mass();
            m_thermo->getPartialMolarEnthalpies(&m_hk(0, j));
            m_kin->getNetProductionRates(&m_wdot(0, j));
        }
    }

    //! Update the transport properties at grid points in the range from `j0`
    //! to `j1`, based on solution `x`. Evaluates the solution at the midpoint
    //! between `j` and `j + 1` to compute the transport properties. For example,
    //! the viscosity at element `j`  is the viscosity evaluated
    //! at the midpoint between `j` and `j + 1`.
    virtual void updateTransport(double* x, size_t j0, size_t j1);

    //! Update the diffusive mass fluxes.
    virtual void updateDiffFluxes(const double* x, size_t j0, size_t j1);

    //! Update the properties (thermo, transport, and diffusion flux).
    //! This function is called in eval after the points which need
    //! to be updated are defined.
    virtual void updateProperties(size_t jg, double* x, size_t jmin, size_t jmax);

    /**
     * Computes the radiative heat loss vector over points jmin to jmax and stores
     * the data in the qdotRadiation variable.
     *
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
    void computeRadiation(double* x, size_t jmin, size_t jmax);

    //! @}

    //! @name Governing Equations
    //! Methods called by eval() to calculate residuals for individual governing
    //! equations.

    /**
     * Evaluate the continuity equation residual.
     *
     * This function calculates the residual of the continuity equation
     * @f[
     *     \frac{d(\rho u)}{dz} + 2\rho V = 0
     * @f]
     *
     * Axisymmetric flame:
     *  The continuity equation propagates information from right-to-left.
     *  The @f$ \rho u @f$ at point 0 is dependent on @f$ \rho u @f$ at point 1,
     *  but not on @f$ \dot{m} @f$ from the inlet.
     *
     * Freely-propagating flame:
     *  The continuity equation propagates information away from a fixed temperature
     *  point that is set in the domain.
     *
     * Unstrained flame:
     *  A specified mass flux; the main example being burner-stabilized flames.
     *
     * The default boundary condition for the continuity equation is
     * (@f$ u = 0 @f$) at the right boundary. Because the equation is a first order
     * equation, only one boundary condition is needed.
     *
     * @param[in] x Local domain state vector, includes variables like temperature,
     *               density, etc.
     * @param[out] rsd Local domain residual vector that stores the continuity
     *                  equation residuals.
     * @param[out] diag Local domain diagonal matrix that controls whether an entry
     *                   has a time-derivative (used by the solver).
     * @param[in] rdt Reciprocal of the timestep.
     * @param[in] jmin The index for the starting point in the local domain grid.
     * @param[in] jmax The index for the ending point in the local domain grid.
     */
    virtual void evalContinuity(double* x, double* rsd, int* diag,
                                double rdt, size_t jmin, size_t jmax);

    /**
     * Evaluate the momentum equation residual.
     *
     * The function calculates the radial momentum equation defined as
     * @f[
     *    \rho u \frac{dV}{dz} + \rho V^2 =
     *    \frac{d}{dz}\left( \mu \frac{dV}{dz} \right) - \Lambda
     * @f]
     *
     * The radial momentum equation is used for axisymmetric flows, and incorporates
     * terms for time and spatial variations of radial velocity (@f$ V @f$). The
     * default boundary condition is zero radial velocity (@f$ V @f$) at the left
     * and right boundary.
     *
     * For argument explanation, see evalContinuity().
     */
    virtual void evalMomentum(double* x, double* rsd, int* diag,
                              double rdt, size_t jmin, size_t jmax);

    /**
     * Evaluate the lambda equation residual.
     *
     * The function calculates the lambda equation as
     * @f[
     *    \frac{d\Lambda}{dz} = 0
     * @f]
     *
     * The lambda equation serves as an eigenvalue that allows the momentum
     * equation and continuity equations to be simultaneously satisfied in
     * axisymmetric flows. The lambda equation propagates information from
     * left-to-right. The default boundary condition is @f$ \Lambda = 0 @f$
     * at the left boundary. The equation is first order and so only one
     * boundary condition is needed.
     *
     * For argument explanation, see evalContinuity().
     */
    virtual void evalLambda(double* x, double* rsd, int* diag,
                            double rdt, size_t jmin, size_t jmax);

    /**
     * Evaluate the energy equation residual.
     *
     * The function calculates the energy equation:
     * @f[
     *   \rho c_p u \frac{dT}{dz} =
     *   \frac{d}{dz}\left( \lambda \frac{dT}{dz} \right)
     *   - \sum_k h_kW_k\dot{\omega}_k
     *   - \sum_k  j_k \frac{dh_k}{dz}
     * @f]
     *
     * The energy equation includes contributions from
     * chemical reactions and diffusion. Default is zero temperature (@f$ T @f$)
     * at the left and right boundaries. These boundary values are updated by the
     * specific boundary object connected to the domain.
     *
     * For argument explanation, see evalContinuity().
     */
    virtual void evalEnergy(double* x, double* rsd, int* diag,
                            double rdt, size_t jmin, size_t jmax);

    /**
     * Evaluate the species equations' residuals.
     *
     * The function calculates the species equations as
     * @f[
     *    \rho u \frac{dY_k}{dz} + \frac{dj_k}{dz} = W_k\dot{\omega}_k
     * @f]
     *
     * The species equations include terms for temporal and spatial variations
     * of species mass fractions (@f$ Y_k @f$). The default boundary condition is zero
     * flux for species at the left and right boundary.
     *
     * For argument explanation, see evalContinuity().
     */
    virtual void evalSpecies(double* x, double* rsd, int* diag,
                             double rdt, size_t jmin, size_t jmax);

    /**
     * Evaluate the electric field equation residual to be zero everywhere.
     *
     * The electric field equation is implemented in the IonFlow class. The default
     * boundary condition is zero electric field (@f$ E @f$) at the boundary,
     * and @f$ E @f$ is zero within the domain.
     *
     * For argument explanation, see evalContinuity().
     */
    virtual void evalElectricField(double* x, double* rsd, int* diag,
                                   double rdt, size_t jmin, size_t jmax);

    //! @} End of Governing Equations

    //! Alternate version of evalContinuity with legacy signature.
    //! Implemented by StFlow; included here to prevent compiler warnings about shadowed
    //! virtual functions.
    //! @deprecated To be removed after %Cantera 3.1.
    virtual void evalContinuity(size_t j, double* x, double* r, int* diag, double rdt);

    /**
     * Evaluate the oxidizer axial velocity equation residual.
     *
     * The function calculates the oxidizer axial velocity equation as
     * @f[
     *    \frac{dU_{o}}{dz} = 0
     * @f]
     *
     * This equation serves as a dummy equation that is used only in the context of
     * two-point flame control, and serves as the way for two interior control points to
     * be specified while maintaining block tridiagonal structure. The default boundary
     * condition is @f$ U_o = 0 @f$ at the right and zero flux at the left boundary.
     *
     * For argument explanation, see evalContinuity().
     */
    virtual void evalUo(double* x, double* rsd, int* diag,
                        double rdt, size_t jmin, size_t jmax);

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

    //! Solution component for oxidizer velocity, @see evalUo
    double Uo(const double* x, size_t j) const {
        return x[index(c_offset_Uo, j)];
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

    //! @name Convective spatial derivatives
    //!
    //! These methods use upwind differencing to calculate spatial derivatives
    //! for velocity, species mass fractions, and temperature. Upwind differencing
    //! is a numerical discretization method that considers the direction of the
    //! flow to improve stability.
    //! @{

    /**
     * This function calculates the spatial derivative of velocity V with respect
     * to z at point j using upwind differencing.
     *
     * @f[
     *   \frac{\partial V}{\partial z} \bigg|_{j} \approx \frac{V(x, j_{\text{loc}}) -
     *   V(x, j_{\text{loc}} - 1)}{m_{\text{dz}}[j_{\text{loc}} - 1]}
     * @f]
     *
     * Where the value of loc is determined by the sign of the axial velocity.
     * If the axial velocity is positive, the value of loc is j. If the axial velocity
     * is negative, the value of loc is j + 1. A positive velocity means that the flow
     * is moving left-to-right.
     *
     * @param[in] x The local domain state vector.
     * @param[in] j The index at which the derivative is computed.
     */
    double dVdz(const double* x, size_t j) const {
        size_t jloc = (u(x, j) > 0.0 ? j : j + 1);
        return (V(x, jloc) - V(x, jloc-1))/m_dz[jloc-1];
    }

    /**
     * This function calculates the spatial derivative of the species mass fraction
     * @f$ Y_k @f$ with respect to z for species k at point j using upwind
     * differencing.
     *
     * @f[
     *   \frac{\partial Y_k}{\partial z} \bigg|_{j} \approx \frac{Y(x, k, j_{\text{loc}})
     *   - Y(x, k, j_{\text{loc}} - 1)}{m_{\text{dz}}[j_{\text{loc}} - 1]}
     * @f]
     *
     * Where the value of loc is determined by the sign of the axial velocity.
     * If the axial velocity is positive, the value of loc is j. If the axial velocity
     * is negative, the value of loc is j + 1. A positive velocity means that the flow
     * is moving left-to-right.
     *
     * @param[in] x The local domain state vector.
     * @param[in] k The species index.
     * @param[in] j The index at which the derivative is computed.
     */
    double dYdz(const double* x, size_t k, size_t j) const {
        size_t jloc = (u(x, j) > 0.0 ? j : j + 1);
        return (Y(x, k, jloc) - Y(x, k, jloc-1))/m_dz[jloc-1];
    }

    /**
     * This function calculates the spatial derivative of temperature T
     * with respect to z at point j using upwind differencing.
     *
     * @f[
     *   \frac{\partial T}{\partial z} \bigg|_{j} \approx \frac{T(x, j_{\text{loc}}) -
     *     T(x, j_{\text{loc}} - 1)}{m_{\text{dz}}[j_{\text{loc}} - 1]}
     * @f]
     *
     * Where the value of loc is determined by the sign of the axial velocity.
     * If the axial velocity is positive, the value of loc is j. If the axial velocity
     * is negative, the value of loc is j + 1. A positive velocity means that the flow
     * is moving left-to-right.
     *
     * @param[in] x The local domain state vector.
     * @param[in] j The index at which the derivative is computed.
     */
    double dTdz(const double* x, size_t j) const {
        size_t jloc = (u(x, j) > 0.0 ? j : j + 1);
        return (T(x, jloc) - T(x, jloc-1))/m_dz[jloc-1];
    }
    //! @}

    /**
     * Compute the shear term from the momentum equation using a central
     * three-point differencing scheme. This term is discretized
     *
     * The term to be discretized is:
     * @f[
     *   \frac{d}{dz}\left(\mu \frac{dV}{dz}\right)
     * @f]
     *
     * Let @f$ A = \mu \frac{dV}{dz} @f$ for simplicity. We'll call this the inner
     * term or inner discretization. In this situation, the inner term is evaluated
     * using a central difference formula, but evaluated at j+1/2 and j-1/2 (halfway
     * between the grid points around point j).
     *
     * The value of @f$ A @f$ at point @f$ j-1/2 @f$ is estimating using a central difference
     * formula:
     *
     * @f[
     *  A_{j-1/2} = \mu_{j-1/2} \frac{V_j - V_{j-1}}{z_j - z_{j-1}}
     * @f]
     *
     * The value of @f$ A @f$ at point @f$ j+1/2 @f$ is:
     *
     * @f[
     *   A_{j+1/2} = \mu_{j+1/2} \frac{V_{j+1} - V_j}{z_{j+1} - z_j}
     * @f]
     *
     * In the implementation, a vector of viscosities are used which actually represents
     * the viscosity evaluated at the midpoints between grid points. In other words,
     * the value of mu[j] is actually @f$ \mu_{j+1/2} @f$. The transport properties
     * between grid points are not the average of the adjacent grid points, but rather,
     * the state of the system at the midpoint is estimated, and then the transport properties
     * are evaluated at that state (T,P,Y) for example.
     *
     * The outer discretization uses a central difference between the j+1/2 and j-1/2
     * locations.
     *
     * @f[
     *  \frac{dA}{dz} \approx \frac{A_{j+1/2} - A_{j-1/2}}{z_{j+1/2} - z_{j-1/2}}
     * @f]
     *
     * Where the values of @f$ z @f$ are: @f$ z_{j+1/2} = z_{j} + 0.5*(z_{j+1} - z_j) = 0.5*(z_{j} + z_{j+1}) @f$ and
     * @f$ z_{j-1/2} = z_j - 0.5*(z_j - z_{j-1}) = 0.5*(z_{j} + z_{j-1}) @f$. The difference between these two values is
     * @f$ z_{j+1/2} - z_{j-1/2} = \frac{z_{j+1} - z_{j-1}}{2} @f$.
     *
     * Substituting these values into the central difference formula gives:
     *
     * @f[
     *   \frac{d}{dz}\left(\mu \frac{dV}{dz}\right) \approx \frac{\mu_{j+1/2} \frac{V_{j+1} - V_j}{z_{j+1} - z_j} -
     *     \mu_{j-1/2} \frac{V_j - V_{j-1}}{z_j - z_{j-1}}}{\frac{z_{j+1} - z_{j-1}}{2}}
     * @f]
     *
     * @param[in] x The local domain state vector.
     * @param[in] j The index at which the derivative is computed.
     */
    double shear(const double* x, size_t j) const {
        double A_left = m_visc[j-1]*(V(x, j) - V(x, j-1)) / (z(j) - z(j-1));
        double A_right = m_visc[j]*(V(x, j+1) - V(x, j)) / (z(j+1) - z(j));
        return 2.0*(A_right - A_left) / (z(j+1) - z(j-1));
    }

    /**
     * Compute the conduction term from the energy equation using a central
     * three-point differencing scheme.
     *
     * The term to be discretized is:
     * @f[
     *   \frac{d}{dz}\left(\lambda \frac{dT}{dz}\right)
     * @f]
     *
     * Let @f$ A = \lambda \frac{dT}{dz} @f$ for simplicity. We'll call this the inner
     * term or inner discretization. In this situation, the inner term is evaluated
     * using a central difference formula, but evaluated at j+1/2 and j-1/2 (halfway
     * between the grid points around point j).
     *
     * The value of @f$ A @f$ at point @f$ j-1/2 @f$ is estimating using a central
     * difference formula:
     *
     * @f[
     *  A_{j-1/2} = \lambda_{j-1/2} \frac{T_j - T_{j-1}}{z_j - z_{j-1}}
     * @f]
     *
     * The value of @f$ A @f$ at point @f$ j+1/2 @f$ is:
     *
     * @f[
     *   A_{j+1/2} = \lambda_{j+1/2} \frac{T_{j+1} - T_j}{z_{j+1} - z_j}
     * @f]
     *
     * In the implementation, a vector of thermal conductivities are used which
     * actually represents the thermal conductivity evaluated at the midpoints between
     * grid points. In other words, the value of lambda[j] is actually
     * @f$ \lambda_{j+1/2} @f$. The transport properties between grid points are not
     * the average of the adjacent grid points, but rather, the state of the system at
     * the midpoint is estimated, and then the transport properties are evaluated at
     * that state (T,P,Y) for example.
     *
     * The outer discretization uses a central difference between the j+1/2 and j-1/2
     * locations.
     *
     * @f[
     *  \frac{dA}{dz} \approx \frac{A_{j+1/2} - A_{j-1/2}}{z_{j+1/2} - z_{j-1/2}}
     * @f]
     *
     * Where the values of @f$ z @f$ are:
     * @f$ z_{j+1/2} = z_{j} + 0.5*(z_{j+1} - z_j) = 0.5*(z_{j} + z_{j+1}) @f$ and
     * @f$ z_{j-1/2} = z_j - 0.5*(z_j - z_{j-1}) = 0.5*(z_{j} + z_{j-1}) @f$. The
     * difference between these two values is
     * @f$ z_{j+1/2} - z_{j-1/2} = \frac{z_{j+1} - z_{j-1}}{2} @f$.
     *
     * Substituting these values into the central difference formula gives:
     *
     * @f[
     *   \frac{d}{dz}\left(\lambda \frac{dT}{dz}\right) \approx \frac{\lambda_{j+1/2}
     *      \frac{T_{j+1} - T_j}{z_{j+1} - z_j} -
     *      \lambda_{j-1/2} \frac{T_j - T_{j-1}}{z_j - z_{j-1}}}{\frac{z_{j+1} - z_{j-1}}{2}}
     * @f]
     *
     * @param[in] x The local domain state vector.
     * @param[in] j The index at which the derivative is computed.
     */
    double conduction(const double* x, size_t j) const {
        double A_left = m_tcon[j-1]*(T(x, j) - T(x, j-1)) / (z(j) - z(j-1));
        double A_right = m_tcon[j]*(T(x, j+1) - T(x, j)) / (z(j+1) - z(j));
        return -2.0*(A_right - A_left) / (z(j+1) - z(j-1));
    }

    /**
     * Array access mapping for a 3D array stored in a 1D vector. Used for
     * accessing data in the m_multidiff member variable.
     *
     *
     * @param[in] k First species index.
     * @param[in] j The grid point index.
     * @param[in] m The second species index.
     */
    size_t mindex(size_t k, size_t j, size_t m) {
        return m*m_nsp*m_nsp + m_nsp*j + k;
    }

    /**
     * Compute the spatial derivative of species specific molar enthalpies using upwind
     * differencing. Updates all species molar enthalpies for all species at point j.
     * Does not return a value, but updates the m_dhk_dz 2D array.
     *
     * This function calculates the spatial derivative of the species specific molar
     * enthalpies @f$ h_k @f$ with respect to z for species k at point j using upwind
     * differencing.
     *
     * @f[
     *   \frac{\partial h_k}{\partial z} \bigg|_{j} \approx \frac{h_k(x, k, j_{\text{loc}})
     *   - h_k(x, k, j_{\text{loc}} - 1)}{m_{\text{dz}}[j_{\text{loc}} - 1]}
     * @f]
     *
     * Where the value of loc is determined by the sign of the axial velocity.
     * If the axial velocity is positive, the value of loc is j. If the axial velocity
     * is negative, the value of loc is j + 1. A positive velocity means that the flow
     * is moving left-to-right.
     *
     * @param[in] k The species index.
     * @param[in] j The index at which the derivative is computed.
     */
    virtual void grad_hk(const double* x, size_t j);

    //---------------------------------------------------------
    //             member data
    //---------------------------------------------------------

    double m_press = -1.0; // pressure

    // Grid spacing. Element j holds the value of grid(j+1) - grid(j)
    vector<double> m_dz;

    // mixture thermo properties
    vector<double> m_rho;  //!< Vector of size #m_nsp to cache densities
    vector<double> m_wtm;  //!< Vector of size #m_nsp to cache mean molecular weights

    // species thermo properties
    vector<double> m_wt;
    vector<double> m_cp;  //!< Vector of size #m_nsp to cache specific heat capacities

    // transport properties
    vector<double> m_visc;
    vector<double> m_tcon;

    //! Vector of size #m_nsp by #m_points for saving density times diffusion
    //! coefficient times species molar mass divided by mean molecular weight
    vector<double> m_diff;

    //! Vector of size #m_nsp by #m_nsp by #m_points for saving multicomponent diffusion
    //! coefficients.
    vector<double> m_multidiff;

    //! Array of size #m_nsp by #m_points for saving thermal diffusion coefficients
    Array2D m_dthermal;

    //! Array of size #m_nsp by #m_points for saving diffusive mass fluxes
    Array2D m_flux;

    //! Array of size #m_nsp by #m_points for saving molar enthalpies
    Array2D m_hk;

    //! Array of size #m_nsp by #m_points-1 for saving enthalpy fluxes
    Array2D m_dhk_dz;

    //! Array of size #m_nsp by #m_points for saving species production rates
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
    ThermoBasis m_fluxGradientBasis = ThermoBasis::molar;
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

    //! Flag for activating two-point flame control
    bool m_twoPointControl = false;

    //! Location of the left control point when two-point control is enabled
    double m_zLeft = Undef;

    //! Temperature of the left control point when two-point control is enabled
    double m_tLeft = Undef;

    //! Location of the right control point when two-point control is enabled
    double m_zRight = Undef;

    //! Temperature of the right control point when two-point control is enabled
    double m_tRight = Undef;

public:
    //! Location of the point where temperature is fixed
    double m_zfixed = Undef;

    //! Temperature at the point used to fix the flame location
    double m_tfixed = -1.0;

private:
    //! Holds the average of the species mass fractions between grid points j and j+1.
    //! Used when building a gas state at the grid midpoints for evaluating transport
    //! properties at the midpoints.
    vector<double> m_ybar;
};

}

#endif
