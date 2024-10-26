/**
 * @file Boundary1D.h
 *
 * Boundary objects for one-dimensional simulations.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BOUNDARY1D_H
#define CT_BOUNDARY1D_H

#include "Domain1D.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "Flow1D.h"

namespace Cantera
{

//! Unique identifier for the left inlet.
const int LeftInlet = 1;

//! Unique identifier for the right inlet.
const int RightInlet = -1;

//! @defgroup bdryGroup Boundaries
//! Boundaries of one-dimensional flow domains.
//! @ingroup onedGroup
//! @{

/**
 * The base class for boundaries between one-dimensional spatial domains. The
 * boundary may have its own internal variables, such as surface species
 * coverages.
 *
 * The boundary types are an inlet, an outlet, a symmetry plane, and a surface.
 *
 * The public methods are all virtual, and the base class implementations throw
 * exceptions.
 */
class Boundary1D : public Domain1D
{
public:
    //! Default constructor
    Boundary1D();

    void init() override {
        _init(1);
    }

    string domainType() const override {
        return "boundary";
    }

    bool isConnector() override {
        return true;
    }

    //! Set the temperature.
    virtual void setTemperature(double t) {
        m_temp = t;
    }

    //! Temperature [K].
    virtual double temperature() {
        return m_temp;
    }

    //! Get the number of species
    virtual size_t nSpecies() {
        return 0;
    }

    //! Set the mole fractions by specifying a string.
    virtual void setMoleFractions(const string& xin) {
        throw NotImplementedError("Boundary1D::setMoleFractions");
    }

    //! Set the mole fractions by specifying an array.
    virtual void setMoleFractions(const double* xin) {
        throw NotImplementedError("Boundary1D::setMoleFractions");
    }

    //! Mass fraction of species k.
    virtual double massFraction(size_t k) {
        throw NotImplementedError("Boundary1D::massFraction");
    }

    //! Set the total mass flow rate [kg/m²/s].
    virtual void setMdot(double mdot) {
        m_mdot = mdot;
    }

    //! Set tangential velocity gradient [1/s] at this boundary.
    virtual void setSpreadRate(double V0) {
        throw NotImplementedError("Boundary1D::setSpreadRate");
    }

    //! Tangential velocity gradient [1/s] at this boundary.
    virtual double spreadRate() {
        throw NotImplementedError("Boundary1D::spreadRate");
    }

    //! The total mass flow rate [kg/m2/s].
    virtual double mdot() {
        return m_mdot;
    }

    void setupGrid(size_t n, const double* z) override {}

    void fromArray(SolutionArray& arr, double* soln) override;

protected:
    //! Initialize member variables based on the adjacent domains.
    //! @param n  Number of state variables associated with the boundary object
    void _init(size_t n);

    Flow1D* m_flow_left = nullptr; //!< Flow domain to the left of this boundary
    Flow1D* m_flow_right = nullptr; //! Flow domain to the right of this boundary
    size_t m_left_nv = 0; //!< Number of state vector components in left flow domain
    size_t m_right_nv = 0; //!< Number of state vector components in right flow domain
    size_t m_left_nsp = 0; //!< Number of species in left flow domain
    size_t m_right_nsp = 0; //!< Number of species in right flow domain
    ThermoPhase* m_phase_left = nullptr; //!< Thermo object used by left flow domain
    ThermoPhase* m_phase_right = nullptr; //!< Thermo object used by right flow domain

    //! Temperature of the boundary.
    double m_temp = 0.0;
    //! Mass flow rate at the boundary.
    double m_mdot = 0.0;
};


/**
 * An inlet.
 * Unstrained flows use an inlet to the left of the flow domain (left-to-right flow).
 * Strained flow configurations may have inlets on the either side of the flow domain.
 */
class Inlet1D : public Boundary1D
{
public:
    //! Default constructor
    Inlet1D();

    //! Constructor with contents
    //! @param solution  Solution representing contents of adjacent flow domain
    //! @param id  Name used to identify this domain
    Inlet1D(shared_ptr<Solution> solution, const string& id="");

    string domainType() const override {
        return "inlet";
    }

    void setSpreadRate(double V0) override;

    double spreadRate() override {
        return m_V0;
    }

    void setTemperature(double T) override;

    void show(const double* x) override;

    size_t nSpecies() override {
        return m_nsp;
    }

    void setMoleFractions(const string& xin) override;
    void setMoleFractions(const double* xin) override;
    double massFraction(size_t k) override {
        return m_yin[k];
    }
    void init() override;
    void eval(size_t jg, double* xg, double* rg, integer* diagg, double rdt) override;
    shared_ptr<SolutionArray> asArray(const double* soln) const override;
    void fromArray(SolutionArray& arr, double* soln) override;

protected:
    //! A marker that indicates whether this is a left inlet or a right inlet.
    int m_ilr;

    //! The spread rate of the inlet [1/s]
    double m_V0 = 0.0;

    size_t m_nsp = 0; //!< Number of species in the adjacent flow domain
    vector<double> m_yin; //!< inlet mass fractions
    string m_xstr; //!< inlet mass fractions. Parsing deferred to init()
    Flow1D* m_flow = nullptr; //!< the adjacent flow domain
};

/**
 * A terminator that does nothing.
 */
class Empty1D : public Boundary1D
{
public:
    //! Default constructor
    Empty1D() = default;

    //! Constructor with contents
    //! @param solution  Solution representing contents
    //! @param id  Name used to identify this domain
    Empty1D(shared_ptr<Solution> solution, const string& id="") : Empty1D() {
        setSolution(solution);
        m_id = id;
    }

    string domainType() const override {
        return "empty";
    }

    void show(const double* x) override {}

    void init() override;

    void eval(size_t jg, double* xg, double* rg, integer* diagg, double rdt) override;

    shared_ptr<SolutionArray> asArray(const double* soln) const override;
};

/**
 * A symmetry plane. The axial velocity u = 0, and all other components have
 * zero axial gradients.
 */
class Symm1D : public Boundary1D
{
public:
    //! Default constructor
    Symm1D() = default;

    //! Constructor with contents
    //! @param solution  Solution representing contents of adjacent flow domain
    //! @param id  Name used to identify this domain
    Symm1D(shared_ptr<Solution> solution, const string& id="") : Symm1D() {
        setSolution(solution);
        m_id = id;
    }

    string domainType() const override {
        return "symmetry-plane";
    }

    void init() override;

    void eval(size_t jg, double* xg, double* rg, integer* diagg, double rdt) override;

    shared_ptr<SolutionArray> asArray(const double* soln) const override;
};


/**
 * An outlet.
 * Flow is assumed to be from left to right.
 */
class Outlet1D : public Boundary1D
{
public:
    //! Default constructor
    Outlet1D() = default;

    //! Constructor with contents
    //! @param solution  Solution representing contents of adjacent flow domain
    //! @param id  Name used to identify this domain
    Outlet1D(shared_ptr<Solution> solution, const string& id="") : Outlet1D() {
        setSolution(solution);
        m_id = id;
    }

    string domainType() const override {
        return "outlet";
    }

    void init() override;

    void eval(size_t jg, double* xg, double* rg, integer* diagg, double rdt) override;

    shared_ptr<SolutionArray> asArray(const double* soln) const override;
};


/**
 * An outlet with specified composition.
 * Flow is assumed to be from left to right.
 */
class OutletRes1D : public Boundary1D
{
public:
    //! Default constructor
    OutletRes1D();

    //! Constructor with contents
    //! @param solution  Solution representing contents of adjacent flow domain
    //! @param id  Name used to identify this domain
    OutletRes1D(shared_ptr<Solution> solution, const string& id="");

    string domainType() const override {
        return "outlet-reservoir";
    }

    void show(const double* x) override {}

    size_t nSpecies() override {
        return m_nsp;
    }

    void setMoleFractions(const string& xin) override;
    void setMoleFractions(const double* xin) override;
    double massFraction(size_t k) override {
        return m_yres[k];
    }

    void init() override;
    void eval(size_t jg, double* xg, double* rg, integer* diagg, double rdt) override;
    shared_ptr<SolutionArray> asArray(const double* soln) const override;
    void fromArray(SolutionArray& arr, double* soln) override;

protected:
    size_t m_nsp = 0; //!< Number of species in the adjacent flow domain
    vector<double> m_yres; //!< Mass fractions in the reservoir
    string m_xstr; //!< Mole fractions in the reservoir
    Flow1D* m_flow = nullptr; //!< The adjacent flow domain
};

/**
 * A non-reacting surface. The axial velocity is zero (impermeable), as is the
 * transverse velocity (no slip). The temperature is specified, and a zero flux
 * condition is imposed for the species.
 */
class Surf1D : public Boundary1D
{
public:
    //! Default constructor
    Surf1D() = default;

    //! Constructor with contents
    //! @param solution  Solution representing contents of adjacent flow domain
    //! @param id  Name used to identify this domain
    Surf1D(shared_ptr<Solution> solution, const string& id="") : Surf1D() {
        setSolution(solution);
        m_id = id;
    }

    string domainType() const override {
        return "surface";
    }

    void init() override;

    void eval(size_t jg, double* xg, double* rg, integer* diagg, double rdt) override;

    shared_ptr<SolutionArray> asArray(const double* soln) const override;
    void fromArray(SolutionArray& arr, double* soln) override;

    //! @deprecated To be removed after Cantera 3.1.
    void show(std::ostream& s, const double* x) override;

    void show(const double* x) override;
};

/**
 * A reacting surface.
 */
class ReactingSurf1D : public Boundary1D
{
public:
    //! Default constructor
    ReactingSurf1D();

    //! Constructor with contents
    //! @param solution  Solution representing contents of adjacent flow domain
    //! @param id  Name used to identify this domain
    ReactingSurf1D(shared_ptr<Solution> solution, const string& id="");

    string domainType() const override {
        return "reacting-surface";
    }

    void setKinetics(shared_ptr<Kinetics> kin) override;

    //! Set whether to solve the equations for the surface species coverages
    void enableCoverageEquations(bool docov) {
        m_enabled = docov;
    }

    //! Indicates whether the equations for the surface species coverages are being
    //! solved
    bool coverageEnabled() {
        return m_enabled;
    }

    string componentName(size_t n) const override;

    void init() override;
    void resetBadValues(double* xg) override;

    void eval(size_t jg, double* xg, double* rg, integer* diagg, double rdt) override;

    shared_ptr<SolutionArray> asArray(const double* soln) const override;
    void fromArray(SolutionArray& arr, double* soln) override;

    void _getInitialSoln(double* x) override {
        m_sphase->getCoverages(x);
    }

    void _finalize(const double* x) override {
        std::copy(x, x+m_nsp, m_fixed_cov.begin());
    }

    void show(const double* x) override;

protected:
    InterfaceKinetics* m_kin = nullptr; //!< surface kinetics mechanism
    SurfPhase* m_sphase = nullptr; //!< phase representing the surface species
    size_t m_nsp = 0; //!< the number of surface phase species
    bool m_enabled = false; //!< True if coverage equations are being solved

    //! temporary vector used to store coverages and production rates. Size is total
    //! number of species in the kinetic mechanism
    vector<double> m_work;

    //! Fixed values of the coverages used when coverage equations are not being solved.
    //! Length is #m_nsp.
    vector<double> m_fixed_cov;
};

//! @} End of bdryGroup

}

#endif
