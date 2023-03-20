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
#include "StFlow.h"

namespace Cantera
{

const int LeftInlet = 1;
const int RightInlet = -1;

/**
 * The base class for boundaries between one-dimensional spatial domains. The
 * boundary may have its own internal variables, such as surface species
 * coverages.
 *
 * The boundary types are an inlet, an outlet, a symmetry plane, and a surface.
 *
 * The public methods are all virtual, and the base class implementations throw
 * exceptions.
 * @ingroup onedim
 */
class Boundary1D : public Domain1D
{
public:
    Boundary1D();

    virtual void init() {
        _init(1);
    }

    virtual string type() const {
        return "boundary";
    }

    virtual bool isConnector() {
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

    virtual size_t nSpecies() {
        return 0;
    }

    //! Set the mole fractions by specifying a std::string.
    virtual void setMoleFractions(const std::string& xin) {
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

    //! Set the total mass flow rate.
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

    virtual void setupGrid(size_t n, const double* z) {}

protected:
    void _init(size_t n);

    StFlow* m_flow_left = nullptr;
    StFlow* m_flow_right = nullptr;
    size_t m_ilr = 0;
    size_t m_left_nv = 0;
    size_t m_right_nv = 0;
    size_t m_left_loc = 0;
    size_t m_right_loc = 0;
    size_t m_left_points = 0;
    size_t m_left_nsp = 0;
    size_t m_right_nsp = 0;
    size_t m_sp_left = 0;
    size_t m_sp_right = 0;
    size_t m_start_left = 0;
    size_t m_start_right = 0;
    ThermoPhase* m_phase_left = nullptr;
    ThermoPhase* m_phase_right = nullptr;
    double m_temp = 0.0;
    double m_mdot = 0.0;
};


/**
 * An inlet.
 * @ingroup onedim
 */
class Inlet1D : public Boundary1D
{
public:
    Inlet1D();

    Inlet1D(shared_ptr<Solution> solution, const string& id="");

    virtual string type() const {
        return "inlet";
    }

    virtual void setSpreadRate(double V0);

    virtual double spreadRate() {
        return m_V0;
    }

    virtual void show(const double* x);

    virtual size_t nSpecies() {
        return m_nsp;
    }

    virtual void setMoleFractions(const std::string& xin);
    virtual void setMoleFractions(const double* xin);
    virtual double massFraction(size_t k) {
        return m_yin[k];
    }
    virtual void init();
    virtual void eval(size_t jg, double* xg, double* rg,
                      integer* diagg, double rdt);
    virtual shared_ptr<SolutionArray> asArray(const double* soln) const;
    virtual void restore(SolutionArray& arr, double* soln);

protected:
    int m_ilr;
    double m_V0 = 0.0;
    size_t m_nsp = 0;
    vector_fp m_yin;
    std::string m_xstr;
    StFlow* m_flow = nullptr;
};

/**
 * A terminator that does nothing.
 * @ingroup onedim
 */
class Empty1D : public Boundary1D
{
public:
    Empty1D() {
        m_type = cEmptyType;
    }

    Empty1D(shared_ptr<Solution> solution, const std::string& id="") : Empty1D() {
        m_solution = solution;
        m_id = id;
    }

    virtual string type() const {
        return "empty";
    }

    virtual void show(const double* x) {}

    virtual void init();

    virtual void eval(size_t jg, double* xg, double* rg,
                      integer* diagg, double rdt);

    virtual shared_ptr<SolutionArray> asArray(const double* soln) const;
    virtual void restore(SolutionArray& arr, double* soln) {}
};

/**
 * A symmetry plane. The axial velocity u = 0, and all other components have
 * zero axial gradients.
 * @ingroup onedim
 */
class Symm1D : public Boundary1D
{
public:
    Symm1D() {
        m_type = cSymmType;
    }

    Symm1D(shared_ptr<Solution> solution, const std::string& id="") : Symm1D() {
        m_solution = solution;
        m_id = id;
    }

    virtual string type() const {
        return "symmetry-plane";
    }

    virtual void init();

    virtual void eval(size_t jg, double* xg, double* rg,
                      integer* diagg, double rdt);

    virtual shared_ptr<SolutionArray> asArray(const double* soln) const;
    virtual void restore(SolutionArray& arr, double* soln) {}
};


/**
 * An outlet.
 * @ingroup onedim
 */
class Outlet1D : public Boundary1D
{
public:
    Outlet1D() {
        m_type = cOutletType;
    }

    Outlet1D(shared_ptr<Solution> solution, const std::string& id="") : Outlet1D() {
        m_solution = solution;
        m_id = id;
    }

    virtual string type() const {
        return "outlet";
    }

    virtual void init();

    virtual void eval(size_t jg, double* xg, double* rg,
                      integer* diagg, double rdt);

    virtual shared_ptr<SolutionArray> asArray(const double* soln) const;
    virtual void restore(SolutionArray& arr, double* soln) {}
};


/**
 * An outlet with specified composition.
 * @ingroup onedim
 */
class OutletRes1D : public Boundary1D
{
public:
    OutletRes1D();

    OutletRes1D(shared_ptr<Solution> solution, const string& id="");

    virtual string type() const {
        return "outlet-reservoir";
    }

    virtual void show(const double* x) {}

    virtual size_t nSpecies() {
        return m_nsp;
    }

    virtual void setMoleFractions(const std::string& xin);
    virtual void setMoleFractions(const double* xin);
    virtual double massFraction(size_t k) {
        return m_yres[k];
    }
    virtual void init();
    virtual void eval(size_t jg, double* xg, double* rg,
                      integer* diagg, double rdt);
    virtual shared_ptr<SolutionArray> asArray(const double* soln) const;
    virtual void restore(SolutionArray& arr, double* soln);

protected:
    size_t m_nsp = 0;
    vector_fp m_yres;
    std::string m_xstr;
    StFlow* m_flow = nullptr;
};

/**
 * A non-reacting surface. The axial velocity is zero (impermeable), as is the
 * transverse velocity (no slip). The temperature is specified, and a zero flux
 * condition is imposed for the species.
 * @ingroup onedim
 */
class Surf1D : public Boundary1D
{
public:
    Surf1D() {
        m_type = cSurfType;
    }

    Surf1D(shared_ptr<Solution> solution, const std::string& id="") : Surf1D() {
        m_solution = solution;
        m_id = id;
    }

    virtual string type() const {
        return "surface";
    }

    virtual void init();

    virtual void eval(size_t jg, double* xg, double* rg,
                      integer* diagg, double rdt);

    virtual shared_ptr<SolutionArray> asArray(const double* soln) const;
    virtual void restore(SolutionArray& arr, double* soln);

    virtual void show(std::ostream& s, const double* x);

    virtual void show(const double* x);
};

/**
 * A reacting surface.
 * @ingroup onedim
 */
class ReactingSurf1D : public Boundary1D
{
public:
    ReactingSurf1D();
    ReactingSurf1D(shared_ptr<Solution> solution, const std::string& id="");

    virtual string type() const {
        return "reacting-surface";
    }

    virtual void setKinetics(shared_ptr<Kinetics> kin);

    //! @deprecated  To be removed after Cantera 3.0; replaced by setKinetics
    void setKineticsMgr(InterfaceKinetics* kin);

    void enableCoverageEquations(bool docov) {
        m_enabled = docov;
    }

    bool coverageEnabled() {
        return m_enabled;
    }

    virtual std::string componentName(size_t n) const;

    virtual void init();
    virtual void resetBadValues(double* xg);

    virtual void eval(size_t jg, double* xg, double* rg,
                      integer* diagg, double rdt);

    virtual shared_ptr<SolutionArray> asArray(const double* soln) const;
    virtual void restore(SolutionArray& arr, double* soln);

    virtual void _getInitialSoln(double* x) {
        m_sphase->getCoverages(x);
    }

    virtual void _finalize(const double* x) {
        std::copy(x, x+m_nsp, m_fixed_cov.begin());
    }

    virtual void show(const double* x);

protected:
    InterfaceKinetics* m_kin = nullptr;
    SurfPhase* m_sphase = nullptr;
    size_t m_surfindex = 0;
    size_t m_nsp = 0;
    bool m_enabled = false;
    vector_fp m_work;
    vector_fp m_fixed_cov;
};

}

#endif
