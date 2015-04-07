/**
 * @file Inlet1D.h
 *
 * Boundary objects for one-dimensional simulations.
 */

/*
 * Copyright 2002-3  California Institute of Technology
 */

#ifndef CT_BDRY1D_H
#define CT_BDRY1D_H

#include "Domain1D.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "StFlow.h"
#include "OneDim.h"
#include "cantera/base/ctml.h"

#include <cstdio>

namespace Cantera
{

const int LeftInlet = 1;
const int RightInlet = -1;

/**
 * The base class for boundaries between one-dimensional spatial
 * domains. The boundary may have its own internal variables, such
 * as surface species coverages.
 *
 * The boundary types are an inlet, an outlet, a symmetry plane,
 * and a surface.
 *
 * The public methods are all virtual, and the base class
 * implementations throw exceptions.
 * @ingroup onedim
 */
class Bdry1D : public Domain1D
{
public:
    Bdry1D();

    virtual void init() {
        _init(1);
    }

    /// Set the temperature.
    virtual void setTemperature(doublereal t) {
        m_temp = t;
    }

    /// Temperature [K].
    virtual doublereal temperature() {
        return m_temp;
    }

    virtual size_t nSpecies() {
        return 0;
    }

    /// Set the mole fractions by specifying a std::string.
    virtual void setMoleFractions(const std::string& xin) {
        err("setMoleFractions");
    }

    /// Set the mole fractions by specifying an array.
    virtual void setMoleFractions(doublereal* xin) {
        err("setMoleFractions");
    }

    /// Mass fraction of species k.
    virtual doublereal massFraction(size_t k) {
        err("massFraction");
        return 0.0;
    }

    /// Set the total mass flow rate.
    virtual void setMdot(doublereal mdot) {
        m_mdot = mdot;
    }

    /// The total mass flow rate [kg/m2/s].
    virtual doublereal mdot() {
        return m_mdot;
    }

    virtual void _getInitialSoln(doublereal* x) {
        writelog("Bdry1D::_getInitialSoln called!\n");
    }

    virtual void setupGrid(size_t n, const doublereal* z) {}

protected:
    void _init(size_t n);

    StFlow* m_flow_left, *m_flow_right;
    size_t m_ilr, m_left_nv, m_right_nv;
    size_t m_left_loc, m_right_loc;
    size_t m_left_points;
    size_t m_nv, m_left_nsp, m_right_nsp;
    size_t m_sp_left, m_sp_right;
    size_t m_start_left, m_start_right;
    ThermoPhase* m_phase_left, *m_phase_right;
    doublereal m_temp, m_mdot;

private:
    void err(const std::string& method) {
        throw CanteraError("Bdry1D::"+method,
                           "attempt to call base class method "+method);
    }
};


/**
 * An inlet.
 * @ingroup onedim
 */
class Inlet1D : public Bdry1D
{
public:
    Inlet1D() : Bdry1D(), m_V0(0.0), m_nsp(0), m_flow(0) {
        m_type = cInletType;
        m_xstr = "";
    }

    /// set spreading rate
    virtual void setSpreadRate(doublereal V0) {
        m_V0 = V0;
        needJacUpdate();
    }

    /// spreading rate
    virtual double spreadRate() {
        return m_V0;
    }

    virtual void showSolution(const doublereal* x) {
        char buf[80];
        sprintf(buf, "    Mass Flux:   %10.4g kg/m^2/s \n", m_mdot);
        writelog(buf);
        sprintf(buf, "    Temperature: %10.4g K \n", m_temp);
        writelog(buf);
        if (m_flow) {
            writelog("    Mass Fractions: \n");
            for (size_t k = 0; k < m_flow->phase().nSpecies(); k++) {
                if (m_yin[k] != 0.0) {
                    sprintf(buf, "        %16s  %10.4g \n",
                            m_flow->phase().speciesName(k).c_str(), m_yin[k]);
                    writelog(buf);
                }
            }
        }
        writelog("\n");
    }

    virtual void _getInitialSoln(doublereal* x) {
        x[0] = m_mdot;
        x[1] = m_temp;
    }

    virtual void _finalize(const doublereal* x) {}

    virtual size_t nSpecies() {
        return m_nsp;
    }

    virtual void setMoleFractions(const std::string& xin);
    virtual void setMoleFractions(doublereal* xin);
    virtual doublereal massFraction(size_t k) {
        return m_yin[k];
    }
    virtual std::string componentName(size_t n) const;
    virtual void init();
    virtual void eval(size_t jg, doublereal* xg, doublereal* rg,
                      integer* diagg, doublereal rdt);
    virtual XML_Node& save(XML_Node& o, const doublereal* const soln);
    virtual void restore(const XML_Node& dom, doublereal* soln, int loglevel);

protected:
    int m_ilr;
    doublereal m_V0;
    size_t m_nsp;
    vector_fp m_yin;
    std::string m_xstr;
    StFlow* m_flow;
};

/**
 * A terminator that does nothing.
 * @ingroup onedim
 */
class Empty1D : public Domain1D
{
public:
    Empty1D() : Domain1D() {
        m_type = cEmptyType;
    }

    virtual std::string componentName(size_t n) const;
    virtual void showSolution(const doublereal* x) {}

    virtual void init();

    virtual void eval(size_t jg, doublereal* xg, doublereal* rg,
                      integer* diagg, doublereal rdt);

    virtual XML_Node& save(XML_Node& o, const doublereal* const soln);
    virtual void restore(const XML_Node& dom, doublereal* soln, int loglevel);
    virtual void _finalize(const doublereal* x) {}
    virtual void _getInitialSoln(doublereal* x) {
        x[0] = 0.0;
    }
};

/**
 * A symmetry plane. The axial velocity u = 0, and all other
 * components have zero axial gradients.
 * @ingroup onedim
 */
class Symm1D : public Bdry1D
{
public:

    Symm1D() : Bdry1D() {
        m_type = cSymmType;
    }

    virtual std::string componentName(size_t n) const;

    virtual void init();

    virtual void eval(size_t jg, doublereal* xg, doublereal* rg,
                      integer* diagg, doublereal rdt);

    virtual XML_Node& save(XML_Node& o, const doublereal* const soln);
    virtual void restore(const XML_Node& dom, doublereal* soln, int loglevel);
    virtual void _finalize(const doublereal* x) {
        ; //m_temp = x[0];
    }
    virtual void _getInitialSoln(doublereal* x) {
        x[0] = m_temp;
    }
};


/**
 *  An outlet.
 */
class Outlet1D : public Bdry1D
{
public:
    Outlet1D() : Bdry1D() {
        m_type = cOutletType;
    }

    virtual std::string componentName(size_t n) const;

    virtual void init();

    virtual void eval(size_t jg, doublereal* xg, doublereal* rg,
                      integer* diagg, doublereal rdt);

    virtual XML_Node& save(XML_Node& o, const doublereal* const soln);
    virtual void restore(const XML_Node& dom, doublereal* soln, int loglevel);
    virtual void _finalize(const doublereal* x) {
        ; //m_temp = x[0];
    }
    virtual void _getInitialSoln(doublereal* x) {
        x[0] = m_temp;
    }
};


/**
 * An outlet with specified composition.
 * @ingroup onedim
 */
class OutletRes1D : public Bdry1D
{
public:
    OutletRes1D() : Bdry1D(), m_nsp(0), m_flow(0) {
        m_type = cOutletResType;
        m_xstr = "";
    }

    virtual void showSolution(const doublereal* x) {}

    virtual void _getInitialSoln(doublereal* x) {
        x[0] = m_temp;
    }

    virtual void _finalize(const doublereal* x) {
        ;
    }

    virtual size_t nSpecies() {
        return m_nsp;
    }

    virtual void setMoleFractions(const std::string& xin);
    virtual void setMoleFractions(doublereal* xin);
    virtual doublereal massFraction(size_t k) {
        return m_yres[k];
    }
    virtual std::string componentName(size_t n) const;
    virtual void init();
    virtual void eval(size_t jg, doublereal* xg, doublereal* rg,
                      integer* diagg, doublereal rdt);
    virtual XML_Node& save(XML_Node& o, const doublereal* const soln);
    virtual void restore(const XML_Node& dom, doublereal* soln, int loglevel);

protected:
    size_t m_nsp;
    vector_fp m_yres;
    std::string m_xstr;
    StFlow* m_flow;
};

/**
 * A non-reacting surface. The axial velocity is zero
 * (impermeable), as is the transverse velocity (no slip). The
 * temperature is specified, and a zero flux condition is imposed
 * for the species.
 */
class Surf1D : public Bdry1D
{
public:
    Surf1D() : Bdry1D() {
        m_type = cSurfType;
    }

    virtual std::string componentName(size_t n) const;

    virtual void init();

    virtual void eval(size_t jg, doublereal* xg, doublereal* rg,
                      integer* diagg, doublereal rdt);

    virtual XML_Node& save(XML_Node& o, const doublereal* const soln);
    virtual void restore(const XML_Node& dom, doublereal* soln, int loglevel);

    virtual void _getInitialSoln(doublereal* x) {
        x[0] = m_temp;
    }

    virtual void _finalize(const doublereal* x) {
        ; //m_temp = x[0];
    }

    virtual void showSolution_s(std::ostream& s, const doublereal* x) {
        s << "-------------------  Surface " << domainIndex() << " ------------------- " << std::endl;
        s << "  temperature: " << m_temp << " K" << "    " << x[0] << std::endl;
    }

    virtual void showSolution(const doublereal* x) {
        char buf[80];
        sprintf(buf, "    Temperature: %10.4g K \n", m_temp);
        writelog(buf);
        writelog("\n");
    }
};

/**
 * A reacting surface.
 * @ingroup onedim
 */
class ReactingSurf1D : public Bdry1D
{
public:
    ReactingSurf1D() : Bdry1D(),
        m_kin(0), m_surfindex(0), m_nsp(0) {
        m_type = cSurfType;
    }

    void setKineticsMgr(InterfaceKinetics* kin) {
        m_kin = kin;
        m_surfindex = kin->surfacePhaseIndex();
        m_sphase = (SurfPhase*)&kin->thermo(m_surfindex);
        m_nsp = m_sphase->nSpecies();
        m_enabled = true;
    }

    void enableCoverageEquations(bool docov) {
        m_enabled = docov;
    }

    virtual std::string componentName(size_t n) const;

    virtual void init();

    virtual void eval(size_t jg, doublereal* xg, doublereal* rg,
                      integer* diagg, doublereal rdt);

    virtual XML_Node& save(XML_Node& o, const doublereal* const soln);
    virtual void restore(const XML_Node& dom, doublereal* soln, int loglevel);

    virtual void _getInitialSoln(doublereal* x) {
        x[0] = m_temp;
        //m_kin->advanceCoverages(1.0);
        m_sphase->getCoverages(x+1);
    }

    virtual void _finalize(const doublereal* x) {
        std::copy(x+1,x+1+m_nsp,m_fixed_cov.begin());
    }

    virtual void showSolution(const doublereal* x) {
        char buf[80];
        sprintf(buf, "    Temperature: %10.4g K \n", x[0]);
        writelog(buf);
        writelog("    Coverages: \n");
        for (size_t k = 0; k < m_nsp; k++) {
            sprintf(buf, "    %20s %10.4g \n", m_sphase->speciesName(k).c_str(),
                    x[k+1]);
            writelog(buf);
        }
        writelog("\n");
    }

protected:
    InterfaceKinetics* m_kin;
    SurfPhase* m_sphase;
    size_t m_surfindex, m_nsp;
    bool m_enabled;
    vector_fp m_work;
    vector_fp m_fixed_cov;
    int dum;
};

}

#endif
