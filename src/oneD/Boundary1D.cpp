//! @file Boundary1D.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/Boundary1D.h"
#include "cantera/oneD/OneDim.h"
#include "cantera/base/ctml.h"
#include "cantera/oneD/StFlow.h"

using namespace std;

namespace Cantera
{

Boundary1D::Boundary1D() : Domain1D(1, 1, 0.0),
    m_flow_left(0), m_flow_right(0),
    m_ilr(0), m_left_nv(0), m_right_nv(0),
    m_left_loc(0), m_right_loc(0),
    m_left_points(0),
    m_left_nsp(0), m_right_nsp(0),
    m_sp_left(0), m_sp_right(0),
    m_start_left(0), m_start_right(0),
    m_phase_left(0), m_phase_right(0), m_temp(0.0), m_mdot(0.0)
{
    m_type = cConnectorType;
}

void Boundary1D::_init(size_t n)
{
    if (m_index == npos) {
        throw CanteraError("Boundary1D::_init",
                           "install in container before calling init.");
    }

    // A boundary object contains only one grid point
    resize(n,1);

    m_left_nsp = 0;
    m_right_nsp = 0;

    // check for left and right flow objects
    if (m_index > 0) {
        Domain1D& r = container().domain(m_index-1);
        if (!r.isConnector()) { // flow domain
            m_flow_left = (StFlow*)&r;
            m_left_nv = m_flow_left->nComponents();
            m_left_points = m_flow_left->nPoints();
            m_left_loc = container().start(m_index-1);
            m_left_nsp = m_left_nv - c_offset_Y;
            m_phase_left = &m_flow_left->phase();
        } else {
            throw CanteraError("Boundary1D::_init",
                "Boundary domains can only be connected on the left to flow "
                "domains, not type {} domains.", r.domainType());
        }
    }

    // if this is not the last domain, see what is connected on the right
    if (m_index + 1 < container().nDomains()) {
        Domain1D& r = container().domain(m_index+1);
        if (!r.isConnector()) { // flow domain
            m_flow_right = (StFlow*)&r;
            m_right_nv = m_flow_right->nComponents();
            m_right_loc = container().start(m_index+1);
            m_right_nsp = m_right_nv - c_offset_Y;
            m_phase_right = &m_flow_right->phase();
        } else {
            throw CanteraError("Boundary1D::_init",
                "Boundary domains can only be connected on the right to flow "
                "domains, not type {} domains.", r.domainType());
        }
    }
}

// ---------------- Inlet1D methods ----------------

Inlet1D::Inlet1D()
    : m_V0(0.0)
    , m_nsp(0)
    , m_flow(0)
{
    m_type = cInletType;
    m_xstr = "";
}

void Inlet1D::showSolution(const double* x)
{
    writelog("    Mass Flux:   {:10.4g} kg/m^2/s \n", m_mdot);
    writelog("    Temperature: {:10.4g} K \n", m_temp);
    if (m_flow) {
        writelog("    Mass Fractions: \n");
        for (size_t k = 0; k < m_flow->phase().nSpecies(); k++) {
            if (m_yin[k] != 0.0) {
                writelog("        {:>16s}  {:10.4g} \n",
                        m_flow->phase().speciesName(k), m_yin[k]);
            }
        }
    }
    writelog("\n");
}

void Inlet1D::setMoleFractions(const std::string& xin)
{
    m_xstr = xin;
    if (m_flow) {
        m_flow->phase().setMoleFractionsByName(xin);
        m_flow->phase().getMassFractions(m_yin.data());
        needJacUpdate();
    }
}

void Inlet1D::setMoleFractions(const double* xin)
{
    if (m_flow) {
        m_flow->phase().setMoleFractions(xin);
        m_flow->phase().getMassFractions(m_yin.data());
        needJacUpdate();
    }
}

void Inlet1D::init()
{
    _init(0);

    // if a flow domain is present on the left, then this must be a right inlet.
    // Note that an inlet object can only be a terminal object - it cannot have
    // flows on both the left and right
    if (m_flow_left) {
        m_ilr = RightInlet;
        m_flow = m_flow_left;
    } else if (m_flow_right) {
        m_ilr = LeftInlet;
        m_flow = m_flow_right;
    } else {
        throw CanteraError("Inlet1D::init","no flow!");
    }

    // components = u, V, T, lambda, + mass fractions
    m_nsp = m_flow->phase().nSpecies();
    m_yin.resize(m_nsp, 0.0);
    if (m_xstr != "") {
        setMoleFractions(m_xstr);
    } else {
        m_yin[0] = 1.0;
    }
}

void Inlet1D::eval(size_t jg, double* xg, double* rg,
                   integer* diagg, double rdt)
{
    if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    if (m_ilr == LeftInlet) {
        // Array elements corresponding to the first point of the flow domain
        double* xb = xg + m_flow->loc();
        double* rb = rg + m_flow->loc();

        // The first flow residual is for u. This, however, is not modified by
        // the inlet, since this is set within the flow domain from the
        // continuity equation.

        // spreading rate. The flow domain sets this to V(0),
        // so for finite spreading rate subtract m_V0.
        rb[c_offset_V] -= m_V0;

        if (m_flow->doEnergy(0)) {
            // The third flow residual is for T, where it is set to T(0).  Subtract
            // the local temperature to hold the flow T to the inlet T.
            rb[c_offset_T] -= m_temp;
        }

        if (m_flow->fixed_mdot()) {
            // The flow domain sets this to -rho*u. Add mdot to specify the mass
            // flow rate.
            rb[c_offset_L] += m_mdot;
        } else {
            // if the flow is a freely-propagating flame, mdot is not specified.
            // Set mdot equal to rho*u, and also set lambda to zero.
            m_mdot = m_flow->density(0)*xb[0];
            rb[c_offset_L] = xb[c_offset_L];
        }

        // add the convective term to the species residual equations
        for (size_t k = 0; k < m_nsp; k++) {
            if (k != m_flow_right->leftExcessSpecies()) {
                rb[c_offset_Y+k] += m_mdot*m_yin[k];
            }
        }

    } else {
        // right inlet
        // Array elements corresponding to the flast point in the flow domain
        double* rb = rg + loc() - m_flow->nComponents();
        rb[c_offset_V] -= m_V0;
        if (m_flow->doEnergy(m_flow->nPoints() - 1)) {
            rb[c_offset_T] -= m_temp; // T
        }
        rb[c_offset_U] += m_mdot; // u
        for (size_t k = 0; k < m_nsp; k++) {
            if (k != m_flow_left->rightExcessSpecies()) {
                rb[c_offset_Y+k] += m_mdot * m_yin[k];
            }
        }
    }
}

XML_Node& Inlet1D::save(XML_Node& o, const double* const soln)
{
    XML_Node& inlt = Domain1D::save(o, soln);
    inlt.addAttribute("type","inlet");
    addFloat(inlt, "temperature", m_temp);
    addFloat(inlt, "mdot", m_mdot);
    for (size_t k=0; k < m_nsp; k++) {
        addFloat(inlt, "massFraction", m_yin[k], "",
                       m_flow->phase().speciesName(k));
    }
    return inlt;
}

void Inlet1D::restore(const XML_Node& dom, double* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    m_mdot = getFloat(dom, "mdot");
    m_temp = getFloat(dom, "temperature");

    m_yin.assign(m_nsp, 0.0);

    for (size_t i = 0; i < dom.nChildren(); i++) {
        const XML_Node& node = dom.child(i);
        if (node.name() == "massFraction") {
            size_t k = m_flow->phase().speciesIndex(node.attrib("type"));
            if (k != npos) {
                m_yin[k] = node.fp_value();
            }
        }
    }
    resize(0, 1);
}

// ------------- Empty1D -------------

void Empty1D::init()
{
    _init(0);
}

void Empty1D::eval(size_t jg, double* xg, double* rg,
     integer* diagg, double rdt)
{
}

XML_Node& Empty1D::save(XML_Node& o, const double* const soln)
{
    XML_Node& symm = Domain1D::save(o, soln);
    symm.addAttribute("type","empty");
    return symm;
}

void Empty1D::restore(const XML_Node& dom, double* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    resize(0, 1);
}

// -------------- Symm1D --------------

void Symm1D::init()
{
    _init(0);
}

void Symm1D::eval(size_t jg, double* xg, double* rg, integer* diagg,
                  double rdt)
{
    if (jg != npos && (jg + 2< firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    double* x = xg + loc();
    double* r = rg + loc();
    integer* diag = diagg + loc();

    if (m_flow_right) {
        size_t nc = m_flow_right->nComponents();
        double* xb = x;
        double* rb = r;
        int* db = diag;
        db[c_offset_V] = 0;
        db[c_offset_T] = 0;
        rb[c_offset_V] = xb[c_offset_V] - xb[c_offset_V + nc]; // zero dV/dz
        if (m_flow_right->doEnergy(0)) {
            rb[c_offset_T] = xb[c_offset_T] - xb[c_offset_T + nc]; // zero dT/dz
        }
    }

    if (m_flow_left) {
        size_t nc = m_flow_left->nComponents();
        double* xb = x - nc;
        double* rb = r - nc;
        int* db = diag - nc;
        db[c_offset_V] = 0;
        db[c_offset_T] = 0;
        rb[c_offset_V] = xb[c_offset_V] - xb[c_offset_V - nc]; // zero dV/dz
        if (m_flow_left->doEnergy(m_flow_left->nPoints() - 1)) {
            rb[c_offset_T] = xb[c_offset_T] - xb[c_offset_T - nc]; // zero dT/dz
        }
    }
}

XML_Node& Symm1D::save(XML_Node& o, const double* const soln)
{
    XML_Node& symm = Domain1D::save(o, soln);
    symm.addAttribute("type","symmetry");
    return symm;
}

void Symm1D::restore(const XML_Node& dom, double* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    resize(0, 1);
}

// -------- Outlet1D --------

OutletRes1D::OutletRes1D()
    : m_nsp(0)
    , m_flow(0)
{
    m_type = cOutletResType;
    m_xstr = "";
}

void Outlet1D::init()
{
    _init(0);

    if (m_flow_right) {
        m_flow_right->setViscosityFlag(false);
    }
    if (m_flow_left) {
        m_flow_left->setViscosityFlag(false);
    }
}

void Outlet1D::eval(size_t jg, double* xg, double* rg, integer* diagg,
                    double rdt)
{
    if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    double* x = xg + loc();
    double* r = rg + loc();
    integer* diag = diagg + loc();

    if (m_flow_right) {
        size_t nc = m_flow_right->nComponents();
        double* xb = x;
        double* rb = r;
        rb[c_offset_U] = xb[c_offset_L];
        if (m_flow_right->doEnergy(0)) {
            rb[c_offset_T] = xb[c_offset_T] - xb[c_offset_T + nc];
        }
        for (size_t k = c_offset_Y; k < nc; k++) {
            rb[k] = xb[k] - xb[k + nc];
        }
    }

    if (m_flow_left) {
        size_t nc = m_flow_left->nComponents();
        double* xb = x - nc;
        double* rb = r - nc;
        int* db = diag - nc;

        // zero Lambda
        if (m_flow_left->fixed_mdot()) {
            rb[c_offset_U] = xb[c_offset_L];
        }

        if (m_flow_left->doEnergy(m_flow_left->nPoints()-1)) {
            rb[c_offset_T] = xb[c_offset_T] - xb[c_offset_T - nc]; // zero T gradient
        }
        size_t kSkip = c_offset_Y + m_flow_left->rightExcessSpecies();
        for (size_t k = c_offset_Y; k < nc; k++) {
            if (k != kSkip) {
                rb[k] = xb[k] - xb[k - nc]; // zero mass fraction gradient
                db[k] = 0;
            }
        }
    }
}

XML_Node& Outlet1D::save(XML_Node& o, const double* const soln)
{
    XML_Node& outlt = Domain1D::save(o, soln);
    outlt.addAttribute("type","outlet");
    return outlt;
}

void Outlet1D::restore(const XML_Node& dom, double* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    resize(0, 1);
}

// -------- OutletRes1D --------

void OutletRes1D::setMoleFractions(const std::string& xres)
{
    m_xstr = xres;
    if (m_flow) {
        m_flow->phase().setMoleFractionsByName(xres);
        m_flow->phase().getMassFractions(m_yres.data());
        needJacUpdate();
    }
}

void OutletRes1D::setMoleFractions(const double* xres)
{
    if (m_flow) {
        m_flow->phase().setMoleFractions(xres);
        m_flow->phase().getMassFractions(m_yres.data());
        needJacUpdate();
    }
}

void OutletRes1D::init()
{
    _init(0);

    if (m_flow_left) {
        m_flow = m_flow_left;
    } else if (m_flow_right) {
        m_flow = m_flow_right;
    } else {
        throw CanteraError("OutletRes1D::init","no flow!");
    }

    m_nsp = m_flow->phase().nSpecies();
    m_yres.resize(m_nsp, 0.0);
    if (m_xstr != "") {
        setMoleFractions(m_xstr);
    } else {
        m_yres[0] = 1.0;
    }
}

void OutletRes1D::eval(size_t jg, double* xg, double* rg,
                       integer* diagg, double rdt)
{
    if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    double* x = xg + loc();
    double* r = rg + loc();
    integer* diag = diagg + loc();

    if (m_flow_right) {
        size_t nc = m_flow_right->nComponents();
        double* xb = x;
        double* rb = r;

        // this seems wrong...
        // zero Lambda
        rb[c_offset_U] = xb[c_offset_L];

        if (m_flow_right->doEnergy(0)) {
            // zero gradient for T
            rb[c_offset_T] = xb[c_offset_T] - xb[c_offset_T + nc];
        }

        // specified mass fractions
        for (size_t k = c_offset_Y; k < nc; k++) {
            rb[k] = xb[k] - m_yres[k-c_offset_Y];
        }
    }

    if (m_flow_left) {
        size_t nc = m_flow_left->nComponents();
        double* xb = x - nc;
        double* rb = r - nc;
        int* db = diag - nc;

        if (!m_flow_left->fixed_mdot()) {
            ;
        } else {
            rb[c_offset_U] = xb[c_offset_L]; // zero Lambda
        }
        if (m_flow_left->doEnergy(m_flow_left->nPoints()-1)) {
            rb[c_offset_T] = xb[c_offset_T] - m_temp; // zero dT/dz
        }
        size_t kSkip = m_flow_left->rightExcessSpecies();
        for (size_t k = c_offset_Y; k < nc; k++) {
            if (k != kSkip) {
                rb[k] = xb[k] - m_yres[k-c_offset_Y]; // fixed Y
                db[k] = 0;
            }
        }
    }
}

XML_Node& OutletRes1D::save(XML_Node& o, const double* const soln)
{
    XML_Node& outlt = Domain1D::save(o, soln);
    outlt.addAttribute("type","outletres");
    addFloat(outlt, "temperature", m_temp, "K");
    for (size_t k=0; k < m_nsp; k++) {
        addFloat(outlt, "massFraction", m_yres[k], "",
                       m_flow->phase().speciesName(k));
    }
    return outlt;
}

void OutletRes1D::restore(const XML_Node& dom, double* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    m_temp = getFloat(dom, "temperature");

    m_yres.assign(m_nsp, 0.0);
    for (size_t i = 0; i < dom.nChildren(); i++) {
        const XML_Node& node = dom.child(i);
        if (node.name() == "massFraction") {
            size_t k = m_flow->phase().speciesIndex(node.attrib("type"));
            if (k != npos) {
                m_yres[k] = node.fp_value();
            }
        }
    }

    resize(0, 1);
}

// -------- Surf1D --------

void Surf1D::init()
{
    _init(0);
}

void Surf1D::eval(size_t jg, double* xg, double* rg,
                  integer* diagg, double rdt)
{
    if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    double* x = xg + loc();
    double* r = rg + loc();

    if (m_flow_right) {
        double* rb = r;
        double* xb = x;
        rb[c_offset_T] = xb[c_offset_T] - m_temp; // specified T
    }

    if (m_flow_left) {
        size_t nc = m_flow_left->nComponents();
        double* rb = r - nc;
        double* xb = x - nc;
        rb[c_offset_T] = xb[c_offset_T] - m_temp; // specified T
    }
}

XML_Node& Surf1D::save(XML_Node& o, const double* const soln)
{
    XML_Node& inlt = Domain1D::save(o, soln);
    inlt.addAttribute("type","surface");
    addFloat(inlt, "temperature", m_temp);
    return inlt;
}

void Surf1D::restore(const XML_Node& dom, double* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    m_temp = getFloat(dom, "temperature");
    resize(0, 1);
}

void Surf1D::showSolution_s(std::ostream& s, const double* x)
{
    s << "-------------------  Surface " << domainIndex() << " ------------------- " << std::endl;
    s << "  temperature: " << m_temp << " K" << std::endl;
}

// -------- ReactingSurf1D --------

ReactingSurf1D::ReactingSurf1D()
    : m_kin(0)
    , m_surfindex(0)
    , m_nsp(0)
{
    m_type = cSurfType;
}

void ReactingSurf1D::setKineticsMgr(InterfaceKinetics* kin)
{
    m_kin = kin;
    m_surfindex = kin->surfacePhaseIndex();
    m_sphase = (SurfPhase*)&kin->thermo(m_surfindex);
    m_nsp = m_sphase->nSpecies();
    m_enabled = true;
}

string ReactingSurf1D::componentName(size_t n) const
{
    if (n < m_nsp) {
        return m_sphase->speciesName(n);
    } else {
        return "<unknown>";
    }
}

void ReactingSurf1D::init()
{
    m_nv = m_nsp;
    _init(m_nsp);
    m_fixed_cov.resize(m_nsp, 0.0);
    m_fixed_cov[0] = 1.0;
    m_work.resize(m_kin->nTotalSpecies(), 0.0);

    for (size_t n = 0; n < m_nsp; n++) {
        setBounds(n, -1.0e-5, 2.0);
    }
}

void ReactingSurf1D::resetBadValues(double* xg) {
    double* x = xg + loc();
    m_sphase->setCoverages(x);
    m_sphase->getCoverages(x);
}

void ReactingSurf1D::eval(size_t jg, double* xg, double* rg,
                          integer* diagg, double rdt)
{
    if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    double* x = xg + loc();
    double* r = rg + loc();
    integer* diag = diagg + loc();

    // set the coverages
    double sum = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        m_work[k] = x[k];
        sum += x[k];
    }
    m_sphase->setTemperature(m_temp);
    m_sphase->setCoveragesNoNorm(m_work.data());

    // set the left gas state to the adjacent point

    size_t leftloc = 0, rightloc = 0;
    size_t pnt = 0;

    if (m_flow_left) {
        leftloc = m_flow_left->loc();
        pnt = m_flow_left->nPoints() - 1;
        m_flow_left->setGas(xg + leftloc, pnt);
    }

    if (m_flow_right) {
        rightloc = m_flow_right->loc();
        m_flow_right->setGas(xg + rightloc, 0);
    }

    m_kin->getNetProductionRates(m_work.data());
    double rs0 = 1.0/m_sphase->siteDensity();
    size_t ioffset = m_kin->kineticsSpeciesIndex(0, m_surfindex);

    if (m_enabled) {
        double maxx = -1.0;
        for (size_t k = 0; k < m_nsp; k++) {
            r[k] = m_work[k + ioffset] * m_sphase->size(k) * rs0;
            r[k] -= rdt*(x[k] - prevSoln(k,0));
            diag[k] = 1;
            maxx = std::max(x[k], maxx);
        }
        r[0] = 1.0 - sum;
        diag[0] = 0;
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            r[k] = x[k] - m_fixed_cov[k];
            diag[k] = 0;
        }
    }

    if (m_flow_right) {
        double* rb = r + m_nsp;
        double* xb = x + m_nsp;
        rb[c_offset_T] = xb[c_offset_T] - m_temp; // specified T
    }
    if (m_flow_left) {
        size_t nc = m_flow_left->nComponents();
        const vector_fp& mwleft = m_phase_left->molecularWeights();
        double* rb = r - nc;
        double* xb = x - nc;
        rb[c_offset_T] = xb[c_offset_T] - m_temp; // specified T
        size_t nSkip = m_flow_left->rightExcessSpecies();
        for (size_t nl = 0; nl < m_left_nsp; nl++) {
            if (nl != nSkip) {
                rb[c_offset_Y+nl] += m_work[nl]*mwleft[nl];
            }
        }
    }
}

XML_Node& ReactingSurf1D::save(XML_Node& o, const double* const soln)
{
    const double* s = soln + loc();
    XML_Node& dom = Domain1D::save(o, soln);
    dom.addAttribute("type","surface");
    addFloat(dom, "temperature", m_temp, "K");
    for (size_t k=0; k < m_nsp; k++) {
        addFloat(dom, "coverage", s[k], "", m_sphase->speciesName(k));
    }
    return dom;
}

void ReactingSurf1D::restore(const XML_Node& dom, double* soln,
                             int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    m_temp = getFloat(dom, "temperature");

    m_fixed_cov.assign(m_nsp, 0.0);
    for (size_t i = 0; i < dom.nChildren(); i++) {
        const XML_Node& node = dom.child(i);
        if (node.name() == "coverage") {
            size_t k = m_sphase->speciesIndex(node.attrib("type"));
            if (k != npos) {
                m_fixed_cov[k] = soln[k] = node.fp_value();
            }
        }
    }
    m_sphase->setCoverages(&m_fixed_cov[0]);

    resize(m_nsp, 1);
}

void ReactingSurf1D::showSolution(const double* x)
{
    writelog("    Temperature: {:10.4g} K \n", m_temp);
    writelog("    Coverages: \n");
    for (size_t k = 0; k < m_nsp; k++) {
        writelog("    {:>20s} {:10.4g} \n", m_sphase->speciesName(k), x[k]);
    }
    writelog("\n");
}
}
