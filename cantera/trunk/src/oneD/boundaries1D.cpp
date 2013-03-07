/**
 * @file boundaries1D.cpp
 */
// Copyright 2002-3  California Institute of Technology

#include "cantera/oneD/Inlet1D.h"

using namespace std;

namespace Cantera
{

Bdry1D::Bdry1D() : Domain1D(1, 1, 0.0),
    m_flow_left(0), m_flow_right(0),
    m_ilr(0), m_left_nv(0), m_right_nv(0),
    m_left_loc(0), m_right_loc(0),
    m_left_points(0), m_nv(0),
    m_left_nsp(0), m_right_nsp(0),
    m_sp_left(0), m_sp_right(0),
    m_start_left(0), m_start_right(0),
    m_phase_left(0), m_phase_right(0), m_temp(0.0), m_mdot(0.0)
{
    m_type = cConnectorType;
}


void Bdry1D::
_init(size_t n)
{
    if (m_index == npos) {
        throw CanteraError("Bdry1D",
                           "install in container before calling init.");
    }

    // A boundary object contains only one grid point
    resize(n,1);

    m_left_nsp = 0;
    m_right_nsp = 0;

    // check for left and right flow objects
    if (m_index > 0) {
        Domain1D& r = container().domain(m_index-1);
        if (r.domainType() == cFlowType) {
            m_flow_left = (StFlow*)&r;
            m_left_nv = m_flow_left->nComponents();
            m_left_points = m_flow_left->nPoints();
            m_left_loc = container().start(m_index-1);
            m_left_nsp = m_left_nv - 4;
            m_phase_left = &m_flow_left->phase();
        } else
            throw CanteraError("Bdry1D::init",
                               "Boundary domains can only be "
                               "connected on the left to flow domains, not type "+int2str(r.domainType())
                               + " domains.");
    }

    // if this is not the last domain, see what is connected on
    // the right
    if (m_index + 1 < container().nDomains()) {
        Domain1D& r = container().domain(m_index+1);
        if (r.domainType() == cFlowType) {
            m_flow_right = (StFlow*)&r;
            m_right_nv = m_flow_right->nComponents();
            m_right_loc = container().start(m_index+1);
            m_right_nsp = m_right_nv - 4;
            m_phase_right = &m_flow_right->phase();
        } else
            throw CanteraError("Bdry1D::init",
                               "Boundary domains can only be "
                               "connected on the right to flow domains, not type "+int2str(r.domainType())
                               + " domains.");
    }
}




//----------------------------------------------------------
//
//   Inlet1D methods
//
//----------------------------------------------------------


void Inlet1D::
setMoleFractions(const std::string& xin)
{
    m_xstr = xin;
    if (m_flow) {
        m_flow->phase().setMoleFractionsByName(xin);
        m_flow->phase().getMassFractions(DATA_PTR(m_yin));
        needJacUpdate();
    }
}

void Inlet1D::
setMoleFractions(doublereal* xin)
{
    if (m_flow) {
        m_flow->phase().setMoleFractions(xin);
        m_flow->phase().getMassFractions(DATA_PTR(m_yin));
        needJacUpdate();
    }
}

string Inlet1D::
componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "mdot";
    case 1:
        return "temperature";
    default:
        break;
    }
    return "unknown";
}

void Inlet1D::
init()
{

    _init(2);

    // set bounds (mdot, T)
    const doublereal lower[2] = {-1.0e5, 200.0};
    const doublereal upper[2] = {1.0e5, 1.e5};
    setBounds(2, lower, 2, upper);

    // set tolerances
    vector_fp rtol(2, 1e-4);
    vector_fp atol(2, 1.e-5);
    setTolerances(2, DATA_PTR(rtol), 2, DATA_PTR(atol));

    // if a flow domain is present on the left, then this must be
    // a right inlet. Note that an inlet object can only be a
    // terminal object - it cannot have flows on both the left and
    // right
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
    m_nsp = m_flow->nComponents() - 4;
    m_yin.resize(m_nsp, 0.0);
    if (m_xstr != "") {
        setMoleFractions(m_xstr);
    } else {
        m_yin[0] = 1.0;
    }
}


void Inlet1D::
eval(size_t jg, doublereal* xg, doublereal* rg,
     integer* diagg, doublereal rdt)
{
    if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* r = rg + loc();
    integer* diag = diagg + loc();
    doublereal* xb, *rb;

    // residual equations for the two local variables

    r[0] = m_mdot - x[0];

    // Temperature
    r[1] = m_temp - x[1];

    // both are algebraic constraints
    diag[0] = 0;
    diag[1] = 0;


    // if it is a left inlet, then the flow solution vector
    // starts 2 to the right in the global solution vector
    if (m_ilr == LeftInlet) {
        xb = x + 2;
        rb = r + 2;

        // The first flow residual is for u. This, however, is not
        // modified by the inlet, since this is set within the flow
        // domain from the continuity equation.

        // spreading rate. The flow domain sets this to V(0),
        // so for finite spreading rate subtract m_V0.
        rb[1] -= m_V0;

        // The third flow residual is for T, where it is set to
        // T(0).  Subtract the local temperature to hold the flow
        // T to the inlet T.
        rb[2] -= x[1];

        // The flow domain sets this to -rho*u. Add mdot to
        // specify the mass flow rate.
        rb[3] += x[0];

        // add the convective term to the species residual equations
        for (size_t k = 1; k < m_nsp; k++) {
            rb[4+k] += x[0]*m_yin[k];
        }

        // if the flow is a freely-propagating flame, mdot is not
        // specified.  Set mdot equal to rho*u, and also set
        // lambda to zero.
        if (!m_flow->fixed_mdot()) {
            r[0] = m_flow->density(0)*xb[0] - x[0];
            rb[3] = xb[3];
        }
    }

    // right inlet.
    else {
        size_t boffset = m_flow->nComponents();
        xb = x - boffset;
        rb = r - boffset;
        rb[1] -= m_V0;
        rb[2] -= x[1]; // T
        rb[0] += x[0]; // u
        for (size_t k = 1; k < m_nsp; k++) {
            rb[4+k] += x[0]*(m_yin[k]);
        }
    }


}

XML_Node& Inlet1D::
save(XML_Node& o, const doublereal* const soln)
{
    const doublereal* s = soln + loc();
    XML_Node& inlt = Domain1D::save(o, soln);
    inlt.addAttribute("type","inlet");
    for (size_t k = 0; k < nComponents(); k++) {
        ctml::addFloat(inlt, componentName(k), s[k], "", "",lowerBound(k), upperBound(k));
    }
    for (size_t k=0; k < m_nsp; k++) {
        ctml::addFloat(inlt, "massFraction", m_yin[k], "",
                       m_flow->phase().speciesName(k));
    }
    return inlt;
}

void Inlet1D::
restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    soln[0] = m_mdot = ctml::getFloat(dom, "mdot", "massflowrate");
    soln[1] = m_temp = ctml::getFloat(dom, "temperature", "temperature");

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
    resize(2,1);
}



//--------------------------------------------------
//      Empty1D
//--------------------------------------------------

string Empty1D::componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "dummy";
    default:
        break;
    }
    return "<unknown>";
}

void Empty1D::
init()   //_init(1);
{
    // set bounds (T)
    const doublereal lower = -1.0;
    const doublereal upper = 1.0;
    setBounds(1, &lower, 1, &upper);

    // set tolerances
    const doublereal rtol = 1e-4;
    const doublereal atol = 1.e-4;
    setTolerances(1, &rtol, 1, &atol);
}

void Empty1D::
eval(size_t jg, doublereal* xg, doublereal* rg,
     integer* diagg, doublereal rdt)
{
    if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* r = rg + loc();
    integer* diag = diagg + loc();
    //        integer *db;

    r[0] = x[0];
    diag[0] = 0;
}

XML_Node& Empty1D::
save(XML_Node& o, const doublereal* const soln)
{
    XML_Node& symm = Domain1D::save(o, soln);
    symm.addAttribute("type","empty");
    return symm;
}

void Empty1D::
restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    resize(1,1);
}



//--------------------------------------------------
//      Symm1D
//--------------------------------------------------

string Symm1D::componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "dummy";
    default:
        break;
    }
    return "<unknown>";
}

void Symm1D::
init()
{
    _init(1);
    // set bounds (T)
    const doublereal lower = -1.0;
    const doublereal upper = 1.0;
    setBounds(1, &lower, 1, &upper);

    // set tolerances
    const doublereal rtol = 1e-4;
    const doublereal atol = 1.e-4;
    setTolerances(1, &rtol, 1, &atol);
}

void Symm1D::
eval(size_t jg, doublereal* xg, doublereal* rg,
     integer* diagg, doublereal rdt)
{
    if (jg != npos && (jg + 2< firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* r = rg + loc();
    integer* diag = diagg + loc();
    doublereal* xb, *rb;
    integer* db;

    r[0] = x[0];
    diag[0] = 0;
    size_t nc;

    if (m_flow_right) {
        nc = m_flow_right->nComponents();
        xb = x + 1;
        rb = r + 1;
        db = diag + 1;
        db[1] = 0;
        db[2] = 0;
        rb[1] = xb[1] - xb[1 + nc];      // zero dV/dz
        rb[2] = xb[2] - xb[2 + nc];      // zero dT/dz
    }

    if (m_flow_left) {
        nc = m_flow_left->nComponents();
        xb = x - nc;
        rb = r - nc;
        db = diag - nc;
        db[1] = 0;
        db[2] = 0;
        rb[1] = xb[1] - xb[1 - nc];      // zero dV/dz
        rb[2] = xb[2] - xb[2 - nc];      // zero dT/dz
    }
}


XML_Node& Symm1D::
save(XML_Node& o, const doublereal* const soln)
{
    XML_Node& symm = Domain1D::save(o, soln);
    symm.addAttribute("type","symmetry");
    return symm;
}

void Symm1D::
restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    resize(1,1);
}


//--------------------------------------------------
//      Outlet1D
//--------------------------------------------------

string Outlet1D::componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "outlet dummy";
    default:
        break;
    }
    return "<unknown>";
}

void Outlet1D::
init()
{
    _init(1);
    // set bounds (T)
    const doublereal lower = -1.0;
    const doublereal upper = 1.0;
    setBounds(1, &lower, 1, &upper);

    // set tolerances
    const doublereal rtol = 1e-4;
    const doublereal atol = 1.e-4;
    setTolerances(1, &rtol, 1, &atol);
    if (m_flow_right) {
        m_flow_right->setViscosityFlag(false);
    }
    if (m_flow_left) {
        m_flow_left->setViscosityFlag(false);
    }
}


void Outlet1D::
eval(size_t jg, doublereal* xg, doublereal* rg,
     integer* diagg, doublereal rdt)
{
    if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* r = rg + loc();
    integer* diag = diagg + loc();
    doublereal* xb, *rb;
    integer* db;

    r[0] = x[0];
    diag[0] = 0;
    size_t nc, k;

    if (m_flow_right) {
        nc = m_flow_right->nComponents();
        xb = x + 1;
        rb = r + 1;
        db = diag + 1;
        rb[0] = xb[3];
        rb[2] = xb[2] - xb[2 + nc];
        for (k = 4; k < nc; k++) {
            //if (m_flow_right->doSpecies(k-4)) {
            rb[k] = xb[k] - xb[k + nc];
            //}
        }
    }

    if (m_flow_left) {
        nc = m_flow_left->nComponents();
        xb = x - nc;
        rb = r - nc;
        db = diag - nc;

        // zero Lambda

        if (!m_flow_left->fixed_mdot()) {
            ;    //                rb[0] = xb[0] - xb[0-nc]; //zero U gradient
        } else {
            rb[0] = xb[3];    // zero Lambda
        }

        rb[2] = xb[2] - xb[2 - nc];  // zero T gradient
        for (k = 5; k < nc; k++) {
            rb[k] = xb[k] - xb[k - nc]; // zero mass fraction gradient
            db[k] = 0;
        }
    }
}


XML_Node& Outlet1D::
save(XML_Node& o, const doublereal* const soln)
{
    XML_Node& outlt = Domain1D::save(o, soln);
    outlt.addAttribute("type","outlet");
    return outlt;
}

void Outlet1D::
restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    resize(1,1);
}




//--------------------------------------------------
//      OutletRes1D
//--------------------------------------------------


void OutletRes1D::
setMoleFractions(const std::string& xres)
{
    m_xstr = xres;
    if (m_flow) {
        m_flow->phase().setMoleFractionsByName(xres);
        m_flow->phase().getMassFractions(DATA_PTR(m_yres));
        needJacUpdate();
    }
}

void OutletRes1D::
setMoleFractions(doublereal* xres)
{
    if (m_flow) {
        m_flow->phase().setMoleFractions(xres);
        m_flow->phase().getMassFractions(DATA_PTR(m_yres));
        needJacUpdate();
    }
}

string OutletRes1D::componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "dummy";
    default:
        break;
    }
    return "<unknown>";
}

void OutletRes1D::
init()
{
    _init(1);
    // set bounds (dummy)
    const doublereal lower = -1.0;
    const doublereal upper = 1.0;
    setBounds(1, &lower, 1, &upper);

    // set tolerances
    const doublereal rtol = 1e-4;
    const doublereal atol = 1.e-4;
    setTolerances(1, &rtol, 1, &atol);

    if (m_flow_left) {
        m_flow = m_flow_left;
    } else if (m_flow_right) {
        m_flow = m_flow_right;
    } else {
        throw CanteraError("OutletRes1D::init","no flow!");
    }

    m_nsp = m_flow->nComponents() - 4;
    m_yres.resize(m_nsp, 0.0);
    if (m_xstr != "") {
        setMoleFractions(m_xstr);
    } else {
        m_yres[0] = 1.0;
    }
}


void OutletRes1D::
eval(size_t jg, doublereal* xg, doublereal* rg,
     integer* diagg, doublereal rdt)
{

    if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* r = rg + loc();
    integer* diag = diagg + loc();
    doublereal* xb, *rb;
    integer* db;

    // drive dummy component to zero
    r[0] = x[0];
    diag[0] = 0;
    size_t nc, k;

    if (m_flow_right) {
        nc = m_flow_right->nComponents();
        xb = x + 1;
        rb = r + 1;
        db = diag + 1;

        // this seems wrong...
        // zero Lambda
        rb[0] = xb[3];

        // zero gradient for T
        rb[2] = xb[2] - xb[2 + nc];

        // specified mass fractions
        for (k = 4; k < nc; k++) {
            rb[k] = xb[k] - m_yres[k-4];
        }
    }

    if (m_flow_left) {

        nc = m_flow_left->nComponents();
        xb = x - nc;
        rb = r - nc;
        db = diag - nc;

        if (!m_flow_left->fixed_mdot()) {
            ;
        } else {
            rb[0] = xb[3];    // zero Lambda
        }
        rb[2] = xb[2] - m_temp; //xb[2] - xb[2 - nc];        // zero dT/dz
        for (k = 5; k < nc; k++) {
            rb[k] = xb[k] - m_yres[k-4];     // fixed Y
            db[k] = 0;
        }
    }
}


XML_Node& OutletRes1D::
save(XML_Node& o, const doublereal* const soln)
{
    XML_Node& outlt = Domain1D::save(o, soln);
    outlt.addAttribute("type","outletres");
    ctml::addFloat(outlt, "temperature", m_temp, "K");
    for (size_t k=0; k < m_nsp; k++) {
        ctml::addFloat(outlt, "massFraction", m_yres[k], "",
                       m_flow->phase().speciesName(k));
    }
    return outlt;
}

void OutletRes1D::
restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    m_temp = ctml::getFloat(dom, "temperature");

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

    resize(1,1);
}


//-----------------------------------------------------------
//
//  Surf1D
//
//-----------------------------------------------------------



string Surf1D::componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "temperature";
    default:
        break;
    }
    return "<unknown>";
}

void Surf1D::
init()
{
    _init(1);
    // set bounds (T)
    const doublereal lower = 200.0;
    const doublereal upper = 1.e5;
    setBounds(1, &lower, 1, &upper);

    // set tolerances
    const doublereal rtol = 1e-4;
    const doublereal atol = 1.e-4;
    setTolerances(1, &rtol, 1, &atol);
}


void Surf1D::
eval(size_t jg, doublereal* xg, doublereal* rg,
     integer* diagg, doublereal rdt)
{
    if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* r = rg + loc();
    integer* diag = diagg + loc();
    doublereal* xb, *rb;

    r[0] = x[0] - m_temp;
    diag[0] = 0;
    size_t nc;

    if (m_flow_right) {
        rb = r + 1;
        xb = x + 1;
        rb[2] = xb[2] - x[0];            // specified T
    }

    if (m_flow_left) {
        nc = m_flow_left->nComponents();
        rb = r - nc;
        xb = x - nc;
        rb[2] = xb[2] - x[0];            // specified T
    }
}

XML_Node& Surf1D::
save(XML_Node& o, const doublereal* const soln)
{
    const doublereal* s = soln + loc();
    //XML_Node& inlt = o.addChild("inlet");
    XML_Node& inlt = Domain1D::save(o, soln);
    inlt.addAttribute("type","surface");
    for (size_t k = 0; k < nComponents(); k++) {
        ctml::addFloat(inlt, componentName(k), s[k], "", "",0.0, 1.0);
    }
    return inlt;
}

void Surf1D::
restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    soln[0] = m_temp = ctml::getFloat(dom, "temperature", "temperature");
    resize(1,1);
}




//-----------------------------------------------------------
//
//  ReactingSurf1D
//
//-----------------------------------------------------------



string ReactingSurf1D::componentName(size_t n) const
{
    if (n == 0) {
        return "temperature";
    } else if (n < m_nsp + 1) {
        return m_sphase->speciesName(n-1);
    } else {
        return "<unknown>";
    }
}

void ReactingSurf1D::
init()
{
    m_nv = m_nsp + 1;
    _init(m_nsp+1);
    m_fixed_cov.resize(m_nsp, 0.0);
    m_fixed_cov[0] = 1.0;
    m_work.resize(m_kin->nTotalSpecies(), 0.0);

    // set bounds
    vector_fp lower(m_nv), upper(m_nv);
    lower[0] = 200.0;
    upper[0] = 1.e5;
    for (size_t n = 0; n < m_nsp; n++) {
        lower[n+1] = -1.0e-5;
        upper[n+1] = 2.0;
    }
    setBounds(m_nv, DATA_PTR(lower), m_nv, DATA_PTR(upper));
    vector_fp rtol(m_nv), atol(m_nv);
    for (size_t n = 0; n < m_nv; n++) {
        rtol[n] = 1.0e-5;
        atol[n] = 1.0e-9;
    }
    atol[0] = 1.0e-4;
    setTolerances(m_nv, DATA_PTR(rtol), m_nv, DATA_PTR(atol));
}


void ReactingSurf1D::
eval(size_t jg, doublereal* xg, doublereal* rg,
     integer* diagg, doublereal rdt)
{
    if (jg != npos && (jg + 2 < firstPoint() || jg > lastPoint() + 2)) {
        return;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* r = rg + loc();
    integer* diag = diagg + loc();
    doublereal* xb, *rb;

    // specified surface temp
    r[0] = x[0] - m_temp;

    // set the coverages
    doublereal sum = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        m_work[k] = x[k+1];
        sum += x[k+1];
    }
    m_sphase->setTemperature(x[0]);
    m_sphase->setCoverages(DATA_PTR(m_work));
    //m_kin->advanceCoverages(1.0);
    //m_sphase->getCoverages(m_fixed_cov.begin());

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

    m_kin->getNetProductionRates(DATA_PTR(m_work));
    doublereal rs0 = 1.0/m_sphase->siteDensity();

    //scale(m_work.begin(), m_work.end(), m_work.begin(), m_mult[0]);

    //        bool enabled = true;
    size_t ioffset = m_kin->kineticsSpeciesIndex(0, m_surfindex);

    if (m_enabled) {
        doublereal maxx = -1.0;
        for (size_t k = 0; k < m_nsp; k++) {
            r[k+1] = m_work[k + ioffset] * m_sphase->size(k) * rs0;
            r[k+1] -= rdt*(x[k+1] - prevSoln(k+1,0));
            diag[k+1] = 1;
            if (x[k+1] > maxx) {
                maxx = x[k+1];
            }
        }
        r[1] = 1.0 - sum;
        diag[1] = 0;
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            r[k+1] = x[k+1] - m_fixed_cov[k];
            diag[k+1] = 0;
        }
    }

    if (m_flow_right) {
        rb = r + 1;
        xb = x + 1;
        rb[2] = xb[2] - x[0];            // specified T
    }
    size_t nc;
    if (m_flow_left) {
        nc = m_flow_left->nComponents();
        const doublereal* mwleft = DATA_PTR(m_phase_left->molecularWeights());
        rb =r - nc;
        xb = x - nc;
        rb[2] = xb[2] - x[0];            // specified T
        for (size_t nl = 1; nl < m_left_nsp; nl++) {
            rb[4+nl] += m_work[nl]*mwleft[nl];
        }
    }
}

XML_Node& ReactingSurf1D::
save(XML_Node& o, const doublereal* const soln)
{
    const doublereal* s = soln + loc();
    XML_Node& dom = Domain1D::save(o, soln);
    dom.addAttribute("type","surface");
    ctml::addFloat(dom, "temperature", s[0], "K");
    for (size_t k=0; k < m_nsp; k++) {
        ctml::addFloat(dom, "coverage", s[k+1], "",
                       m_sphase->speciesName(k));
    }
    return dom;
}

void ReactingSurf1D::
restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    soln[0] = m_temp = ctml::getFloat(dom, "temperature");

    m_fixed_cov.assign(m_nsp, 0.0);
    for (size_t i = 0; i < dom.nChildren(); i++) {
        const XML_Node& node = dom.child(i);
        if (node.name() == "coverage") {
            size_t k = m_sphase->speciesIndex(node.attrib("type"));
            if (k != npos) {
                m_fixed_cov[k] = soln[k+1] = node.fp_value();
            }
        }
    }
    m_sphase->setCoverages(&m_fixed_cov[0]);

    resize(m_nsp+1,1);
}
}

