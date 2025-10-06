//! @file Boundary1D.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/SolutionArray.h"
#include "cantera/oneD/Boundary1D.h"
#include "cantera/oneD/OneDim.h"
#include "cantera/oneD/Flow1D.h"

using namespace std;

namespace Cantera
{

Boundary1D::Boundary1D() : Domain1D(1, 1, 0.0)
{
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
        if (!r.isConnector()) { // multi-point domain
            m_left_nv = r.nComponents();
            if (m_left_nv > c_offset_Y) {
                m_left_nsp = m_left_nv - c_offset_Y;
            } else {
                m_left_nsp = 0;
            }
            m_flow_left = dynamic_cast<Flow1D*>(&r);
            if (m_flow_left != nullptr) {
                m_phase_left = &m_flow_left->phase();
            }
        } else {
            throw CanteraError("Boundary1D::_init",
                "Boundary domains can only be connected on the left to flow "
                "domains, not '{}' domains.", r.type());
        }
    }

    // if this is not the last domain, see what is connected on the right
    if (m_index + 1 < container().nDomains()) {
        Domain1D& r = container().domain(m_index+1);
        if (!r.isConnector()) { // multi-point domain
            m_right_nv = r.nComponents();
            if (m_right_nv > c_offset_Y) {
                m_right_nsp = m_right_nv - c_offset_Y;
            } else {
                m_right_nsp = 0;
            }
            m_flow_right = dynamic_cast<Flow1D*>(&r);
            if (m_flow_right != nullptr) {
                m_phase_right = &m_flow_right->phase();
            }
        } else {
            throw CanteraError("Boundary1D::_init",
                "Boundary domains can only be connected on the right to flow "
                "domains, not '{}' domains.", r.type());
        }
    }
}

void Boundary1D::fromArray(const shared_ptr<SolutionArray>& arr)
{
    setMeta(arr->meta());
}

// ---------------- Inlet1D methods ----------------

Inlet1D::Inlet1D()
{
}

Inlet1D::Inlet1D(shared_ptr<Solution> solution, const string& id)
    : Inlet1D()
{
    setSolution(solution);
    m_id = id;
}


//! set spreading rate
void Inlet1D::setSpreadRate(double V0)
{
    m_V0 = V0;
    needJacUpdate();
}

void Inlet1D::setTemperature(double T)
{
    Boundary1D::setTemperature(T);
    // Adjust flow domain temperature bounds based on inlet temperature
    if (m_flow != nullptr && m_flow->lowerBound(c_offset_T) >= m_temp) {
        m_flow->setBounds(c_offset_T, m_temp - 5.0, m_flow->upperBound(c_offset_T));
    }
}

void Inlet1D::show(const double* x)
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

void Inlet1D::setMoleFractions(const string& xin)
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
    if (m_flow_left && !m_flow_right) {
        if (!m_flow_left->isStrained()) {
            throw CanteraError("Inlet1D::init",
                "Right inlets with right-to-left flow are only supported for "
                "strained flow configurations.");
        }
        m_ilr = RightInlet;
        m_flow = m_flow_left;
    } else if (m_flow_right) {
        m_ilr = LeftInlet;
        m_flow = m_flow_right;
    } else {
        throw CanteraError("Inlet1D::init", "Inlet1D is not properly connected.");
    }

    // components = u, V, T, Lambda, + mass fractions
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

        if (m_flow->doEnergy(0)) {
            // The third flow residual is for T, where it is set to T(0).  Subtract
            // the local temperature to hold the flow T to the inlet T.
            rb[c_offset_T] -= m_temp;
        } else {
            rb[c_offset_T] -= m_flow->T_fixed(0);
        }

        if (m_flow->isFree()) {
            // if the flow is a freely-propagating flame, mdot is not specified.
            // Set mdot equal to rho*u.
            m_mdot = m_flow->density(0) * xb[c_offset_U];
        } else if (m_flow->isStrained()) { // axisymmetric flow
            if (m_flow->twoPointControlEnabled()) {
                // When using two-point control, the mass flow rate at the left inlet is
                // not specified. Instead, the mass flow rate is dictated by the
                // velocity at the left inlet, which comes from the U variable. The
                // default boundary condition specified in the Flow1D.cpp file already
                // handles this case. We only need to update the stored value of m_mdot
                // so that other equations that use the quantity are consistent.
                m_mdot = m_flow->density(0)*xb[c_offset_U];
            } else {
                // The flow domain sets this to -rho*u. Add mdot to specify the mass
                // flow rate
                rb[c_offset_L] += m_mdot;
            }

            // spreading rate. The flow domain sets this to V(0),
            // so for finite spreading rate subtract m_V0.
            rb[c_offset_V] -= m_V0;
        } else { // unstrained flow
            rb[c_offset_U] = m_flow->density(0) * xb[c_offset_U] - m_mdot;
        }

        // add the convective term to the species residual equations
        for (size_t k = 0; k < m_nsp; k++) {
            if (k != m_flow_right->leftExcessSpecies()) {
                rb[c_offset_Y+k] += m_mdot*m_yin[k];
            }
        }

    } else {
        // right inlet (should only be used for counter-flow flames)
        // Array elements corresponding to the last point in the flow domain
        double* rb = rg + loc() - m_flow->nComponents();
        double* xb = xg + loc() - m_flow->nComponents();
        size_t last_index = m_flow->nPoints() - 1;

        rb[c_offset_V] -= m_V0;
        if (m_flow->doEnergy(m_flow->nPoints() - 1)) {
            rb[c_offset_T] -= m_temp; // T
        } else {
            rb[c_offset_T] -= m_flow->T_fixed(m_flow->nPoints() - 1);
        }

        if (m_flow->twoPointControlEnabled()) { // For point control adjustments
            // At the right boundary, the mdot is dictated by the velocity at the right
            // boundary, which comes from the Uo variable. The variable Uo is the
            // left-moving velocity and has a negative value, so the mass flow has to be
            // negated to give a positive value when using Uo.
            m_mdot = -m_flow->density(last_index) * xb[c_offset_Uo];
        }
        rb[c_offset_U] += m_mdot;

        for (size_t k = 0; k < m_nsp; k++) {
            if (k != m_flow_left->rightExcessSpecies()) {
                rb[c_offset_Y+k] += m_mdot * m_yin[k];
            }
        }
    }
}

shared_ptr<SolutionArray> Inlet1D::toArray(bool normalize) const
{
    AnyMap meta = Boundary1D::getMeta();
    meta["mass-flux"] = m_mdot;
    auto arr = SolutionArray::create(m_solution, 1, meta);

    // set gas state (using pressure from adjacent domain)
    double pressure = m_flow->phase().pressure();
    auto phase = m_solution->thermo();
    phase->setState_TPY(m_temp, pressure, m_yin.data());
    vector<double> data(phase->stateSize());
    phase->saveState(data);

    arr->setState(0, data);
    if (normalize) {
        arr->normalize();
    }
    return arr;
}

void Inlet1D::fromArray(const shared_ptr<SolutionArray>& arr)
{
    Boundary1D::setMeta(arr->meta());
    arr->setLoc(0);
    auto phase = arr->thermo();
    auto meta = arr->meta();
    m_temp = phase->temperature();
    if (meta.hasKey("mass-flux")) {
        m_mdot = meta.at("mass-flux").asDouble();
    } else {
        // convert data format used by Python h5py export (Cantera < 3.0)
        auto aux = arr->getAuxiliary(0);
        m_mdot = phase->density() * aux.at("velocity").as<double>();
    }
    phase->getMassFractions(m_yin.data());
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

shared_ptr<SolutionArray> Empty1D::toArray(bool normalize) const
{
    AnyMap meta = Boundary1D::getMeta();
    return SolutionArray::create(m_solution, 0, meta);
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

shared_ptr<SolutionArray> Symm1D::toArray(bool normalize) const
{
    AnyMap meta = Boundary1D::getMeta();
    return SolutionArray::create(m_solution, 0, meta);
}

// -------- Outlet1D --------

OutletRes1D::OutletRes1D()
{
}

OutletRes1D::OutletRes1D(shared_ptr<Solution> solution, const string& id)
    : OutletRes1D()
{
    setSolution(solution);
    m_id = id;
}

void Outlet1D::init()
{
    _init(0);

    if (m_flow_right) {
        throw CanteraError("Outlet1D::init",
            "Left outlets with right-to-left flow are not supported.");
    }
    if (m_flow_left) {
        m_flow_left->setViscosityFlag(false);
    } else {
        throw CanteraError("Outlet1D::init", "Outlet1D is not connected.");
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

    // flow is left-to-right
    size_t nc = m_flow_left->nComponents();
    double* xb = x - nc;
    double* rb = r - nc;
    int* db = diag - nc;

    size_t last = m_flow_left->nPoints() - 1;
    if (m_flow_left->doEnergy(last)) {
        rb[c_offset_T] = xb[c_offset_T] - xb[c_offset_T - nc]; // zero T gradient
    } else {
        rb[c_offset_T] = xb[c_offset_T] - m_flow_left->T_fixed(last);
    }
    size_t kSkip = c_offset_Y + m_flow_left->rightExcessSpecies();
    for (size_t k = c_offset_Y; k < nc; k++) {
        if (k != kSkip) {
            rb[k] = xb[k] - xb[k - nc]; // zero mass fraction gradient
            db[k] = 0;
        }
    }
}

shared_ptr<SolutionArray> Outlet1D::toArray(bool normalize) const
{
    AnyMap meta = Boundary1D::getMeta();
    return SolutionArray::create(m_solution, 0, meta);
}

// -------- OutletRes1D --------

void OutletRes1D::setMoleFractions(const string& xres)
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

    if (m_flow_right) {
        throw CanteraError("OutletRes1D::init",
            "Left outlets with right-to-left flow are not supported.");
    }
    if (m_flow_left) {
        m_flow = m_flow_left;
    } else {
        throw CanteraError("OutletRes1D::init", "no flow!");
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

    size_t nc = m_flow_left->nComponents();
    double* xb = x - nc;
    double* rb = r - nc;
    int* db = diag - nc;

    size_t last = m_flow_left->nPoints() - 1;
    if (m_flow_left->doEnergy(last)) {
        rb[c_offset_T] = xb[c_offset_T] - xb[c_offset_T - nc]; // zero T gradient
    } else {
        rb[c_offset_T] = xb[c_offset_T] - m_flow_left->T_fixed(last);
    }
    size_t kSkip = m_flow_left->rightExcessSpecies();
    for (size_t k = c_offset_Y; k < nc; k++) {
        if (k != kSkip) {
            rb[k] = xb[k] - m_yres[k-c_offset_Y]; // fixed Y
            db[k] = 0;
        }
    }
}

shared_ptr<SolutionArray> OutletRes1D::toArray(bool normalize) const
{
    AnyMap meta = Boundary1D::getMeta();
    meta["temperature"] = m_temp;
    auto arr = SolutionArray::create(m_solution, 1, meta);

    // set gas state (using pressure from adjacent domain)
    double pressure = m_flow->phase().pressure();
    auto phase = m_solution->thermo();
    phase->setState_TPY(m_temp, pressure, &m_yres[0]);
    vector<double> data(phase->stateSize());
    phase->saveState(data);

    arr->setState(0, data);
    if (normalize) {
        arr->normalize();
    }
    return arr;
}

void OutletRes1D::fromArray(const shared_ptr<SolutionArray>& arr)
{
    Boundary1D::setMeta(arr->meta());
    arr->setLoc(0);
    auto phase = arr->thermo();
    m_temp = phase->temperature();
    auto Y = phase->massFractions();
    std::copy(Y, Y + m_nsp, &m_yres[0]);
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

shared_ptr<SolutionArray> Surf1D::toArray(bool normalize) const
{
    AnyMap meta = Boundary1D::getMeta();
    meta["temperature"] = m_temp;
    return SolutionArray::create(m_solution, 0, meta);
}

void Surf1D::fromArray(const shared_ptr<SolutionArray>& arr)
{
    auto meta = arr->meta();
    m_temp = meta["temperature"].asDouble();
    meta.erase("temperature");
    Boundary1D::setMeta(meta);
}

void Surf1D::show(const double* x)
{
    writelog("    Temperature: {:10.4g} K \n\n", m_temp);
}

// -------- ReactingSurf1D --------

ReactingSurf1D::ReactingSurf1D()
    : m_kin(0)
    , m_nsp(0)
{
}

ReactingSurf1D::ReactingSurf1D(shared_ptr<Solution> solution, const string& id)
{
    auto phase = std::dynamic_pointer_cast<SurfPhase>(solution->thermo());
    if (!phase) {
        throw CanteraError("ReactingSurf1D::ReactingSurf1D",
            "Detected incompatible ThermoPhase type '{}'", solution->thermo()->type());
    }
    auto kin = std::dynamic_pointer_cast<InterfaceKinetics>(solution->kinetics());
    if (!kin) {
        throw CanteraError("ReactingSurf1D::ReactingSurf1D",
            "Detected incompatible kinetics type '{}'",
            solution->kinetics()->kineticsType());
    }
    setSolution(solution);
    m_id = id;
    m_kin = kin.get();
    m_sphase = phase.get();
    m_nsp = m_sphase->nSpecies();
    m_enabled = true;
}

void ReactingSurf1D::setKinetics(shared_ptr<Kinetics> kin)
{
    auto sol = Solution::create();
    sol->setThermo(kin->reactionPhase());
    sol->setKinetics(kin);
    sol->setTransportModel("none");
    setSolution(sol);
    m_kin = dynamic_pointer_cast<InterfaceKinetics>(kin).get();
    m_sphase = dynamic_pointer_cast<SurfPhase>(kin->reactionPhase()).get();
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

    if (m_enabled) {
        for (size_t k = 0; k < m_nsp; k++) {
            r[k] = m_work[k] * m_sphase->size(k) * rs0;
            r[k] -= rdt*(x[k] - prevSoln(k,0));
            diag[k] = 1;
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
        const vector<double>& mwleft = m_phase_left->molecularWeights();
        double* rb = r - nc;
        double* xb = x - nc;
        rb[c_offset_T] = xb[c_offset_T] - m_temp; // specified T
        size_t nSkip = m_flow_left->rightExcessSpecies();
        size_t l_offset = 0;
        ThermoPhase* left_thermo = &m_flow_left->phase();
        for (size_t nth = 0; nth < m_kin->nPhases(); nth++) {
            if (&m_kin->thermo(nth) == left_thermo) {
                l_offset = m_kin->kineticsSpeciesIndex(0, nth);
                break;
            }
        }
        for (size_t nl = 0; nl < m_left_nsp; nl++) {
            if (nl != nSkip) {
                rb[c_offset_Y+nl] += m_work[nl + l_offset]*mwleft[nl];
            }
        }
    }
}

shared_ptr<SolutionArray> ReactingSurf1D::toArray(bool normalize) const
{
    if (!m_state) {
        throw CanteraError("ReactingSurf1D::toArray",
            "Domain needs to be installed in a container before calling toArray.");
    }
    double* soln = m_state->data() + m_iloc;
    AnyMap meta = Boundary1D::getMeta();
    meta["temperature"] = m_temp;
    meta["phase"]["name"] = m_sphase->name();
    AnyValue source = m_sphase->input().getMetadata("filename");
    meta["phase"]["source"] = source.empty() ? "<unknown>" : source.asString();

    // set state of surface phase
    m_sphase->setState_TP(m_temp, m_sphase->pressure());
    m_sphase->setCoverages(soln);
    vector<double> data(m_sphase->stateSize());
    m_sphase->saveState(data.size(), &data[0]);

    auto arr = SolutionArray::create(m_solution, 1, meta);
    arr->setState(0, data);
    if (normalize) {
        arr->normalize();
    }
    return arr;
}

void ReactingSurf1D::fromArray(const shared_ptr<SolutionArray>& arr)
{
    if (!m_state) {
        throw CanteraError("Domain1D::fromArray",
            "Domain needs to be installed in a container before calling fromArray.");
    }
    resize(nComponents(), arr->size());
    m_container->resize();
    double* soln = m_state->data() + m_iloc;

    Boundary1D::setMeta(arr->meta());
    arr->setLoc(0);
    auto surf = std::dynamic_pointer_cast<SurfPhase>(arr->thermo());
    if (!surf) {
        throw CanteraError("ReactingSurf1D::fromArray",
            "Restoring of coverages requires surface phase");
    }
    m_temp = surf->temperature();
    surf->getCoverages(soln);
    _finalize(soln);
}

void ReactingSurf1D::show(const double* x)
{
    writelog("    Temperature: {:10.4g} K \n", m_temp);
    writelog("    Coverages: \n");
    for (size_t k = 0; k < m_nsp; k++) {
        writelog("    {:>20s} {:10.4g} \n", m_sphase->speciesName(k), x[k]);
    }
    writelog("\n");
}
}
