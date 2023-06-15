//! @file Boundary1D.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/SolutionArray.h"
#include "cantera/oneD/Boundary1D.h"
#include "cantera/oneD/OneDim.h"
#include "cantera/oneD/StFlow.h"

using namespace std;

namespace Cantera
{

Boundary1D::Boundary1D() : Domain1D(1, 1, 0.0)
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
        if (!r.isConnector()) { // multi-point domain
            m_left_nv = r.nComponents();
            if (m_left_nv > c_offset_Y) {
                m_left_nsp = m_left_nv - c_offset_Y;
            } else {
                m_left_nsp = 0;
            }
            m_left_loc = container().start(m_index-1);
            m_left_points = r.nPoints();
            m_flow_left = dynamic_cast<StFlow*>(&r);
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
            m_right_loc = container().start(m_index+1);
            m_flow_right = dynamic_cast<StFlow*>(&r);
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

// ---------------- Inlet1D methods ----------------

Inlet1D::Inlet1D()
{
    m_type = cInletType;
}

Inlet1D::Inlet1D(shared_ptr<Solution> solution, const string& id)
    : Inlet1D()
{
    m_solution = solution;
    m_id = id;
}


//! set spreading rate
void Inlet1D::setSpreadRate(double V0)
{
    m_V0 = V0;
    needJacUpdate();
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
        // Array elements corresponding to the last point in the flow domain
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

shared_ptr<SolutionArray> Inlet1D::asArray(const double* soln) const
{
    AnyMap meta = Boundary1D::getMeta();
    meta["mass-flux"] = m_mdot;
    auto arr = SolutionArray::create(m_solution, 1, meta);

    // set gas state (using pressure from adjacent domain)
    double pressure = m_flow->phase().pressure();
    auto phase = m_solution->thermo();
    phase->setState_TPY(m_temp, pressure, m_yin.data());
    vector_fp data(phase->stateSize());
    phase->saveState(data);

    arr->setState(0, data);
    return arr;
}

void Inlet1D::fromArray(SolutionArray& arr, double* soln)
{
    Boundary1D::setMeta(arr.meta());
    arr.setLoc(0);
    auto phase = arr.thermo();
    auto meta = arr.meta();
    m_temp = phase->temperature();
    if (meta.hasKey("mass-flux")) {
        m_mdot = meta.at("mass-flux").asDouble();
    } else {
        // convert data format used by Python h5py export (Cantera < 3.0)
        auto aux = arr.getAuxiliary(0);
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

shared_ptr<SolutionArray> Empty1D::asArray(const double* soln) const
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

shared_ptr<SolutionArray> Symm1D::asArray(const double* soln) const
{
    AnyMap meta = Boundary1D::getMeta();
    return SolutionArray::create(m_solution, 0, meta);
}

// -------- Outlet1D --------

OutletRes1D::OutletRes1D()
{
    m_type = cOutletResType;
}

OutletRes1D::OutletRes1D(shared_ptr<Solution> solution, const string& id)
    : OutletRes1D()
{
    m_solution = solution;
    m_id = id;
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

shared_ptr<SolutionArray> Outlet1D::asArray(const double* soln) const
{
    AnyMap meta = Boundary1D::getMeta();
    return SolutionArray::create(m_solution, 0, meta);
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

shared_ptr<SolutionArray> OutletRes1D::asArray(const double* soln) const
{
    AnyMap meta = Boundary1D::getMeta();
    meta["temperature"] = m_temp;
    auto arr = SolutionArray::create(m_solution, 1, meta);

    // set gas state (using pressure from adjacent domain)
    double pressure = m_flow->phase().pressure();
    auto phase = m_solution->thermo();
    phase->setState_TPY(m_temp, pressure, &m_yres[0]);
    vector_fp data(phase->stateSize());
    phase->saveState(data);

    arr->setState(0, data);
    return arr;
}

void OutletRes1D::fromArray(SolutionArray& arr, double* soln)
{
    Boundary1D::setMeta(arr.meta());
    arr.setLoc(0);
    auto phase = arr.thermo();
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

shared_ptr<SolutionArray> Surf1D::asArray(const double* soln) const
{
    AnyMap meta = Boundary1D::getMeta();
    meta["temperature"] = m_temp;
    return SolutionArray::create(m_solution, 0, meta);
}

void Surf1D::fromArray(SolutionArray& arr, double* soln)
{
    Boundary1D::setMeta(arr.meta());
    arr.setLoc(0);
    m_temp = arr.thermo()->temperature();
}

void Surf1D::show(std::ostream& s, const double* x)
{
    s << "-------------------  Surface " << domainIndex() << " ------------------- " << std::endl;
    s << "  temperature: " << m_temp << " K" << std::endl;
}

void Surf1D::show(const double* x)
{
    writelog("    Temperature: {:10.4g} K \n\n", m_temp);
}

// -------- ReactingSurf1D --------

ReactingSurf1D::ReactingSurf1D()
    : m_kin(0)
    , m_surfindex(0)
    , m_nsp(0)
{
    m_type = cSurfType;
}

ReactingSurf1D::ReactingSurf1D(shared_ptr<Solution> solution, const std::string& id)
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
    m_solution = solution;
    m_id = id;
    m_kin = kin.get();
    m_sphase = phase.get();

    m_surfindex = m_kin->reactionPhaseIndex();
    m_nsp = m_sphase->nSpecies();
    m_enabled = true;
}

void ReactingSurf1D::setKinetics(shared_ptr<Kinetics> kin)
{
    m_solution = Solution::create();
    m_solution->setThermo(kin->reactionPhase());
    m_solution->setKinetics(kin);
    m_solution->setTransportModel("none");
    m_kin = dynamic_pointer_cast<InterfaceKinetics>(kin).get();
    m_surfindex = kin->reactionPhaseIndex();
    m_sphase = dynamic_pointer_cast<SurfPhase>(kin->reactionPhase()).get();
    m_nsp = m_sphase->nSpecies();
    m_enabled = true;
}

void ReactingSurf1D::setKineticsMgr(InterfaceKinetics* kin)
{
    warn_deprecated("ReactingSurf1D::setKineticsMgr",
        "To be removed after Cantera 3.0. Replaced by Domain1D::setKinetics.");
    m_kin = kin;
    m_surfindex = kin->reactionPhaseIndex();
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
        for (size_t k = 0; k < m_nsp; k++) {
            r[k] = m_work[k + ioffset] * m_sphase->size(k) * rs0;
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
        const vector_fp& mwleft = m_phase_left->molecularWeights();
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

shared_ptr<SolutionArray> ReactingSurf1D::asArray(const double* soln) const
{
    AnyMap meta = Boundary1D::getMeta();
    meta["temperature"] = m_temp;
    meta["phase"]["name"] = m_sphase->name();
    AnyValue source = m_sphase->input().getMetadata("filename");
    meta["phase"]["source"] = source.empty() ? "<unknown>" : source.asString();

    // set state of surface phase
    m_sphase->setState_TP(m_temp, m_sphase->pressure());
    m_sphase->setCoverages(soln);
    vector_fp data(m_sphase->stateSize());
    m_sphase->saveState(data.size(), &data[0]);

    auto arr = SolutionArray::create(m_solution, 1, meta);
    arr->setState(0, data);
    return arr;
}

void ReactingSurf1D::fromArray(SolutionArray& arr, double* soln)
{
    Boundary1D::setMeta(arr.meta());
    arr.setLoc(0);
    auto surf = std::dynamic_pointer_cast<SurfPhase>(arr.thermo());
    if (!surf) {
        throw CanteraError("ReactingSurf1D::fromArray",
            "Restoring of coverages requires surface phase");
    }
    m_temp = surf->temperature();
    surf->getCoverages(soln);
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
