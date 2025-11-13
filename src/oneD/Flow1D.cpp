//! @file Flow1D.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/SolutionArray.h"
#include "cantera/oneD/Flow1D.h"
#include "cantera/oneD/refine.h"
#include "cantera/transport/Transport.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/numerics/funcs.h"
#include "cantera/base/global.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/AnyMap.h"


using namespace std;

namespace Cantera
{

namespace { // restrict scope of auxiliary variables to local translation unit

const map<string, size_t> componentMap = {
    {"velocity", c_offset_U}, // axial velocity [m/s]
    {"spreadRate", c_offset_V}, // strain rate
    {"T", c_offset_T}, // temperature [kelvin]
    {"Lambda", c_offset_L}, // (1/r)dP/dr
    {"eField", c_offset_E}, // electric field
    {"Uo", c_offset_Uo}, // oxidizer axial velocity [m/s]
};

} // end unnamed namespace

Flow1D::Flow1D(ThermoPhase* ph, size_t nsp, size_t points) :
    Domain1D(nsp+c_offset_Y, points),
    m_nsp(nsp)
{
    warn_deprecated("Flow1D::Flow1D(ThermoPhase*, size_t, size_t)",
        "To be removed after Cantera 3.2. Use constructor using Solution instead.");
    _init(ph, nsp, points);
}

void Flow1D::_init(ThermoPhase* ph, size_t nsp, size_t points)
{
    m_points = points;

    if (ph == 0) {
        return; // used to create a dummy object
    }
    m_thermo = ph;

    size_t nsp2 = m_thermo->nSpecies();
    if (nsp2 != m_nsp) {
        m_nsp = nsp2;
        Domain1D::resize(m_nsp+c_offset_Y, points);
    }

    // make a local copy of the species molecular weight vector
    m_wt = m_thermo->molecularWeights();

    // set pressure based on associated thermo object
    setPressure(m_thermo->pressure());

    // the species mass fractions are the last components in the solution
    // vector, so the total number of components is the number of species
    // plus the offset of the first mass fraction.
    m_nv = c_offset_Y + m_nsp;

    // Turn off the energy equation at all points
    m_do_energy.resize(m_points,false);

    m_diff.resize(m_nsp*m_points);
    m_multidiff.resize(m_nsp*m_nsp*m_points);
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_hk.resize(m_nsp, m_points, 0.0);
    m_dhk_dz.resize(m_nsp, m_points - 1, 0.0);
    m_ybar.resize(m_nsp);
    m_qdotRadiation.resize(m_points, 0.0);

    //-------------- default solution bounds --------------------
    setBounds(c_offset_U, -1e20, 1e20); // no bounds on u
    setBounds(c_offset_V, -1e20, 1e20); // no bounds on V
    setBounds(c_offset_T, 200.0, 2*m_thermo->maxTemp()); // temperature bounds
    setBounds(c_offset_L, -1e20, 1e20); // Lambda should be negative
    setBounds(c_offset_E, -1e20, 1e20); // no bounds on electric field
    setBounds(c_offset_Uo, -1e20, 1e20); // no bounds on Uo
    // mass fraction bounds
    for (size_t k = 0; k < m_nsp; k++) {
        setBounds(c_offset_Y+k, -1.0e-7, 1.0e5);
    }

    //-------------------- grid refinement -------------------------
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    m_refiner->setActive(c_offset_L, false);
    m_refiner->setActive(c_offset_Uo, false);

    vector<double> gr;
    for (size_t ng = 0; ng < m_points; ng++) {
        gr.push_back(1.0*ng/m_points);
    }
    setupGrid(m_points, gr.data());

    // Initialize the radiation object (hardcoded for now)
    std::string propertyModel = "TabularPlanckMean";
    std::string solverModel = "OpticallyThin";

    // Define lambdas for T(x, j) and X(x, k, j)
    auto Tfunc = [this](const double* x, size_t j) {
        return this->T(x, j);
    };
    auto Xfunc = [this](const double* x, size_t k, size_t j) {
        return this->X(x, k, j);
    };

    double emissivityLeft = 0.0;
    double emissivityRight = 0.0;

    // Optionally override from 'radiation-parameters.yaml'
    try {
        AnyMap radcfg = AnyMap::fromYamlFile("radiation-parameters.yaml");
        if (radcfg.hasKey("Radiation")) {
            auto& rad = radcfg["Radiation"].as<AnyMap>();
            if (rad.hasKey("property")) {
                propertyModel = rad["property"].asString();
            }
            if (rad.hasKey("solver")) {
                solverModel = rad["solver"].asString();
            }
            if (rad.hasKey("emissivity")) {
                auto& e = rad["emissivity"].as<AnyMap>();
                if (e.hasKey("left")) {
                    emissivityLeft = e["left"].asDouble();
                }
                if (e.hasKey("right")) {
                    emissivityRight = e["right"].asDouble();
                }
            }
        }
    } catch (CanteraError&) {
        // Optional file; fall back to defaults silently
    }
    bool usedCustomRadlib = false;
    // If a RadLib property is selected, allow optional RadLib-specific parameters
    // to be passed from the 'Radiation.radlib' block in radiation-parameters.yaml
    if (propertyModel.rfind("RadLib.", 0) == 0 || propertyModel.rfind("radlib-", 0) == 0) {
        try {
            AnyMap radcfg = AnyMap::fromYamlFile("radiation-parameters.yaml");
            if (radcfg.hasKey("Radiation")) {
                auto& rad = radcfg["Radiation"].as<AnyMap>();
                if (rad.hasKey("radlib")) {
                    auto& ropt = rad["radlib"].as<AnyMap>();
                    double fvsoot = ropt.hasKey("fvsoot") ? ropt["fvsoot"].asDouble() : 0.0;
                    int nGray = ropt.hasKey("nGray") ? static_cast<int>(ropt["nGray"].asInt()) : 25;
                    double Tref = ropt.hasKey("Tref") ? ropt["Tref"].asDouble() : 1500.0;
                    double Pref = ropt.hasKey("Pref") ? ropt["Pref"].asDouble() : m_thermo->pressure();

                    std::unique_ptr<RadiationPropertyCalculator> props;
                    if (propertyModel == "RadLib.PlanckMean" || propertyModel == "radlib-pm") {
                        props = makeRadLibProps("RadLib.PlanckMean", m_thermo, fvsoot);
                    } else if (propertyModel == "RadLib.WSGG" || propertyModel == "radlib-wsgg") {
                        props = makeRadLibProps("RadLib.WSGG", m_thermo, fvsoot);
                    } else if (propertyModel == "RadLib.RCSLW" || propertyModel == "radlib-rcslw") {
                        props = makeRadLibProps("RadLib.RCSLW", m_thermo, fvsoot, nGray, Tref, Pref);
                    }

                    if (props) {
                        std::unique_ptr<RadiationSolver> solver;
                        if (solverModel == "OpticallyThin") {
                            solver = std::make_unique<OpticallyThinSolver>();
                        } else {
                            throw CanteraError("Flow1D::Flow1D",
                                "Unknown radiation solver model: {}", solverModel);
                        }
                        m_radiation = std::make_unique<Radiation1D>(
                            m_thermo, m_thermo->pressure(), m_points,
                            Tfunc, Xfunc, std::move(props), std::move(solver));
                        m_radiation->setBoundaryEmissivities(emissivityLeft, emissivityRight);
                        usedCustomRadlib = true;
                    }
                }
            }
        } catch (CanteraError&) {
            // optional config; ignore and fall back to defaults handled below
        }
    }

    if (!usedCustomRadlib) {
        m_radiation = createRadiation1D(
            propertyModel,
            solverModel,
            m_thermo,
            m_thermo->pressure(),
            m_points,
            Tfunc,
            Xfunc,
            emissivityLeft,
            emissivityRight
        );
    }
}

Flow1D::Flow1D(shared_ptr<ThermoPhase> th, size_t nsp, size_t points)
    : Flow1D(th.get(), nsp, points)
{
    auto sol = Solution::create();
    sol->setThermo(th);
    setSolution(sol);
}

Flow1D::Flow1D(shared_ptr<Solution> phase, const string& id, size_t points)
    : Domain1D(phase->thermo()->nSpecies()+c_offset_Y, points)
    , m_nsp(phase->thermo()->nSpecies())
{
    _init(phase->thermo().get(), phase->thermo()->nSpecies(), points);

    m_solution = phase;
    m_solution->thermo()->addSpeciesLock();
    m_id = id;
    m_kin = m_solution->kinetics().get();
    m_trans = m_solution->transport().get();
    if (m_trans->transportModel() == "none") {
        throw CanteraError("Flow1D::Flow1D",
            "An appropriate transport model\nshould be set when instantiating the "
            "Solution ('gas') object.");
    }
    m_solution->registerChangedCallback(this, [this]() {
        _setKinetics(m_solution->kinetics());
        _setTransport(m_solution->transport());
    });
}

Flow1D::~Flow1D()
{
    if (m_solution) {
        m_solution->removeChangedCallback(this);
    }
}

double Flow1D::lambda(const double* x, size_t j) const {
    warn_deprecated("Flow1D::lambda",
        "To be removed after Cantera 3.2. Component 'lambda' is renamed to 'Lambda'.");
    return x[index(c_offset_L, j)];
}

string Flow1D::domainType() const {
    if (m_isFree) {
        return "free-flow";
    }
    if (m_usesLambda) {
        return "axisymmetric-flow";
    }
    return "unstrained-flow";
}

void Flow1D::setKinetics(shared_ptr<Kinetics> kin)
{
    warn_deprecated("Flow1D::setKinetics",
        "After Cantera 3.2, a change of domain contents after instantiation "
        "will be disabled.");
    _setKinetics(kin);
}

void Flow1D::_setKinetics(shared_ptr<Kinetics> kin)
{
    m_kin = kin.get();
    m_solution->setKinetics(kin);
}

void Flow1D::setTransport(shared_ptr<Transport> trans)
{
    warn_deprecated("Flow1D::setTransport",
        "After Cantera 3.2, a change of domain contents after instantiation "
        "will be disabled.");
    _setTransport(trans);
}

void Flow1D::_setTransport(shared_ptr<Transport> trans)
{
    if (!trans) {
        throw CanteraError("Flow1D::_setTransport", "Unable to set empty transport.");
    }
    m_trans = trans.get();
    if (m_trans->transportModel() == "none") {
        throw CanteraError("Flow1D::_setTransport", "Invalid Transport model 'none'.");
    }
    m_do_multicomponent = (m_trans->transportModel() == "multicomponent" ||
        m_trans->transportModel() == "multicomponent-CK");

    m_diff.resize(m_nsp * m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp * m_nsp*m_points);
    }
    m_dthermal.resize(m_nsp, m_points, 0.0);
    m_solution->setTransport(trans);
}

void Flow1D::resize(size_t ncomponents, size_t points)
{
    Domain1D::resize(ncomponents, points);
    m_rho.resize(m_points, 0.0);
    m_wtm.resize(m_points, 0.0);
    m_cp.resize(m_points, 0.0);
    m_visc.resize(m_points, 0.0);
    m_tcon.resize(m_points, 0.0);

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
    }
    m_dthermal.resize(m_nsp, m_points, 0.0);
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_hk.resize(m_nsp, m_points, 0.0);
    m_dhk_dz.resize(m_nsp, m_points - 1, 0.0);
    m_do_energy.resize(m_points,false);
    m_qdotRadiation.resize(m_points, 0.0);
    m_fixedtemp.resize(m_points);

    m_dz.resize(m_points-1);
    m_z.resize(m_points);
}

void Flow1D::setupGrid(size_t n, const double* z)
{
    resize(m_nv, n);

    m_z[0] = z[0];
    for (size_t j = 1; j < m_points; j++) {
        if (z[j] <= z[j-1]) {
            throw CanteraError("Flow1D::setupGrid",
                               "grid points must be monotonically increasing");
        }
        m_z[j] = z[j];
        m_dz[j-1] = m_z[j] - m_z[j-1];
    }
}

void Flow1D::resetBadValues(double* xg)
{
    double* x = xg + loc();
    for (size_t j = 0; j < m_points; j++) {
        double* Y = x + m_nv*j + c_offset_Y;
        m_thermo->setMassFractions(Y);
        m_thermo->getMassFractions(Y);
    }
}

void Flow1D::setTransportModel(const string& model)
{
    if (model == "none") {
        throw CanteraError("Flow1D::setTransportModel",
            "Invalid Transport model 'none'.");
    }
    m_solution->setTransportModel(model);
    Flow1D::_setTransport(m_solution->transport());
}

string Flow1D::transportModel() const {
    return m_trans->transportModel();
}

void Flow1D::setFluxGradientBasis(ThermoBasis fluxGradientBasis) {
    m_fluxGradientBasis = fluxGradientBasis;
    if (transportModel() != "mixture-averaged-CK" &&
        transportModel() != "mixture-averaged") {
        warn_user("Flow1D::setFluxGradientBasis",
                  "Setting fluxGradientBasis only affects "
                  "the mixture-averaged diffusion model.");
    }
}

void Flow1D::_getInitialSoln(double* x)
{
    for (size_t j = 0; j < m_points; j++) {
        T(x,j) = m_thermo->temperature();
        m_thermo->getMassFractions(&Y(x, 0, j));
        m_rho[j] = m_thermo->density();
    }
}

void Flow1D::setGas(const double* x, size_t j)
{
    m_thermo->setTemperature(T(x,j));
    const double* yy = x + m_nv*j + c_offset_Y;
    m_thermo->setMassFractions_NoNorm(yy);
    m_thermo->setPressure(m_press);
}

void Flow1D::setGasAtMidpoint(const double* x, size_t j)
{
    m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
    const double* yy_j = x + m_nv*j + c_offset_Y;
    const double* yy_j_plus1 = x + m_nv*(j+1) + c_offset_Y;
    for (size_t k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yy_j[k] + yy_j_plus1[k]);
    }
    m_thermo->setMassFractions_NoNorm(m_ybar.data());
    m_thermo->setPressure(m_press);
}

void Flow1D::_finalize(const double* x)
{
    if (!(m_do_multicomponent ||
          m_trans->transportModel() == "mixture-averaged" ||
          m_trans->transportModel() == "mixture-averaged-CK")
        && m_do_soret) {

        throw CanteraError("Flow1D::_finalize",
            "Thermal diffusion (the Soret effect) is enabled, but it "
            "only ompatible with the mixture-averaged and multicomponent "
            "transport models.");
    }

    size_t nz = m_zfix.size();
    bool e = m_do_energy[0];
    for (size_t j = 0; j < m_points; j++) {
        if (e || nz == 0) {
            m_fixedtemp[j] = T(x, j);
        } else {
            double zz = (z(j) - z(0))/(z(m_points - 1) - z(0));
            double tt = linearInterp(zz, m_zfix, m_tfix);
            m_fixedtemp[j] = tt;
        }
    }
    if (e) {
        solveEnergyEqn();
    }

    if (m_isFree) {
        // If the domain contains the temperature fixed point, make sure that it
        // is correctly set. This may be necessary when the grid has been modified
        // externally.
        if (m_tfixed != Undef) {
            for (size_t j = 0; j < m_points; j++) {
                if (z(j) == m_zfixed) {
                    return; // fixed point is already set correctly
                }
            }

            for (size_t j = 0; j < m_points - 1; j++) {
                // Find where the temperature profile crosses the current
                // fixed temperature.
                if ((T(x, j) - m_tfixed) * (T(x, j+1) - m_tfixed) <= 0.0) {
                    m_tfixed = T(x, j+1);
                    m_zfixed = z(j+1);
                    return;
                }
            }
        }
    }
}

void Flow1D::eval(size_t jGlobal, double* xGlobal, double* rsdGlobal,
                  integer* diagGlobal, double rdt)
{
    // If evaluating a Jacobian, and the global point is outside the domain of
    // influence for this domain, then skip evaluating the residual
    if (jGlobal != npos && (jGlobal + 1 < firstPoint() || jGlobal > lastPoint() + 1)) {
        return;
    }

    // start of local part of global arrays
    double* x = xGlobal + loc();
    double* rsd = rsdGlobal + loc();
    integer* diag = diagGlobal + loc();

    size_t jmin, jmax;
    if (jGlobal == npos) { // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jGlobal == 0) ? 0 : jGlobal - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    updateProperties(jGlobal, x, jmin, jmax);

    if (m_do_radiation) { // Calculation of qdotRadiation
        computeRadiation(x, jmin, jmax);
    }

    evalContinuity(x, rsd, diag, rdt, jmin, jmax);
    evalMomentum(x, rsd, diag, rdt, jmin, jmax);
    evalEnergy(x, rsd, diag, rdt, jmin, jmax);
    evalLambda(x, rsd, diag, rdt, jmin, jmax);
    evalElectricField(x, rsd, diag, rdt, jmin, jmax);
    evalUo(x, rsd, diag, rdt, jmin, jmax);
    evalSpecies(x, rsd, diag, rdt, jmin, jmax);
}

void Flow1D::updateProperties(size_t jg, double* x, size_t jmin, size_t jmax)
{
    // properties are computed for grid points from j0 to j1
    size_t j0 = std::max<size_t>(jmin, 1) - 1;
    size_t j1 = std::min(jmax+1,m_points-1);

    updateThermo(x, j0, j1);
    if (jg == npos || m_force_full_update) {
        // update transport properties only if a Jacobian is not being
        // evaluated, or if specifically requested
        updateTransport(x, j0, j1);
    }
    if (jg == npos) {
        double* Yleft = x + index(c_offset_Y, jmin);
        m_kExcessLeft = distance(Yleft, max_element(Yleft, Yleft + m_nsp));
        double* Yright = x + index(c_offset_Y, jmax);
        m_kExcessRight = distance(Yright, max_element(Yright, Yright + m_nsp));
    }

    // update the species diffusive mass fluxes whether or not a
    // Jacobian is being evaluated
    updateDiffFluxes(x, j0, j1);
}

void Flow1D::updateTransport(double* x, size_t j0, size_t j1)
{
     if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            double wtm = m_thermo->meanMolecularWeight();
            double rho = m_thermo->density();
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMultiDiffCoeffs(m_nsp, &m_multidiff[mindex(0,0,j)]);

            // Use m_diff as storage for the factor outside the summation
            for (size_t k = 0; k < m_nsp; k++) {
                m_diff[k+j*m_nsp] = m_wt[k] * rho / (wtm*wtm);
            }

            m_tcon[j] = m_trans->thermalConductivity();
            if (m_do_soret) {
                m_trans->getThermalDiffCoeffs(m_dthermal.ptrColumn(0) + j*m_nsp);
            }
        }
    } else { // mixture averaged transport
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);

            if (m_fluxGradientBasis == ThermoBasis::molar) {
                m_trans->getMixDiffCoeffs(&m_diff[j*m_nsp]);
            } else {
                m_trans->getMixDiffCoeffsMass(&m_diff[j*m_nsp]);
            }

            double rho = m_thermo->density();

            if (m_fluxGradientBasis == ThermoBasis::molar) {
                double wtm = m_thermo->meanMolecularWeight();
                for (size_t k=0; k < m_nsp; k++) {
                    m_diff[k+j*m_nsp] *= m_wt[k] * rho / wtm;
                }
            } else {
                for (size_t k=0; k < m_nsp; k++) {
                    m_diff[k+j*m_nsp] *= rho;
                }
            }
            m_tcon[j] = m_trans->thermalConductivity();
            if (m_do_soret) {
                m_trans->getThermalDiffCoeffs(m_dthermal.ptrColumn(0) + j*m_nsp);
            }
        }
    }
}

void Flow1D::updateDiffFluxes(const double* x, size_t j0, size_t j1)
{
    if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            double dz = z(j+1) - z(j);
            for (size_t k = 0; k < m_nsp; k++) {
                double sum = 0.0;
                for (size_t m = 0; m < m_nsp; m++) {
                    sum += m_wt[m] * m_multidiff[mindex(k,m,j)] * (X(x,m,j+1)-X(x,m,j));
                }
                m_flux(k,j) = sum * m_diff[k+j*m_nsp] / dz;
            }
        }
    } else {
        for (size_t j = j0; j < j1; j++) {
            double sum = 0.0;
            double dz = z(j+1) - z(j);
            if (m_fluxGradientBasis == ThermoBasis::molar) {
                for (size_t k = 0; k < m_nsp; k++) {
                    m_flux(k,j) = m_diff[k+m_nsp*j] * (X(x,k,j) - X(x,k,j+1))/dz;
                    sum -= m_flux(k,j);
                }
            } else {
                for (size_t k = 0; k < m_nsp; k++) {
                    m_flux(k,j) = m_diff[k+m_nsp*j] * (Y(x,k,j) - Y(x,k,j+1))/dz;
                    sum -= m_flux(k,j);
                }
            }
            // correction flux to ensure that \sum_k Y_k V_k = 0.
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,j) += sum*Y(x,k,j);
            }
        }
    }

    if (m_do_soret) {
        for (size_t m = j0; m < j1; m++) {
            double gradlogT = 2.0 * (T(x,m+1) - T(x,m)) /
                              ((T(x,m+1) + T(x,m)) * (z(m+1) - z(m)));
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,m) -= m_dthermal(k,m)*gradlogT;
            }
        }
    }
}

void Flow1D::computeRadiation(double* x, size_t jmin, size_t jmax)
{
    m_radiation->computeRadiation(x, jmin, jmax, m_qdotRadiation);
}

void Flow1D::setRadiationModels(const std::string& propertyModel,
                                const std::string& solverModel)
{
    // Rebuild the Radiation1D object with the requested models, preserving
    // current boundary emissivities
    double emissivityLeft = 0.0;
    double emissivityRight = 0.0;
    if (m_radiation) {
        emissivityLeft = m_radiation->leftEmissivity();
        emissivityRight = m_radiation->rightEmissivity();
    }
    m_radiation = createRadiation1D(propertyModel, solverModel,
        m_thermo,
        m_press,
        m_points,
        // point accessors you already use:
        [this](const double* x, size_t j){ return this->T(x, j); },
        [this](const double* x, size_t k, size_t j){ return this->X(x, k, j); },
        emissivityLeft,
        emissivityRight
    );
}

void Flow1D::evalContinuity(double* x, double* rsd, int* diag,
                            double rdt, size_t jmin, size_t jmax)
{
    // The left boundary has the same form for all cases.
    if (jmin == 0) { // left boundary
        rsd[index(c_offset_U, jmin)] = -(rho_u(x, jmin+1) - rho_u(x, jmin))/m_dz[jmin]
                                       -(density(jmin+1)*V(x, jmin+1)
                                       + density(jmin)*V(x, jmin));
        diag[index(c_offset_U, jmin)] = 0; // Algebraic constraint
    }

    if (jmax == m_points - 1) { // right boundary
        if (m_usesLambda) { // zero mass flux
            rsd[index(c_offset_U, jmax)] = rho_u(x, jmax);
        } else { // zero gradient, same for unstrained or free-flow
            rsd[index(c_offset_U, jmax)] = rho_u(x, jmax) - rho_u(x, jmax-1);
        }
        diag[index(c_offset_U, jmax)] = 0; // Algebraic constraint
    }

    // j0 and j1 are constrained to only interior points
    size_t j0 = std::max<size_t>(jmin, 1);
    size_t j1 = std::min(jmax, m_points-2);
    if (m_usesLambda) { // "axisymmetric-flow"
        for (size_t j = j0; j <= j1; j++) { // interior points
            // For "axisymmetric-flow", the continuity equation  propagates the
            // mass flow rate information to the left (j+1 -> j) from the value
            // specified at the right boundary. The Lambda information propagates
            // in the opposite direction.
            rsd[index(c_offset_U, j)] = -(rho_u(x, j+1) - rho_u(x, j))/m_dz[j]
                                        -(density(j+1)*V(x, j+1) + density(j)*V(x, j));
            diag[index(c_offset_U, j)] = 0; // Algebraic constraint
        }
    } else if (m_isFree) { // "free-flow"
        for (size_t j = j0; j <= j1; j++) {
            // terms involving V are zero as V=0 by definition
            if (z(j) > m_zfixed) {
                rsd[index(c_offset_U, j)] = -(rho_u(x, j) - rho_u(x, j-1))/m_dz[j-1];
            } else if (z(j) == m_zfixed) {
                if (m_do_energy[j]) {
                    rsd[index(c_offset_U, j)] = (T(x, j) - m_tfixed);
                } else {
                    rsd[index(c_offset_U, j)] = (rho_u(x, j) - m_rho[0]*0.3); // why 0.3?
                }
            } else { // z(j) < m_zfixed
                rsd[index(c_offset_U, j)] = -(rho_u(x, j+1) - rho_u(x, j))/m_dz[j];
            }
            diag[index(c_offset_U, j)] = 0; // Algebraic constraint
        }
    } else { // "unstrained-flow" (fixed mass flow rate)
        for (size_t j = j0; j <= j1; j++) {
            rsd[index(c_offset_U, j)] = rho_u(x, j) - rho_u(x, j-1);
            diag[index(c_offset_U, j)] = 0; // Algebraic constraint
        }
    }
}

void Flow1D::evalMomentum(double* x, double* rsd, int* diag,
                          double rdt, size_t jmin, size_t jmax)
{
    if (!m_usesLambda) { //disable this equation
        for (size_t j = jmin; j <= jmax; j++) {
            rsd[index(c_offset_V, j)] = V(x, j);
            diag[index(c_offset_V, j)] = 0;
        }
        return;
    }

    if (jmin == 0) { // left boundary
        rsd[index(c_offset_V, jmin)] = V(x, jmin);
    }

    if (jmax == m_points - 1) { // right boundary
        rsd[index(c_offset_V, jmax)] = V(x, jmax);
        diag[index(c_offset_V, jmax)] = 0;
    }

    // j0 and j1 are constrained to only interior points
    size_t j0 = std::max<size_t>(jmin, 1);
    size_t j1 = std::min(jmax, m_points-2);
    for (size_t j = j0; j <= j1; j++) { // interior points
        rsd[index(c_offset_V, j)] = (shear(x, j) - Lambda(x, j)
                                     - rho_u(x, j) * dVdz(x, j)
                                     - m_rho[j] * V(x, j) * V(x, j)) / m_rho[j];
        if (!m_twoPointControl) {
            rsd[index(c_offset_V, j)] -= rdt * (V(x, j) - V_prev(j));
            diag[index(c_offset_V, j)] = 1;
        } else {
            diag[index(c_offset_V, j)] = 0;
        }
    }
}

void Flow1D::evalLambda(double* x, double* rsd, int* diag,
                        double rdt, size_t jmin, size_t jmax)
{
    if (!m_usesLambda) { // disable this equation
        for (size_t j = jmin; j <= jmax; j++) {
            rsd[index(c_offset_L, j)] = Lambda(x, j);
            diag[index(c_offset_L, j)] = 0;
        }
        return;
    }

    if (jmin == 0) { // left boundary
        if (m_twoPointControl) {
            rsd[index(c_offset_L, jmin)] = Lambda(x, jmin+1) - Lambda(x, jmin);
        } else {
            rsd[index(c_offset_L, jmin)] = -rho_u(x, jmin);
        }
    }

    if (jmax == m_points - 1) { // right boundary
        rsd[index(c_offset_L, jmax)] = Lambda(x, jmax) - Lambda(x, jmax-1);
        diag[index(c_offset_L, jmax)] = 0;
    }

    // j0 and j1 are constrained to only interior points
    size_t j0 = std::max<size_t>(jmin, 1);
    size_t j1 = std::min(jmax, m_points-2);
    for (size_t j = j0; j <= j1; j++) { // interior points
        if (m_twoPointControl) {
            if (z(j) == m_zLeft) {
                rsd[index(c_offset_L, j)] = T(x,j) - m_tLeft;
            } else if (z(j) > m_zLeft) {
                rsd[index(c_offset_L, j)] = Lambda(x, j) - Lambda(x, j-1);
            } else if (z(j) < m_zLeft) {
                rsd[index(c_offset_L, j)] = Lambda(x, j) - Lambda(x, j+1);
            }
        } else {
            rsd[index(c_offset_L, j)] = Lambda(x, j) - Lambda(x, j-1);
        }
        diag[index(c_offset_L, j)] = 0;
    }
}

void Flow1D::evalEnergy(double* x, double* rsd, int* diag,
                        double rdt, size_t jmin, size_t jmax)
{
    if (jmin == 0) { // left boundary
        rsd[index(c_offset_T, jmin)] = T(x, jmin);
    }

    if (jmax == m_points - 1) { // right boundary
        rsd[index(c_offset_T, jmax)] = T(x, jmax);
    }

    // j0 and j1 are constrained to only interior points
    size_t j0 = std::max<size_t>(jmin, 1);
    size_t j1 = std::min(jmax, m_points-2);
    for (size_t j = j0; j <= j1; j++) {
        if (m_do_energy[j]) {
            grad_hk(x, j);
            double sum = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                double flxk = 0.5*(m_flux(k, j-1) + m_flux(k, j));
                sum += m_wdot(k, j)*m_hk(k, j);
                sum += flxk * m_dhk_dz(k, j) / m_wt[k];
            }

            rsd[index(c_offset_T, j)] = - m_cp[j]*rho_u(x, j)*dTdz(x, j)
                                        - conduction(x, j) - sum;
            rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);
            rsd[index(c_offset_T, j)] -= (m_qdotRadiation[j] / (m_rho[j] * m_cp[j]));
            if (!m_twoPointControl || (m_z[j] != m_tLeft && m_z[j] != m_tRight)) {
                rsd[index(c_offset_T, j)] -= rdt*(T(x, j) - T_prev(j));
                diag[index(c_offset_T, j)] = 1;
            } else {
                diag[index(c_offset_T, j)] = 0;
            }
        } else {
            // residual equations if the energy equation is disabled
            rsd[index(c_offset_T, j)] = T(x, j) - T_fixed(j);
            diag[index(c_offset_T, j)] = 0;
        }
    }
}

void Flow1D::evalUo(double* x, double* rsd, int* diag,
                    double rdt, size_t jmin, size_t jmax)
{
    if (!m_twoPointControl) { // disable this equation
        for (size_t j = jmin; j <= jmax; j++) {
            rsd[index(c_offset_Uo, j)] = Uo(x, j);
            diag[index(c_offset_Uo, j)] = 0;
        }
        return;
    }

    if (jmin == 0) { // left boundary
        rsd[index(c_offset_Uo, jmin)] = Uo(x, jmin+1) - Uo(x, jmin);
        diag[index(c_offset_Uo, jmin)] = 0;
    }

    if (jmax == m_points - 1) { // right boundary
        if(m_twoPointControl) {
            rsd[index(c_offset_Uo, jmax)] = Uo(x, jmax) - Uo(x, jmax-1);
        }
        diag[index(c_offset_Uo, jmax)] = 0;
    }

    // j0 and j1 are constrained to only interior points
    size_t j0 = std::max<size_t>(jmin, 1);
    size_t j1 = std::min(jmax, m_points-2);
    for (size_t j = j0; j <= j1; j++) { // interior points
        if (m_twoPointControl) {
            if (z(j) == m_zRight) {
                rsd[index(c_offset_Uo, j)] = T(x, j) - m_tRight;
            } else if (z(j) > m_zRight) {
                rsd[index(c_offset_Uo, j)] = Uo(x, j) - Uo(x, j-1);
            } else if (z(j) < m_zRight) {
                rsd[index(c_offset_Uo, j)] = Uo(x, j) - Uo(x, j+1);
            }
        }
        diag[index(c_offset_Uo, j)] = 0;
    }
}

void Flow1D::evalSpecies(double* x, double* rsd, int* diag,
                         double rdt, size_t jmin, size_t jmax)
{
    if (jmin == 0) { // left boundary
        double sum = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            sum += Y(x,k,jmin);
            rsd[index(c_offset_Y+k, jmin)] = -(m_flux(k, jmin) +
                                                rho_u(x, jmin) * Y(x, k, jmin));
        }
        rsd[index(c_offset_Y + leftExcessSpecies(), jmin)] = 1.0 - sum;
        diag[index(c_offset_Y + leftExcessSpecies(), jmin)] = 0;
    }

    if (jmax == m_points - 1) { // right boundary
        double sum = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            sum += Y(x,k,jmax);
            rsd[index(k+c_offset_Y, jmax)] = m_flux(k, jmax-1) +
                                             rho_u(x, jmax)*Y(x, k, jmax);
        }
        rsd[index(c_offset_Y + rightExcessSpecies(), jmax)] = 1.0 - sum;
        diag[index(c_offset_Y + rightExcessSpecies(), jmax)] = 0;
    }

    // j0 and j1 are constrained to only interior points
    size_t j0 = std::max<size_t>(jmin, 1);
    size_t j1 = std::min(jmax, m_points-2);
    for (size_t j = j0; j <= j1; j++) {
        for (size_t k = 0; k < m_nsp; k++) {
            double convec = rho_u(x, j)*dYdz(x, k, j);
            double diffus = 2*(m_flux(k, j) - m_flux(k, j-1)) / (z(j+1) - z(j-1));
            rsd[index(c_offset_Y + k, j)] = (m_wt[k]*m_wdot(k, j)
                                              - convec - diffus) / m_rho[j]
                                            - rdt*(Y(x, k, j) - Y_prev(k, j));
            diag[index(c_offset_Y + k, j)] = 1;
        }
    }
}

void Flow1D::evalElectricField(double* x, double* rsd, int* diag,
                               double rdt, size_t jmin, size_t jmax)
{
    for (size_t j = jmin; j <= jmax; j++) {
        // The same value is used for left/right/interior points
        rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
    }
}

void Flow1D::show(const double* x)
{
    writelog("    Pressure:  {:10.4g} Pa\n", m_press);

    Domain1D::show(x);

    if (m_do_radiation) {
        writeline('-', 79, false, true);
        writelog("\n          z      radiative heat loss");
        writeline('-', 79, false, true);
        for (size_t j = 0; j < m_points; j++) {
            writelog("\n {:10.4g}        {:10.4g}", m_z[j], m_qdotRadiation[j]);
        }
        writelog("\n");
    }
}

string Flow1D::componentName(size_t n) const
{
    switch (n) {
    case c_offset_U:
        return "velocity";
    case c_offset_V:
        return "spreadRate";
    case c_offset_T:
        return "T";
    case c_offset_L:
        return "Lambda";
    case c_offset_E:
        return "eField";
    case c_offset_Uo:
        return "Uo";
    default:
        if (n >= c_offset_Y && n < (c_offset_Y + m_nsp)) {
            return m_thermo->speciesName(n - c_offset_Y);
        } else {
            throw IndexError("Flow1D::componentName", "component", n, m_nv);
        }
    }
}

size_t Flow1D::componentIndex(const string& name, bool checkAlias) const
{
    if (componentMap.count(name)) {
        return componentMap.at(name);
    }
    for (size_t n = c_offset_Y; n < m_nsp + c_offset_Y; n++) {
        if (componentName(n) == name) {
            return n;
        }
    }
    if (checkAlias) {
        const auto& aliasMap = _componentAliasMap();
        if (aliasMap.count(name)) {
            return componentIndex(aliasMap.at(name), false);
        }
    }
    throw CanteraError("Flow1D::componentIndex",
                       "Component '{}' not found", name);
}

bool Flow1D::hasComponent(const string& name, bool checkAlias) const
{
    if (componentMap.count(name)) {
        return true;
    }
    for (size_t n = c_offset_Y; n < m_nsp + c_offset_Y; n++) {
        if (componentName(n) == name) {
            return true;
        }
    }
    if (checkAlias && _componentAliasMap().count(name)) {
        return true;
    }
    return false;
}

bool Flow1D::componentActive(size_t n) const
{
    switch (n) {
    case c_offset_V: // spreadRate
        return m_usesLambda;
    case c_offset_L: // Lambda
        return m_usesLambda;
    case c_offset_E: // eField
        return false;
    case c_offset_Uo: // oxidizer velocity for two-point control
        return twoPointControlEnabled();
    default:
        return true;
    }
}
AnyMap Flow1D::getMeta() const
{
    AnyMap state = Domain1D::getMeta();
    state["transport-model"] = m_trans->transportModel();

    state["phase"]["name"] = m_thermo->name();
    AnyValue source = m_thermo->input().getMetadata("filename");
    state["phase"]["source"] = source.empty() ? "<unknown>" : source.asString();

    state["radiation-enabled"] = m_do_radiation;
    if (m_do_radiation) {
        state["emissivity-left"] = (m_radiation ? m_radiation->leftEmissivity() : m_epsilon_left);
        state["emissivity-right"] = (m_radiation ? m_radiation->rightEmissivity() : m_epsilon_right);
    }

    set<bool> energy_flags(m_do_energy.begin(), m_do_energy.end());
    if (energy_flags.size() == 1) {
        state["energy-enabled"] = m_do_energy[0];
    } else {
        state["energy-enabled"] = m_do_energy;
    }

    state["Soret-enabled"] = m_do_soret;

    state["flux-gradient-basis"] = static_cast<long int>(m_fluxGradientBasis);

    state["refine-criteria"]["ratio"] = m_refiner->maxRatio();
    state["refine-criteria"]["slope"] = m_refiner->maxDelta();
    state["refine-criteria"]["curve"] = m_refiner->maxSlope();
    state["refine-criteria"]["prune"] = m_refiner->prune();
    state["refine-criteria"]["grid-min"] = m_refiner->gridMin();
    state["refine-criteria"]["max-points"] =
        static_cast<long int>(m_refiner->maxPoints());

    if (m_zfixed != Undef) {
        state["fixed-point"]["location"] = m_zfixed;
        state["fixed-point"]["temperature"] = m_tfixed;
    }

    // Two-point control meta data
    if (m_twoPointControl) {
        state["continuation-method"]["type"] = "two-point";
        state["continuation-method"]["left-location"] = m_zLeft;
        state["continuation-method"]["right-location"] = m_zRight;
        state["continuation-method"]["left-temperature"] = m_tLeft;
        state["continuation-method"]["right-temperature"] = m_tRight;
    }

    return state;
}

void Flow1D::updateState(size_t loc)
{
    if (!m_state) {
        throw CanteraError("Flow1D::updateState",
            "Domain needs to be installed in a container.");
    }
    if (!m_solution) {
        throw CanteraError("Flow1D::updateState",
            "Domain does not have associated Solution object.");
    }
    const double* soln = m_state->data() + m_iloc + m_nv * loc;
    m_solution->thermo()->setMassFractions_NoNorm(soln + c_offset_Y);
    m_solution->thermo()->setState_TP(*(soln + c_offset_T), m_press);
}

void Flow1D::getValues(const string& component, vector<double>& values) const
{
    if (!m_state) {
        throw CanteraError("Flow1D::getValues",
            "Domain needs to be installed in a container.");
    }
    if (values.size() != nPoints()) {
        throw ArraySizeError("Flow1D::getValues", values.size(), nPoints());
    }
    auto i = componentIndex(component);
    if (!componentActive(i)) {
        warn_user(
            "Flow1D::getValues", "Component '{}' is not used by '{}'.",
            component, domainType());
    }
    const double* soln = m_state->data() + m_iloc;
    for (size_t j = 0; j < nPoints(); j++) {
        values[j] = soln[index(i,j)];
    }
}

void Flow1D::setValues(const string& component, const vector<double>& values)
{
    if (!m_state) {
        throw CanteraError("Flow1D::setValues",
            "Domain needs to be installed in a container.");
    }
    if (values.size() != nPoints()) {
        throw ArraySizeError("Flow1D::_etValues", values.size(), nPoints());
    }
    auto i = componentIndex(component);
    if (!componentActive(i)) {
        throw CanteraError(
            "Flow1D::setValues", "Component '{}' is not used by '{}'.",
            component, domainType());
    }
    double* soln = m_state->data() + m_iloc;
    for (size_t j = 0; j < nPoints(); j++) {
        soln[index(i,j)] = values[j];
    }
}

void Flow1D::getResiduals(const string& component, vector<double>& values) const
{
    if (!m_state) {
        throw CanteraError("Flow1D::getResiduals",
            "Domain needs to be installed in a container.");
    }
    if (values.size() != nPoints()) {
        throw ArraySizeError("Flow1D::getResiduals", values.size(), nPoints());
    }
    auto i = componentIndex(component);
    if (!componentActive(i)) {
        warn_user(
            "Flow1D::getResiduals", "Component '{}' is not used by '{}'.",
            component, domainType());
    }
    const double* soln = m_container->_workVector().data() + m_iloc;
    for (size_t j = 0; j < nPoints(); j++) {
        values[j] = soln[index(i,j)];
    }
}

void Flow1D::setProfile(const string& component,
                        const vector<double>& pos, const vector<double>& values)
{
    if (!m_state) {
        throw CanteraError("Flow1D::setProfile",
            "Domain needs to be installed in a container.");
    }
    if (pos.size() != values.size()) {
        throw CanteraError(
            "Flow1D::setProfile", "Vectors for positions and values must have same "
            "size.\nSizes are {} and {}, respectively.", pos.size(), values.size());
    }
    if (pos.front() != 0.0 || pos.back() != 1.0) {
        throw CanteraError("Flow1D::setProfile",
            "'pos' vector must span the range [0, 1]. Got a vector spanning "
            "[{}, {}] instead.", pos.front(), pos.back());
    }
    auto i = componentIndex(component);
    if (!componentActive(i)) {
        throw CanteraError(
            "Flow1D::setProfile", "Component '{}' is not used by '{}'.",
            component, domainType());
    }
    double* soln = m_state->data() + m_iloc;
    double z0 = zmin();
    double zDelta = zmax() - z0;
    for (size_t j = 0; j < nPoints(); j++) {
        double frac = (z(j) - z0)/zDelta;
        double v = linearInterp(frac, pos, values);
        soln[index(i,j)] = v;
    }
}

void Flow1D::setFlatProfile(const string& component, double value)
{
    if (!m_state) {
        throw CanteraError("Flow1D::setFlatProfile",
            "Domain needs to be installed in a container.");
    }
    auto i = componentIndex(component);
    if (!componentActive(i)) {
        throw CanteraError(
            "Flow1D::setFlatProfile", "Component '{}' is not used by '{}'.",
            component, domainType());
    }
    double* soln = m_state->data() + m_iloc;
    for (size_t j = 0; j < nPoints(); j++) {
        soln[index(i,j)] = value;
    }
}

shared_ptr<SolutionArray> Flow1D::toArray(bool normalize)
{
    if (!m_state) {
        throw CanteraError("Flow1D::toArray",
            "Domain needs to be installed in a container before calling toArray.");
    }
    double* soln = m_state->data() + m_iloc;
    auto arr = SolutionArray::create(
        m_solution, static_cast<int>(nPoints()), getMeta());
    arr->addExtra("grid", false); // leading entry
    AnyValue value;
    value = m_z;
    arr->setComponent("grid", value);
    vector<double> data(nPoints());
    for (size_t i = 0; i < nComponents(); i++) {
        if (componentActive(i)) {
            auto name = componentName(i);
            for (size_t j = 0; j < nPoints(); j++) {
                data[j] = soln[index(i, j)];
            }
            if (!arr->hasComponent(name)) {
                arr->addExtra(name, componentIndex(name) > c_offset_Y);
            }
            value = data;
            arr->setComponent(name, value);
        }
    }
    updateThermo(soln, 0, m_points-1);
    value = m_rho;
    arr->setComponent("D", value); // use density rather than pressure

    if (m_do_radiation) {
        arr->addExtra("radiative-heat-loss", true); // add at end
        value = m_qdotRadiation;
        arr->setComponent("radiative-heat-loss", value);
    }

    if (normalize) {
        arr->normalize();
    }
    return arr;
}

void Flow1D::fromArray(const shared_ptr<SolutionArray>& arr)
{
    if (!m_state) {
        throw CanteraError("Domain1D::fromArray",
            "Domain needs to be installed in a container before calling fromArray.");
    }
    resize(nComponents(), arr->size());
    m_container->resize();
    double* soln = m_state->data() + m_iloc;

    Domain1D::setMeta(arr->meta());
    arr->setLoc(0);
    auto phase = arr->thermo();
    setPressure(phase->pressure());

    const auto grid = arr->getComponent("grid").as<vector<double>>();
    setupGrid(nPoints(), &grid[0]);
    setMeta(arr->meta()); // can affect which components are active

    for (size_t i = 0; i < nComponents(); i++) {
        if (!componentActive(i)) {
            continue;
        }
        string name = componentName(i);
        if (arr->hasComponent(name)) {
            const vector<double> data = arr->getComponent(name).as<vector<double>>();
            for (size_t j = 0; j < nPoints(); j++) {
                soln[index(i,j)] = data[j];
            }
        } else if (name == "Lambda" && arr->hasComponent("lambda")) {
            // edge case: 'lambda' is renamed to 'Lambda' in Cantera 3.2
            const auto data = arr->getComponent("lambda").as<vector<double>>();
            for (size_t j = 0; j < nPoints(); j++) {
                soln[index(i,j)] = data[j];
            }
        } else {
            warn_user("Flow1D::fromArray", "Saved state does not contain values for "
                "component '{}' in domain '{}'.", name, id());
        }
    }

    updateProperties(npos, soln, 0, m_points - 1);
    _finalize(soln);
}

void Flow1D::setMeta(const AnyMap& state)
{
    if (state.hasKey("energy-enabled")) {
        const AnyValue& ee = state["energy-enabled"];
        if (ee.isScalar()) {
            m_do_energy.assign(nPoints(), ee.asBool());
        } else {
            m_do_energy = ee.asVector<bool>(nPoints());
        }
    }

    if (state.hasKey("transport-model")) {
        setTransportModel(state["transport-model"].asString());
    }

    if (state.hasKey("Soret-enabled")) {
        m_do_soret = state["Soret-enabled"].asBool();
    }

    if (state.hasKey("flux-gradient-basis")) {
        m_fluxGradientBasis = static_cast<ThermoBasis>(
                state["flux-gradient-basis"].asInt());
    }

    if (state.hasKey("radiation-enabled")) {
        m_do_radiation = state["radiation-enabled"].asBool();
        if (m_do_radiation) {
            m_epsilon_left = state["emissivity-left"].asDouble();
            m_epsilon_right = state["emissivity-right"].asDouble();
            if (m_radiation) { m_radiation->setBoundaryEmissivities(m_epsilon_left, m_epsilon_right);}
        }
    }

    if (state.hasKey("refine-criteria")) {
        const AnyMap& criteria = state["refine-criteria"].as<AnyMap>();
        double ratio = criteria.getDouble("ratio", m_refiner->maxRatio());
        double slope = criteria.getDouble("slope", m_refiner->maxDelta());
        double curve = criteria.getDouble("curve", m_refiner->maxSlope());
        double prune = criteria.getDouble("prune", m_refiner->prune());
        m_refiner->setCriteria(ratio, slope, curve, prune);

        if (criteria.hasKey("grid-min")) {
            m_refiner->setGridMin(criteria["grid-min"].asDouble());
        }
        if (criteria.hasKey("max-points")) {
            m_refiner->setMaxPoints(criteria["max-points"].asInt());
        }
    }

    if (state.hasKey("fixed-point")) {
        m_zfixed = state["fixed-point"]["location"].asDouble();
        m_tfixed = state["fixed-point"]["temperature"].asDouble();
    }

    // Two-point control meta data
    if (state.hasKey("continuation-method")) {
        const AnyMap& cm = state["continuation-method"].as<AnyMap>();
        if (cm["type"] == "two-point") {
            m_twoPointControl = true;
            m_zLeft = cm["left-location"].asDouble();
            m_zRight = cm["right-location"].asDouble();
            m_tLeft = cm["left-temperature"].asDouble();
            m_tRight = cm["right-temperature"].asDouble();
        } else {
            warn_user("Flow1D::setMeta", "Unknown continuation method '{}'.",
                cm["type"].asString());
        }
    }
}

void Flow1D::solveEnergyEqn(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = true;
        }
    } else {
        if (!m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = true;
    }
    m_refiner->setActive(c_offset_U, true);
    m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    if (changed) {
        needJacUpdate();
    }
}

size_t Flow1D::getSolvingStage() const
{
    throw NotImplementedError("Flow1D::getSolvingStage",
        "Not used by '{}' objects.", type());
}

void Flow1D::setSolvingStage(const size_t stage)
{
    throw NotImplementedError("Flow1D::setSolvingStage",
        "Not used by '{}' objects.", type());
}

void Flow1D::solveElectricField(size_t j)
{
    throw NotImplementedError("Flow1D::solveElectricField",
        "Not used by '{}' objects.", type());
}

void Flow1D::fixElectricField(size_t j)
{
    throw NotImplementedError("Flow1D::fixElectricField",
        "Not used by '{}' objects.", type());
}

bool Flow1D::doElectricField(size_t j) const
{
    throw NotImplementedError("Flow1D::doElectricField",
        "Not used by '{}' objects.", type());
}

void Flow1D::setBoundaryEmissivities(double e_left, double e_right)
{
    m_radiation->setBoundaryEmissivities(e_left, e_right);
}

void Flow1D::fixTemperature(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = false;
        }
    } else {
        if (m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = false;
    }
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    if (changed) {
        needJacUpdate();
    }
}

void Flow1D::grad_hk(const double* x, size_t j)
{
    size_t jloc = (u(x, j) > 0.0 ? j : j + 1);
    for(size_t k = 0; k < m_nsp; k++) {
        m_dhk_dz(k, j) =  (m_hk(k, jloc) - m_hk(k, jloc-1))/m_dz[jloc-1];
    }
}

// Two-point control functions
double Flow1D::leftControlPointTemperature() const
{
    if (m_twoPointControl) {
        if (m_zLeft != Undef) {
            return m_tLeft;
        } else {
            throw CanteraError("Flow1D::leftControlPointTemperature",
                "Invalid operation: left control point location is not set");
        }
    } else {
        throw CanteraError("Flow1D::leftControlPointTemperature",
                "Invalid operation: two-point control is not enabled.");
    }
}

double Flow1D::leftControlPointCoordinate() const
{
    if (m_twoPointControl) {
        if (m_zLeft != Undef) {
            return m_zLeft;
        } else {
            throw CanteraError("Flow1D::leftControlPointCoordinate",
                "Invalid operation: left control point location is not set");
        }
    } else {
        throw CanteraError("Flow1D::leftControlPointCoordinate",
                "Invalid operation: two-point control is not enabled.");
    }
}

void Flow1D::setLeftControlPointTemperature(double temperature)
{
    if (m_twoPointControl) {
        if (m_zLeft != Undef) {
            m_tLeft = temperature;
        } else {
            throw CanteraError("Flow1D::setLeftControlPointTemperature",
                "Invalid operation: left control point location is not set");
        }
    } else {
        throw CanteraError("Flow1D::setLeftControlPointTemperature",
                "Invalid operation: two-point control is not enabled.");
    }
}

void Flow1D::setLeftControlPointCoordinate(double z_left)
{
    if (m_twoPointControl) {
            m_zLeft = z_left;
    } else {
        throw CanteraError("Flow1D::setLeftControlPointCoordinate",
                "Invalid operation: two-point control is not enabled.");
    }
}

double Flow1D::rightControlPointTemperature() const
{
    if (m_twoPointControl) {
        if (m_zRight != Undef) {
            return m_tRight;
        } else {
            throw CanteraError("Flow1D::rightControlPointTemperature",
                "Invalid operation: right control point location is not set");
        }
    } else {
        throw CanteraError("Flow1D::rightControlPointTemperature",
                "Invalid operation: two-point control is not enabled.");
    }
}

double Flow1D::rightControlPointCoordinate() const
{
    if (m_twoPointControl) {
        if (m_zRight != Undef) {
            return m_zRight;
        } else {
            throw CanteraError("Flow1D::rightControlPointCoordinate",
                "Invalid operation: right control point location is not set");
        }
    } else {
        throw CanteraError("Flow1D::rightControlPointCoordinate",
                "Invalid operation: two-point control is not enabled.");
    }
}

void Flow1D::setRightControlPointTemperature(double temperature)
{
    if (m_twoPointControl) {
        if (m_zRight != Undef) {
            m_tRight = temperature;
        } else {
            throw CanteraError("Flow1D::setRightControlPointTemperature",
                "Invalid operation: right control point location is not set");
        }
    } else {
        throw CanteraError("Flow1D::setRightControlPointTemperature",
                "Invalid operation: two-point control is not enabled.");
    }
}

void Flow1D::setRightControlPointCoordinate(double z_right)
{
    if (m_twoPointControl) {
            m_zRight = z_right;
    } else {
        throw CanteraError("Flow1D::setRightControlPointCoordinate",
                "Invalid operation: two-point control is not enabled.");
    }
}

void Flow1D::enableTwoPointControl(bool twoPointControl)
{
    if (isStrained()) {
        m_twoPointControl = twoPointControl;
        // Prevent finding spurious solutions with negative velocity (outflow) at either
        // inlet.
        setBounds(c_offset_V, -1e-5, 1e20);
    } else {
        throw CanteraError("Flow1D::enableTwoPointControl",
            "Invalid operation: two-point control can only be used"
            "with axisymmetric flames.");
    }
}

} // namespace
