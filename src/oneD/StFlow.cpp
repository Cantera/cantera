//! @file StFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/StFlow.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"

using namespace std;

namespace Cantera
{

StFlow::StFlow(ThermoPhase* ph, size_t nsp, size_t points) :
    Domain1D(nsp+c_offset_Y, points),
    m_press(-1.0),
    m_nsp(nsp),
    m_thermo(0),
    m_kin(0),
    m_trans(0),
    m_epsilon_left(0.0),
    m_epsilon_right(0.0),
    m_do_soret(false),
    m_do_multicomponent(false),
    m_do_radiation(false),
    m_kExcessLeft(0),
    m_kExcessRight(0),
    m_zfixed(Undef),
    m_tfixed(-1.)
{
    if (ph->type() == "IdealGas") {
        m_thermo = static_cast<IdealGasPhase*>(ph);
    } else {
        throw CanteraError("StFlow::StFlow",
                           "Unsupported phase type: need 'IdealGasPhase'");
    }
    m_type = cFlowType;
    m_points = points;

    if (ph == 0) {
        return; // used to create a dummy object
    }

    size_t nsp2 = m_thermo->nSpecies();
    if (nsp2 != m_nsp) {
        m_nsp = nsp2;
        Domain1D::resize(m_nsp+c_offset_Y, points);
    }

    // make a local copy of the species molecular weight vector
    m_wt = m_thermo->molecularWeights();

    // the species mass fractions are the last components in the solution
    // vector, so the total number of components is the number of species
    // plus the offset of the first mass fraction.
    m_nv = c_offset_Y + m_nsp;

    // enable all species equations by default
    m_do_species.resize(m_nsp, true);

    // but turn off the energy equation at all points
    m_do_energy.resize(m_points,false);

    m_diff.resize(m_nsp*m_points);
    m_multidiff.resize(m_nsp*m_nsp*m_points);
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_ybar.resize(m_nsp);
    m_qdotRadiation.resize(m_points, 0.0);

    //-------------- default solution bounds --------------------
    setBounds(0, -1e20, 1e20); // no bounds on u
    setBounds(1, -1e20, 1e20); // V
    setBounds(2, 200.0, 2*m_thermo->maxTemp()); // temperature bounds
    setBounds(3, -1e20, 1e20); // lambda should be negative

    // mass fraction bounds
    for (size_t k = 0; k < m_nsp; k++) {
        setBounds(c_offset_Y+k, -1.0e-7, 1.0e5);
    }

    //-------------------- grid refinement -------------------------
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    m_refiner->setActive(c_offset_L, false);

    vector_fp gr;
    for (size_t ng = 0; ng < m_points; ng++) {
        gr.push_back(1.0*ng/m_points);
    }
    setupGrid(m_points, gr.data());

    // Find indices for radiating species
    AbsorptionSpeciesList = {"CO2","H2O","CO","CH4","OH","N2O","NH3",
        "H2O2","H2CO","NO2","NO","HCN","C2H2","C2H4","C2H6","C2N2",
        "C4H2","HCOOH","HNO3","O3"};
    //AbsorptionSpeciesList = {"H2O"};

    for(std::string spIter : AbsorptionSpeciesList) {
        AbsorptionSpeciesMap.insert({spIter, m_thermo->speciesIndex(spIter)});
    }
}

void StFlow::resize(size_t ncomponents, size_t points)
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
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_do_energy.resize(m_points,false);
    m_qdotRadiation.resize(m_points, 0.0);
    m_fixedtemp.resize(m_points);

    m_dz.resize(m_points-1);
    m_z.resize(m_points);
}

void StFlow::setupGrid(size_t n, const doublereal* z)
{
    resize(m_nv, n);

    m_z[0] = z[0];
    for (size_t j = 1; j < m_points; j++) {
        if (z[j] <= z[j-1]) {
            throw CanteraError("StFlow::setupGrid",
                               "grid points must be monotonically increasing");
        }
        m_z[j] = z[j];
        m_dz[j-1] = m_z[j] - m_z[j-1];
    }
}

void StFlow::resetBadValues(double* xg)
{
    double* x = xg + loc();
    for (size_t j = 0; j < m_points; j++) {
        double* Y = x + m_nv*j + c_offset_Y;
        m_thermo->setMassFractions(Y);
        m_thermo->getMassFractions(Y);
    }
}

void StFlow::setTransport(Transport& trans)
{
    m_trans = &trans;
    m_do_multicomponent = (m_trans->transportType() == "Multi");

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
}

void StFlow::_getInitialSoln(double* x)
{
    for (size_t j = 0; j < m_points; j++) {
        T(x,j) = m_thermo->temperature();
        m_thermo->getMassFractions(&Y(x, 0, j));
    }
}

void StFlow::setGas(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(T(x,j));
    const doublereal* yy = x + m_nv*j + c_offset_Y;
    m_thermo->setMassFractions_NoNorm(yy);
    m_thermo->setPressure(m_press);
}

void StFlow::setGasAtMidpoint(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
    const doublereal* yyj = x + m_nv*j + c_offset_Y;
    const doublereal* yyjp = x + m_nv*(j+1) + c_offset_Y;
    for (size_t k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
    }
    m_thermo->setMassFractions_NoNorm(m_ybar.data());
    m_thermo->setPressure(m_press);
}

void StFlow::_finalize(const doublereal* x)
{
    if (!m_do_multicomponent && m_do_soret) {
        throw CanteraError("StFlow::_finalize",
            "Thermal diffusion (the Soret effect) is enabled, and requires "
            "using a multicomponent transport model.");
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

    if (domainType() == cFreeFlow) {
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

void StFlow::eval(size_t jg, doublereal* xg,
                  doublereal* rg, integer* diagg, doublereal rdt)
{
    // if evaluating a Jacobian, and the global point is outside the domain of
    // influence for this domain, then skip evaluating the residual
    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;
    if (jg == npos) { // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    updateProperties(jg, x, jmin, jmax);
    evalResidual(x, rsd, diag, rdt, jmin, jmax);
}

void StFlow::updateProperties(size_t jg, double* x, size_t jmin, size_t jmax)
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

void StFlow::evalResidual(double* x, double* rsd, int* diag,
                          double rdt, size_t jmin, size_t jmax)
{
    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    // calculation of qdotRadiation

    // The simple radiation model used was established by Y. Liu and B. Rogg [Y.
    // Liu and B. Rogg, Modelling of thermally radiating diffusion flames with
    // detailed chemistry and transport, EUROTHERM Seminars, 17:114-127, 1991].
    // This model uses the optically thin limit and the gray-gas approximation
    // to simply calculate a volume specified heat flux out of the Planck
    // absorption coefficients, the boundary emissivities and the temperature.
    // The model considers only CO2 and H2O as radiating species. Polynomial
    // lines calculate the species Planck coefficients for H2O and CO2. The data
    // for the lines is taken from the RADCAL program [Grosshandler, W. L.,
    // RADCAL: A Narrow-Band Model for Radiation Calculations in a Combustion
    // Environment, NIST technical note 1402, 1993]. The coefficients for the
    // polynomials are taken from [http://www.sandia.gov/TNF/radiation.html].

    if (m_do_radiation) {
        // variable definitions for the Planck absorption coefficient and the
        // radiation calculation:
        doublereal k_P_ref = 1.0*OneAtm;
                             
        // Temperatures for Planck optical path length evaluation, K
        
        std::vector<double> TemperatureOPL = {200.0, 300.0, 400.0, 500.0, 600.0,
            700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0,
            1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0,
            2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0};
            
        const int OPL_table_size = TemperatureOPL.size();
        
        // Planck optical path length, m        
        std::vector<double> PlanckOPL_CO2 = {0.028819231725166233, 0.03851088876725997,
            0.03899126980906671, 0.032980228203593136, 0.030267348946179023, 0.030743052568319492,
            0.03349552664432657, 0.03815025367862556, 0.044695166534542856, 0.0533451125340091,
            0.06448815851162609, 0.07867018394357354, 0.09659983965382253, 0.11916430710589461,
            0.14745252150882412, 0.18278437065235484, 0.22674394556072777, 0.2812154755687562,
            0.34842470656989205, 0.43098095501152517, 0.5319261720083194, 0.6547800464563425,
            0.8035960044763693, 0.983014582040391, 1.1983144211620795, 1.4554774859486679,
            1.7612497326430492, 2.123189910047587, 2.549740894200347, 3.050301729813575,
            3.6352339058247307, 4.316015725482628, 5.1051815487411565, 6.016472443851096};
        std::vector<double> PlanckOPL_H2O = {0.007180563211295323, 0.019365261372768307,
            0.034867966945125076, 0.05230634537613443, 0.07276548153611018, 0.09699925744256371,
            0.12530329568293222, 0.15801661948972606, 0.19580389577973367, 0.23972350297492273,
            0.2912064962199706, 0.35202916571833975, 0.4243127637866092, 0.510551491576356,
            0.6136629518326753, 0.7370547168162068, 0.8847050123750548, 1.061254758533107,
            1.2721110105236084, 1.5235621242734583, 1.822906018849437, 2.178591373741365,
            2.6003745205523052, 3.099489796041002, 3.688832242804873, 4.383157080523857,
            5.1993025902907215, 6.156408948413962, 7.276171359454505, 8.583098608019881,
            10.104757796299838, 11.87209392358693, 13.919670784850705, 16.28600423253945};
        std::vector<double> PlanckOPL_CO = {30.32780456023482, 1.5025772744944488,
            0.48568677330878235, 0.31675988274374905, 0.2809201570439709, 0.28946895791957333,
            0.32239055951412515, 0.37417537060525935, 0.4437961211963581, 0.5320760473126184,
            0.6408278564965408, 0.7725511001710472, 0.9303212767894896, 1.1177046866345075,
            1.3387967943387578, 1.598138988644575, 1.9007086226166436, 2.252003287440127,
            2.6578585623777644, 3.1245488139145126, 3.658747887146305, 4.2674749977544435,
            4.958101084203171, 5.7383843347079315, 6.616293182716349, 7.600244912979516,
            8.698845819200764, 9.920965846773354, 11.276052909388223, 12.77335136618167,
            14.42279395291501, 16.234355340729095, 18.218430714508237, 20.3856192883815};
        std::vector<double> PlanckOPL_CH4 = {0.6939169829570981, 0.2276100632203682,
            0.18937636253724896, 0.19815582280865973, 0.2200025240118599, 0.2534843471195877,
            0.3046841472491959, 0.3823986600184083, 0.4987623999795097, 0.6713411345124516,
            0.9259289806598495, 1.3002850424277999, 1.849189501228091, 2.6564473100322914,
            3.8304102615717746, 5.5344722927178935, 7.994950093215335, 11.526582452228748,
            16.563118973976213, 23.69700078531631, 33.73017627335858, 47.73920940736953,
            67.1570036554854, 93.87559477159049, 130.37319123976783, 179.87209418262805,
            246.53023232205777, 335.6763859293513, 454.0906008652992, 610.3469544483086,
            815.2123010446278, 1082.1303102619831, 1427.7830218942922, 1872.756285795192};
        std::vector<double> PlanckOPL_OH = {0.009691578817706287, 0.030049601179558346,
            0.06765464326386915, 0.12708598866615728, 0.21115912062446549, 0.3200686268532878,
            0.4522763199038495, 0.6063045700105023, 0.7818952037702664, 0.980162180031524,
            1.2032232614189733, 1.453779407069445, 1.7348075847988176, 2.0493997120339555,
            2.4006745319830265, 2.7917401665539523, 3.2256754466853086, 3.7054691023958872,
            4.23395822768527, 4.813594141021041, 5.446035279732838, 6.131334910130047,
            6.866656642167146, 7.644231944896166, 8.44868549777347, 9.253875110774946,
            10.020328484366413, 10.694938827176772, 11.215426325696, 11.520789003674768,
            11.56668421121504, 11.33994851879643, 10.864608936023622, 10.19585170664347};
        std::vector<double> PlanckOPL_N2O = {0.05604186150530222, 0.04991781435183458,
            0.04001827500393808, 0.036150506950346904, 0.03710719640348422, 0.04180788677648103,
            0.05018622510260605, 0.06285657185631298, 0.08094410959148458, 0.10604926499576005,
            0.1402734389982863, 0.18627028306484095, 0.24730302890664432, 0.32730217346266893,
            0.430908689506439, 0.5635079200362725, 0.7312449464709093, 0.9410221482388516,
            1.2004777818389483, 1.5179552157303737, 1.9024578243125905, 2.363587182232489,
            2.911482794626652, 3.5567676931747707, 4.310436165387461, 5.18386757038981,
            6.188658945204697, 7.336641004525838, 8.639751162834212, 10.110065918784443,
            11.759665166149045, 13.600597157732668, 15.644893691074058, 17.9045229031322};
        std::vector<double> PlanckOPL_NH3 = {0.015853113133010365, 0.01926676385432445,
            0.026415824676239527, 0.03897633796851235, 0.05888157907336428, 0.08915105723552613,
            0.13355327539124195, 0.1960730549285638, 0.28048642863977286, 0.3904687215819763,
            0.5303141578042857, 0.7058788022739142, 0.925329453596419, 1.199577110754361,
            1.5425833020806663, 1.9716075385766434, 2.507579090345429, 3.175515361108895,
            4.005017660292852, 5.030732259316466, 6.292943371621803, 7.838191292724922,
            9.719463627047077, 11.997068455181829, 14.738974683444567, 18.02055752855426,
            21.925832857763382, 26.548007401632894, 31.987628299146124, 38.35527359161105,
            45.770040671919816, 54.35989812846817, 64.26225033961265, 75.62383376369114};
        std::vector<double> PlanckOPL_H2O2 = {0.0061240643487373726, 0.014358654955073104,
            0.02372207214832294, 0.03509853861756539, 0.05041287435495022, 0.07158406668799877,
            0.10079925820290134, 0.14072600135394597, 0.1946390235350878, 0.26650037055313247,
            0.3610028583328736, 0.48358138449031224, 0.6403937124726201, 0.8382784165969821,
            1.0846948121183866, 1.3876522056917315, 1.7556322223573986, 2.1975119501537628,
            2.722491030278221, 3.340011592913925, 4.059702755928071, 4.891311219712732,
            5.844656956185561, 6.929574259055882, 8.155886163115925, 9.533355173169655,
            11.071661643874863, 12.780377074549552, 14.668934187405775, 16.746625569074368,
            19.022581230831573, 21.50576621761172, 24.204937981243575, 27.128677559801766};
        std::vector<double> PlanckOPL_H2CO = {1.3851378337328886, 0.29966707757879113,
            0.14782184113680305, 0.11321165888878547, 0.10667714260745761, 0.11262176302472546,
            0.12740720489674168, 0.15069851692677394, 0.18357419050287804, 0.22805478242524305,
            0.2870105594836183, 0.3642229104940912, 0.464518116525414, 0.5939600794802649,
            0.7600869230903206, 0.9722007453399172, 1.2417019725072342, 1.5824881073722192,
            2.011408496739065, 2.5487901402831663, 3.2190405187317457, 4.051345642825396,
            5.080436220225772, 6.347502018651657, 7.901169358594316, 9.798657874500417,
            12.107045287048898, 14.904679352835634, 18.28279597617536, 22.347272505969286,
            27.220616304368882, 33.044168742639215, 39.9804690379451, 48.216066577742545};            
        std::vector<double> PlanckOPL_NO2 = {0.18967185355199417, 0.04428846239231481,
            0.028416106886409966, 0.028275414773864766, 0.03422377575827348, 0.04560041986270607,
            0.06381585795362235, 0.09150569222640514, 0.13253556724406976, 0.1922144503571911,
            0.27758651879037544, 0.39778260471688276, 0.564395391544332, 0.791908677828003,
            1.0981494235689724, 1.5047895936423843, 2.037869062204293, 2.728368349256181,
            3.6128356924408402, 4.7340233488511325, 6.141594179781417, 7.892863536806141,
            10.053585199881873, 12.698770218856987, 15.91363338637731, 19.79442325583844,
            24.449499989156376, 30.000309241474096, 36.58245963378199, 44.34696368197642,
            53.461298650892346, 64.1107214629735, 76.49956856745997, 90.85268289351858};
        std::vector<double> PlanckOPL_NO = {12.49908285125374, 1.3459019533042522,
            0.6011825839900713, 0.4748113234371179, 0.47742862256124935, 0.5370442220062525,
            0.6376193053372756, 0.7764870585324466, 0.9555723091275765, 1.1788559185177765,
            1.4514629725987682, 1.7793270399233352, 2.1689573586466553, 2.627448906362707,
            3.1623084968864976, 3.781539187660205, 4.49349909489309, 5.307064481132861,
            6.231329265582485, 7.275885774773834, 8.450611521315984, 9.765632286526694,
            11.231430006248385, 12.858500531681175, 14.657751390671171, 16.639987005055776,
            18.816075851493107, 21.19670280326081, 23.793015918675906, 26.61491214156301,
            29.673133893532462, 32.97754031606446, 36.53746836667566, 40.3622966223194};
        std::vector<double> PlanckOPL_HCN = {0.028635211450903347, 0.03757526481800404,
            0.061354128074583666, 0.099041934820521, 0.14986680872525904, 0.2112860122528754,
            0.2821606517643699, 0.3648298945235665, 0.4645287980339424, 0.5924143164459074,
            0.7542345476291972, 0.9604408381860212, 1.2237054355039148, 1.5594639696461987,
            1.986552912551171, 2.527919503080509, 3.211397414785762, 4.070563895597286,
            5.145712796528703, 6.484886092047459, 8.144985631727018, 10.193065855766827,
            12.707593159716264, 15.779896881489458, 19.51568964842916, 24.03665133695894,
            29.482042907652733, 36.01052705210453, 43.802207627836076, 53.06005801030469,
            64.01228918330375, 76.91425579492252, 92.05099992084168, 109.73826951057255};
        std::vector<double> PlanckOPL_C2H2 = {0.010929525364339873, 0.01392974120512336,
            0.02405443205882381, 0.04386508025067529, 0.07855266240755755, 0.1345672180286686,
            0.21996737264216254, 0.3461990874422214, 0.5305814978967278, 0.798986833061833,
            1.18893005089752, 1.7534946898618935, 2.5663306054757027, 3.727852933362788,
            5.372529800880729, 7.677251759511393, 10.870517832651087, 15.242285828331596,
            21.154028397137477, 29.048837769380295, 39.46135154206141, 53.02686745889668,
            70.48958296253277, 92.71027377239228, 120.6719646353374, 155.48532979799032,
            198.39202277505726, 250.76660329085342, 314.1168538786097, 390.087262270045,
            480.45331608083166, 587.1190933644476, 712.117974644975, 857.6066079444056};
        std::vector<double> PlanckOPL_C2H4 = {0.01654680957667988, 0.02492716585471849,
            0.0460983767688889, 0.09024057897411016, 0.17647152655342777, 0.33699482855742513,
            0.6253345063344479, 1.130036554810777, 1.9963077428931275, 3.4594563059971244,
            5.895932142906171, 9.90032557003133, 16.40006735361465, 26.82443744854234,
            43.35046571234703, 69.25728060994958, 109.43037257639607, 171.07357189711118,
            264.7028534520616, 405.52279353849264, 615.3126797518005, 924.9981677988787,
            1378.1207657380173, 2035.4777754072459, 2981.2965059567937, 4331.43337255344,
            6244.0321979790415, 8933.516689165475, 12688.676814602117, 17895.963748387712,
            25069.504472771623, 34888.98940246357, 48247.80480294268, 66314.77809464224};
        std::vector<double> PlanckOPL_C2H6 = {0.38479552555827473, 0.330335973067736,
            0.4412642688261878, 0.7234904395389986, 1.3047333754840624, 2.4582817867678797,
            4.728954047788753, 9.185141033669707, 17.899543021230603, 34.850795560924084,
            67.58893234939407, 130.26844542962587, 249.09122926085138, 471.9318246172849,
            885.130031885148, 1642.3336783887491, 3013.480869817436, 5466.6396914389925,
            9803.371388627636, 17379.457798083575, 30460.7202083096, 52790.33601327037,
            90480.336797945, 153407.43112853833, 257355.273339689, 427303.61580306507,
            702366.9359563864, 1143267.5642510268, 1843354.9491235327, 2944846.812649131,
            4662659.170473707, 7318745.5983123155, 11391692.329199707, 17587602.579628713};
        std::vector<double> PlanckOPL_C2N2 = {0.005068498559075421, 0.00794066986454325,
            0.011326690006467697, 0.0162428641418049, 0.023999391745476977, 0.03635127840260863,
            0.05584914808886189, 0.08625527503938223, 0.1330651839347766, 0.20418883262957285,
            0.3108436682130422, 0.4687112261455473, 0.699423764955656, 1.0324479058510616,
            1.507453218412811, 2.1772618202142606, 3.1114775402859705, 4.400946829376023,
            6.163150049802422, 8.54873395424021, 11.749306098379188, 16.006768024926618,
            21.624282462568466, 28.979328225842597, 38.538809411076464, 50.8768807730018,
            66.6955409114437, 86.84833705611524, 112.36801734231409, 144.49770280025078,
            184.7272696900555, 234.83366678173417, 296.92842179122425, 373.5097598037753};
        std::vector<double> PlanckOPL_C4H2 = {0.008654385959387364, 0.013221048546382841,
            0.023444997915175448, 0.04492768707054329, 0.09001749431190405, 0.1838640708259934,
            0.3765391487893722, 0.7654115015689522, 1.5353750432524051, 3.0294321701285036,
            5.869712810223213, 11.160509846374135, 20.82201029603419, 38.1284600107405,
            68.56051866929708, 121.13277007136199, 210.4288561712099, 359.67571007269487,
            605.3226740532338, 1003.7744131679052, 1641.1702181707305, 2647.42124377185,
            4216.140762151766, 6632.661018999894, 10313.01303955769, 15857.676154221808,
            24125.02967067159, 36330.98667841451, 54182.78820638309, 80057.91631401845,
            117240.68658474933, 170234.22688919335, 245167.68464599433, 350326.5016591533};
        std::vector<double> PlanckOPL_HCOOH = {0.008500553910125858, 0.009363497122291741,
            0.01283858517065757, 0.01895258140461621, 0.02911801074174551, 0.045793592077955574,
            0.07290068695530225, 0.11656063070671369, 0.18618867990418764, 0.2960655295828281,
            0.46756168130170694, 0.7322445505013475, 1.1361637959006956, 1.745696069866044,
            2.655423918456038, 3.9986624722081285, 5.9613689634328955, 8.800366773348818,
            12.86702911179098, 18.637779834813458, 26.753121652310554, 38.06717802336468,
            53.710058801743436, 75.16625089593722, 104.37193238145707, 143.8355434073111,
            196.78635835589836, 267.3562198557535, 360.8011265516046, 483.76938983787875,
            644.6270384585974, 853.8465812437338, 1124.47299892604, 1472.682051768596};
        std::vector<double> PlanckOPL_HNO3 = {0.006384561215332105, 0.007945818865473717,
            0.011369606874221595, 0.018185332096860546, 0.030917439657766183, 0.05402995712473861,
            0.0952035718822803, 0.16731021826266373, 0.29144824024627225, 0.501507552422749,
            0.8509050408969643, 1.4223515155012039, 2.341789738404136, 3.797990929638042,
            6.069741193008762, 9.563079645542468, 14.861694346418044, 22.794395726651224,
            34.52449765491417, 51.66720541762567, 76.44215125950534, 111.87016843316668,
            162.0251932450192, 232.3537440607607, 330.07832708074693, 464.7023518143656,
            648.6382771969264, 897.9856543786034, 1233.4873962465554, 1681.7015375997942,
            2276.4250569677397, 3060.4203470895163, 4087.4985253227364, 5425.027223873235};
        std::vector<double> PlanckOPL_O3 = {0.06090151715462474, 0.04094792168395479,
            0.0493930255271001, 0.06968082884526133, 0.10301858705107333, 0.1546948114590627,
            0.23334115457989626, 0.35145488255968, 0.5262795613218485, 0.7808920244469918,
            1.1454338063269254, 1.658468204302404, 2.368506909036788, 3.335619739243394,
            4.633232474776676, 6.350081219500141, 8.592252521288646, 11.485497395583318,
            15.17760366231008, 19.84106558507588, 25.67573909058048, 32.91200786023972,
            41.81367130298549, 52.68162087288378, 65.85701739068608, 81.72543958358092,
            100.7206904850768, 123.3288052772987, 150.09273327833398, 181.61685811709106,
            218.5722288990807, 261.7005632383693, 311.82119282475406, 369.8355675948166};        
            
        std::map<std::string, std::vector<double>> PlanckOPLMap;

        PlanckOPLMap.insert({"CO2", PlanckOPL_CO2});
        PlanckOPLMap.insert({"H2O", PlanckOPL_H2O});
        PlanckOPLMap.insert({"CO", PlanckOPL_CO});
        PlanckOPLMap.insert({"CH4", PlanckOPL_CH4});
        PlanckOPLMap.insert({"OH", PlanckOPL_OH});
        PlanckOPLMap.insert({"N2O", PlanckOPL_N2O});
        PlanckOPLMap.insert({"NH3", PlanckOPL_NH3});
        PlanckOPLMap.insert({"H2O2", PlanckOPL_H2O2});
        PlanckOPLMap.insert({"H2CO", PlanckOPL_H2CO});
        PlanckOPLMap.insert({"NO2", PlanckOPL_NO2});
        PlanckOPLMap.insert({"NO", PlanckOPL_NO});
        PlanckOPLMap.insert({"HCN", PlanckOPL_HCN});

        PlanckOPLMap.insert({"C2H2", PlanckOPL_C2H2});
        PlanckOPLMap.insert({"C2H4", PlanckOPL_C2H4});
        PlanckOPLMap.insert({"C2H6", PlanckOPL_C2H6});
        PlanckOPLMap.insert({"C2N2", PlanckOPL_C2N2});
        PlanckOPLMap.insert({"C4H2", PlanckOPL_C4H2});
        PlanckOPLMap.insert({"HCOOH", PlanckOPL_HCOOH});
        PlanckOPLMap.insert({"HNO3", PlanckOPL_HNO3});
        PlanckOPLMap.insert({"O3", PlanckOPL_O3});

        // natural logarithms of reversed optical path lengths
        // for linear interpolation in logarithm scale
        
        // calculation of the two boundary values
        double boundary_Rad_left = m_epsilon_left * StefanBoltz * pow(T(x, 0), 4);
        double boundary_Rad_right = m_epsilon_right * StefanBoltz * pow(T(x, m_points - 1), 4);

        double coef = 0.0;

        // loop over all grid points
        for (size_t j = jmin; j < jmax; j++) {
            // helping variable for the calculation
            double radiative_heat_loss = 0;
            // temperature table interval search  
            int T_index = 0;
            
            for (int k = 0; k < OPL_table_size; k++) {
                if (T(x, j) < TemperatureOPL[k]) {
                    if (T(x, j) < TemperatureOPL[0]) {
                        T_index = 0;//lower table limit
                    }
                    else {
                        T_index = k;
                    }
                    break;
                }
                else {
                    T_index=OPL_table_size-1;//upper table limit
                }
            }
            // calculation of the mean Planck absorption coefficient
            double k_P = 0.0;
            
            for(std::string sp_name : AbsorptionSpeciesList) {
                // absorption coefficient for specie
                if (AbsorptionSpeciesMap[sp_name] != npos) {
                    double k_P_specie = 0.0;
                    if ((T_index == 0) || (T_index == OPL_table_size-1)) {
                        coef=log(1.0/PlanckOPLMap[sp_name][T_index]);
                    }
                    else {
                        coef=log(1.0/PlanckOPLMap[sp_name][T_index-1])+
                        (log(1.0/PlanckOPLMap[sp_name][T_index])-log(1.0/PlanckOPLMap[sp_name][T_index-1]))*
                        (T(x, j)-TemperatureOPL[T_index-1])/(TemperatureOPL[T_index]-TemperatureOPL[T_index-1]);                         
                    }
                    k_P_specie = exp(coef);

                    k_P_specie /= k_P_ref;
                    k_P += m_press * X(x, AbsorptionSpeciesMap[sp_name], j) * k_P_specie;
                }            
            }
            
            // calculation of the radiative heat loss term
            radiative_heat_loss = 2 * k_P *(2 * StefanBoltz * pow(T(x, j), 4)
            - boundary_Rad_left - boundary_Rad_right);

            // set the radiative heat loss vector
            m_qdotRadiation[j] = radiative_heat_loss;
        }
    }

    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {
            // these may be modified by a boundary object

            // Continuity. This propagates information right-to-left, since
            // rho_u at point 0 is dependent on rho_u at point 1, but not on
            // mdot from the inlet.
            rsd[index(c_offset_U,0)] =
                -(rho_u(x,1) - rho_u(x,0))/m_dz[0]
                -(density(1)*V(x,1) + density(0)*V(x,0));

            // the inlet (or other) object connected to this one will modify
            // these equations by subtracting its values for V, T, and mdot. As
            // a result, these residual equations will force the solution
            // variables to the values for the boundary object
            rsd[index(c_offset_V,0)] = V(x,0);
            if (doEnergy(0)) {
                rsd[index(c_offset_T,0)] = T(x,0);
            } else {
                rsd[index(c_offset_T,0)] = T(x,0) - T_fixed(0);
            }
            rsd[index(c_offset_L,0)] = -rho_u(x,0);

            // The default boundary condition for species is zero flux. However,
            // the boundary object may modify this.
            double sum = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                sum += Y(x,k,0);
                rsd[index(c_offset_Y + k, 0)] =
                    -(m_flux(k,0) + rho_u(x,0)* Y(x,k,0));
            }
            rsd[index(c_offset_Y + leftExcessSpecies(), 0)] = 1.0 - sum;

            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, 0)] = x[index(c_offset_E, j)];
        } else if (j == m_points - 1) {
            evalRightBoundary(x, rsd, diag, rdt);
            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
        } else { // interior points
            evalContinuity(j, x, rsd, diag, rdt);
            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];

            //------------------------------------------------
            //    Radial momentum equation
            //
            //    \rho dV/dt + \rho u dV/dz + \rho V^2
            //       = d(\mu dV/dz)/dz - lambda
            //-------------------------------------------------
            rsd[index(c_offset_V,j)]
            = (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j)
               - m_rho[j]*V(x,j)*V(x,j))/m_rho[j]
              - rdt*(V(x,j) - V_prev(j));
            diag[index(c_offset_V, j)] = 1;

            //-------------------------------------------------
            //    Species equations
            //
            //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
            //   = M_k\omega_k
            //-------------------------------------------------
            getWdot(x,j);
            for (size_t k = 0; k < m_nsp; k++) {
                double convec = rho_u(x,j)*dYdz(x,k,j);
                double diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))
                                / (z(j+1) - z(j-1));
                rsd[index(c_offset_Y + k, j)]
                = (m_wt[k]*(wdot(k,j))
                   - convec - diffus)/m_rho[j]
                  - rdt*(Y(x,k,j) - Y_prev(k,j));
                diag[index(c_offset_Y + k, j)] = 1;
            }

            //-----------------------------------------------
            //    energy equation
            //
            //    \rho c_p dT/dt + \rho c_p u dT/dz
            //    = d(k dT/dz)/dz
            //      - sum_k(\omega_k h_k_ref)
            //      - sum_k(J_k c_p_k / M_k) dT/dz
            //-----------------------------------------------
            if (m_do_energy[j]) {
                setGas(x,j);

                // heat release term
                const vector_fp& h_RT = m_thermo->enthalpy_RT_ref();
                const vector_fp& cp_R = m_thermo->cp_R_ref();
                double sum = 0.0;
                double sum2 = 0.0;
                for (size_t k = 0; k < m_nsp; k++) {
                    double flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));
                    sum += wdot(k,j)*h_RT[k];
                    sum2 += flxk*cp_R[k]/m_wt[k];
                }
                sum *= GasConstant * T(x,j);
                double dtdzj = dTdz(x,j);
                sum2 *= GasConstant * dtdzj;

                rsd[index(c_offset_T, j)] = - m_cp[j]*rho_u(x,j)*dtdzj
                                            - divHeatFlux(x,j) - sum - sum2;
                rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);
                rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
                rsd[index(c_offset_T, j)] -= (m_qdotRadiation[j] / (m_rho[j] * m_cp[j]));
                diag[index(c_offset_T, j)] = 1;
            } else {
                // residual equations if the energy equation is disabled
                rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
                diag[index(c_offset_T, j)] = 0;
            }

            rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
            diag[index(c_offset_L, j)] = 0;
        }
    }
}

void StFlow::updateTransport(doublereal* x, size_t j0, size_t j1)
{
     if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            doublereal wtm = m_thermo->meanMolecularWeight();
            doublereal rho = m_thermo->density();
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
            m_trans->getMixDiffCoeffs(&m_diff[j*m_nsp]);
            m_tcon[j] = m_trans->thermalConductivity();
        }
    }
}

void StFlow::showSolution(const doublereal* x)
{
    writelog("    Pressure:  {:10.4g} Pa\n", m_press);

    Domain1D::showSolution(x);

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

void StFlow::updateDiffFluxes(const doublereal* x, size_t j0, size_t j1)
{
    if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            double dz = z(j+1) - z(j);
            for (size_t k = 0; k < m_nsp; k++) {
                doublereal sum = 0.0;
                for (size_t m = 0; m < m_nsp; m++) {
                    sum += m_wt[m] * m_multidiff[mindex(k,m,j)] * (X(x,m,j+1)-X(x,m,j));
                }
                m_flux(k,j) = sum * m_diff[k+j*m_nsp] / dz;
            }
        }
    } else {
        for (size_t j = j0; j < j1; j++) {
            double sum = 0.0;
            double wtm = m_wtm[j];
            double rho = density(j);
            double dz = z(j+1) - z(j);
            for (size_t k = 0; k < m_nsp; k++) {
                m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
                m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
                sum -= m_flux(k,j);
            }
            // correction flux to insure that \sum_k Y_k V_k = 0.
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

string StFlow::componentName(size_t n) const
{
    switch (n) {
    case 0:
        return "velocity";
    case 1:
        return "spread_rate";
    case 2:
        return "T";
    case 3:
        return "lambda";
    case 4:
        return "eField";
    default:
        if (n >= c_offset_Y && n < (c_offset_Y + m_nsp)) {
            return m_thermo->speciesName(n - c_offset_Y);
        } else {
            return "<unknown>";
        }
    }
}

size_t StFlow::componentIndex(const std::string& name) const
{
    if (name=="u") {
        warn_deprecated("StFlow::componentIndex",
                        "To be changed after Cantera 2.5. "
                        "Solution component 'u' renamed to 'velocity'");
        return 0;
    } else if (name=="velocity") {
        return 0;
    } else if (name=="V") {
        warn_deprecated("StFlow::componentIndex",
                        "To be changed after Cantera 2.5. "
                        "Solution component 'V' renamed to 'spread_rate'");
        return 1;
    } else if (name=="spread_rate") {
        return 1;
    } else if (name=="T") {
        return 2;
    } else if (name=="lambda") {
        return 3;
    } else if (name == "eField") {
        return 4;
    } else {
        for (size_t n=c_offset_Y; n<m_nsp+c_offset_Y; n++) {
            if (componentName(n)==name) {
                return n;
            }
        }
        throw CanteraError("StFlow1D::componentIndex",
                           "no component named " + name);
    }
}

void StFlow::restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    Domain1D::restore(dom, soln, loglevel);
    vector<string> ignored;
    size_t nsp = m_thermo->nSpecies();
    vector_int did_species(nsp, 0);

    vector<XML_Node*> str = dom.getChildren("string");
    for (size_t istr = 0; istr < str.size(); istr++) {
        const XML_Node& nd = *str[istr];
        writelog(nd["title"]+": "+nd.value()+"\n");
    }

    double pp = getFloat(dom, "pressure", "pressure");
    setPressure(pp);
    vector<XML_Node*> d = dom.child("grid_data").getChildren("floatArray");
    vector_fp x;
    size_t np = 0;
    bool readgrid = false, wrote_header = false;
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        if (nm == "z") {
            getFloatArray(fa,x,false);
            np = x.size();
            if (loglevel >= 2) {
                writelog("Grid contains {} points.\n", np);
            }
            readgrid = true;
            setupGrid(np, x.data());
        }
    }
    if (!readgrid) {
        throw CanteraError("StFlow::restore",
                           "domain contains no grid points.");
    }

    debuglog("Importing datasets:\n", loglevel >= 2);
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        getFloatArray(fa,x,false);
        if (nm == "u") {
            debuglog("axial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "axial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_U,j)] = x[j];
            }
        } else if (nm == "z") {
            ; // already read grid
        } else if (nm == "V") {
            debuglog("radial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "radial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_V,j)] = x[j];
            }
        } else if (nm == "T") {
            debuglog("temperature   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "temperature array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_T,j)] = x[j];
            }

            // For fixed-temperature simulations, use the imported temperature
            // profile by default.  If this is not desired, call
            // setFixedTempProfile *after* restoring the solution.
            vector_fp zz(np);
            for (size_t jj = 0; jj < np; jj++) {
                zz[jj] = (grid(jj) - zmin())/(zmax() - zmin());
            }
            setFixedTempProfile(zz, x);
        } else if (nm == "L") {
            debuglog("lambda   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("StFlow::restore",
                                   "lambda arary size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_L,j)] = x[j];
            }
        } else if (m_thermo->speciesIndex(nm) != npos) {
            debuglog(nm+"   ", loglevel >= 2);
            if (x.size() == np) {
                size_t k = m_thermo->speciesIndex(nm);
                did_species[k] = 1;
                for (size_t j = 0; j < np; j++) {
                    soln[index(k+c_offset_Y,j)] = x[j];
                }
            }
        } else {
            ignored.push_back(nm);
        }
    }

    if (loglevel >=2 && !ignored.empty()) {
        writelog("\n\n");
        writelog("Ignoring datasets:\n");
        size_t nn = ignored.size();
        for (size_t n = 0; n < nn; n++) {
            writelog(ignored[n]+"   ");
        }
    }

    if (loglevel >= 1) {
        for (size_t ks = 0; ks < nsp; ks++) {
            if (did_species[ks] == 0) {
                if (!wrote_header) {
                    writelog("Missing data for species:\n");
                    wrote_header = true;
                }
                writelog(m_thermo->speciesName(ks)+" ");
            }
        }
    }

    if (dom.hasChild("energy_enabled")) {
        getFloatArray(dom, x, false, "", "energy_enabled");
        if (x.size() == nPoints()) {
            for (size_t i = 0; i < x.size(); i++) {
                m_do_energy[i] = (x[i] != 0);
            }
        } else if (!x.empty()) {
            throw CanteraError("StFlow::restore", "energy_enabled is length {}"
                               "but should be length {}", x.size(), nPoints());
        }
    }

    if (dom.hasChild("species_enabled")) {
        getFloatArray(dom, x, false, "", "species_enabled");
        if (x.size() == m_nsp) {
            for (size_t i = 0; i < x.size(); i++) {
                m_do_species[i] = (x[i] != 0);
            }
        } else if (!x.empty()) {
            // This may occur when restoring from a mechanism with a different
            // number of species.
            if (loglevel > 0) {
                warn_user("StFlow::restore", "species_enabled is "
                    "length {} but should be length {}. Enabling all species "
                    "equations by default.", x.size(), m_nsp);
            }
            m_do_species.assign(m_nsp, true);
        }
    }

    if (dom.hasChild("refine_criteria")) {
        XML_Node& ref = dom.child("refine_criteria");
        refiner().setCriteria(getFloat(ref, "ratio"), getFloat(ref, "slope"),
                              getFloat(ref, "curve"), getFloat(ref, "prune"));
        refiner().setGridMin(getFloat(ref, "grid_min"));
    }

    if (domainType() == cFreeFlow) {
        getOptionalFloat(dom, "t_fixed", m_tfixed);
        getOptionalFloat(dom, "z_fixed", m_zfixed);
    }
}

XML_Node& StFlow::save(XML_Node& o, const doublereal* const sol)
{
    Array2D soln(m_nv, m_points, sol + loc());
    XML_Node& flow = Domain1D::save(o, sol);
    flow.addAttribute("type",flowType());

    if (m_desc != "") {
        addString(flow,"description",m_desc);
    }
    XML_Node& gv = flow.addChild("grid_data");
    addFloat(flow, "pressure", m_press, "Pa", "pressure");

    addFloatArray(gv,"z",m_z.size(), m_z.data(),
                  "m","length");
    vector_fp x(soln.nColumns());

    soln.getRow(c_offset_U, x.data());
    addFloatArray(gv,"u",x.size(),x.data(),"m/s","velocity");

    soln.getRow(c_offset_V, x.data());
    addFloatArray(gv,"V",x.size(),x.data(),"1/s","rate");

    soln.getRow(c_offset_T, x.data());
    addFloatArray(gv,"T",x.size(),x.data(),"K","temperature");

    soln.getRow(c_offset_L, x.data());
    addFloatArray(gv,"L",x.size(),x.data(),"N/m^4");

    for (size_t k = 0; k < m_nsp; k++) {
        soln.getRow(c_offset_Y+k, x.data());
        addFloatArray(gv,m_thermo->speciesName(k),
                      x.size(),x.data(),"","massFraction");
    }
    if (m_do_radiation) {
        addFloatArray(gv, "radiative_heat_loss", m_z.size(),
            m_qdotRadiation.data(), "W/m^3", "specificPower");
    }
    vector_fp values(nPoints());
    for (size_t i = 0; i < nPoints(); i++) {
        values[i] = m_do_energy[i];
    }
    addNamedFloatArray(flow, "energy_enabled", nPoints(), &values[0]);

    values.resize(m_nsp);
    for (size_t i = 0; i < m_nsp; i++) {
        values[i] = m_do_species[i];
    }
    addNamedFloatArray(flow, "species_enabled", m_nsp, &values[0]);

    XML_Node& ref = flow.addChild("refine_criteria");
    addFloat(ref, "ratio", refiner().maxRatio());
    addFloat(ref, "slope", refiner().maxDelta());
    addFloat(ref, "curve", refiner().maxSlope());
    addFloat(ref, "prune", refiner().prune());
    addFloat(ref, "grid_min", refiner().gridMin());
    if (m_zfixed != Undef) {
        addFloat(flow, "z_fixed", m_zfixed, "m");
        addFloat(flow, "t_fixed", m_tfixed, "K");
    }
    return flow;
}

void StFlow::solveEnergyEqn(size_t j)
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

void StFlow::setBoundaryEmissivities(doublereal e_left, doublereal e_right)
{
    if (e_left < 0 || e_left > 1) {
        throw CanteraError("StFlow::setBoundaryEmissivities",
            "The left boundary emissivity must be between 0.0 and 1.0!");
    } else if (e_right < 0 || e_right > 1) {
        throw CanteraError("StFlow::setBoundaryEmissivities",
            "The right boundary emissivity must be between 0.0 and 1.0!");
    } else {
        m_epsilon_left = e_left;
        m_epsilon_right = e_right;
    }
}

void StFlow::fixTemperature(size_t j)
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

void StFlow::evalRightBoundary(double* x, double* rsd, int* diag, double rdt)
{
    size_t j = m_points - 1;

    // the boundary object connected to the right of this one may modify or
    // replace these equations. The default boundary conditions are zero u, V,
    // and T, and zero diffusive flux for all species.

    rsd[index(c_offset_V,j)] = V(x,j);
    doublereal sum = 0.0;
    rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
    diag[index(c_offset_L, j)] = 0;
    for (size_t k = 0; k < m_nsp; k++) {
        sum += Y(x,k,j);
        rsd[index(k+c_offset_Y,j)] = m_flux(k,j-1) + rho_u(x,j)*Y(x,k,j);
    }
    rsd[index(c_offset_Y + rightExcessSpecies(), j)] = 1.0 - sum;
    diag[index(c_offset_Y + rightExcessSpecies(), j)] = 0;
    if (domainType() == cAxisymmetricStagnationFlow) {
        rsd[index(c_offset_U,j)] = rho_u(x,j);
        if (m_do_energy[j]) {
            rsd[index(c_offset_T,j)] = T(x,j);
        } else {
            rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
        }
    } else if (domainType() == cFreeFlow) {
        rsd[index(c_offset_U,j)] = rho_u(x,j) - rho_u(x,j-1);
        rsd[index(c_offset_T,j)] = T(x,j) - T(x,j-1);
    }
}

void StFlow::evalContinuity(size_t j, double* x, double* rsd, int* diag, double rdt)
{
    //algebraic constraint
    diag[index(c_offset_U, j)] = 0;
    //----------------------------------------------
    //    Continuity equation
    //
    //    d(\rho u)/dz + 2\rho V = 0
    //----------------------------------------------
    if (domainType() == cAxisymmetricStagnationFlow) {
        // Note that this propagates the mass flow rate information to the left
        // (j+1 -> j) from the value specified at the right boundary. The
        // lambda information propagates in the opposite direction.
        rsd[index(c_offset_U,j)] =
            -(rho_u(x,j+1) - rho_u(x,j))/m_dz[j]
            -(density(j+1)*V(x,j+1) + density(j)*V(x,j));
    } else if (domainType() == cFreeFlow) {
        if (grid(j) > m_zfixed) {
            rsd[index(c_offset_U,j)] =
                - (rho_u(x,j) - rho_u(x,j-1))/m_dz[j-1]
                - (density(j-1)*V(x,j-1) + density(j)*V(x,j));
        } else if (grid(j) == m_zfixed) {
            if (m_do_energy[j]) {
                rsd[index(c_offset_U,j)] = (T(x,j) - m_tfixed);
            } else {
                rsd[index(c_offset_U,j)] = (rho_u(x,j)
                                            - m_rho[0]*0.3);
            }
        } else if (grid(j) < m_zfixed) {
            rsd[index(c_offset_U,j)] =
                - (rho_u(x,j+1) - rho_u(x,j))/m_dz[j]
                - (density(j+1)*V(x,j+1) + density(j)*V(x,j));
        }
    }
}

} // namespace
