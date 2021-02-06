//! @file Sub.cpp
/*
 * The Substance class
 * D. Goodwin, Caltech Nov. 1996
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/tpx/Sub.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

using std::string;
using namespace Cantera;

namespace {
// these correspond to ordering withing propertyFlag::type
std::string propertySymbols[] = {"H", "S", "U", "V", "P", "T"};
}

namespace tpx
{
Substance::Substance() :
    T(Undef),
    Rho(Undef),
    Tslast(Undef),
    Rhf(Undef),
    Rhv(Undef),
    Pst(Undef),
    m_energy_offset(0.0),
    m_entropy_offset(0.0),
    kbr(0)
{
}

void Substance::setStdState(double h0, double s0, double t0, double p0)
{
    Set(PropertyPair::TP, t0, p0);
    double hh = h();
    double ss = s();
    double hoff = h0 - hh;
    double soff = s0 - ss;
    m_entropy_offset += soff;
    m_energy_offset += hoff;
}

double Substance::P()
{
    return TwoPhase() ? Ps() : Pp();
}

const double DeltaT = 0.000001;

double Substance::cv()
{
    if (TwoPhase(true)) {
        // While cv can be calculated for the two-phase region (the state can
        // always be continuously varied along an isochor on a T-v diagram),
        // this calculation is currently not implemented
        return std::numeric_limits<double>::quiet_NaN();
    }

    double Tsave = T, dt = 1.e-4 * T;
    double x0 = x();
    double T1 = std::max(Tmin(), Tsave - dt);
    double T2 = std::min(Tmax(), Tsave + dt);

    set_T(T1);
    double x1 = x();
    if ((x0 == 1.0 || x0 == 0.0) && x1 != x0) {
        // If the initial state was pure liquid or pure vapor, and the state at
        // T-dT is not, just take a one-sided difference
        T1 = Tsave;
        set_T(T1);
    }
    double s1 = s();

    set_T(T2);
    double x2 = x();
    if ((x0 == 1.0 || x0 == 0.0) && x2 != x0) {
        // If the initial state was pure liquid or pure vapor, and the state at
        // T+dT is not, just take a one-sided difference
        T2 = Tsave;
        set_T(T2);
    }
    double s2 = s();

    set_T(Tsave);
    return T*(s2 - s1)/(T2-T1);
}

double Substance::cp()
{
    if (TwoPhase(true)) {
        // In the two-phase region, cp is infinite
        return std::numeric_limits<double>::infinity();
    }

    double Tsave = T, dt = 1.e-4 * T;
    double RhoSave = Rho;
    double T1, T2, s1, s2;
    double p0 = P();

    if (Rho >= Rhf) {
        // initial state is pure liquid
        T1 = std::max(Tmin(), Tsave - dt);
        Set(PropertyPair::TP, T1, p0);
        s1 = s();
        try {
            T2 = std::min(Tsat(p0), Tsave + dt);
        } catch (CanteraError&) {
            // Tsat does not exist beyond two-phase region
            T2 = Tsave + dt;
        }
        if (T2 < Tsave + dt) {
            Set(PropertyPair::TX, T2, 0.);
        } else {
            Set(PropertyPair::TP, T2, p0);
        }
        s2 = s();
    } else {
        // initial state is pure vapor
        try {
            T1 = std::max(Tsat(p0), Tsave - dt);
        } catch (CanteraError&) {
            // Tsat does not exist beyond two-phase region
            T1 = Tsave - dt;
        }
        if (T1 > Tsave - dt) {
            Set(PropertyPair::TX, T1, 1.);
        } else {
            Set(PropertyPair::TP, T1, p0);
        }
        s1 = s();
        T2 = std::min(Tmax(), Tsave + dt);
        Set(PropertyPair::TP, T2, p0);
        s2 = s();
    }

    Set(PropertyPair::TV, Tsave, 1.0 / RhoSave);
    return T * (s2 - s1) / (T2 - T1);
}

double Substance::thermalExpansionCoeff()
{
    if (TwoPhase(true)) {
        // In the two-phase region, the thermal expansion coefficient is
        // infinite
        return std::numeric_limits<double>::infinity();
    }

    double Tsave = T, dt = 1.e-4 * T;
    double RhoSave = Rho;
    double T1, T2, v1, v2;
    double p0 = P();

    if (Rho >= Rhf) {
        // initial state is pure liquid
        T1 = std::max(Tmin(), Tsave - dt);
        Set(PropertyPair::TP, T1, p0);
        v1 = v();
        try {
            T2 = std::min(Tsat(p0), Tsave + dt);
        } catch (CanteraError&) {
            // Tsat does not exist beyond two-phase region
            T2 = Tsave + dt;
        }
        if (T2 < Tsave + dt) {
            Set(PropertyPair::TX, T2, 0.);
        } else {
            Set(PropertyPair::TP, T2, p0);
        }
        v2 = v();
    } else {
        // initial state is pure vapor
        try {
            T1 = std::max(Tsat(p0), Tsave - dt);
        } catch (CanteraError&) {
            // Tsat does not exist beyond two-phase region
            T1 = Tsave - dt;
        }
        if (T1 > Tsave - dt) {
            Set(PropertyPair::TX, T1, 1.);
        } else {
            Set(PropertyPair::TP, T1, p0);
        }
        v1 = v();
        T2 = std::min(Tmax(), Tsave + dt);
        Set(PropertyPair::TP, T2, p0);
        v2 = v();
    }

    Set(PropertyPair::TV, Tsave, 1.0 / RhoSave);
    return 2.0 * (v2 - v1) / ((v2 + v1) * (T2 - T1));
}

double Substance::isothermalCompressibility()
{
    if (TwoPhase(true)) {
        // In the two-phase region, the isothermal compressibility is infinite
        return std::numeric_limits<double>::infinity();
    }

    double Psave = P(), dp = 1.e-4 * Psave;
    double RhoSave = Rho;
    double P1, P2, v1, v2;
    double v0 = v();

    if (Rho >= Rhf) {
        // initial state is pure liquid
        P1 = Psave - dp;
        if (T > Tmin() && T <= Tcrit()) {
            P1 = std::max(Ps(), P1);
        }
        if (P1 > Psave - dp) {
            Set(PropertyPair::PX, P1, 0.);
        } else {
            Set(PropertyPair::TP, T, P1);
        }
        v1 = v();
        P2 = Psave + dp;
        Set(PropertyPair::TP, T, P2);
        v2 = v();
    } else {
        // initial state is pure vapor
        P1 = std::max(1.e-7, Psave - dp);
        Set(PropertyPair::TP, T, P1);
        v1 = v();
        P2 = Psave + dp;
        if (T > Tmin() && T <= Tcrit()) {
            P2 = std::min(Ps(), P2);
        }
        if (P2 < Psave + dp) {
            Set(PropertyPair::PX, P2, 1.);
        } else {
            Set(PropertyPair::TP, T, P2);
        }
        v2 = v();
    }

    Set(PropertyPair::TV, T, 1.0 / RhoSave);
    return -(v2 - v1) / (v0 * (P2 - P1));
}

double Substance::dPsdT()
{
    double tsave = T;
    double ps1 = Ps();
    double dpdt;
    if (T + DeltaT < Tcrit()) {
        T += DeltaT;
        dpdt = (Ps() - ps1)/DeltaT;
    } else {
        T -= DeltaT;
        // Ps() will fail for illegal temperature
        dpdt = (ps1 - Ps())/DeltaT;
    }
    T = tsave;
    return dpdt;
}

int Substance::TwoPhase(bool strict)
{
    if (T >= Tcrit()) {
        return 0;
    }
    update_sat();
    if (strict) {
        return ((Rho > Rhv) && (Rho < Rhf) ? 1 : 0);
    }
    return ((Rho >= Rhv) && (Rho <= Rhf) ? 1 : 0);
}

double Substance::x()
{
    if (T >= Tcrit()) {
        return (1.0/Rho < Vcrit() ? 0.0 : 1.0);
    } else {
        update_sat();
        if (Rho <= Rhv) {
            return 1.0;
        } else if (Rho >= Rhf) {
            return 0.0;
        } else {
            double vv = 1.0/Rhv;
            double vl = 1.0/Rhf;
            return (1.0/Rho - vl)/(vv - vl);
        }
    }
}

double Substance::Tsat(double p)
{
    double Tsave = T;
    if (p <= 0. || p > Pcrit()) {
        // @TODO: this check should also fail below the triple point pressure
        // (calculated as Ps() for Tmin() as it is not available as a numerical
        // parameter). A bug causing low-temperature TP setters to fail for
        // Heptane (GitHub issue #605) prevents this at the moment.
        // Without this check, a calculation attempt fails with the less
        // meaningful error "No convergence" below.
        throw CanteraError("Substance::Tsat",
                           "Illegal pressure value: {}", p);
    }
    T = Tsave;
    Tsave = T;

    int LoopCount = 0;
    double tol = 1.e-6*p;
    if (T < Tmin() || T > Tcrit()) {
        T = 0.5*(Tcrit() - Tmin());
    }
    double dp = 10*tol;
    while (fabs(dp) > tol) {
        T = std::min(T, Tcrit());
        T = std::max(T, Tmin());
        dp = p - Ps();
        double dt = dp/dPsdT();
        double dta = fabs(dt);
        double dtm = 0.1*T;
        if (dta > dtm) {
            dt = dt*dtm/dta;
        }
        T += dt;
        LoopCount++;
        if (LoopCount > 100) {
            T = Tsave;
            throw CanteraError("Substance::Tsat", "No convergence: p = {}", p);
        }
    }
    double tsat = T;
    T = Tsave;
    return tsat;
}

// property tolerances
static const double TolAbsH = 1.e-4; // J/kg
static const double TolAbsU = 1.e-4;
static const double TolAbsS = 1.e-7;
static const double TolAbsP = 0.0; // Pa, this is supposed to be zero
static const double TolAbsV = 1.e-8;
static const double TolAbsT = 1.e-3;
static const double TolRel = 1.e-8;

void Substance::Set(PropertyPair::type XY, double x0, double y0)
{
    double temp;

    /* if inverted (PT) switch order and change sign of XY (TP = -PT) */
    if (XY < 0) {
        std::swap(x0, y0);
        XY = static_cast<PropertyPair::type>(-XY);
    }

    switch (XY) {
    case PropertyPair::TV:
        set_T(x0);
        set_v(y0);
        break;
    case PropertyPair::HP:
        if (Lever(Pgiven, y0, x0, propertyFlag::H)) {
            return;
        }
        set_xy(propertyFlag::H, propertyFlag::P,
               x0, y0, TolAbsH, TolAbsP, TolRel, TolRel);
        break;
    case PropertyPair::SP:
        if (Lever(Pgiven, y0, x0, propertyFlag::S)) {
            return;
        }
        set_xy(propertyFlag::S, propertyFlag::P,
               x0, y0, TolAbsS, TolAbsP, TolRel, TolRel);
        break;
    case PropertyPair::PV:
        if (Lever(Pgiven, x0, y0, propertyFlag::V)) {
            return;
        }
        set_xy(propertyFlag::P, propertyFlag::V,
               x0, y0, TolAbsP, TolAbsV, TolRel, TolRel);
        break;
    case PropertyPair::TP:
        set_T(x0);
        if (x0 < Tcrit()) {
            if (fabs(y0 - Ps()) / y0 < TolRel) {
                throw CanteraError("Substance::Set",
                                   "Saturated mixture detected: use vapor "
                                   "fraction to specify state instead");
            } else if (y0 < Ps()) {
                Set(PropertyPair::TX, x0, 1.0);
            } else {
                Set(PropertyPair::TX, x0, 0.0);
            }
        }
        set_xy(propertyFlag::T, propertyFlag::P,
               x0, y0, TolAbsT, TolAbsP, TolRel, TolRel);
        break;
    case PropertyPair::UV:
        set_xy(propertyFlag::U, propertyFlag::V,
               x0, y0, TolAbsU, TolAbsV, TolRel, TolRel);
        break;
    case PropertyPair::ST:
        if (Lever(Tgiven, y0, x0, propertyFlag::S)) {
            return;
        }
        set_xy(propertyFlag::S, propertyFlag::T,
               x0, y0, TolAbsS, TolAbsT, TolRel, TolRel);
        break;
    case PropertyPair::SV:
        set_xy(propertyFlag::S, propertyFlag::V,
               x0, y0, TolAbsS, TolAbsV, TolRel, TolRel);
        break;
    case PropertyPair::UP:
        if (Lever(Pgiven, y0, x0, propertyFlag::U)) {
            return;
        }
        set_xy(propertyFlag::U, propertyFlag::P,
               x0, y0, TolAbsU, TolAbsP, TolRel, TolRel);
        break;
    case PropertyPair::VH:
        set_xy(propertyFlag::V, propertyFlag::H,
               x0, y0, TolAbsV, TolAbsH, TolRel, TolRel);
        break;
    case PropertyPair::TH:
        set_xy(propertyFlag::T, propertyFlag::H,
               x0, y0, TolAbsT, TolAbsH, TolRel, TolRel);
        break;
    case PropertyPair::SH:
        set_xy(propertyFlag::S, propertyFlag::H,
               x0, y0, TolAbsS, TolAbsH, TolRel, TolRel);
        break;
    case PropertyPair::PX:
        temp = Tmin();
        set_T(temp);
        if (y0 > 1.0 || y0 < 0.0) {
            throw CanteraError("Substance::Set",
                               "Invalid vapor fraction, {}", y0);
        }
        if (x0 < Ps() || x0 > Pcrit()) {
            throw CanteraError("Substance::Set",
                               "Illegal pressure value: {} (supercritical "
                               "or below triple point)", x0);
        }
        temp = Tsat(x0);
        set_T(temp);
        update_sat();
        Rho = 1.0/((1.0 - y0)/Rhf + y0/Rhv);
        break;
    case PropertyPair::TX:
        if (y0 > 1.0 || y0 < 0.0) {
            throw CanteraError("Substance::Set",
                               "Invalid vapor fraction, {}", y0);
        }
        if (x0 < Tmin() || x0 > Tcrit()) {
            throw CanteraError("Substance::Set",
                               "Illegal temperature value: {} "
                               "(supercritical or below triple point)", x0);
        }
        set_T(x0);
        update_sat();
        Rho = 1.0/((1.0 - y0)/Rhf + y0/Rhv);
        break;
    default:
        throw CanteraError("Substance::Set", "Invalid input.");
    }
}

//------------------ Protected and Private Functions -------------------

void Substance::set_Rho(double r0)
{
    if (r0 > 0.0) {
        Rho = r0;
    } else {
        throw CanteraError("Substance::set_Rho", "Invalid density: {}", r0);
    }
}

void Substance::set_T(double t0)
{
    if ((t0 >= Tmin()) && (t0 <= Tmax())) {
        T = t0;
    } else {
        throw CanteraError("Substance::set_T", "Illegal temperature: {}", t0);
    }
}

void Substance::set_v(double v0)
{
    if (v0 > 0) {
        Rho = 1.0/v0;
    } else {
        throw CanteraError("Substance::set_v",
                           "Negative specific volume: {}", v0);
    }
}

double Substance::Ps()
{
    if (T < Tmin() || T > Tcrit()) {
        throw CanteraError("Substance::Ps",
                           "Illegal temperature value: {}", T);
    }
    update_sat();
    return Pst;
}

void Substance::update_sat()
{
    if (T != Tslast) {
        double Rho_save = Rho;
        double pp = Psat(); // illegal temperatures are caught by Psat()
        double lps = log(pp);
        // trial value = Psat from correlation
        int i;
        for (i = 0; i<20; i++) {
            if (i==0) {
                Rho = ldens(); // trial value = liquid density
            } else {
                Rho = Rhf;
            }
            set_TPp(T,pp);
            Rhf = Rho; // sat liquid density

            double gf = hp() - T*sp();
            if (i==0) {
                Rho = pp*MolWt()/(GasConstant*T); // trial value = ideal gas
            } else {
                Rho = Rhv;
            }
            set_TPp(T,pp);

            Rhv = Rho; // sat vapor density
            double gv = hp() - T*sp();
            double dg = gv - gf;
            if (Rhv > Rhf) {
                std::swap(Rhv, Rhf);
                dg = - dg;
            }

            if (fabs(dg) < 0.001) {
                break;
            }
            double dp = dg/(1.0/Rhv - 1.0/Rhf);
            double psold = pp;
            if (fabs(dp) > pp) {
                lps -= dg/(pp*(1.0/Rhv - 1.0/Rhf));
                pp = exp(lps);
            } else {
                pp -= dp;
                lps = log(pp);
            }
            if (pp > Pcrit()) {
                pp = psold + 0.5*(Pcrit() - psold);
                lps = log(pp);
            } else if (pp < 0.0) {
                pp = psold/2.0;
                lps = log(pp);
            }
        }

        if (i >= 20) {
            throw CanteraError("Substance::update_sat", "No convergence");
        } else {
            Pst = pp;
            Tslast = T;
        }
        Rho = Rho_save;
    }
}

double Substance::vprop(propertyFlag::type ijob)
{
    switch (ijob) {
    case propertyFlag::H:
        return hp();
    case propertyFlag::S:
        return sp();
    case propertyFlag::U:
        return up();
    case propertyFlag::V:
        return vp();
    case propertyFlag::P:
        return Pp();
    default:
        throw CanteraError("Substance::vprop", "Invalid job index");
    }
}

int Substance::Lever(int itp, double sat, double val, propertyFlag::type ifunc)
{
    double psat;
    double Tsave = T;
    double Rhosave = Rho;
    if (itp == Tgiven) {
        if (sat >= Tcrit()) {
            return 0;
        }
        T = sat;
        psat = Ps();
    } else if (itp == Pgiven) {
        if (sat >= Pcrit()) {
            return 0;
        }
        psat = sat;
        try {
            T = Tsat(psat);
        } catch (CanteraError&) {
            // Failure to converge here is not an error
            T = Tsave;
            Rho = Rhosave;
            return 0;
        }
    } else {
        throw CanteraError("Substance::Lever", "General error");
    }
    Set(PropertyPair::TX, T, 1.0);
    double Valg = vprop(ifunc);
    Set(PropertyPair::TX, T, 0.0);
    double Valf = vprop(ifunc);
    if (val >= Valf && val <= Valg) {
        double xx = (val - Valf)/(Valg - Valf);
        double vv = (1.0 - xx)/Rhf + xx/Rhv;
        set_v(vv);
        return 1;
    } else {
        T = Tsave;
        Rho = Rhosave;
        return 0;
    }
}

void Substance::set_xy(propertyFlag::type ifx, propertyFlag::type ify,
                       double X, double Y,
                       double atx, double aty,
                       double rtx, double rty)
{
    double v_here, t_here;
    double dvs1 = 2.0*Vcrit();
    double dvs2 = 0.7*Vcrit();
    int LoopCount = 0;

    double v_save = 1.0/Rho;
    double t_save = T;

    if ((T == Undef) && (Rho == Undef)) {
        // new object, try to pick a "reasonable" starting point
        Set(PropertyPair::TV,Tcrit()*1.1,Vcrit()*1.1);
        t_here = T;
        v_here = 1.0/Rho;
    } else if (Rho == Undef) {
        // new object, try to pick a "reasonable" starting point
        Set(PropertyPair::TV,T,Vcrit()*1.1);
        t_here = T;
        v_here = 1.0/Rho;
    } else {
        v_here = v_save;
        t_here = t_save;
    }

    double Xa = fabs(X);
    double Ya = fabs(Y);
    while (true) {
        double x_here = prop(ifx);
        double y_here = prop(ify);
        double err_x = fabs(X - x_here);
        double err_y = fabs(Y - y_here);

        if ((err_x < atx + rtx*Xa) && (err_y < aty + rty*Ya)) {
            break;
        }

        /* perturb t */
        double dt = 0.001*t_here;
        if (t_here + dt > Tmax()) {
            dt *= -1.0;
        }

        /* perturb v */
        double dv = 0.001*v_here;
        if (v_here <= Vcrit()) {
            dv *= -1.0;
        }

        /* derivatives with respect to T */
        Set(PropertyPair::TV, t_here + dt, v_here);
        double dxdt = (prop(ifx) - x_here)/dt;
        double dydt = (prop(ify) - y_here)/dt;

        /* derivatives with respect to v */
        Set(PropertyPair::TV, t_here, v_here + dv);
        double dxdv = (prop(ifx) - x_here)/dv;
        double dydv = (prop(ify) - y_here)/dv;

        double det = dxdt*dydv - dydt*dxdv;
        dt = ((X - x_here)*dydv - (Y - y_here)*dxdv)/det;
        dv = ((Y - y_here)*dxdt - (X - x_here)*dydt)/det;

        double dvm = 0.2*v_here;
        if (v_here < dvs1) {
            dvm *= 0.5;
        }
        if (v_here < dvs2) {
            dvm *= 0.5;
        }
        double dtm = 0.1*t_here;
        double dva = fabs(dv);
        double dta = fabs(dt);
        if (dva > dvm) {
            dv *= dvm/dva;
        }
        if (dta > dtm) {
            dt *= dtm/dta;
        }
        v_here += dv;
        t_here += dt;
        t_here = clip(t_here, Tmin(), Tmax());
        if (v_here <= 0.0) {
            v_here = 0.0001;
        }
        Set(PropertyPair::TV, t_here, v_here);
        LoopCount++;
        if (LoopCount > 200) {
            std::string msg = fmt::format("No convergence. {} = {}, {} = {}",
                propertySymbols[ifx], X, propertySymbols[ify], Y);
            if (t_here == Tmin()) {
                msg += fmt::format("\nAt temperature limit (Tmin = {})", Tmin());
            } else if (t_here == Tmax()) {
                msg += fmt::format("\nAt temperature limit (Tmax = {})", Tmax());
            }
            throw CanteraError("Substance::set_xy", msg);
        }
    }
}

double Substance::prop(propertyFlag::type ijob)
{
    if (ijob == propertyFlag::P) {
        return P();
    }
    if (ijob == propertyFlag::T) {
        return T;
    }
    double xx = x();
    if ((xx > 0.0) && (xx < 1.0)) {
        double Rho_save = Rho;
        Rho = Rhv;
        double vp = vprop(ijob);
        Rho = Rhf;
        double lp = vprop(ijob);
        double pp = (1.0 - xx)*lp + xx*vp;
        Rho = Rho_save;
        return pp;
    } else {
        return vprop(ijob);
    }
}

static const double ErrP = 1.e-7;
static const double Big = 1.e30;

void Substance::BracketSlope(double Pressure)
{
    if (kbr == 0) {
        dv = (v_here < Vcrit() ? -0.05*v_here : 0.2*v_here);
        if (Vmin > 0.0) {
            dv = 0.2*v_here;
        }
        if (Vmax < Big) {
            dv = -0.05*v_here;
        }
    } else {
        double dpdv = (Pmax - Pmin)/(Vmax - Vmin);
        v_here = Vmax;
        P_here = Pmax;
        dv = dvbf*(Pressure - P_here)/dpdv;
        dvbf = 0.5*dvbf;
    }
}

void Substance::set_TPp(double Temp, double Pressure)
{
    kbr = 0;
    dvbf = 1.0;
    Vmin = 0.0;
    Vmax = Big;
    Pmin = Big;
    Pmax = 0.0;
    double dvs1 = 2.0*Vcrit();
    double dvs2 = 0.7*Vcrit();
    int LoopCount = 0;

    double v_save = 1.0/Rho;
    T = Temp;
    v_here = vp();

    // loop
    while (P_here = Pp(),
            fabs(Pressure - P_here) >= ErrP* Pressure || LoopCount == 0) {
        if (P_here < 0.0) {
            BracketSlope(Pressure);
        } else {
            dv = 0.001*v_here;
            if (v_here <= Vcrit()) {
                dv *= -1.0;
            }
            Set(PropertyPair::TV, Temp, v_here+dv);
            double dpdv = (Pp() - P_here)/dv;
            if (dpdv > 0.0) {
                BracketSlope(Pressure);
            } else {
                if ((P_here > Pressure) && (v_here > Vmin)) {
                    Vmin = v_here;
                } else if ((P_here < Pressure) && (v_here < Vmax)) {
                    Vmax = v_here;
                }
                if (v_here == Vmin) {
                    Pmin = P_here;
                }
                if (v_here == Vmax) {
                    Pmax = P_here;
                }
                if (Vmin >= Vmax) {
                    throw CanteraError("Substance::set_TPp", "Vmin >= Vmax");
                } else if ((Vmin > 0.0) && (Vmax < Big)) {
                    kbr = 1;
                }
                dvbf = 1.0;
                if (dpdv == 0.0) {
                    dvbf = 0.5;
                    BracketSlope(Pressure);
                } else {
                    dv = (Pressure - P_here)/dpdv;
                }
            }
        }
        double dvm = 0.2*v_here;
        if (v_here < dvs1) {
            dvm *= 0.5;
        }
        if (v_here < dvs2) {
            dvm *= 0.5;
        }
        if (kbr != 0) {
            double vt = v_here + dv;
            if ((vt < Vmin) || (vt > Vmax)) {
                dv = Vmin + (Pressure - Pmin)*(Vmax - Vmin)/(Pmax - Pmin) - v_here;
            }
        }
        double dva = fabs(dv);
        if (dva > dvm) {
            dv *= dvm/dva;
        }
        v_here += dv;
        if (dv == 0.0) {
            throw CanteraError("Substance::set_TPp", "dv = 0 and no convergence");
        }
        Set(PropertyPair::TV, Temp, v_here);
        LoopCount++;
        if (LoopCount > 200) {
            Set(PropertyPair::TV, Temp, v_save);
            throw CanteraError("Substance::set_TPp",
                               "No convergence for P = {}, T = {}\n"
                               "(P* = {}, V* = {})", Pressure, Temp,
                               Pressure/Pcrit(), v_save/Vcrit());
        }
    }
    Set(PropertyPair::TV, Temp,v_here);
}
}
