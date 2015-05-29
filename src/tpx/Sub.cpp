//! @file Sub.cpp
/*
 * The Substance class
 * D. Goodwin, Caltech Nov. 1996
 */
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

double Substance::P()
{
    return TwoPhase() ? Ps() : Pp();
}

const double DeltaT = 0.000001;

double Substance::dPsdT()
{
    double tsave = T;
    double ps1 = Ps();
    T = T + DeltaT;
    double dpdt = (Ps() - ps1)/DeltaT;
    T = tsave;
    return dpdt;
}

int Substance::TwoPhase()
{
    if (T >= Tcrit()) {
        return 0;
    }
    update_sat();
    return ((Rho < Rhf) && (Rho > Rhv) ? 1 : 0);
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
    if (p <= 0.0 || p > Pcrit()) {
        throw TPX_Error("Substance::Tsat", "illegal pressure value");
    }
    int LoopCount = 0;
    double tol = 1.e-6*p;
    double Tsave = T;
    if (T < Tmin()) {
        T = 0.5*(Tcrit() - Tmin());
    }
    if (T >= Tcrit()) {
        T = 0.5*(Tcrit() - Tmin());
    }
    double dp;
    do {
        if (T > Tcrit()) {
            T = Tcrit() - 0.001;
        }
        if (T < Tmin()) {
            T = Tmin() + 0.001;
        }
        dp = p - Ps();
        double dt = dp/dPsdT();
        double dta = fabs(dt);
        double dtm = 0.1*T;
        if (dta > dtm) {
            dt = dt*dtm/dta;
        }
        T = T + dt;
        LoopCount++;
        if (LoopCount > 100) {
            T = Tsave;
            throw TPX_Error("Substance::Tsat", "No convergence");
        }
    } while (fabs(dp) > tol);
    double tsat = T;
    T = Tsave;
    return tsat;
}

// absolute tolerances
static const double TolAbsH = 0.0001; // J/kg
static const double TolAbsU = 0.0001;
static const double TolAbsS = 1.e-7;
static const double TolAbsP = 0.000; // Pa
static const double TolAbsV = 1.e-8;
static const double TolAbsT = 1.e-3;
static const double TolRel = 3.e-8;

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
        if (x0 < Tcrit()) {
            set_T(x0);
            if (y0 < Ps()) {
                Set(PropertyPair::TX, x0, 1.0);
            } else {
                Set(PropertyPair::TX, x0, 0.0);
            }
        } else {
            set_T(x0);
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
        temp = Tsat(x0);
        if (y0 > 1.0 || y0 < 0.0) {
            throw TPX_Error("Substance::Set",
                            "Invalid vapor fraction, " + fp2str(y0));
        } else if (temp >= Tcrit()) {
            throw TPX_Error("Substance::Set",
                            "Can't set vapor fraction above the critical point");
        } else {
            set_T(temp);
            update_sat();
            Rho = 1.0/((1.0 - y0)/Rhf + y0/Rhv);
        }
        break;

    case PropertyPair::TX:
        if (y0 > 1.0 || y0 < 0.0) {
            throw TPX_Error("Substance::Set",
                            "Invalid vapor fraction, " + fp2str(y0));
        } else if (x0 >= Tcrit()) {
            throw TPX_Error("Substance::Set",
                            "Can't set vapor fraction above the critical point");
        } else {
            set_T(x0);
            update_sat();
            Rho = 1.0/((1.0 - y0)/Rhf + y0/Rhv);
        }
        break;

    default:
        throw TPX_Error("Substance::Set", "Invalid input.");
    }
}

//------------------ Protected and Private Functions -------------------

void Substance::set_Rho(double r0)
{
    if (r0 > 0.0) {
        Rho = r0;
    } else {
        throw TPX_Error("Substance::set_Rho", "Invalid density: " + fp2str(r0));
    }
}

void Substance::set_T(double t0)
{
    if ((t0 >= Tmin()) && (t0 <= Tmax())) {
        T = t0;
    } else {
        throw TPX_Error("Substance::set_T", "illegal temperature: " + fp2str(t0));
    }
}

void Substance::set_v(double v0)
{
    if (v0 > 0) {
        Rho = 1.0/v0;
    } else {
        throw TPX_Error("Substance::set_v",
                        "negative specific volume: "+fp2str(v0));
    }
}

double Substance::Ps()
{
    if (T < Tmin() || T > Tcrit()) {
        throw TPX_Error("Substance::Ps",
                        "illegal temperature value "+fp2str(T));
    }
    update_sat();
    return Pst;
}

void Substance::update_sat()
{
    if ((T != Tslast) && (T < Tcrit())) {
        double Rho_save = Rho;

        double pp = Psat();
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
                Rho = pp*MolWt()/(8314.0*T); // trial value = ideal gas
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

            if (fabs(dg) < 0.001 && Rhf > Rhv) {
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
        if (Rhf <= Rhv) {
            throw TPX_Error("Substance::update_sat",
                            "wrong root found for sat. liquid or vapor at P = "+fp2str(pp));
        }

        if (i >= 20) {
            throw TPX_Error("substance::update_sat","no convergence");
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
        throw TPX_Error("Substance::vprop", "invalid job index");
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
        } catch (TPX_Error&) {
            // Failure to converge here is not an error
            T = Tsave;
            Rho = Rhosave;
            return 0;
        }
    } else {
        throw TPX_Error("Substance::Lever","general error");
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
            std::string msg = "No convergence. " +
                propertySymbols[ifx] + " = " + fp2str(X) + ", " +
                propertySymbols[ify] + " = " + fp2str(Y);
            if (t_here == Tmin()) {
                msg += "\nAt temperature limit (Tmin = " + fp2str(Tmin()) + ")";
            } else if (t_here == Tmax()) {
                msg += "\nAt temperature limit (Tmax = " + fp2str(Tmax()) + ")";
            }
            throw TPX_Error("Substance::set_xy", msg);
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
                    throw TPX_Error("Substance::set_TPp","Vmin >= Vmax");
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
            throw TPX_Error("Substance::set_TPp","dv = 0 and no convergence");
        }
        Set(PropertyPair::TV, Temp, v_here);
        LoopCount++;
        if (LoopCount > 100) {
            Set(PropertyPair::TV, Temp, v_save);
            throw TPX_Error("Substance::set_TPp",string("no convergence for ")
                            +"P* = "+fp2str(Pressure/Pcrit())+". V* = "
                            +fp2str(v_save/Vcrit()));
        }
    }
    Set(PropertyPair::TV, Temp,v_here);
}
}
