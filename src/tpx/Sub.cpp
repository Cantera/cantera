/*
 * The Substance class
 * D. Goodwin, Caltech Nov. 1996
 */
#include "cantera/tpx/Sub.h"
#include <math.h>
#include <fstream>
#include <stdio.h>

using std::string;

namespace tpx
{

static string fp2str(double x, string fmt="%g")
{
    char buf[30];
    sprintf(buf, fmt.c_str(), x);
    return string(buf);
}

/*
    static string int2str(int n, string fmt="%d") {
        char buf[30];
        sprintf(buf, fmt.c_str(), n);
        return string(buf);
    }
*/

string TPX_Error::ErrorMessage = "";
string TPX_Error::ErrorProcedure = "";

string errorMsg(int flag)
{
    switch (flag) {
    case NoConverge:
        return "no convergence";
    case GenError:
        return "general error";
    case InvalidInput:
        return "invalid input";
    case TempError:
        return "temperature error";
    case PresError:
        return "pressure error";
    default:
        return "(unknown error)";
    }
}

//-------------- Public Member Functions --------------

Substance::Substance() :
    T(Undef),
    Rho(Undef),
    Tslast(Undef),
    Rhf(Undef),
    Rhv(Undef),
    Pst(Undef),
    Err(0),
    m_energy_offset(0.0),
    m_entropy_offset(0.0),
    kbr(0)
{
}

/// Pressure [Pa]. If two phases are present, return the
/// saturation pressure; otherwise return the pressure
/// computed directly from the underlying eos.
double Substance::P()
{
    double ppp = (TwoPhase() ? Ps() : Pp());
    return (Err ? Undef : ppp);
}

const double DeltaT = 0.000001;

/// The derivative of the saturation pressure
/// with respect to temperature.
double Substance::dPsdT()
{
    double ps1, tsave, dpdt;
    tsave = T;
    ps1 = Ps();
    set_T(T + DeltaT);
    dpdt = (Ps() - ps1)/DeltaT;
    set_T(tsave);
    return (Err ? Undef : dpdt);
}

/// true if a liquid/vapor mixture, false otherwise
int Substance::TwoPhase()
{
    if (T >= Tcrit()) {
        return 0;
    }
    update_sat();
    return ((Rho < Rhf) && (Rho > Rhv) ? 1 : 0);
}

/// Vapor fraction.
/// If T >= Tcrit, 0 is returned for v < Vcrit, and 1 is
/// returned if v > Vcrit.
double Substance::x()
{
    double xx, vv, vl;
    if (T >= Tcrit()) {
        return (1.0/Rho < Vcrit() ? 0.0 : 1.0);
    } else {
        update_sat();
        if (Rho <= Rhv) {
            xx = 1.0;
        } else if (Rho >= Rhf) {
            xx = 0.0;
        } else {
            vv = 1.0/Rhv;
            vl = 1.0/Rhf;
            xx = (1.0/Rho - vl)/(vv - vl);
        }
        return (Err ? Undef : xx);
    }
}

/// Saturation temperature at pressure p.
double Substance::Tsat(double p)
{
    double Tsave, p_here, dp, dt, dpdt, dta,
           dtm, tsat;
    if (Err || (p <= 0.0) || (p > Pcrit())) {
        throw TPX_Error("Substance::Tsat","illegal pressure value");
        set_Err(PresError);
        return Undef;
    }
    int LoopCount = 0;
    double tol = 1.e-6*p;
    Tsave = T;
    if (T < Tmin()) {
        T = 0.5*(Tcrit() - Tmin());
    }
    if (T >= Tcrit()) {
        T = 0.5*(Tcrit() - Tmin());
    }
    do {
        if (Err) {
            break;
        }
        if (T > Tcrit()) {
            T = Tcrit() - 0.001;
        }
        if (T < Tmin()) {
            T = Tmin() + 0.001;
        }
        p_here = Ps();
        dpdt = dPsdT();
        dp = p - p_here;
        dt = dp/dpdt;
        dta = fabs(dt);
        dtm = 0.1*T;
        if (dta > dtm) {
            dt = dt*dtm/dta;
        }
        T = T + dt;
        LoopCount++;
        if (LoopCount > 100) {
            T = Tsave;
            set_Err(NoConverge);
            return Undef;
        }
    } while (fabs(dp) > tol);
    tsat = T;
    T = Tsave;
    return (Err ? Undef : tsat);
}


// absolute tolerances

static const double TolAbsH = 0.0001; // J/kg
static const double TolAbsU = 0.0001;
static const double TolAbsS = 1.e-7;
static const double TolAbsP = 0.000; // Pa
static const double TolAbsV = 1.e-8;
static const double TolAbsT = 1.e-3;
static const double TolRel = 3.e-8;

void Substance::Set(int XY, double x0, double y0)
{
    double temp;

    clear_Err();  // clear error flag

    /* if inverted (PT) switch order and change sign of XY (TP = -PT) */
    if (XY < 0) {
        double tmp = x0;
        x0 = y0;
        y0 = tmp;
        XY *= -1;
    }

    switch (XY) {
    case TV:
        set_T(x0);
        set_v(y0);
        break;

    case HP:
        if (Lever(Pgiven, y0, x0, EvalH)) {
            return;
        }
        set_xy(EvalH, EvalP, x0, y0, TolAbsH, TolAbsP, TolRel, TolRel);
        break;

    case SP:
        if (Lever(Pgiven, y0, x0, EvalS)) {
            return;
        }
        set_xy(EvalS, EvalP, x0, y0, TolAbsS, TolAbsP, TolRel, TolRel);
        break;

    case PV:
        if (Lever(Pgiven, x0, y0, EvalV)) {
            return;
        }
        set_xy(EvalP, EvalV, x0, y0, TolAbsP, TolAbsV, TolRel, TolRel);
        break;

    case TP:
        if (x0 < Tcrit()) {
            set_T(x0);
            if (y0 < Ps()) {
                Set(TX, x0, Vapor);
            } else {
                Set(TX, x0, Liquid);
            }
        } else {
            set_T(x0);
        }
        set_xy(EvalT, EvalP, x0, y0, TolAbsT, TolAbsP, TolRel, TolRel);
        break;

    case UV:
        set_xy(EvalU, EvalV, x0, y0, TolAbsU, TolAbsV, TolRel, TolRel);
        break;

    case ST:
        if (Lever(Tgiven, y0, x0, EvalS)) {
            return;
        }
        set_xy(EvalS, EvalT, x0, y0, TolAbsS, TolAbsT, TolRel, TolRel);
        break;

    case SV:
        set_xy(EvalS, EvalV, x0, y0, TolAbsS, TolAbsV, TolRel, TolRel);
        break;

    case UP:
        if (Lever(Pgiven, y0, x0, EvalU)) {
            return;
        }
        set_xy(EvalU, EvalP, x0, y0, TolAbsU, TolAbsP, TolRel, TolRel);
        break;

    case VH:
        set_xy(EvalV, EvalH, x0, y0, TolAbsV, TolAbsH, TolRel, TolRel);
        break;

    case TH:
        set_xy(EvalT, EvalH, x0, y0, TolAbsT, TolAbsH, TolRel, TolRel);
        break;

    case SH:
        set_xy(EvalS, EvalH, x0, y0, TolAbsS, TolAbsH, TolRel, TolRel);
        break;

    case PX:
        temp = Tsat(x0);
        if ((y0 >= 0.0) && (y0 <= 1.0) && (temp < Tcrit())) {
            set_T(temp);
            update_sat();
            Rho = 1.0/((1.0 - y0)/Rhf + y0/Rhv);
        } else {
            set_Err(InvalidInput);
        }
        break;

    case TX:
        if ((y0 >= 0.0) && (y0 <= 1.0) && (x0 < Tcrit())) {
            set_T(x0);
            update_sat();
            Rho = 1.0/((1.0 - y0)/Rhf + y0/Rhv);
        } else {
            set_Err(InvalidInput);
        }
        break;

    default:
        set_Err(InvalidInput);
    }
    if (Err) {
        T = Undef;
        Rho = Undef;
        Tslast = Undef;
        Rhf = Undef;
        Rhv = Undef;
    }
}

//------------------ Protected and Private Functions -------------------

void Substance::set_Rho(double r0)
{
    if (r0 > 0.0) {
        Rho = r0;
    } else {
        set_Err(InvalidInput);
    }
}

void Substance::set_T(double t0)
{
    if ((t0 >= Tmin()) && (t0 <= Tmax())) {
        T = t0;
    } else {
        throw TPX_Error("Substance::set_T",
                        "illegal temperature value "+fp2str(t0));
        set_Err(TempError);
    }
}

void Substance::set_v(double v0)
{
    if (v0 > 0) {
        Rho = 1.0/v0;
    } else {
        throw TPX_Error("Substance::set_v",
                        "negative specific volume: "+fp2str(v0));
        set_Err(InvalidInput);
    }
}

double Substance::Ps()
{
    if (T < Tmin() || T > Tcrit()) {
        throw TPX_Error("Substance::Ps",
                        "illegal temperature value "+fp2str(T));
        set_Err(TempError);
        return Undef;
    }
    update_sat();
    return Pst;
}

// update saturated liquid and vapor densities and saturation pressure

void Substance::Set_meta(double phase, double pp)
{
    if (phase == Liquid) {
        Rho = ldens();    // trial value = liquid dens
    } else {
        Rho = pp*MolWt()/(8314.0*T);    // trial value = ideal gas
    }
    set_TPp(T, pp);
}

void Substance::update_sat()
{
    if ((T != Tslast) && (T < Tcrit())) {
        double Rho_save = Rho;
        double gf, gv, dg, dp, dlp, psold;

        double pp = Psat();
        double lps = log(pp);
        // trial value = Psat from correlation
        int i;

        for (i = 0; i<20; i++) {
            if (i==0) {
                Rho = ldens();                // trial value = liquid density
            } else {
                Rho = Rhf;
            }
            set_TPp(T,pp);
            Rhf = Rho;                    // sat liquid density

            gf = hp() - T*sp();
            if (i==0) {
                Rho = pp*MolWt()/(8314.0*T);  // trial value = ideal gas
            } else {
                Rho = Rhv;
            }
            set_TPp(T,pp);

            Rhv = Rho;                    // sat vapor density
            gv = hp() - T*sp();
            dg = gv - gf;

            if (Rhv > Rhf) {
                std::swap(Rhv, Rhf);
                dg = - dg;
            }

            if (fabs(dg) < 0.001 && Rhf > Rhv) {
                break;
            }
            dp = dg/(1.0/Rhv - 1.0/Rhf);
            psold = pp;

            if (fabs(dp) > pp) {
                dlp = dg/(pp*(1.0/Rhv - 1.0/Rhf));
                lps -= dlp;
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
            Pst = Undef;
            Rhv = Undef;
            Rhf = Undef;
            Tslast = Undef;
            throw TPX_Error("substance::update_sat","no convergence");
            set_Err(NoConverge);
        } else {
            Pst = pp;
            Tslast = T;
        }
        Rho = Rho_save;
    }
}

double Substance::vprop(int ijob)
{
    switch (ijob) {
    case EvalH:
        return hp();
    case EvalS:
        return sp();
    case EvalU:
        return up();
    case EvalV:
        return vp();
    case EvalP:
        return Pp();
    default:
        return Undef;
    }
}

int Substance::Lever(int itp, double sat, double val, int ifunc)
{
    /*
     * uses lever rule to set state in the dome. Returns 1 if in dome,
     * 0 if not, in which case state not set.
     */
    double Valf, Valg, Tsave, Rhosave, xx, vv, psat;
    Tsave = T;
    Rhosave = Rho;
    if (itp == Tgiven) {
        if (sat >= Tcrit()) {
            return 0;
        }
        set_T(sat);
        psat = Ps();
    } else if (itp == Pgiven) {
        if (sat >= Pcrit()) {
            return 0;
        }
        psat = sat;
        T = Tsat(psat);
        if (T == Undef) {
            Err = 0;
            T = Tsave;
            Rho = Rhosave;
            return 0;
        }
    } else {
        throw TPX_Error("Substance::Lever","general error");
        set_Err(GenError);
        return GenError;
    }
    Set(TX, T, Vapor);
    Valg = vprop(ifunc);
    Set(TX, T, Liquid);
    Valf = vprop(ifunc);
    if (Err) {
        return Err;
    } else if ((val >= Valf) && (val <= Valg)) {
        xx = (val - Valf)/(Valg - Valf);
        vv = (1.0 - xx)/Rhf + xx/Rhv;
        set_v(vv);
        return 1;
    } else {
        T = Tsave;
        Rho = Rhosave;
        return 0;
    }
}


void Substance::set_xy(int ifx, int ify, double X, double Y,
                       double atx, double aty,
                       double rtx, double rty)
{

    double v_here, t_here, dv, dt, dxdt, dydt, dxdv, dydv,
           det, x_here, y_here, dvm, dtm, dva, dta;
    double Xa, Ya, err_x, err_y;
    double dvs1 = 2.0*Vcrit();
    double dvs2 = 0.7*Vcrit();
    int LoopCount = 0;

    double v_save = 1.0/Rho;
    double t_save = T;

    if (Err) {
        return;
    }

    if ((T == Undef) && (Rho == Undef)) {  // new object, try to pick
        Set(TV,Tcrit()*1.1,Vcrit()*1.1);   // "reasonable" starting point
        t_here = T;
        v_here = 1.0/Rho;
    } else if (Rho == Undef) { // new object, try to pick
        Set(TV,T,Vcrit()*1.1);   // "reasonable" starting point
        t_here = T;
        v_here = 1.0/Rho;
    } else {
        v_here = v_save;
        t_here = t_save;
    }

    Xa = fabs(X);
    Ya = fabs(Y);


    // loop
    do {
        if (Err) {
            break;
        }
        x_here = prop(ifx);
        y_here = prop(ify);
        err_x = fabs(X - x_here);
        err_y = fabs(Y - y_here);

        if ((err_x < atx + rtx*Xa) && (err_y < aty + rty*Ya)) {
            break;
        }

        /* perturb t */
        dt = 0.001*t_here;
        if (t_here + dt > Tmax()) {
            dt *= -1.0;
        }

        /* perturb v */
        dv = 0.001*v_here;
        if (v_here <= Vcrit()) {
            dv *= -1.0;
        }

        /* derivatives with respect to T */
        Set(TV, t_here + dt, v_here);
        dxdt = (prop(ifx) - x_here)/dt;
        dydt = (prop(ify) - y_here)/dt;

        /* derivatives with respect to v */
        Set(TV, t_here, v_here + dv);
        dxdv = (prop(ifx) - x_here)/dv;
        dydv = (prop(ify) - y_here)/dv;

        det = dxdt*dydv - dydt*dxdv;
        dt = ((X - x_here)*dydv - (Y - y_here)*dxdv)/det;
        dv = ((Y - y_here)*dxdt - (X - x_here)*dydt)/det;

        dvm = 0.2*v_here;
        if (v_here < dvs1) {
            dvm *= 0.5;
        }
        if (v_here < dvs2) {
            dvm *= 0.5;
        }
        dtm = 0.1*t_here;
        dva = fabs(dv);
        dta = fabs(dt);
        if (dva > dvm) {
            dv *= dvm/dva;
        }
        if (dta > dtm) {
            dt *= dtm/dta;
        }
        v_here += dv;
        t_here += dt;
        if (t_here >= Tmax()) {
            t_here = Tmax() - 0.001;
        } else if (t_here <= Tmin()) {
            t_here = Tmin() + 0.001;
        }
        if (v_here <= 0.0) {
            v_here = 0.0001;
        }
        Set(TV, t_here, v_here);
        LoopCount++;
        if (LoopCount > 200) {
            throw TPX_Error("Substance::set_xy","no convergence");
            set_Err(NoConverge);
            break;
        }
    } while (1);
}


double Substance::prop(int ijob)
{
    double xx, pp, lp, vp, Rho_save;
    if (ijob == EvalP) {
        return P();
    }
    if (ijob == EvalT) {
        return T;
    }
    xx = x();
    if ((xx > 0.0) && (xx < 1.0)) {
        Rho_save = Rho;
        Rho = Rhv;
        vp = vprop(ijob);
        Rho = Rhf;
        lp = vprop(ijob);
        pp = (1.0 - xx)*lp + xx*vp;
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
    set_T(Temp);
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
            Set(TV, Temp, v_here+dv);
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
                    set_Err(GenError);
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
            set_Err(NoConverge);
            return;
        }
        Set(TV, Temp, v_here);
        LoopCount++;
        if (LoopCount > 100) {
            Set(TV, Temp, v_save);
            throw TPX_Error("Substance::set_TPp",string("no convergence for ")
                            +"P* = "+fp2str(Pressure/Pcrit())+". V* = "
                            +fp2str(v_save/Vcrit()));
            set_Err(NoConverge);
            return;
        }
    }
    Set(TV, Temp,v_here);
}
}
