/**
 *  @file converters/writelog.cpp
 *
 */

// Copyright 2001  California Institute of Technology

#include <cstdio>
#include "writelog.h"

using namespace std;

namespace ckr
{


/// format a string for the log file message printed when starting a new task
string newTask(std::string msg)
{
    string s = "\n";
    s += msg + "...\n";
    return s;
}


/// print the falloff parameters for a pressure-dependent reaction
bool writeFalloff(int type, const vector_fp& c, std::ostream& log)
{

    log.precision(6);
    log.width(0);
    log.flags(ios::uppercase);

    //    bool ok = true;
    switch (type) {

    case Lindemann:
        log << "   Lindemann falloff function" << endl;
        return true;

    case Troe:
        log << "   Troe falloff function: " << endl;
        if (c.size() == 3) {
            log << "      alpha, T***, T* = (" << c[0] << ", " << c[1]
                << ", " << c[2] << ")" << endl;
        } else if (c.size() == 4) {
            log << "      alpha, T***, T*, T** = (" << c[0] << ", " << c[1]
                << ", " << c[2] << ", " << c[3] << ")" << endl;
        } else {
            for (size_t n = 0; n < c.size(); n++) {
                log << c[n] << ", ";
                log << endl;
            }
            log << "###### ERROR #####   incorrect number of parameters" << endl;
            return false;
        }
        return true;

    case SRI:
        log << "   SRI falloff function: " << endl;
        if (c.size() == 3) {
            log << "      a, b, c = (" << c[0] << ", " << c[1]
                << ", " << c[2] << ")" << endl;
        } else if (c.size() == 5) {
            log << "      a, b, c, d, e = (" << c[0] << ", " << c[1]
                << ", " << c[2] << ", " << c[3] << ", " << c[4]
                << ")" << endl;
        } else {
            for (size_t n = 0; n < c.size(); n++) {
                log << c[n] << ", ";
                log << endl;
            }
            log << "##### ERROR #####  incorrect number of parameters" << endl;
            return false;
        }
        return true;

    default:
        log << "unknown falloff type: " << type << endl;
        return false;
    }
}


/// print the rate coefficient parameters
bool writeRateCoeff(const RateCoeff& k, std::ostream& log)
{
    log.precision(10);
    log.width(0);
    log.flags(ios::uppercase);

    bool ok = true;
    size_t nb;

    switch (k.type) {

    case Arrhenius:
        log <<" A, n, E = (" << k.A << ", " << k.n
            << ", " << k.E << ")" << endl;
        break;

    case LandauTeller:
        log << "A, n, E, B, C = (" << k.A << ", " << k.n
            << ", " << k.E << ", " << k.B << ", " << k.C
            << ") *** Landau-Teller ***" << endl;
        break;

    case Jan:
        log <<" A, n, E = (" << k.A << ", " << k.n
            << ", " << k.E << ") *** JAN *** " << endl;
        nb = k.b.size();
        for (size_t n = 0; n < nb; n++) {
            log << "   b" << n+1 << ": " << k.b[n] << endl;
        }
        if (nb != 9) log
                    << "warning: number of b coefficients should be 9."
                    << endl;
        break;

    case Fit1:
        log <<" A, n, E = (" << k.A << ", " << k.n
            << ", " << k.E << ") *** FIT1 *** " << endl;
        nb = k.b.size();
        for (size_t n = 0; n < nb; n++) {
            log << "   b" << n+1 << ": " << k.b[n] << endl;
        }
        if (nb != 9) log
                    << "warning: number of b coefficients should be 4."
                    << endl;
        break;

    default:
        log << "unknown rate coefficient type: " << k.type << endl;
        ok = false;
    }
    return ok;
}

/**
 * Write onto an output stream the chemical equation for a reaction.
 */
void printReactionEquation(std::ostream& f, const Reaction& r)
{
    //        r.write(f);
    f << reactionEquation(r);
}


/**
 * Write to a string the chemical equation for a reaction.
 */
string reactionEquation(const Reaction& r)
{
    string s = "";
    int nr = static_cast<int>(r.reactants.size());
    int np = static_cast<int>(r.products.size());
    int k;
    double m;
    char buf[20];

    for (k = 0; k < nr; k++) {
        m =  r.reactants[k].number;
        if (m != 1.0) {
            sprintf(buf,"%g",m);
            s += string(buf);
            s += " ";
        }
        s +=  r.reactants[k].name;
        if (k < nr - 1) {
            s += " + ";
        }
    }

    if (r.isFalloffRxn) {
        s += " (+ " + r.thirdBody + ")";
    } else if (r.isThreeBodyRxn) {
        s += " + " + r.thirdBody;
    }
    if (r.isReversible) {
        s += " <=> ";
    } else {
        s += " => ";
    }

    for (k = 0; k < np; k++) {
        m =  r.products[k].number;
        if (m != 1.0) {
            sprintf(buf,"%g",m);
            s += string(buf);
            s += " ";
        }
        s += r.products[k].name;
        if (k < np - 1) {
            s += " + ";
        }
    }
    if (r.isFalloffRxn) {
        s += " (+ " + r.thirdBody + ")";
    } else if (r.isThreeBodyRxn) {
        s += " + " + r.thirdBody;
    }
    return s;
}



/**
 *
 * Write a summary of the properties of one species to the log file.
 *  @param log     log file output stream
 *  @param spec    instance of Species class
 */

void writeSpeciesData(std::ostream& log, const Species& spec)
{

    if (!spec.id.empty()) {
        log << endl << "   id/date: " << spec.id << endl;
    } else {
        log << " ... " << endl;
    }

    log << "   phase: "
        << spec.phase << endl
        << "   composition: (";

    for (size_t ie = 0; ie < spec.elements.size(); ie++) {
        if (!spec.elements[ie].name.empty()) {
            log.flags(ios::fixed);
            log.precision(0);
            if (ie > 0) {
                log << ", ";
            }
            log << spec.elements[ie].number << " "
                << spec.elements[ie].name;
        }
    }
    log << ")";

    if (spec.thermoFormatType == 0) {
        log.flags(ios::showpoint | ios::fixed);
        log.precision(2);
        log << endl << "   Tlow, Tmid, Thigh: (" << spec.tlow << ", " <<
            spec.tmid << ", " << spec.thigh << ")" << endl << endl;
        log << "   coefficients (low, high):" << endl;
        log.flags(ios::scientific | ios::uppercase | ios::internal);
        log.precision(8);
        for (int j = 0; j < 7; j++) {
            log << "   a" << j + 1;
            log.setf(ios::showpos);
            log << "  \t" << spec.lowCoeffs[j]
                << "  \t" << spec.highCoeffs[j] << endl;
            log.unsetf(ios::showpos);
        }
        log << endl;
    } else if (spec.thermoFormatType == 1) {
        log.flags(ios::showpoint | ios::fixed);
        log.precision(2);
        log << endl;
        log << "Number of temp regions = " << spec.nTempRegions << endl;
        for (int i = 0; i < spec.nTempRegions; i++) {
            log << "   Tlow, Thigh: (" << spec.minTemps[i] << ", "
                << spec.maxTemps[i] << ")" << endl << endl;
            log << "   coefficients :" << endl;
            log.flags(ios::scientific | ios::uppercase | ios::internal);
            log.precision(8);
            vector_fp& cc = *spec.region_coeffs[i];
            for (int j = 0; j < 9; j++) {
                log << "   a" << j + 1;
                log.setf(ios::showpos);
                log << "  \t" << cc[j]  << endl;
                log.unsetf(ios::showpos);
            }
            log << endl;
        }
    }
}


}


