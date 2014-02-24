/**
 * @file Nuclei.h Provides class Nucleus
 */

#ifndef CT_NUCL_H
#define CT_NUCL_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"

namespace Cantera
{

/**
 * Represents atomic nuclei. These classes only provide minimal
 * information, and are designed only to handle nuclear statistics
 * effects on spectra.
 * @ingroup spectroscopy
 * @deprecated incomplete / abandoned
 */
class Nucleus
{
public:
    Nucleus(const std::string& symbol,
            int nP, int nN, doublereal spin) : m_np(nP),
        m_nn(nN), m_spin(spin),
        m_sym(symbol) {
            warn_deprecated("class Nucleus");
        }
    virtual ~Nucleus() {}
    int nProtons() {
        return m_np;
    }
    int mNeutrons() {
        return m_nn;
    }
    doublereal spin() {
        return m_spin;
    }
    int multiplicity() {
        return (int)(2*m_spin) + 1;
    }
    int atomicNumber() {
        return m_np;
    }
    std::string symbol() {
        return m_sym;
    }

    bool operator==(Nucleus& b) {
        if (m_np == b.m_np && m_nn == b.m_nn) {
            return true;
        } else {
            return false;
        }
    }

    bool operator!=(Nucleus& b) {
        return !(*this == b);
    }

    bool isBoson() {
        return (m_spin - std::floor(m_spin) < 0.001);
    }


protected:

    int m_np;       //< Number of protons
    int m_nn;       //< Number of electrons
    doublereal m_spin;     //< Spin.
    std::string m_sym;  //< Symbol.
};

inline Nucleus* HydrogenNucleus()
{
    return new Nucleus("H", 1, 0, 0.5);
}
inline Nucleus* DeuteriumNucleus()
{
    return new Nucleus("D", 1, 1, 1.0);
}
inline Nucleus* TritiumNucleus()
{
    return new Nucleus("T", 1, 2, 0.5);
}
inline Nucleus* He3Nucleus()
{
    return new Nucleus("He3", 2, 1, 0.5);
}
inline Nucleus* He4Nucleus()
{
    return new Nucleus("He3", 2, 2, 0.0);
}
inline Nucleus* C12nucleus()
{
    return new Nucleus("C12", 6, 6, 0.0);
}
inline Nucleus* C13nucleus()
{
    return new Nucleus("C13", 6, 7, 0.5);
}
inline Nucleus* N14nucleus()
{
    return new Nucleus("N14", 7, 7, 1.0);
}
inline Nucleus* N15nucleus()
{
    return new Nucleus("N15", 7, 8, 0.5);
}
inline Nucleus* O16nucleus()
{
    return new Nucleus("O16", 8, 8, 0.0);
}
inline Nucleus* O17nucleus()
{
    return new Nucleus("O17", 8, 9, 2.5);
}
inline Nucleus* O18nucleus()
{
    return new Nucleus("O18", 8, 10, 0.0);
}
inline Nucleus* F19nucleus()
{
    return new Nucleus("F19", 9, 10, 0.5);
}


} // Cantera

#endif
