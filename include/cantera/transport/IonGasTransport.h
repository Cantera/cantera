/**
 * @file IonGasTransport.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ION_GAS_TRANSPORT_H
#define CT_ION_GAS_TRANSPORT_H

#include "MixTransport.h"

namespace Cantera
{
//! Class IonGasTransport implements Stockmayer-(n,6,4) model for transport of ions.
/*!
 * As implemented here, only binary transport between neutrals and ions is considered
 * for calculating mixture-average diffusion coefficients and mobilities. When
 * polarizability is not provide for an ion, LJ model is used instead of n64 model.
 * Only neutral species are considered for thermal conductivity and viscosity.
 *
 * References for Stockmayer-(n,6,4) model: Selle and Riedel @cite selle1999,
 * @cite selle2000; Han et al. @cite han2015; Chiflikian @cite chiflikian1995; and
 * Viehland et al. @cite viehland1975.
 *
 * Stockmayer-(n,6,4) model is not suitable for collision between O2/O2-
 * due to resonant charge transfer. Therefore, an experimental collision
 * data is used instead.
 *
 * Data taken from @cite prager2005.
 *
 * @ingroup tranprops
 */
class IonGasTransport : public MixTransport
{
public:
    IonGasTransport() = default;

    virtual std::string transportModel() const {
        return "ionized-gas";
    }

    virtual void init(ThermoPhase* thermo, int mode, int log_level);

    //! Viscosity of the mixture  (kg/m/s).
    //! Only Neutral species contribute to Viscosity.
    virtual double viscosity();

    //! Returns the mixture thermal conductivity (W/m/K).
    //! Only Neutral species contribute to thermal conductivity.
    virtual double thermalConductivity();

    //! The mobilities for ions in gas.
    //! The ion mobilities are calculated by Blanc's law.
    virtual void getMobilities(double* const mobi);

    //! The mixture transport for ionized gas.
    //! The binary transport between two charged species is neglected.
    virtual void getMixDiffCoeffs(double* const d);

    /*! The electrical conductivity (Siemens/m).
     * \f[
     *     \sigma = \sum_k{\left|C_k\right| \mu_k \frac{X_k P}{k_b T}}
     * \f]
     */
    virtual double electricalConductivity();

protected:
    //! setup parameters for n64 model
    void setupN64();

    //! Generate polynomial fits to the binary diffusion coefficients.
    //! Use Stockmayer-(n,6,4) model for collision between charged and neutral species.
    virtual void fitDiffCoeffs(MMCollisionInt& integrals);

    /*!
     * Collision integral of omega11 of n64 collision model.
     * The collision integral was fitted by Han et al. using the table
     * by Viehlan et al.
     * Note: Han release the range to 1000, but Selle suggested that
     * a high temperature model is needed for T* > 10.
     */
    double omega11_n64(const double tstar, const double gamma);

    //! electrical properties
    vector_fp m_speciesCharge;

    //! index of ions (exclude electron.)
    std::vector<size_t> m_kIon;

    //! index of neutral species
    std::vector<size_t> m_kNeutral;

    //! index of electron
    size_t m_kElectron = npos;

    //! parameter of omega11 of n64
    DenseMatrix m_gamma;

    //! polynomial of the collision integral for O2/O2-
    vector_fp m_om11_O2;
};

}

#endif
