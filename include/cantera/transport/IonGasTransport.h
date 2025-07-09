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

    string transportModel() const override {
        return "ionized-gas";
    }

    void init(ThermoPhase* thermo, int mode, int log_level=-7) override;

    //! Viscosity of the mixture  (kg/m/s).
    //! Only Neutral species contribute to Viscosity.
    double viscosity() override;

    //! Returns the mixture thermal conductivity (W/m/K).
    //! Only Neutral species contribute to thermal conductivity.
    double thermalConductivity() override;

    //! The mobilities for ions in gas.
    //! The ion mobilities are calculated by Blanc's law.
    void getMobilities(double* const mobi) override;

    //! The mixture transport for ionized gas.
    //! The binary transport between two charged species is neglected.
    void getMixDiffCoeffs(double* const d) override;

    /**
     * The electrical conductivity (Siemens/m).
     * @f[
     *     \sigma = \sum_k{\left|C_k\right| \mu_k \frac{X_k P}{k_b T}}
     * @f]
     */
    double electricalConductivity() override;

protected:
    //! setup parameters for n64 model
    void setupN64();

    //! Generate polynomial fits to the binary diffusion coefficients.
    //! Use Stockmayer-(n,6,4) model for collision between charged and neutral species.
    void fitDiffCoeffs(MMCollisionInt& integrals) override;

    /**
     * Collision integral of omega11 of n64 collision model.
     * The collision integral was fitted by Han et al. using the table
     * by Viehlan et al.
     * Note: Han release the range to 1000, but Selle suggested that
     * a high temperature model is needed for T* > 10.
     */
    double omega11_n64(const double tstar, const double gamma);

    //! electrical properties
    vector<double> m_speciesCharge;

    //! index of ions (exclude electron.)
    vector<size_t> m_kIon;

    //! index of neutral species
    vector<size_t> m_kNeutral;

    //! index of electron
    size_t m_kElectron = npos;

    //! parameter of omega11 of n64
    DenseMatrix m_gamma;

    //! polynomial of the collision integral for O2/O2-
    vector<double> m_om11_O2;
};

}

#endif
