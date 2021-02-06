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
 * References for Stockmayer-(n,6,4) model:
 *
 * 1. Selle, Stefan, and Uwe Riedel. "Transport properties of ionized species."
 *    Annals of the New York Academy of Sciences 891.1 (1999): 72-80.
 * 2. Selle, Stefan, and Uwe Riedel. "Transport coefficients of reacting air at
 *    high temperatures." 38th Aerospace Sciences Meeting and Exhibit. 1999.
 * 3. Han, Jie, et al. "Numerical modelling of ion transport in flames."
 *    Combustion Theory and Modelling 19.6 (2015): 744-772.
 *    DOI: 10.1080/13647830.2015.1090018
 * 4. Chiflikian, R. V. "The analog of Blanc’s law for drift velocities
 *    of electrons in gas mixtures in weakly ionized plasma."
 *    Physics of Plasmas 2.10 (1995): 3902-3909.
 * 5. Viehland, L. A., et al. "Tables of transport collision integrals for
 *    (n, 6, 4) ion-neutral potentials." Atomic Data and Nuclear Data Tables
 *    16.6 (1975): 495-514.
 *
 * Stockmayer-(n,6,4) model is not suitable for collision between O2/O2-
 * due to resonant charge transfer. Therefore, an experimental collision
 * data is used instead.
 *
 * Data taken from:
 *
 * Prager, Jens. Modeling and simulation of charged species in
 * lean methane-oxygen flames. Diss. 2005. Page 104.
 *
 * @ingroup tranprops
 */
class IonGasTransport : public MixTransport
{
public:
    IonGasTransport();

    virtual std::string transportType() const {
        return "Ion";
    }

    virtual void init(thermo_t* thermo, int mode, int log_level);

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
    size_t m_kElectron;

    //! parameter of omega11 of n64
    DenseMatrix m_gamma;

    //! polynomial of the collision integral for O2/O2-
    vector_fp m_om11_O2;
};

}

#endif
