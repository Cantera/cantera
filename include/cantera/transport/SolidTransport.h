/**
 *  @file SolidTransport.h
 *   Header file for defining the class SolidTransport, which handles transport
 *   of ions within solid phases
 *  (see \ref tranprops and \link Cantera::SolidTransport SolidTransport \endlink).
 */

// Copyright 2003  California Institute of Technology

#ifndef CT_SOLIDTRAN_H
#define CT_SOLIDTRAN_H

#include "LTPspecies.h"
#include "TransportBase.h"
#include "cantera/numerics/DenseMatrix.h"

namespace Cantera
{
//! Class SolidTransport implements transport properties for solids.
//! @ingroup tranprops
class SolidTransport : public Transport
{
public:
    SolidTransport();
    SolidTransport(const SolidTransport& right);
    SolidTransport&  operator=(const SolidTransport& right);
    virtual Transport* duplMyselfAsTransport() const;

    virtual int model() const {
        return cSolidTransport;
    }

    //! Returns the ionic conductivity of the phase
    /*!
     *  The thermo phase needs to be updated (temperature) prior to calling this.
     *  The ionConductivity calculation is handled by subclasses of
     *  LTPspecies as specified in the input file.
     */
    virtual doublereal ionConductivity() ;

    //! Returns the thermal conductivity of the phase
    /*!
     *  The thermo phase needs to be updated (temperature) prior to calling this.
     *  The thermalConductivity calculation is handled by subclasses of
     *  LTPspecies as specified in the input file.
     *
     *  There is also a legacy method to evaluate
     * \f[
     * \lambda = A T^n \exp(-E/RT)
     * \f]
     */
    virtual doublereal thermalConductivity();

    //! Returns the electron conductivity of the phase
    /*!
     * The thermo phase needs to be updated (temperature) prior to calling
     * this. The ionConductivity calculation is handled by subclasses of
     * LTPspecies as specified in the input file.
     *
     * There is also a legacy multicomponent diffusion approach to electrical
     * conductivity.
     */
    virtual doublereal electricalConductivity();

    /*!
     * The diffusivity of defects in the solid (m^2/s). The thermo phase needs
     *  to be updated (temperature) prior to calling this. The
     *  defectDiffusivity calculation is handled by subclasses of LTPspecies
     *  as specified in the input file.
     */
    virtual doublereal defectDiffusivity();

    /**
     * The activity of defects in the solid.
     * At some point this should be variable and the diffusion coefficient should depend on it.
     *
     * The thermo phase needs to be updated (temperature) prior to calling this.
     * The defectActivity calculation is handled by subclasses of
     * LTPspecies as specified in the input file.
     */
    virtual doublereal defectActivity();

    /*
     * The diffusion coefficients are computed from
     *
     * \f[
     * D_k = A_k T^{n_k} \exp(-E_k/RT).
     * \f]
     *
     * The diffusion coefficients are only non-zero for species for which
     * parameters have been specified using method setParameters.
     *  @todo HEWSON WONDERS IF THE FOLLOWING ARE RELEVANT??
     */
    virtual void getMixDiffCoeffs(doublereal* const d);

    virtual void getMobilities(doublereal* const mobil);

    //! @deprecated
    virtual void setParameters(const int n, const int k, const doublereal* const p);

    friend class TransportFactory;

protected:
    //! Initialize the transport object
    /*!
     * Here we change all of the internal dimensions to be sufficient. We get
     * the object ready to do property evaluations. A lot of the input
     * required to do property evaluations is contained in the
     * SolidTransportParams class that is filled in TransportFactory.
     *
     * @param tr  Transport parameters for all of the species
     *            in the phase.
     */
    virtual bool initSolid(SolidTransportData& tr);

private:
    //! Model type for the ionic conductivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* m_ionConductivity;

    //! Model type for the thermal conductivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* m_thermalConductivity;

    //! Model type for the electrical conductivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* m_electConductivity;

    //! Model type for the defectDiffusivity -- or more like a defect diffusivity in the context of the solid phase.
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* m_defectDiffusivity;

    //! Model type for the defectActivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* m_defectActivity;

    //! number of mobile species
    size_t m_nmobile;

    //! Coefficient for the diffusivity of species within a solid
    /*!
     *  This is with respect to the lattice
     *  units = m**2 / s
     *  vector of length m_nmobile
     */
    vector_fp m_Adiff;

    //! Temperature power coefficient for the diffusivity of species in a solid
    /*!
     *  vector of length m_nmobile
     */
    vector_fp m_Ndiff;

    //! Arrhenius factor for the species diffusivities of a solid
    /*!
     *  units = temperature
     *  vector of length m_nmobile
     */
    vector_fp m_Ediff;

    //! Index of mobile species to global species
    /*!
     *  vector of length m_nmobile
     */
    vector_int m_sp;

    //! Coefficient for the thermal conductivity of a solid
    /*!
     *  units = kg m / s3 /K   = W/m/K
     */
    doublereal m_Alam;

    //! Temperature power coefficient for the thermal conductivity of a solid
    doublereal m_Nlam;

    //! Arrhenius factor for the thermal conductivity of a solid
    /*!
     *  units = temperature
     */
    doublereal m_Elam;

    //! extra fp array of length nSpecies()
    vector_fp m_work;
};
}
#endif
