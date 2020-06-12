/**
 * @file BinarySolutionTabulatedThermo.h
 * Header file for an binary solution model with tabulated standard state
 * thermodynamic data (see \ref thermoprops and class
 * \link Cantera::BinarySolutionTabulatedThermo BinarySolutionTabulatedThermo\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BINARYSOLUTIONTABULATEDTHERMO_H
#define CT_BINARYSOLUTIONTABULATEDTHERMO_H

#include "IdealSolidSolnPhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

//! Overloads the virtual methods of class IdealSolidSolnPhase to implement
//! tabulated standard state thermodynamics for one species in a binary
//! solution.
/**
 *
 *  BinarySolutionTabulatedThermo is derived from IdealSolidSolnPhase, but
 *  overwrites the standard state thermodynamic data using tabulated data,
 *  as provided by the user in the input file.  This ends up being useful for
 *  certain non-ideal / non-dilute species where the interaction potentials, as
 *  a function of composition / solute mole fraction, are not easily represented
 *  by any closed-form equation of state.
 *
 *  A good example of this type of phase is intercalation-based lithium storage
 *  materials used for lithium-ion battery electrodes.  Measuring the open
 *  circuit voltage \f$ E_eq \f$, relative to a reference electrode, as a
 *  function of lithium mole fraction and as a function of temperature, provides
 *  a means to evaluate the gibbs free energy of reaction:
 *
 *  \f[
 *  \Delta g_{\rm rxn} = -\frac{E_eq}{nF}
 *  \f]
 *
 *  where \f$ n\f$ is the charge number transferred to the phase, via the
 *  reaction, and \f$ F \f$ is Faraday's constant.  The gibbs energy of
 *  reaction, in turn, can be separated into enthalpy and entropy of reaction
 *  components:
 *
 *  \f[
 *  \Delta g_{\rm rxn} = \Delta h_{\rm rxn} - T\Delta s_{\rm rxn}
 *  \f]
 *  \f[
 *  \frac{d\Delta g_{\rm rxn}}{dT} =  - \Delta s_{\rm rxn}
 *  \f]
 *
 *  For the tabulated binary phase, the user identifies a 'tabulated' species,
 *  while the other is considered the 'reference' species.  The standard state
 *  thermo variables for the tabulated species therefore incorporate any and all
 *  excess energy contributions, and are calculated according to the reaction
 *  energy terms:
 *
 *  \f[
 *  \Delta h_{\rm rxn} = \sum_k \nu_k h^{\rm o}_k
 *  \f]
 *  \f[
 *  \Delta s_{\rm rxn} = \sum_k \nu_k s^{\rm o}_k + RT\ln\left(\prod_k\left(\frac{c_k}{c^{\rm o}_k} \right)^{\nu_k}\right)
 *  \f]
 *
 *  Where the 'reference' species is automatically assigned standard state
 *  thermo variables \f$ h^{\rm o} = 0\f$ and \f$ s^{\rm o} = 0\f$, and standard
 *  state thermo variables for species in any other phases are calculated
 *  according to the rules specified in that phase definition.
 *
 *  The present model is intended for modeling non-ideal, tabulated
 *  thermodynamics for binary solutions where the tabulated species is
 *  incorporated via an electrochemical reaction, such that the open circuit
 *  voltage can be measured, relative to a counter electrode species with
 *  standard state thermo properties \f$ h^{\rm o} = 0\f$.
 *  It is possible that this can be generalized such that this assumption about
 *  the counter-electrode is not required. At present, this is left as future
 *  work.
 *
 *  The user therefore provides a table of three equally-sized vectors of
 *  tabulated data:
 *
 *  - \f$ x_{\rm tab}\f$ = array of mole fractions for the tabulated species
 *                         at which measurements were conducted and thermo
 *                         data are provided.
 *  - \f$ h_{\rm tab}\f$ = \f$ F\left(-E_{\rm eq}\left(x,T^{\rm o} \right) + T^{\rm o} \frac{dE_{\rm eq}\left(x,T^{\rm o} \right)}{dT}\right) \f$
 *  - \f$ s_{\rm tab}\f$ = \f$ F \left(\frac{dE_{\rm eq}\left(x,T^{\rm o} \right)}{dT} + s_{\rm counter}^{\rm o} \right) \f$
 *
 *  where \f$ E_{\rm eq}\left(x,T^{\rm o} \right) \f$ and \f$ \frac{dE_{\rm eq}\left(x,T^{\rm o} \right)}{dT} \f$
 *  are the experimentally-measured open circuit voltage and derivative in
 *  open circuit voltage with respect to temperature, respectively, both
 *  measured as a mole fraction of \f$ x \f$ for the tabulated species and at a
 *  temperature of \f$ T^{\rm o} \f$.  The arrays \f$ h_{\rm tab}\f$ and
 *  \f$ s_{\rm tab}\f$ must be the same length as the \f$ x_{\rm tab}\f$ array.
 *
 *  From these tabulated inputs, the standard state thermodynamic properties
 *  for the tabulated species (subscript \f$ k\f$, tab) are calculated as:
 *
 *  \f[
 *   h^{\rm o}_{k,\,{\rm tab}} =  h_{\rm tab}
 *  \f]
 *  \f[
 *   s^{\rm o}_{k,\,{\rm tab}} =  s_{\rm tab} + R\ln\frac{x_{k,\,{\rm tab}}}{1-x_{k,\,{\rm tab}}} + \frac{R}{F} \ln\left(\frac{c^{\rm o}_{k,\,{\rm ref}}}{c^{\rm o}_{k,\,{\rm tab}}}\right)
 *  \f]
 *
 *  Now, whenever the composition has changed, the lookup/interpolation of the
 *  tabulated thermo data is performed to update the standard state
 *  thermodynamic data for the tabulated species.
 *
 * @ingroup thermoprops
 */
class BinarySolutionTabulatedThermo : public IdealSolidSolnPhase
{
public:
    //! Default constructor for BinarySolutionTabulatedThermo
    BinarySolutionTabulatedThermo();

    //! Construct and initialize an BinarySolutionTabulatedThermo ThermoPhase object
    //! directly from an input file
    /*!
     * This constructor will also fully initialize the object.
     *
     * @param infile File name for the input file containing information
     *               for this phase
     * @param id     The name of this phase. This is used to look up
     *               the phase in the input file.
     */
    BinarySolutionTabulatedThermo(const std::string& infile, const std::string& id="");

    //! Construct and initialize an BinarySolutionTabulatedThermo ThermoPhase object
    //! directly from an XML database
    /*!
     * @param root   XML tree containing a description of the phase.
     *               The tree must be positioned at the XML element
     *               named phase with id, "id", on input to this routine.
     * @param id     The name of this phase. This is used to look up
     *               the phase in the XML datafile.
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    BinarySolutionTabulatedThermo(XML_Node& root, const std::string& id="");

    virtual std::string type() const {
        return "BinarySolutionTabulatedThermo";
    }

    virtual void initThermo();
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id_);

protected:
    //! If the compositions have changed, update the tabulated thermo lookup
    virtual void compositionChanged();

    //! Species thermodynamics interpolation functions
    std::pair<double,double> interpolate(double x) const;

    //! Current tabulated species index
    size_t m_kk_tab;

    //! Current tabulated species mole fraction
    mutable double m_xlast;

    //! Tabulated contribution to h0[m_kk_tab] at the current composition
    mutable double m_h0_tab;

    //! Tabulated contribution to s0[m_kk_tab] at the current composition
    mutable double m_s0_tab;

    //! Vector for storing tabulated thermo
    vector_fp m_molefrac_tab;
    vector_fp m_enthalpy_tab;
    vector_fp m_entropy_tab;

private:
    virtual void _updateThermo() const;
};
}

#endif
