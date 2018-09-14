/**
 * @file BinarySolutionTabulatedThermo.h
 * Header file for an binary solution model with tabulated standard state
 * thermodynamic data (see \ref thermoprops and class
 * \link Cantera::BinarySolutionTabulatedThermo BinarySolutionTabulatedThermo\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_BINARYSOLUTIONTABULATEDTHERMO_H
#define CT_BINARYSOLUTIONTABULATEDTHERMO_H

#include "IdealSolidSolnPhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

//!   Overloads the virtual methods of class IdealSolidSolnPhase to implement the
//!   tabulated thermodynamics for one species.
/**
 *
 *
 * @ingroup thermoprops
 */
class BinarySolutionTabulatedThermo : public IdealSolidSolnPhase
{
public:
    /**
     * Constructor for BinarySolutionTabulatedThermo.
     * The generalized concentrations can have three different forms
     * depending on the value of the member attribute #m_formGC, which
     * is supplied in the constructor or read from the XML data file.
     *
     * @param formCG This parameter initializes the #m_formGC variable.
     */
    BinarySolutionTabulatedThermo(int formCG=0);

    //! Construct and initialize an BinarySolutionTabulatedThermo ThermoPhase object
    //! directly from an ASCII input file
    /*!
     * This constructor will also fully initialize the object.
     * The generalized concentrations can have three different forms
     * depending on the value of the member attribute #m_formGC, which
     * is supplied in the constructor or read from the XML data file.
     *
     * @param infile File name for the XML datafile containing information
     *               for this phase
     * @param id     The name of this phase. This is used to look up
     *               the phase in the XML datafile.
     * @param formCG This parameter initializes the #m_formGC variable.
     */
    BinarySolutionTabulatedThermo(const std::string& infile, const std::string& id="", int formCG=0);

    //! Construct and initialize an BinarySolutionTabulatedThermo ThermoPhase object
    //! directly from an XML database
    /*!
     * The generalized concentrations can have three different forms
     * depending on the value of the member attribute #m_formGC, which
     * is supplied in the constructor and/or read from the data file.
     *
     * @param root   XML tree containing a description of the phase.
     *               The tree must be positioned at the XML element
     *               named phase with id, "id", on input to this routine.
     * @param id     The name of this phase. This is used to look up
     *               the phase in the XML datafile.
     * @param formCG This parameter initializes the #m_formGC variable.
     */
    BinarySolutionTabulatedThermo(XML_Node& root, const std::string& id="", int formCG=0);

    virtual std::string type() const {
        return "BinarySolutionTabulatedThermo";
    }

    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id_);

protected:

    int m_formGC;

    double m_Pref;

    double m_Pcurrent;

    vector_fp m_speciesMolarVolume;

    //! If the compositions have changed, update the tabulated thermo lookup
    virtual void compositionChanged();

    //! Species thermodynamics interpolation functions
    double* interpolate(double x) const;

    //! Current tabulated species index
    size_t m_kk_tab;

    //! Current tabulated species mole fraction
    double m_xlast;

    //! Vector for storing tabulated thermo
    vector_fp m_molefrac_tab;
    vector_fp m_enthalpy_tab;
    vector_fp m_entropy_tab;

private:
    void _updateThermo();
};
}

#endif
