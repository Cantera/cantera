//! @file IdealGasMix.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CXX_IDEALGASMIX
#define CXX_IDEALGASMIX

#pragma message("warning: IdealGasMix.h is deprecated and will be removed after Cantera 2.5.")

#include "thermo/IdealGasPhase.h"
#include "kinetics/GasKinetics.h"
#include "kinetics/importKinetics.h"
#include "base/stringUtils.h"

namespace Cantera
{

//! Convenience class which inherits from both IdealGasPhase and GasKinetics
/*!
 * @deprecated To be removed after Cantera 2.5.
 *             Replaceable with Solution and/or IdealGasPhase/GasKinetics.
 */
class IdealGasMix :
    public IdealGasPhase,
    public GasKinetics
{
public:
    IdealGasMix() : m_ok(false), m_r(0) {}

    IdealGasMix(const std::string& infile, std::string id_="") :
        m_ok(false), m_r(0)
    {
        warn_deprecated("class IdealGasMix",
            "To be removed after Cantera 2.5. "
            "Replaceable with Solution and/or IdealGasPhase/GasKinetics.");
        m_r = get_XML_File(infile);
        m_id = id_;
        if (id_ == "-") {
            id_ = "";
        }
        m_ok = buildSolutionFromXML(*m_r,
                                    m_id, "phase", this, this);
        if (!m_ok) throw CanteraError("IdealGasMix",
                                          "Cantera::buildSolutionFromXML returned false");
    }

    IdealGasMix(XML_Node& root,
                std::string id_) : m_ok(false), m_r(&root), m_id(id_) {
        warn_deprecated("class IdealGasMix",
            "To be removed after Cantera 2.5. "
            "Replaceable with Solution and/or IdealGasPhase/GasKinetics.");
        m_ok = buildSolutionFromXML(root, id_, "phase", this, this);
    }

    IdealGasMix(const IdealGasMix& other) : m_ok(false),
        m_r(other.m_r),
        m_id(other.m_id) {
        warn_deprecated("class IdealGasMix",
            "To be removed after Cantera 2.5. "
            "Replaceable with Solution and/or IdealGasPhase/GasKinetics.");
        m_ok = buildSolutionFromXML(*m_r, m_id, "phase", this, this);
    }

    bool operator!() {
        return !m_ok;
    }
    bool ready() const {
        return m_ok;
    }
    friend std::ostream& operator<<(std::ostream& s, IdealGasMix& mix) {
        std::string r = mix.report(true);
        s << r;
        return s;
    }

protected:
    bool m_ok;
    XML_Node* m_r;
    std::string m_id;
};
}

#endif
