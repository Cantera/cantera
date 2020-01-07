//! @file YamlWriter.h Declaration for class Cantera::YamlWriter.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_YAMLWRITER_H
#define CT_YAMLWRITER_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class Solution;

//! A class for generating full YAML input files from multiple data sources
class YamlWriter
{
public:
    YamlWriter();

    //! Include a phase definition for the specified Solution object
    void addPhase(shared_ptr<Solution> soln);

    std::string toYamlString() const;
    void toYamlFile(const std::string& filename) const;

    //! For output floating point values, set the maximum number of digits to
    //! the right of the decimal point. The default is 15 digits.
    void setPrecision(long int n) {
        m_float_precision = n;
    }

    //! By default user-defined data present in the input is preserved on
    //! output. This method can be used to skip output of user-defined data
    //! fields which are not directly used by Cantera.
    void skipUserDefined(bool skip=true) {
        m_skip_user_defined = skip;
    }

protected:
    std::vector<shared_ptr<Solution>> m_phases;

    //! @see setPrecision()
    long int m_float_precision;

    //! @see skipUserDefined()
    bool m_skip_user_defined;
};

}

#endif
