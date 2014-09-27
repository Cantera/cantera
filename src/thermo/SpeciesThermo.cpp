#include "cantera/thermo/SpeciesThermo.h"

namespace Cantera {

bool SpeciesThermo::ready(size_t nSpecies) {
    if (m_installed.size() < nSpecies) {
        return false;
    }
    for (size_t k = 0; k < nSpecies; k++) {
        if (!m_installed[k]) {
            return false;
        }
    }
    return true;
}

void SpeciesThermo::markInstalled(size_t k) {
    if (k >= m_installed.size()) {
        m_installed.resize(k+1, false);
    }
    m_installed[k] = true;
}

}
