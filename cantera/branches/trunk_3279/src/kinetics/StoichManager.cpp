#include "cantera/kinetics/importKinetics.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

#include <cstring>

using namespace ctml;
using namespace std;

namespace Cantera
{

    //! Add a single reaction to the list of reactions that this
    //! stoichiometric manager object handles.
    /*!
     * @param rxn  Reaction index of the current reaction. This is used
     *             as an index into vectors which have length n_total_rxn.
     * @param k    This is a vector of integer values specifying the
     *             species indices. The length of this vector species
     *             the number of different species in the description.
     *             The value of the entries are the species indices.
     *             These are used as indexes into vectors which have
     *             length n_total_species.
     *  @param order This is a vector of the same length as vector k.
     *         The order is used for the routine power(), which produces
     *         a power law expression involving the species vector.
     *  @param stoich  This is used to handle fractional stoichiometric coefficients
     *                 on the product side of irreversible reactions.
     */
    void StoichManagerN::add(size_t rxn, const std::vector<size_t>& k, const vector_fp& order,
             const vector_fp& stoich) {
	//printf ("add called\n");
        if (order.size() != k.size()) {
           throw CanteraError("StoichManagerN::add()", "size of order and species arrays differ");    
        }
        if (stoich.size() != k.size()) {
           throw CanteraError("StoichManagerN::add()", "size of stoich and species arrays differ");    
        }
        bool frac = false;
        for (size_t n = 0; n < stoich.size(); n++) {
            if (fmod(stoich[n], 1.0) || fmod(order[n], 1.0)) {
                frac = true;
                break;
            }
        }
        if (frac || k.size() > 3) {
            m_cn_list.push_back(C_AnyN(rxn, k, order, stoich));
        } else {
            // Try to express the reaction with unity stoichiometric
            // coefficients (by repeating species when necessary) so that the
            // simpler 'multiply' function can be used to compute the rate
            // instead of 'power'.
            std::vector<size_t> kRep;
            for (size_t n = 0; n < k.size(); n++) {
                for (size_t i = 0; i < stoich[n]; i++)
                    kRep.push_back(k[n]);
            }

            switch (kRep.size()) {
            case 1:
                m_c1_list.push_back(C1(rxn, kRep[0]));
                break;
            case 2:
                m_c2_list.push_back(C2(rxn, kRep[0], kRep[1]));
                break;
            case 3:
                m_c3_list.push_back(C3(rxn, kRep[0], kRep[1], kRep[2]));
                break;
            default:
                m_cn_list.push_back(C_AnyN(rxn, k, order, stoich));
            }
        }
    }
}

