/**
 *  @file Group.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_RXNPATH_GROUP
#define CT_RXNPATH_GROUP

#include "cantera/base/ct_defs.h"

namespace Cantera
{

/**
 * Class Group is an internal class used by class ReactionPath. It
 * represents some subset of the atoms of a molecule.
 */
class Group
{
public:
    Group() : m_sign(-999) { }
    Group(size_t n) : m_sign(0) {
        m_comp.resize(n,0);
    }
    Group(const vector_int& elnumbers) :
        m_comp(elnumbers), m_sign(0) {
        validate();
    }
    Group(const std::vector<size_t>& elnumbers) :
        m_comp(elnumbers.size()), m_sign(0) {
        for (size_t i = 0; i < elnumbers.size(); i++) {
            m_comp[i] = int(elnumbers[i]);
        }
        validate();
    }
    Group(const Group& g) :
        m_comp(g.m_comp), m_sign(g.m_sign) { }
    Group& operator=(const Group& g) {
        if (&g != this) {
            m_comp = g.m_comp;
            m_sign = g.m_sign;
        }
        return *this;
    }
    virtual ~Group() {}

    /**
     * Decrement the atom numbers by those in group 'other'.
     */
    void operator-=(const Group& other) {
        for (size_t m = 0; m < m_comp.size(); m++) {
            m_comp[m] -= other.m_comp[m];
        }
        validate();
    }
    void operator+=(const Group& other) {
        for (size_t m = 0; m < m_comp.size(); m++) {
            m_comp[m] += other.m_comp[m];
        }
        validate();
    }
    void operator*=(int a) {
        for (size_t m = 0; m < m_comp.size(); m++) {
            m_comp[m] *= a;
        }
        validate();
    }
    bool operator==(const Group& other) const {
        for (size_t m = 0; m < m_comp.size(); m++) {
            if (m_comp[m] != other.m_comp[m]) {
                return false;
            }
        }
        return true;
    }
    friend Group operator-(const Group& g1, const Group& g2) {
        Group diff(g1);
        diff -= g2;
        return diff;
    }
    friend Group operator+(const Group& g1, const Group& g2) {
        Group sum(g1);
        sum += g2;
        return sum;
    }
    /*!
     * A group is 'valid' if all of its nonzero atom numbers have
     * the same sign, either positive or negative. This method
     * checks for this, and if the group is not valid it sets
     * m_sign to -999, and sets all atom numbers to zero.
     */
    void validate();

    /**
     * True if all non-zero atom numbers have the same sign.
     */
    bool valid() const {
        return (m_sign != -999);
    }
    bool operator!() const {
        return (m_sign == -999);
    }
    int sign() const {
        return m_sign;
    }
    size_t size() const {
        return m_comp.size();
    }

    /// Number of atoms in the group (>= 0)
    int nAtoms() const {
        int sum = 0;
        for (size_t m = 0; m < m_comp.size(); m++) {
            sum += std::abs(m_comp[m]);
        }
        return sum;
    }
    /// Number of atoms of element m (positive or negative)
    int nAtoms(size_t m) const {
        if (m_comp.empty()) {
            return 0;
        }
        return m_comp[m];
    }

    std::ostream& fmt(std::ostream& s, const std::vector<std::string>& esymbols) const;

    friend std::ostream& operator<<(std::ostream& s,
                                    const Group& g);

private:
    vector_int m_comp;
    int m_sign;
};

}

#endif
