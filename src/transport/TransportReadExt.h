#ifndef TRANSPORTREADEXT_H
#define TRANSPORTREADEXT_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
#include "StringFunct.h"


using namespace std;

namespace Cantera {


class CollisionPair
{
public:

    /**
     * Constructs a collision pair given the collision pair name.  Assumes the
     * net charge is neutral.
     */
    CollisionPair(const std::string& pair)
    {
        std::vector<std::string> tokens;
        std::string str = pair;

        String::eraseAll(str, ".");
        String::tokenize(str, tokens, "-");

        initialize(tokens[0], tokens[1]);

    }

    /**
     * Constructs a collision pair given two species objects.
     */
    CollisionPair(
        const string s1,
        const string s2)
    {
        initialize(s1, s2);
    }

    /**
     * Equality operator.
     */
    bool operator == (const CollisionPair& right) const {
        return m_collision == right.m_collision;
    }

    /**
     * Less than operator used to sort collision pairs.
     */
    bool operator < (const CollisionPair& right) const {
        return m_collision < right.m_collision;
    }

    /**
     * Returns the name of the collision pair.
     */
    const std::string& name() const {
        return m_collision;
    }

    /**
     * Returns the charge of this collision pair. Negative values are attractive
     * positive values are attractive and zero represents neutral collisions.
     */
 /*   int charge() const {
        return m_charge;
    }
*/

    /**
     * Returns the name of the first species in the pair.
     */
    const std::string& speciesName1() const {
        return m_s1;
    }

    /**
     * Returns the name of the second species in the pair.
     */
    const std::string& speciesName2() const {
        return m_s2;
    }

private:

    /**
     * Some common initialization stuff.
     */
    void initialize(const std::string &name1, const std::string &name2)
    {
        m_s1 = name1;
        m_s2 = name2;

        if (m_s1 < m_s2)
            m_collision = "." + m_s1 + "-" + m_s2 + ".";
        else
            m_collision = "." + m_s2 + "-" + m_s1 + ".";
    }

private:

    std::string m_collision;
    std::string m_s1;
    std::string m_s2;

}; // CollisionPair



// Read coefficients for transport properties from external files
class TransportReadExt
{
    public:
        TransportReadExt();
	void getCoefficients_omega11(string species1, string species2, vector<double> &coefficients);
	void getCoefficients_omega22(string species1, string species2, vector<double> &coefficients);
	bool isValidCollisionString(const std::string &str);
	void load_coeff();
	void getCoefficients_bstar(string species1, string species2, vector<double> &coefficients);

	vector<double> C_omega11;
	vector<double> C_omega22;
	vector<string> Pairs;
	vector<double> C_bstar;

};

};

#endif // TRANSPORTREADEXT_H


