// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TEMPERATURE_TABLE_H
#define TEMPERATURE_TABLE_H

#include <vector>
#include <algorithm>

using std::vector;

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/**
 *   This Class constructs a vector of temperature from which to make
 *   a table.
 */
class TemperatureTable
{

public:
    size_t NPoints;
    bool Include298;
    double Tlow;                 //!<   Min temperature for thermo data fit
    double Thigh;                //!<   Max temperature for thermo table
    double DeltaT;
    vector<double> T;
    size_t numAddedTs;
    vector<double> AddedTempVector;
public:
    /*
     * Default constructor for TemperatureTable()
     */
    TemperatureTable(const size_t nPts = 14,
                     const bool inc298 = true,
                     const double tlow = 300.,
                     const double deltaT = 100.,
                     const size_t numAdded = 0,
                     const double* addedTempVector = 0) :
        NPoints(nPts),
        Include298(inc298),
        Tlow(tlow),
        DeltaT(deltaT),
        T(0),
        numAddedTs(numAdded) {
        /****************************/
        AddedTempVector.resize(numAdded, 0.0);
        for (size_t i = 0; i < numAdded; i++) {
            AddedTempVector[i] = addedTempVector[i];
        }
        T.resize(NPoints, 0.0);
        double TCurrent = Tlow;
        for (size_t i = 0; i < NPoints; i++) {
            T[i] = TCurrent;
            TCurrent += DeltaT;
        }
        if (Include298) {
            T.push_back(298.15);
            NPoints++;
        }
        if (numAdded > 0) {
            T.resize(NPoints+numAdded, 0.0);
            for (size_t i = 0; i < numAdded; i++) {
                T[i+NPoints] = addedTempVector[i];
            }
            NPoints += numAdded;
        }
        std::sort(T.begin(), T.end());
    }
    /***********************************************************************/
    /***********************************************************************/
    /***********************************************************************/
    /*
     * Destructor
     */
    ~TemperatureTable() {}

    /***********************************************************************/
    /***********************************************************************/
    /***********************************************************************/
    /*
     *  Overloaded operator[]
     *
     *       return the array value in the vector
     */
    double operator[](const int i) {
        return T[i];
    }
    /***********************************************************************/
    /***********************************************************************/
    /***********************************************************************/
    /*
     *  size()
     */
    size_t size() {
        return NPoints;
    }
    /***********************************************************************/
    /***********************************************************************/
    /***********************************************************************/
    /*
     * Block assignment and copy constructors: not needed.
     */
private:
    TemperatureTable(const TemperatureTable&);
    TemperatureTable& operator=(const TemperatureTable&);
};
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
#endif
