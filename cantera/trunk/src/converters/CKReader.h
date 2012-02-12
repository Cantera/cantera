/**
 *  @file CKReader.h
 *
 */

// Copyright 2001  California Institute of Technology

#ifndef CKR_CKRREADER_H
#define CKR_CKRREADER_H

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "CKParser.h"

#include <string>
#include <vector>

namespace ckr
{

class Group
{
public:

    /// Construct a new empty Group object
    Group() : name("<empty>"), index(-1) {}

    Group(const std::string& nm) : name(nm), index(-1) {}

    /// Destructor
    ~Group() {}

    std::string name;                 //!<  name
    int index;                   //!<  index number
    std::map<std::string, double> comp;    //!<  elemental composition

    bool operator==(const Group& g) const {
        return (name == g.name);
    }
    bool operator!=(const Group& g) const {
        return !(*this == g);
    }
};

/// a list (vector) of Groups
typedef std::vector<Group>      groupList;


/**
 *  Chemkin file reader class. Class CKReader parses and validates a file
 *  containing a description of a chemical reaction mechanism in Chemkin
 *  format. See the Examples section for examples of how CKReader is
 *  used in user programs.
 */

class CKReader
{
public:

    /**
     * Constructor. Construct a new CKReader instance. By default,
     * validation is enabled, as well as verbose output to the log file.
     */
    CKReader() : verbose(true), validate(true), debug(false) {}

    /// Destructor. Does nothing.
    ~CKReader() {}

    elementList   elements;     ///<  a list of Element objects
    speciesList   species;      ///<  a list of Species objects
    reactionList  reactions;    ///<  a list of Reaction objects
    groupList     groups;       ///<  a list of Groups
    speciesTable  speciesData;  ///<  a map from species names to Species objects
    ReactionUnits units;        ///<  reaction units

    /**
     * Read and optionally validate a Chemkin input file.
     * @param inputFile  path to the input file.
     * @param thermoDatabase  path to the species thermodynamic property database.
     * If no database is required, enter a null string.
     * @param logFile file to write logging and error messages to.
     * @return true if no errors encountered, false otherwise.
     */
    bool read(const std::string& inputFile,
              const std::string& thermoDatabase, const std::string& logFile);

    void write(string outputFile);  ///< not implemented.

    bool verbose;         ///<  print detailed messages to log file
    bool validate;        ///<  validate elements, species, and reaction
    bool debug;           ///<  enable debugging output

private:

    //    void validateElements(ostream& log);
    bool validateSpecies(ostream& log);    ///< validate the species.
    bool validateReactions(ostream& log);  ///< validate the reactions.
    bool writeReactions(ostream& log);
};


bool checkBalance(ostream& f, speciesTable& speciesData, reactionList& r,
                  vector<int>& unbalanced, double tolerance=1.0e-3);
bool checkThermo(ostream& f, speciesList& species, double tol);

bool filter(const string& infile, const string& database,
            const string& outfile, const vector<int>& species, const vector<int>& reactions);

}


#endif



