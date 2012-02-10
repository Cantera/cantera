/**
 *  @file CKParser.h
 *
 */

// Copyright 2001  California Institute of Technology

#ifndef CKR_CKPARSER_H
#define CKR_CKPARSER_H

#include <fstream>
#include <string>
#include <iostream>

#include "ckr_defs.h"
#include "Element.h"
#include "Species.h"
#include "Reaction.h"

namespace ckr
{


// typedefs

/// readability constants
//@{
const int NoThermoDatabase = 10;
const int HasTempRange = 11;
//@}


/// Exception class for syntax errors.
class CK_SyntaxError : public CK_Exception
{
public:
    CK_SyntaxError(std::ostream& f, const std::string& s, int linenum = -1);
    std::ostream& m_out;
};

/**
 * Chemkin mechanism file parser. For internal use by class CKReader.
 */
class CKParser
{

public:

    CKParser(std::string ckfile, std::ostream* log);
    CKParser(std::istream* infile, const std::string& fname,
             std::ostream* log);

    /// Destructor.
    ~CKParser() {}

    bool readElementSection(elementList& elements);
    bool readSpeciesSection(speciesList& species);
    bool readThermoSection(std::vector<std::string>& names,
                           speciesTable& speciesData, vector_fp& temp,
                           int& optionFlag, std::ostream& log);
    bool readReactionSection(const std::vector<std::string>& speciesNames,
                             std::vector<std::string>& elementNames,
                             reactionList& reactions, ReactionUnits& units);
    bool advanceToKeyword(const std::string& kw, const std::string& stop);
    bool verbose;
    bool debug;

    bool readNASA9ThermoSection(std::vector<std::string>& names,
                                speciesTable& species, vector_fp& temp,
                                int& optionFlag, std::ostream& log);

    void readNASA9ThermoRecord(Species& sp);

private:

    //! Local value of the line number being read
    /*!
     *  This is used for debug IO printout purposes
     */
    int m_line;

    std::string m_buf;
    std::string m_comment;

    //! This is the input file that is read
    /*!
     *  It's an istream
     */
    std::istream* m_ckfile;

    std::string m_ckfilename;

    //! Pointer to the ostream for writing debugging output log info
    std::ostream* m_log;

    bool m_nasafmt;

    //! Boolean indicating new NASA input file format
    /*!
     *  If this is true, a completely different input file parser is
     *  used.
     */
    bool m_nasa9fmt;

    char m_last_eol;
    void readThermoRecord(Species& sp);

    //!    Get a line from the input file, and return it in string s.
    /*!
     *  If the line contains a comment character (!), then return only the
     *  portion preceding it.  Non-printing characters are replaced by
     *  spaces.
     *
     *  The input file is m_ckfile, an istream.
     *
     *  @param s        On return, s contains the line read from the
     *                  input file.
     *  @param comment  On return, comment contains the text following the
     *                  comment character on the line, if any.
     */
    void getCKLine(std::string& s, std::string& comment);

    void putCKLine(std::string& s, std::string& comment);
    void missingAuxData(const std::string& kw);

    void checkSpeciesName(std::string spname);
};
}


#endif
