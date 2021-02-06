/**
 *  @file plots.h Contains declarations for utility functions for outputing to
 *       plotting programs.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLOTS_H
#define CT_PLOTS_H

#include "cantera/base/Array.h"
#include "cantera/base/ctexceptions.h"

#include <fstream>

namespace Cantera
{
//! Write a Plotting file
/*!
 * @param fname      Output file name
 * @param fmt        Either TEC or XL or CSV
 * @param plotTitle  Title of the plot
 * @param names      vector of variable names
 * @param data       N x M data array.
 *                     data(n,m) is the m^th value of the n^th variable.
 */
void writePlotFile(const std::string& fname, const std::string& fmt,
                   const std::string& plotTitle, const std::vector<std::string> &names,
                   const Array2D& data);

//! Write a Tecplot data file.
/*!
 * @param s        output stream
 * @param title    plot title
 * @param names    vector of variable names
 * @param data      N x M data array.
 *                 data(n,m) is the m^th value of the n^th variable.
 */
void outputTEC(std::ostream& s, const std::string& title,
               const std::vector<std::string>& names,
               const Array2D& data);

//! Write an Excel spreadsheet in 'csv' form.
/*!
 * @param s          output stream
 * @param title      plot title
 * @param names      vector of variable names
 * @param data       N x M data array.
 *                        data(n,m) is the m^th value of the n^th variable.
 */
void outputExcel(std::ostream& s, const std::string& title,
                 const std::vector<std::string>& names,
                 const Array2D& data);
}

#endif
