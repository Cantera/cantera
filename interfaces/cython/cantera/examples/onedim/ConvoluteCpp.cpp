/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: ConvoluteCpp.cpp
//
// Description: This file contains implementation of five classes: FpArray3D, FpArray2D, ConvoluteCpp,
// CreateOutputFile, and CreateDiffusionFlameletsFile.  All of these classes implement capability originally
// written in python and converted to C++ for faster execution.  These C++ classes are accessed from
// python via Cython functions.
//
// Creation Date: 23-Feb-2017
//
// Revision History:
// SWCR            Date       Reason
// ------------------------------------------------------------------------------------------------------------------
//                                                                                                                   
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017 The University Of Dayton. All Rights Reserved.
//
// No part of this program may be photocopied, transferred, or otherwise reproduced in machine or
// human readable form without the prior written consent of The University Of Dayton.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES,
// INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// UNIVERSITY OF DAYTON OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
// OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  THE UNIVERSITY OF DAYTON
// HAS NO OBLIGATION TO SUPPORT THE SOFTWARE.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "ConvoluteCPP.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class FpArray3D
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// FpArrray3d constructor that initializes member variables that define array dimensions and array memory location.
/// </summary>
///
/// <param name="zlen">length of z dimension (major)</param>
/// <param name="ylen">length of y dimension</param>
/// <param name="xlen">length of x dimension (minor)</param>
/// <param name="data">pointer to contiguous python memory array</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FpArray3D::FpArray3D(int zlen, int ylen, int xlen, double *data):
    m_Zlen(zlen), m_Ylen(ylen), m_Xlen(xlen), m_data(data) {}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Gets an array data element via 3 dimensional indices (i.e. x,y,z).
/// </summary>
///
/// <param name="z">index of z dimension (major)</param>
/// <param name="y">index of y dimension</param>
/// <param name="x">index of x dimension (minor)</param>
///
/// <returns>array element</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FpArray3D::GetData(int z, int y, int x)
{
    return m_data[(z*m_Ylen*m_Xlen)+(y*m_Xlen)+x];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Sets an array data element via 3 dimensional indices (i.e. x,y,z).
/// </summary>
///
/// <param name="z">index of z dimension (major)</param>
/// <param name="y">index of y dimension</param>
/// <param name="x">index of x dimension (minor)</param>
/// <param name="value">new value of data element</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FpArray3D::SetData(int z, int y, int x, double value)
{
    m_data[(z*m_Ylen*m_Xlen) + (y*m_Xlen) + x] = value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Returns length of z dimension.
/// </summary>
///
/// <returns>length of z dimension</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FpArray3D::GetZlen()
{
    return m_Zlen;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Returns length of y dimension.
/// </summary>
///
/// <returns>length of y dimension</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FpArray3D::GetYlen()
{
    return m_Ylen;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Returns length of x dimension.
/// </summary>
///
/// <returns>length of x dimension</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FpArray3D::GetXlen()
{
    return m_Xlen;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class FpArray2D
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// FpArrray2d constructor that initializes member variables that define array dimensions and array memory location.
/// </summary>
///
/// <param name="ylen">length of y dimension</param>
/// <param name="xlen">length of x dimension (minor)</param>
/// <param name="data">pointer to contiguous python memory array</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FpArray2D::FpArray2D(int ylen, int xlen, double *data) :
   m_Ylen(ylen),  m_Xlen(xlen), m_data(data) {}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Gets an array data element via 2 dimensional indices (i.e. x,y).
/// </summary>
///
/// <param name="y">index of y dimension</param>
/// <param name="x">index of x dimension (minor)</param>
///
/// <returns>array element</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FpArray2D::GetData(int y, int x)
{
    return m_data[(y*m_Xlen) + x];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Sets an array data element via 2 dimensional indices (i.e. x,y).
/// </summary>
///
/// <param name="y">index of y dimension</param>
/// <param name="x">index of x dimension (minor)</param>
/// <param name="value">new value of data element</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FpArray2D::SetData(int y, int x, double value)
{
    m_data[(y*m_Xlen) + x] = value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Returns length of y dimension.
/// </summary>
///
/// <returns>length of y dimension</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FpArray2D::GetYlen()
{
    return m_Ylen; 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Returns length of x dimension.
/// </summary>
///
/// <returns>length of x dimension</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FpArray2D::GetXlen()
{ 
    return m_Xlen;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class ConvoluteCpp
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Default ConvoluteCpp constructor that initializes member variables that define pointers to array locations to 
/// nullptr.
/// </summary>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ConvoluteCpp::ConvoluteCpp():
    m_Betas(nullptr), m_Betas2(nullptr), m_DictEntry(nullptr), m_Zz(nullptr), m_rho(false)
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Alternate ConvoluteCpp constructor that initializes member variables that define pointers to array locations to 
/// nullptr.  Sets rho flag to true if specified name matches.
/// </summary>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ConvoluteCpp::ConvoluteCpp(char *name):
    m_Betas(nullptr), m_Betas2(nullptr), m_DictEntry(nullptr), m_Zz(nullptr)
{
    m_rho = (strcmp(name, "rho") == 0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Creates a 3 dimensional array mapping to python Betas contiguous memory pointer.
/// </summary>
///
/// <param name="z">length of z dimension (major)</param>
/// <param name="y">length of y dimension</param>
/// <param name="x">length of x dimension (minor)</param>
/// <param name="data">pointer to contiguous python memory array</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConvoluteCpp::SetBetas(int zlen, int ylen, int xlen, double *data)
{
    if (m_Betas != nullptr)
    {
        delete m_Betas;
    }
    m_Betas = new FpArray3D(zlen, ylen, xlen, data);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Creates a 3 dimensional array mapping to python Betas2 contiguous memory pointer.
/// </summary>
///
/// <param name="z">length of z dimension (major)</param>
/// <param name="y">length of y dimension</param>
/// <param name="x">length of x dimension (minor)</param>
/// <param name="data">pointer to contiguous python memory array</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConvoluteCpp::SetBetas2(int zlen, int ylen, int xlen, double *data)
{
    if (m_Betas2 != nullptr)
    {
        delete m_Betas2;
    }
    m_Betas2 = new FpArray3D(zlen, ylen, xlen, data);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Creates a 2 dimensional array mapping to python dictionary entry contiguous memory pointer.
/// </summary>
///
/// <param name="y">length of y dimension</param>
/// <param name="x">length of x dimension (minor)</param>
/// <param name="data">pointer to contiguous python memory array</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConvoluteCpp::SetDictEntry(int ylen, int xlen, double *data)
{
    if (m_DictEntry != nullptr)
    {
        delete m_DictEntry;
    }
    m_DictEntry = new FpArray2D(ylen, xlen, data);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Creates a 2 dimensional array mapping to python Zz contiguous memory pointer.
/// </summary>
///
/// <param name="y">length of y dimension</param>
/// <param name="x">length of x dimension (minor)</param>
/// <param name="data">pointer to contiguous python memory array</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConvoluteCpp::SetZz(int ylen, int xlen, double *data)
{
    if (m_Zz != nullptr)
    {
        delete m_Zz;
    }
    m_Zz = new FpArray2D(ylen, xlen, data);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Saves the length and a pointer to the python Zx contiguous memory pointer.
/// </summary>
///
/// <param name="x">length of x dimension (minor)</param>
/// <param name="data">pointer to contiguous python memory array</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConvoluteCpp::SetZx(int xlen, double *data)
{
    m_Zx = data;
    m_Xlen = xlen;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Saves the length and a pointer to the python Zy contiguous memory pointer.
/// </summary>
///
/// <param name="y">length of y dimension (minor)</param>
/// <param name="data">pointer to contiguous python memory array</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConvoluteCpp::SetZy(int ylen, double *data)
{
    m_Zy = data;
    m_Ylen = ylen;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Computes Zx, Zy, and Zz for rho using Betas, Betas2, and dictEntry.
/// </summary>
///
/// <param name="cvar">value of c variance</param>
/// <param name="zvar">value of z variance</param>
/// <param name="zmeanIndex">index into zMean</param>
/// <param name="zvarIndex">index into zVar</param>
/// <param name="cmeanIndex">index into cMean</param>
/// <param name="cvarIndex">index into cVar</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConvoluteCpp::EvalRho(double cvar, double zvar, int zmeanIndex, int zvarIndex, int cmeanIndex, int cvarIndex)
{
    for (int i = 0; i < m_Xlen; i++)
    {
        double beta = m_Betas->GetData(i, zmeanIndex, zvarIndex);
        if (cvar == 0)
        {
            // Model "A": treats C as a Dirac delta
            m_Zx[i] = beta / m_DictEntry->GetData(i, cmeanIndex);
        }
        else
        {
            for (int j = 0; j < m_Ylen; j++)
            {
                // Attempt to construct the integrand
                double beta2 = m_Betas2->GetData(j, cmeanIndex, cvarIndex);
                if (zvar == 0)
                {
                    // Model "A": treats Z as a Dirac delta
                    m_Zy[j] = beta2 / m_DictEntry->GetData(zmeanIndex, j);
                }
                else
                {
                    // Model "A": treats Z as a Dirac delta
                    double temp = beta * beta2 / m_DictEntry->GetData(i, j);
                    m_Zz->SetData(i, j, temp);
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Computes Zx, Zy, and Zz for non-rho entries using Betas, Betas2, and dictEntry.
/// </summary>
///
/// <param name="cvar">value of c variance</param>
/// <param name="zvar">value of z variance</param>
/// <param name="zmeanIndex">index into zMean</param>
/// <param name="zvarIndex">index into zVar</param>
/// <param name="cmeanIndex">index into cMean</param>
/// <param name="cvarIndex">index into cVar</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConvoluteCpp::EvalNormal(double cvar, double zvar, int zmeanIndex, int zvarIndex, int cmeanIndex, int cvarIndex)
{
    for (int i = 0; i < m_Xlen; i++)
    {
        // Attempt to construct the integrand
        double beta = m_Betas->GetData(i, zmeanIndex, zvarIndex);
        if (cvar == 0.0)
        {
            // Model "A": treats C as a Dirac delta
            m_Zx[i] = beta * m_DictEntry->GetData(i, cmeanIndex);
        }
        else
        {
            for (int j = 0; j < m_Ylen; j++)
            {
                double beta2 = m_Betas2->GetData(j, cmeanIndex, cvarIndex);
                if (zvar == 0.0)
                {
                    // Model "A": treats Z as a Dirac delta
                    m_Zy[j] = beta2 * m_DictEntry->GetData(zmeanIndex, j);
                }
                else
                {
                    // Model "A": treats Z as a Dirac delta
                    double temp = beta * beta2 * m_DictEntry->GetData(i, j);
                    m_Zz->SetData(i, j, temp);
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Computes Zx, Zy, and Zz using Betas, Betas2, and dictEntry.
/// </summary>
///
/// <param name="cvar">value of c variance</param>
/// <param name="zvar">value of z variance</param>
/// <param name="zmeanIndex">index into zMean</param>
/// <param name="zvarIndex">index into zVar</param>
/// <param name="cmeanIndex">index into cMean</param>
/// <param name="cvarIndex">index into cVar</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConvoluteCpp::Eval(double cvar, double zvar, int zmeanIndex, int zvarIndex, int cmeanIndex, int cvarIndex)
{
    if (m_rho == true)
    {
        EvalRho(cvar, zvar, zmeanIndex, zvarIndex, cmeanIndex, cvarIndex);
    }
    else
    {
        EvalNormal(cvar, zvar, zmeanIndex, zvarIndex, cmeanIndex, cvarIndex);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Test function to insure Zz array storage mapped correctly.
/// </summary>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void  ConvoluteCpp::EvalZzArray()
{
    for (int i = 0; i < m_Zz->GetYlen(); i++)
    {
        for (int j = 0; j < m_Zz->GetXlen(); j++)
        {
            m_Zz->SetData(i, j, i * 1000 + j);
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Test function to insure Zx array storage mapped correctly.
/// </summary>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void  ConvoluteCpp::EvalZxArray()
{
    for (int i = 0; i < m_Xlen; i++)
    {
        m_Zx[i] += i;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Test function to insure Betas array storage mapped correctly.
/// </summary>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void  ConvoluteCpp::EvalBetaArray()
{    
    for (int k = 0; k < m_Betas->GetZlen(); k++)
    {
        for (int i = 0; i < m_Betas->GetYlen(); i++)
        {
            for (int j = 0; j < m_Betas->GetXlen(); j++)
            {
                m_Betas->SetData(k, i, j, k*10000 + i * 100 + j);
            }
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class CreateOutputFile
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Default CreateOutputFile constructor that initializes member variables to nullptr.
/// </summary>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CreateOutputFile::CreateOutputFile()
{
    m_File = nullptr;
    m_OutputData = nullptr;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Creates/Opens the output file.  Creates the outputData array mapping.  Initializes the Zmean array mapping.
/// Initializes local variable indicating if the file will be using a compressed format.
/// </summary>
///
/// <param name="filename">output file name</param>
/// <param name="ylen">length of y dimension of OutputData mapping</param>
/// <param name="xlen">length of x dimension of OutputData mapping</param>
/// <param name="data">pointer to python OutputData contiguous array</param>
/// <param name="zmeanlen">length of zMean array</param>
/// <param name="ZMean">pointer to python zMean contiguous array</param>
/// <param name="compressed">flag indicating if file should be created using a compressed format</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CreateOutputFile::CreateOutputFile(char *filename, int ylen, int xlen, double *data, int zmeanlen, double *ZMean, int compressed)
{
    m_OutputData = new FpArray2D(ylen, xlen, data);
    m_ZMean = ZMean;
    m_Zmeanlen = zmeanlen;
    fopen_s(&m_File, filename, "a");
    m_Compressed = compressed;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Writes a line of formatted text to the output file.
/// </summary>
///
/// <param name="cvar">value of c variance</param>
/// <param name="zvar">value of z variance</param>
/// <param name="cmean">value of c mean</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateOutputFile::Write(double cvar, double zvar, double cmean)
{   
    for (int j = 0; j < m_OutputData->GetXlen(); j++)
    {
        if (m_Compressed == 0)
        {
            fprintf(m_File, "%10.6e\t", m_ZMean[j]);
            fprintf(m_File, "%10.6e\t", cmean);
            fprintf(m_File, "%10.6e\t", zvar);
            fprintf(m_File, "%10.6e\t", cvar);
        }
        for (int i = 0; i < m_OutputData->GetYlen(); i++)
        {
            fprintf(m_File, "%10.6e\t", m_OutputData->GetData(i,j));
        }
        fprintf(m_File, "\n");
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Closes the output file.
/// </summary>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateOutputFile::Close()
{
    fclose(m_File);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class CreateDiffusionFlameletsFile
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Default CreateDiffusionFlameletsFile constructor that initializes  member variable to nullptr.
/// </summary>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CreateDiffusionFlameletsFile::CreateDiffusionFlameletsFile()
{
    m_File = nullptr;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Creates/Opens the output file.
/// </summary>
///
/// <param name="filename">output file name</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CreateDiffusionFlameletsFile::CreateDiffusionFlameletsFile(char *filename)
{
    fopen_s(&m_File, filename, "w");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Writes the output file header.
/// </summary>
///
/// <param name="maxC">max C</param>
/// <param name="numSpecies">number of species</param>
/// <param name="gridpoints">number of gridpoints</param>
/// <param name="pressure">value of pressure</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateDiffusionFlameletsFile::WriteHeader(double maxC, int numSpecies, int gridpoints, double pressure)
{
    fprintf(m_File, "HEADER\n");
    fprintf(m_File, "PREMIX_CHI 0.0000\n");
    fprintf(m_File, "C  %6.6f\n", maxC);
    fprintf(m_File, "NUMOFSPECIES %3d\n", numSpecies);
    fprintf(m_File, "GRIDPOINTS	%3d\n", gridpoints);
    fprintf(m_File, "PRESSURE %6.2f\n", pressure);
    fprintf(m_File, "BODY");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Writes a section of text for the specified species.
/// </summary>
///
/// <param name="len">number of data items to write out</param>
/// <param name="data">pointer to data item array</param>
/// <param name="name">species name</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateDiffusionFlameletsFile::WriteData(int len, double *data, char *name)
{
    fprintf(m_File, "\n%s\n", name);
    for (int j = 0; j < len; j++)
    {
        fprintf(m_File, "%15.9e ", data[j]);
        if (((j + 1) % 5 == 0) && ((j + 1) != len))
        {
            fprintf(m_File, "\n");
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Writes a section of text for the specified species.  
/// Makes sure data is greater then minimum value or outputs 0.0.
/// </summary>
///
/// <param name="len">number of data items to write out</param>
/// <param name="data">pointer to data item array</param>
/// <param name="name">species name</param>
/// <param name="minValue">minimum data value allowed</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateDiffusionFlameletsFile::WriteDataMin(int len, double *data, char *name, double minValue)
{
    fprintf(m_File, "\n%s\n", name);
    for (int j = 0; j < len; j++)
    {
        double outdata = data[j];
        if (outdata < minValue)
        {
            outdata = 0.0;
        }
        fprintf(m_File, "%15.9e ", outdata);
        if (((j + 1) % 5 == 0) && ((j + 1) != len))
        {
            fprintf(m_File, "\n");
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Writes a PREMIX-CDOT section of text when flow is laminar.  
/// </summary>
///
/// <param name="ylen">ylen for mapping into net production rate array</param>
/// <param name="xlen">xlen for mapping into net production rate array</param>
/// <param name="nprData">pointer to net production rate contiguous python array</param>
/// <param name="mwlen">length of molecular weight python array</param>
/// <param name="molecularWeight">pointer to molecular weight contiguous python array</param>
/// <param name="ialen">length of index array</param>
/// <param name="indexArray">pointer to index contiguous python array</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateDiffusionFlameletsFile::WriteLaminarDataReaction(int ylen, int xlen, double *nprData,
                                                     int mwlen, double *molecularWeight, 
                                                     int ialen, int *indexArray)
{
    // create 2D array object
    FpArray2D netProductionRates(ylen, xlen, nprData);
    fprintf(m_File, "\n%s\n", "PREMIX_CDOT");
    fprintf(m_File, "%15.9e ", 0.0);

    for (int j = 1; j < xlen-1; j++)
    {
        double source = 0.0;
        for (int i = 0; i < ialen; i++)
        {
            source += netProductionRates.GetData(indexArray[i], j) * molecularWeight[indexArray[i]];
        }
        fprintf(m_File, "%15.9e ", source);
        if (((j + 1) % 5 == 0) && ((j + 1) != xlen))
        {
            fprintf(m_File, "\n");
        }
    }
    fprintf(m_File, "%15.9e ", 0.0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Writes a PREMIX-CDOT section of text when flow is turbulent.  
/// </summary>
///
/// <param name="ylen">ylen for mapping into net production rate array</param>
/// <param name="xlen">xlen for mapping into net production rate array</param>
/// <param name="nprData">pointer to net production rate contiguous python array</param>
/// <param name="mwlen">length of temperature python array</param>
/// <param name="molecularWeight">pointer to temperature contiguous python array</param>
/// <param name="ialen">length of index array</param>
/// <param name="indexArray">pointer to index contiguous python array</param>
/// <param name="pressure">value of pressure</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateDiffusionFlameletsFile::WriteTurbulentDataReaction(int ylen, int xlen, double *nprData,
    int tlen, double *temperature,
    int ialen, int *indexArray, double pressure)
{
    // create 2D array object
    FpArray2D netProductionRates(ylen, xlen, nprData);
    fprintf(m_File, "\n%s\n", "PREMIX_CDOT");
    fprintf(m_File, "%15.9e ", 0.0);

    const double UniversalGasConstant = 8314.5;
    for (int j = 1; j < xlen - 1; j++)
    {
        double source = 0.0;
        for (int i = 0; i < ialen; i++)
        {
            source += netProductionRates.GetData(indexArray[i], j);
        }
        double result = source / (pressure / (temperature[j] * UniversalGasConstant));
        fprintf(m_File, "%15.9e ", result);
        if (((j + 1) % 5 == 0) && ((j + 1) != xlen))
        {
            fprintf(m_File, "\n");
        }
    }
    fprintf(m_File, "%15.9e ", 0.0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Writes a data producation section of text using the given species name.  
/// </summary>
///
/// <param name="len">len for mapping into net production rate array</param>
/// <param name="nprData">pointer to net production rate contiguous python array</param>
/// <param name="molecularWeight">value of molecular weight</param>
/// <param name="name">species name</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateDiffusionFlameletsFile::WriteDataProduction(int len, double *nprData, double molecularWeight, char *name)
{
    // create 2D array object
    fprintf(m_File, "\n%s\n", name);
    for (int j = 0; j < len; j++)
    {
        double source = nprData[j] * molecularWeight;
        fprintf(m_File, "%15.9e ", source);
        if (((j + 1) % 5 == 0) && ((j + 1) != len))
        {
            fprintf(m_File, "\n");
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// <summary>
/// Writes  two carriage returns.  
/// </summary>
///
/// <param name="len">len for mapping into net production rate array</param>
/// <param name="nprData">pointer to net production rate contiguous python array</param>
/// <param name="molecularWeight">value of molecular weight</param>
/// <param name="name">species name</param>
///
/// <returns>none</returns>
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CreateDiffusionFlameletsFile::WriteCR()
{
    fprintf(m_File, "\n\n");
}

void CreateDiffusionFlameletsFile::Close()
{
    fclose(m_File);
}
