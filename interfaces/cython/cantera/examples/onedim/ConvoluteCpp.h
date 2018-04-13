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

#pragma once

#include <vector>
#include <stdio.h>

//! Vector of doubles.
typedef std::vector<double> vector_fp;
//! Vector of ints
typedef std::vector<int> vector_int;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// 
/// <version>1.0.0</version>
/// <author>Robert B. Olding</author>
///
/// <summary>
/// The FpArray3D class provides access to a three dimensional array mapped to a contiguous one dimensional 
/// array of python memory.
/// </summary>
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class FpArray3D
{
public:
    FpArray3D(int zlen, int ylen, int xlen, double *data);
    inline double GetData(int z, int y, int x);
    inline void SetData(int z, int y, int x, double value);
    inline int GetZlen();
    inline int GetYlen();
    inline int GetXlen();

private:
    const int m_Xlen;
    const int m_Ylen;
    const int m_Zlen;
    double *m_data;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// 
/// <version>1.0.0</version>
/// <author>Robert B. Olding</author>
///
/// <summary>
/// The FpArray2D class provides access to a two dimensional array mapped to a contiguous one dimensional 
/// array of python memory.
/// </summary>
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class FpArray2D
{
public:
    FpArray2D(int ylen, int xlen, double *data);
    inline double GetData(int y, int x);
    inline void SetData(int y, int x, double value);
    inline int GetYlen();
    inline int GetXlen();

private:
    const int m_Xlen;
    const int m_Ylen;
    double *m_data;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// 
/// <version>1.0.0</version>
/// <author>Robert B. Olding</author>
///
/// <summary>
/// The ConvoluteCpp class implements some of the slower python functions implemented in the Convolution routine.
/// </summary>
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ConvoluteCpp
{

public:
    ConvoluteCpp();
    ConvoluteCpp(char *name);

    void SetBetas(int zlen, int ylen, int xlen, double *data);
    void SetBetas2(int zlen, int ylen, int xlen, double *data);
    void SetDictEntry(int ylen, int xlen, double *data);
    void SetZz(int ylen, int xlen, double *data);
    void SetZx(int xlen, double *data);
    void SetZy(int ylen, double *data);

    void EvalNormal(double cvar, double zvar, int zmeanIndex, int zvarIndex, int cmeanIndex, int cvarIndex);
    void EvalRho(double cvar, double zvar, int zmeanIndex, int zvarIndex, int cmeanIndex, int cvarIndex);
    void Eval(double cvar, double zvar, int zmeanIndex, int zvarIndex, int cmeanIndex, int cvarIndex);

    void EvalZzArray();
    void EvalZxArray();
    void EvalBetaArray();

private:
    FpArray3D *m_Betas;
    FpArray3D *m_Betas2;
    FpArray2D *m_DictEntry;
    FpArray2D *m_Zz;
    double *m_Zx;
    double *m_Zy;
    int m_Xlen;
    int m_Ylen;
    bool m_rho;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// 
/// <version>1.0.0</version>
/// <author>Robert B. Olding</author>
///
/// <summary>
/// The CreateOutputFile class implements fast file output via fprintf for the PdfScript python routine.
/// </summary>
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class CreateOutputFile
{
public:
    CreateOutputFile();
    CreateOutputFile(char *filename, int ylen, int xlen, double *data, int zmeanlen, double *ZMean, int m_Compressed);
    void Write(double cvar, double zvar, double cmean);
    void Close();

private:
    FpArray2D *m_OutputData;
    FILE *m_File;
    int m_Zmeanlen;
    double* m_ZMean;
    int m_Compressed;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// 
/// <version>1.0.0</version>
/// <author>Robert B. Olding</author>
///
/// <summary>
/// The CreateDiffusionFlameletsFile class implements fast file output via fprintf for the DiffusionFlaemelts
/// python routine.
/// </summary>
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class CreateDiffusionFlameletsFile
{
public:
    CreateDiffusionFlameletsFile();
    CreateDiffusionFlameletsFile(char *filename);
    void WriteHeader(double maxC, int numSpecies, int gridpoints, double pressuree);
    void WriteData(int len, double *data, char *nam);
    void WriteDataMin(int len, double *data, char *name, double minValue);
    void WriteLaminarDataReaction(int ylen, int xlen, double *nprData,
        int mwlen, double *molecularWeight,
        int ialen, int *indexArray);
    void WriteTurbulentDataReaction(int ylen, int xlen, double *nprData,
        int tlen, double *temperature,
        int ialen, int *indexArray, double pressure);
    void WriteDataProduction(int len, double *nprData, double molecularWeight, char *name);
    void WriteCR();
    void Close();

private:
    FILE *m_File;
};
