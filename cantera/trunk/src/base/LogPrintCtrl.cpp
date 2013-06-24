/**
 * @file LogPrintCtrl.cpp
 *    Declarations for a simple class that augments the logfile printing capabilities
 *   (see \ref Cantera::LogPrintCtrl).
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "LogPrintCtrl.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera
{

LogPrintCtrl::LogPrintCtrl(int Ndec) :
    m_ffss(0),
    m_pc(0)
{
    warn_deprecated("class LogPrintCtrl");
    m_ffss = new std::ostream(m_os.rdbuf());
    m_pc = new PrintCtrl(*m_ffss, Ndec);
}

LogPrintCtrl::~LogPrintCtrl()
{
    delete m_pc;
    delete m_ffss;
}

void LogPrintCtrl::pr_de_c10(const double din, int p, const int wMin,
                             const int wMax)
{
    m_pc->pr_de_c10(din, p, wMin, wMax);
    writelog(m_os.str());
    m_os.str("");
}

void LogPrintCtrl::pr_de(const double d, int sigDigIn, const int wMinIn,
                         const int wMaxIn)
{
    m_pc->pr_de(d, sigDigIn, wMinIn, wMaxIn);
    writelog(m_os.str());
    m_os.str("");
}

double LogPrintCtrl::cropAbs10(const double d, int Ndec) const
{
    return m_pc->cropAbs10(d, Ndec);
}

double LogPrintCtrl::cropSigDigits(const double d, int nSig) const
{
    return m_pc->cropSigDigits(d, nSig);
}

int LogPrintCtrl::setNdec(int Ndec)
{
    return m_pc->setNdec(Ndec);
}

int LogPrintCtrl::setSigDigits(int nSigDigits)
{
    return m_pc->setSigDigits(nSigDigits);
}

int LogPrintCtrl::setWmin(int wmin)
{
    return m_pc->setWmin(wmin);
}

int LogPrintCtrl::setWmax(int wmax)
{
    return m_pc->setWmax(wmax);
}

}
