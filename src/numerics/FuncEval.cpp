#include "cantera/numerics/FuncEval.h"
#include <sstream>

namespace Cantera
{

FuncEval::FuncEval()
    : m_suppress_errors(false)
{
}

int FuncEval::eval_nothrow(double t, double* y, double* ydot)
{
    try {
        eval(t, y, ydot, m_sens_params.data());
    } catch (CanteraError& err) {
        if (suppressErrors()) {
            m_errors.push_back(err.what());
        } else {
            writelog(err.what());
        }
        return 1; // possibly recoverable error
    } catch (std::exception& err) {
        if (suppressErrors()) {
            m_errors.push_back(err.what());
        } else {
            writelog("FuncEval::eval_nothrow: unhandled exception:\n");
            writelog(err.what());
            writelogendl();
        }
        return -1; // unrecoverable error
    } catch (...) {
        std::string msg = "FuncEval::eval_nothrow: unhandled exception"
            " of unknown type\n";
        if (suppressErrors()) {
            m_errors.push_back(msg);
        } else {
            writelog(msg);
        }
        return -1; // unrecoverable error
    }
    return 0; // successful evaluation
}

int FuncEval::eval_nothrow(double t, double* y, double* ydot, double* r)
{
    try {
        eval(t, y, ydot, m_sens_params.data(), r);
    } catch (CanteraError& err) {
        if (suppressErrors()) {
            m_errors.push_back(err.what());
        } else {
            writelog(err.what());
        }
        return 1; // possibly recoverable error
    } catch (std::exception& err) {
        if (suppressErrors()) {
            m_errors.push_back(err.what());
        } else {
            writelog("FuncEval::eval_nothrow: unhandled exception:\n");
            writelog(err.what());
            writelogendl();
        }
        return -1; // unrecoverable error
    } catch (...) {
        std::string msg = "FuncEval::eval_nothrow: unhandled exception"
            " of unknown type\n";
        if (suppressErrors()) {
            m_errors.push_back(msg);
        } else {
            writelog(msg);
        }
        return -1; // unrecoverable error
    }
    return 0; // successful evaluation
}

std::string FuncEval::getErrors() const {
    std::stringstream errs;
    for (const auto& err : m_errors) {
        errs << err;
        errs << "\n";
    }
    return errs.str();
}

void FuncEval::eval(double t, double* y, double* ydot, double* p)
{
    throw CanteraError("FuncEval::eval", "Not implemented!");
}

void FuncEval::eval(double t, double* y, double* ydot, double* p, double* residual)
{
    throw CanteraError("FuncEval::eval", "Not implemented!");
}

}
