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

std::string FuncEval::getErrors() const {
    std::stringstream errs;
    for (const auto& err : m_errors) {
        errs << err;
        errs << "\n";
    }
    return errs.str();
}

int FuncEval::preconditioner_setup_nothrow(double t, double* y, double gamma)
{
    try {
        preconditionerSetup(t, y, gamma);
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
            writelog("FuncEval::preconditioner_setup_nothrow: unhandled exception:\n");
            writelog(err.what());
            writelogendl();
        }
        return -1; // unrecoverable error
    } catch (...) {
        std::string msg = "FuncEval::preconditioner_setup_nothrow: unhandled exception"
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

int FuncEval::preconditioner_solve_nothrow(double* rhs, double* output)
{
    try {
        preconditionerSolve(rhs, output); // perform preconditioner solve
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
            writelog("FuncEval::preconditioner_solve_nothrow: unhandled exception:\n");
            writelog(err.what());
            writelogendl();
        }
        return -1; // unrecoverable error
    } catch (...) {
        std::string msg = "FuncEval::preconditioner_solve_nothrow: unhandled exception"
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

}
