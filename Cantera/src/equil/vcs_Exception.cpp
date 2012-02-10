
#include "vcs_Exception.h"

namespace VCSnonideal
{

vcsError::vcsError(std::string proc, std::string msg, int errorCode) :
    m_proc(proc),
    m_msg(msg),
    m_errorCode(errorCode)
{

}

}
