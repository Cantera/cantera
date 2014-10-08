#include "cantera/base/logger.h"

#include <fstream>

class fileLog: public Cantera::Logger
{
public:
    explicit fileLog(const std::string& fName) {
        m_fName = fName;
        m_fs.open(fName.c_str(), std::ios::out);
    }

    virtual void write(const std::string& msg) {
        m_fs << msg;
    }

    virtual void writeendl() {
        m_fs << std::endl;
    }

    virtual ~fileLog() {
        m_fs.close();
    }

    std::string m_fName;
    std::fstream m_fs;
};
