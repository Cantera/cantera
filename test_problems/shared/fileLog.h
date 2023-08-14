#include "cantera/base/logger.h"

#include <fstream>

class fileLog: public Cantera::Logger
{
public:
    explicit fileLog(const std::string& fName) {
        m_fName = fName;
        m_fs.open(fName, std::ios::out);
    }

    void write(const std::string& msg) override {
        m_fs << msg;
    }

    void writeendl() override {
        m_fs << std::endl;
    }

    std::string m_fName;
    std::fstream m_fs;
};
