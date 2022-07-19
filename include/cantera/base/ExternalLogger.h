// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EXTERNAL_LOGGER_H
#define CT_EXTERNAL_LOGGER_H

#include "logger.h"

namespace Cantera {

class ExternalLogger : public Logger
{
public:
    explicit ExternalLogger(Writer writer) {
        if (writer == nullptr) {
            throw CanteraError("ExternalLogger::ExternalLogger",
                "Argument “writer” must not be null!");
        }

        m_writer = writer;
    }

    void write(const std::string& msg) override {
        m_writeBuffer << msg;
    }

    void writeendl() override {
        m_writer(LogLevel::INFO, "Info", m_writeBuffer.str().c_str());

        m_writeBuffer.clear();
    }

    void warn(const std::string& warning, const std::string& msg) override {
        m_writer(LogLevel::WARN, warning.c_str(), msg.c_str());
    }

    void error(const std::string& msg) override {
        m_writer(LogLevel::ERROR, "Error", msg.c_str());
    }

private:
    std::stringstream m_writeBuffer;

    Writer m_writer = nullptr;
};

}

#endif
