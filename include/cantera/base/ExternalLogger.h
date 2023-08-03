// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EXTERNAL_LOGGER_H
#define CT_EXTERNAL_LOGGER_H

#include "logger.h"

namespace Cantera {

//! Logger that delegates to an external source via a callback to produce log output.
//! @ingroup logGroup
class ExternalLogger : public Logger
{
public:
    explicit ExternalLogger(LogCallback writer) {
        if (writer == nullptr) {
            throw CanteraError("ExternalLogger::ExternalLogger",
                "Argument “writer” must not be null!");
        }

        m_writer = writer;
    }

    void write(const string& msg) override {
        m_writeBuffer.append(msg);

        if (!m_writeBuffer.empty() && m_writeBuffer.back() == '\n') {
            // This is a bit strange, but the terminal new line is interpreted to mean
            // “end of message”, so we want to pop it from the message itself.
            // The other side of the logger will be in charge of deciding whether
            // “messages” will have a terminal new line or not.
            m_writeBuffer.pop_back();

            m_writer(LogLevel::INFO, "Info", m_writeBuffer.c_str());

            m_writeBuffer.erase();
        }
    }

    void writeendl() override {
        m_writer(LogLevel::INFO, "Info", m_writeBuffer.c_str());

        m_writeBuffer.erase();
    }

    void warn(const string& warning, const string& msg) override {
        m_writer(LogLevel::WARN, warning.c_str(), msg.c_str());
    }

    void error(const string& msg) override {
        m_writer(LogLevel::ERROR, "Error", msg.c_str());
    }

private:
    string m_writeBuffer;

    LogCallback m_writer = nullptr;
};

}

#endif
