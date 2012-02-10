/**
 *  @file PID_Controller.h
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_PID_H
#define CT_PID_H

namespace Cantera
{

class PID_Controller
{

public:

    /// Default constructor.
    PID_Controller() : m_v0(Undef), m_p(Undef), m_i(Undef), m_d(Undef),
        m_setpoint(Undef), m_last(Undef), m_time(Undef),
        m_xint(Undef), m_out(Undef), m_dt(Undef) {}

    /**
     * Copy constructor. Gains and setpoint are copied, but not
     * the internal parameters defining the state of the
     * controller.  Method 'reset' must be called for the copy
     * before using it.
     */
    PID_Controller(const PID_Controller& pid)
        : m_v0(pid.m_v0), m_p(pid.m_p), m_i(pid.m_i), m_d(pid.m_d),
          m_setpoint(pid.m_setpoint),
          m_last(Undef), m_time(Undef), m_xint(Undef) {}

    /**
     * Assignment operator. @see Copy constructor.
     */
    PID_Controller& operator=(const PID_Controller& pid) {
        if (this == &pid) {
            return *this;
        }
        m_v0 = pid.m_v0;
        m_p = pid.m_p;
        m_i = pid.m_i;
        m_d = pid.m_d;
        m_setpoint = pid.m_setpoint;
        m_last = Undef;
        m_time = Undef;
        m_xint = Undef;
        return *this;
    }

    /**
     * Reset the start time to time, and the current value of
     * the input to input. Sets the integrated error signal to zero.
     */
    void reset(doublereal time = 0.0, doublereal input = 0.0) {
        m_time = time;
        m_last = input;
        m_xint = 0.0;
        m_out = m_v0;
        m_dt = 1.0;
        m_maxerr = 0.0;
    }

    doublereal setpoint(doublereal y = Undef) {
        if (y != Undef) {
            m_setpoint = y;
        }
        return m_setpoint;
    }

    bool getGains(vector_fp& gains) {
        gains.resize(4);
        return getGains(4, gains.begin());
    }

    bool getGains(int n, doublereal* gains) {
        if (n < 4) {
            return false;
        }
        gains[0] = m_v0;
        gains[1] = m_p;
        gains[2] = m_i;
        gains[3] = m_d;
        return true;
    }

    bool setGains(const vector_fp& gains) {
        return setGains(int(gains.size()), gains.begin());
    }

    bool setGains(int n, const doublereal* gains) {
        if (n < 4) {
            return false;
        }
        m_v0 = gains[0];
        m_p = gains[1];
        m_i = gains[2];
        m_d = gains[3];
        if (m_p < 0.0 || m_i < 0.0 || m_d < 0.0) {
            return false;
        }
        return true;
    }

    void update(doublereal time, doublereal input) {
        if (time <= m_time) {
            return;
        }
        doublereal err = input - m_setpoint;
        if (fabs(err) > m_maxerr) {
            m_maxerr = fabs(err);
        }
        m_dt = time - m_time;
        m_xint += (0.5*(input + m_last) - m_setpoint) * m_dt;
        m_last = input;
        m_time = time;
        doublereal xdot = (input - m_last)/m_dt;
        m_out = m_v0 - m_p*(input - m_setpoint) - m_i*m_xint
                - m_d*xdot;
    }

    doublereal output(doublereal input) {
        return fmaxx(0.0,
                     m_out - (m_p + m_d/m_dt + 0.5*m_i*m_dt)*(input - m_last));
    }

    doublereal maxError() {
        return m_maxerr;
    }

    bool ready() {
        return (m_time != Undef
                && m_setpoint != Undef
                && m_v0 != Undef);
    }

protected:

    doublereal
    m_v0,
    m_p,
    m_i,
    m_d,
    m_setpoint,
    m_maxerr,
    m_last,
    m_time,
    m_xint,
    m_out,
    m_dt;
};

}

#endif

