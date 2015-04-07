#ifndef CT_XML_WRITER
#define CT_XML_WRITER

// Note: this class is only used by TransportFactory, and is likely to
// go away a a future date.

#include <iostream>

namespace Cantera
{

//////////////////////  XML_Writer  //////////////////////////

class XML_Writer
{
public:
    XML_Writer(std::ostream& output_) :
        m_s(output_), _indent("   "), _level(0) {}
    virtual ~XML_Writer() {}
    std::ostream& m_s;

    std::string _indent;
    int _level;

    std::ostream& output() {
        return m_s;
    }

    inline std::string XML_filter(const std::string& name) {
        int ns = static_cast<int>(name.size());
        std::string nm(name);
        for (int m = 0; m < ns; m++)
            if (name[m] == ' '
                    || name[m] == '('
                    || name[m] == ')') {
                nm[m] = '_';
            }
        return nm;
    }

    /**
    * XML_comment()
     *
     *  Add a comment element to the current XML output file
     *  Comment elements start with <!-- and end with -->
     *  Comments are indented according to the current lvl,
     *  _level
     *
     *  input
     * ---------
     *    s : Output stream containing the XML file
     *    comment : Reference to a string containing the comment
     */
    inline void XML_comment(std::ostream& s, const std::string& comment) {
        for (int n = 0; n < _level; n++) {
            s << _indent;
        }
        s << "<!--" << comment << "-->" << std::endl;
    }

    inline void XML_open(std::ostream& s, const std::string& tag, const std::string& p = "") {
        for (int n = 0; n < _level; n++) {
            s << _indent;
        }
        _level++;
        s << "<" << XML_filter(tag) << p << ">" << std::endl;
    }

    inline void XML_close(std::ostream& s, const std::string& tag) {
        _level--;
        for (int n = 0; n < _level; n++) {
            s << _indent;
        }
        s << "</" << XML_filter(tag) << ">" << std::endl;
    }

    template<class T>
    void XML_item(std::ostream& s, const std::string& tag, T value) {
        for (int n = 0; n < _level; n++) {
            s << _indent;
        }
        s << "<" << XML_filter(tag) << ">"
          << value << "</" << tag << ">" << std::endl;
    }

    template<class iter>
    void XML_writeVector(std::ostream& s, const std::string& indent,
                         const std::string& name, int vsize, iter v) {
        int ni;
        for (ni = 0; ni < _level; ni++) {
            s << _indent;
        }
        s << "<" << XML_filter(name) << "> ";

        int n = vsize;
        int n5 = n/5;
        int i, j, k = 0;
        for (j = 0; j < n5; j++) {
            for (i = 0; i < 5; i++) {
                s << v[k] << (k < n - 1 ? ", " : "");
                k++;
            }
            if (j < n5-1) {
                s << std::endl;
                for (ni = 0; ni < _level; ni++) {
                    s << _indent;
                }
            }
        }
        for (i = k; i < n; i++) {
            s << v[k] << (k < n - 1 ? ", " : "");
            k++;
        }

        s << "</" << XML_filter(name) << ">" << std::endl;
    }
};

}

#endif
