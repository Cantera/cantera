/**
 * @file xml.cpp
 * Classes providing support for XML data files. These classes
 * implement only those aspects of XML required to read, write, and
 * manipulate CTML data files.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/xml.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"
#include "cantera/base/utilities.h"

#include <sstream>
#include <fstream>

using namespace std;

namespace Cantera
{
////////////////////// exceptions ////////////////////////////

//! Class representing a generic XML error condition
class XML_Error : public CanteraError
{
protected:
    //! Constructor
    /*!
     * Note, we don't actually post the error in this class. Therefore, this
     * class can't be used externally. Therefore, it's a protected constructor.
     *
     * @param file Name of the XML file being processed
     * @param line Number number where the error occurred.
     */
    XML_Error(const std::string& file, int line) {
        m_msg = fmt::format("Error in XML file '{}' at line {}.\n", file, line);
    }

    virtual std::string getMessage() const {
        return m_msg;
    }

    //! destructor
    virtual ~XML_Error() throw() {
    }

    //! String message for the error
    std::string m_msg;
};

//! Class representing a specific type of XML file formatting error
/*!
 * An XML tag is not matched
 */
class XML_TagMismatch : public XML_Error
{
public:
    //! Constructor
    /*!
     * An XML element must have the same opening and closing name.
     *
     * @param  opentag    String representing the opening of the XML bracket
     * @param  closetag   String representing the closing of the XML bracket
     * @param  filename   Name of the XML file being processed
     * @param  line       Line number where the error occurred.
     */
    XML_TagMismatch(const std::string& opentag, const std::string& closetag,
                    const std::string& filename, int line) :
        XML_Error(filename, line) {
        m_msg += fmt::format("<{}> paired with </{}>.\n", opentag, closetag);
    }

    virtual std::string getClass() const {
        return "XML_TagMismatch";
    }
};

//! Class representing a specific type of XML file formatting error
/*!
 * An XML_Node doesn't have a required child node
 */
class XML_NoChild : public XML_Error
{
public:
    //! Constructor
    /*!
     * An XML element doesn't have the required child node
     *
     * @param  p       XML_Node to write a string error message
     * @param  parent  Name of the parent node
     * @param  child   Name of the required child node
     * @param  filename Name of the XML file being processed
     * @param  line    Line number where the error occurred.
     */
    XML_NoChild(const XML_Node* p, const std::string& parent,
                std::string child, const std::string& filename, int line) :
        XML_Error(filename, line) {
        m_msg += fmt::format("The XML Node <{}> does not contain a required "
            "child node named <{}>.\nExisting children are named:\n",
            parent, child);
        for (auto cnode : p->children()) {
            m_msg += fmt::format("    <{}>\n", cnode->name());
        }
    }

    virtual std::string getClass() const {
        return "XML_NoChild";
    }
};

//////////////////// XML_Reader methods ///////////////////////

XML_Reader::XML_Reader(std::istream& input) :
    m_s(input),
    m_line(0)
{
}

void XML_Reader::getchr(char& ch)
{
    m_s.get(ch);
    if (ch == '\n') {
        m_line++;
    }
}

//! Find the first position of a character, q, in string, s, which is not
//! immediately preceded by the backslash character
/*!
 * @param s        Input string
 * @param q        Search for this character
 * @param istart   Defaults to 0
 */
static string::size_type findUnbackslashed(const std::string& s, const char q,
        std::string::size_type istart = 0)
{
    size_t icurrent = istart;
    while (true) {
        size_t iloc = s.find(q, icurrent);
        if (iloc == string::npos || iloc == 0) {
            return iloc;
        }
        char cm1 = s[iloc-1];
        if (cm1 == '\\') {
            if (iloc >= (s.size() -1)) {
                return string::npos;
            }
            icurrent = iloc + 1;
        } else {
            return iloc;
        }
    }
}

int XML_Reader::findQuotedString(const std::string& s, std::string& rstring) const
{
    const char q1 = '\'';
    const char q2 = '"';
    rstring = "";
    char qtype = ' ';
    string::size_type iloc1, iloc2, ilocStart = 0;
    iloc1 = findUnbackslashed(s, q1);
    iloc2 = findUnbackslashed(s, q2);
    if (iloc2 != string::npos) {
        ilocStart = iloc2;
        qtype = q2;
    }
    if (iloc1 != string::npos && iloc1 < ilocStart) {
        ilocStart = iloc1;
        qtype = q1;
    }
    if (qtype == ' ') {
        return 0;
    }

    iloc1 = findUnbackslashed(s, qtype, ilocStart+1);

    if (iloc1 == string::npos) {
        return 0;
    }

    // Define the return string by the two endpoints. Strip the surrounding
    // quotes as well
    rstring = s.substr(ilocStart + 1, iloc1 - 1);

    // Return the first character position past the quotes
    return static_cast<int>(iloc1)+1;
}

void XML_Reader::parseTag(const std::string& tag, std::string& name,
                          std::map<std::string, std::string>& attribs) const
{
    string s = trimCopy(tag);
    size_t iloc = s.find(' ');
    if (iloc != string::npos) {
        name = s.substr(0, iloc);
        s = trimCopy(s.substr(iloc+1,s.size()));
        if (s[s.size()-1] == '/') {
            name += "/";
        }

        // get attributes
        while (true) {
            iloc = s.find('=');
            if (iloc == string::npos) {
                break;
            }
            string attr = trimCopy(s.substr(0,iloc));
            if (attr == "") {
                break;
            }
            s = trimCopy(s.substr(iloc+1,s.size()));
            string val;
            iloc = findQuotedString(s, val);
            attribs[attr] = val;
            if (iloc != string::npos) {
                if (iloc < s.size()) {
                    s = trimCopy(s.substr(iloc,s.size()));
                } else {
                    break;
                }
            }
        }
    } else {
        name = s;
    }
}

std::string XML_Reader::readTag(std::map<std::string, std::string>& attribs)
{
    string tag = "";
    bool incomment = false;
    char ch = '-';
    while (true) {
        if (m_s.eof() || (getchr(ch), ch == '<')) {
            break;
        }
    }
    char ch1 = ' ', ch2 = ' ';
    while (true) {
        if (m_s.eof()) {
            tag = "EOF";
            break;
        }
        ch2 = ch1;
        ch1 = ch;
        getchr(ch);
        if (ch == '-') {
            if (ch1 == '-' && ch2 == '!') {
                incomment = true;
                tag = "-";
            }
        } else if (ch == '>') {
            if (incomment) {
                if (ch1 == '-' && ch2 == '-') {
                    break;
                }
            } else {
                break;
            }
        }
        if (isprint(ch)) {
            tag += ch;
        }
    }
    if (incomment) {
        attribs.clear();
        return tag;
    } else {
        string name;
        parseTag(tag, name, attribs);
        return name;
    }
}

std::string XML_Reader::readValue()
{
    string tag = "";
    char ch = '\n';
    bool front = true;
    while (true) {
        if (m_s.eof()) {
            break;
        }
        char lastch = ch;
        getchr(ch);
        if (ch == '\n') {
            front = true;
        } else if (ch != ' ') {
            front = false;
        }
        if (ch == '<') {
            m_s.putback(ch);
            break;
        }
        if (front && lastch == ' ' && ch == ' ') {
            ;
        } else {
            tag += ch;
        }
    }
    return trimCopy(tag);
}

//////////////////////////  XML_Node  /////////////////////////////////

XML_Node::XML_Node(const std::string& nm, XML_Node* const parent_) :
    m_name(nm),
    m_parent(parent_),
    m_root(0),
    m_locked(false),
    m_iscomment(false),
    m_linenum(0)
{
    if (!parent_) {
        m_root = this;
    } else {
        m_root = &parent_->root();
    }
}

XML_Node::XML_Node(const XML_Node& right) :
    m_parent(0),
    m_root(0),
    m_locked(false),
    m_iscomment(right.m_iscomment),
    m_linenum(right.m_linenum)
{
    m_root = this;
    m_name = right.m_name;
    m_value = right.m_value;
    right.copy(this);
}

XML_Node& XML_Node::operator=(const XML_Node& right)
{
    if (&right != this) {
        for (size_t i = 0; i < m_children.size(); i++) {
            if (m_children[i] && m_children[i]->parent() == this) {
                delete m_children[i];
                m_children[i] = 0;
            }
        }
        m_children.resize(0);
        right.copy(this);
    }
    return *this;
}

XML_Node::~XML_Node()
{
    if (m_locked) {
        writelog("XML_Node::~XML_Node: deleted a locked XML_Node: "+name());
    }
    for (size_t i = 0; i < m_children.size(); i++) {
        if (m_children[i] && m_children[i]->parent() == this) {
            delete m_children[i];
            m_children[i] = 0;
        }
    }
}

void XML_Node::clear()
{
    for (size_t i = 0; i < m_children.size(); i++) {
        if (m_children[i] && m_children[i]->parent() == this) {
            delete m_children[i];
            m_children[i] = 0;
        }
    }
    m_value.clear();
    m_childindex.clear();
    m_attribs.clear();
    m_children.clear();

    m_iscomment = false;
    m_linenum = 0;
}

XML_Node& XML_Node::mergeAsChild(XML_Node& node)
{
    m_children.push_back(&node);
    m_childindex.insert({node.name(), m_children.back()});
    node.setRoot(root());
    node.setParent(this);
    return *m_children.back();
}

XML_Node& XML_Node::addChild(const XML_Node& node)
{
    return mergeAsChild(*(new XML_Node(node)));
}

XML_Node& XML_Node::addChild(const std::string& sname)
{
    return mergeAsChild(*(new XML_Node(sname, this)));
}

XML_Node& XML_Node::addChild(const std::string& name, const std::string& value)
{
    XML_Node& c = addChild(name);
    c.addValue(value);
    return c;
}

XML_Node& XML_Node::addChild(const std::string& name, const doublereal value,
                             const std::string& fmt)
{
    XML_Node& c = addChild(name);
    c.addValue(value, fmt);
    return c;
}

void XML_Node::removeChild(const XML_Node* const node)
{
    auto i = find(m_children.begin(), m_children.end(), node);
    m_children.erase(i);
    m_childindex.erase(node->name());
}

void XML_Node::addComment(const std::string& comment)
{
    addChild("comment", comment);
}

void XML_Node::addValue(const std::string& val)
{
    m_value = val;
    if (m_name == "comment") {
        m_iscomment = true;
    }
}

void XML_Node::addValue(const doublereal val, const std::string& fmt)
{
    m_value = trimCopy(fmt::sprintf(fmt, val));
}

std::string XML_Node::value() const
{
    return m_value;
}

std::string XML_Node::value(const std::string& cname) const
{
    return child(cname).value();
}

std::string XML_Node::operator()(const std::string& cname) const
{
    return value(cname);
}

doublereal XML_Node::fp_value() const
{
    return fpValueCheck(m_value);
}

integer XML_Node::int_value() const
{
    return std::atoi(m_value.c_str());
}

void XML_Node::addAttribute(const std::string& attrib, const std::string& value)
{
    m_attribs[attrib] = value;
}

void XML_Node::addAttribute(const std::string& attrib,
                            const doublereal vvalue, const std::string& fmt)
{
    m_attribs[attrib] = fmt::sprintf(fmt, vvalue);
}

void XML_Node::addAttribute(const std::string& aattrib, const int vvalue)
{
    m_attribs[aattrib] = fmt::format("{}", vvalue);
}

void XML_Node::addAttribute(const std::string& aattrib, const size_t vvalue)
{
    m_attribs[aattrib] = fmt::format("{}", vvalue);
}

std::string XML_Node::operator[](const std::string& attr) const
{
    return attrib(attr);
}

std::string XML_Node::attrib(const std::string& attr) const
{
    return getValue<string,string>(m_attribs, attr, "");
}

std::map<std::string,std::string>& XML_Node::attribs()
{
    return m_attribs;
}

const std::map<std::string,std::string>& XML_Node::attribsConst() const
{
    return m_attribs;
}

void XML_Node::setLineNumber(const int n)
{
    m_linenum = n;
}

int XML_Node::lineNumber() const
{
    return m_linenum;
}

XML_Node* XML_Node::parent() const
{
    return m_parent;
}

XML_Node* XML_Node::setParent(XML_Node* const p)
{
    m_parent = p;
    return p;
}

bool XML_Node::hasChild(const std::string& ch) const
{
    return (m_childindex.find(ch) != m_childindex.end());
}

bool XML_Node::hasAttrib(const std::string& a) const
{
    return (m_attribs.find(a) != m_attribs.end());
}

std::string XML_Node::id() const
{
    if (hasAttrib("id")) {
        return attrib("id");
    }
    return "";
}

XML_Node& XML_Node::child(const size_t n) const
{
    return *m_children[n];
}

const std::vector<XML_Node*>& XML_Node::children() const
{
    return m_children;
}

size_t XML_Node::nChildren(const bool discardComments) const
{
    if (discardComments) {
        size_t count = 0;
        for (size_t i = 0; i < m_children.size(); i++) {
            XML_Node* xc = m_children[i];
            if (!(xc->isComment())) {
                count++;
            }
        }
        return count;
    }
    return m_children.size();
}

bool XML_Node::isComment() const
{
    return m_iscomment;
}

void XML_Node::_require(const std::string& a, const std::string& v) const
{
    if (hasAttrib(a) && attrib(a) == v) {
        return;
    }
    string msg="XML_Node "+name()+" is required to have an attribute named " + a +
               " with the value \"" + v +"\", but instead the value is \"" + attrib(a);
    throw CanteraError("XML_Node::_require", msg);
}

XML_Node* XML_Node::findNameID(const std::string& nameTarget,
                               const std::string& idTarget) const
{
    XML_Node* scResult = 0;
    std::string idattrib = id();
    if (name() == nameTarget && (idTarget == "" || idTarget == idattrib)) {
        return const_cast<XML_Node*>(this);
    }
    for (size_t n = 0; n < m_children.size(); n++) {
        XML_Node* sc = m_children[n];
        if (sc->name() == nameTarget) {
            if (idTarget == "") {
                return sc;
            }
            idattrib = sc->id();
            if (idTarget == idattrib) {
                return sc;
            }
        }
    }
    for (size_t n = 0; n < m_children.size(); n++) {
        XML_Node* sc = m_children[n];
        scResult = sc->findNameID(nameTarget, idTarget);
        if (scResult) {
            return scResult;
        }
    }
    return scResult;
}

XML_Node* XML_Node::findNameIDIndex(const std::string& nameTarget,
                                    const std::string& idTarget, const int index_i) const
{
    XML_Node* scResult = 0;
    std::string idattrib = id();
    std::string ii = attrib("index");
    std::string index_s = fmt::format("{}", index_i);
    int iMax = -1000000;
    if (name() == nameTarget && (idTarget == "" || idTarget == idattrib) && index_s == ii) {
        return const_cast<XML_Node*>(this);
    }
    for (size_t n = 0; n < m_children.size(); n++) {
        XML_Node* sc = m_children[n];
        if (sc->name() == nameTarget) {
            ii = sc->attrib("index");
            int indexR = atoi(ii.c_str());
            idattrib = sc->id();
            if ((idTarget == idattrib || idTarget == "") && index_s == ii) {
                return sc;
            }
            if (indexR > iMax) {
                scResult = sc;
                iMax = indexR;
            }
        }
    }
    return scResult;
}

XML_Node* XML_Node::findID(const std::string& id_, const int depth) const
{
    if (hasAttrib("id") && attrib("id") == id_) {
        return const_cast<XML_Node*>(this);
    }
    if (depth > 0) {
        for (size_t i = 0; i < nChildren(); i++) {
            XML_Node* r = m_children[i]->findID(id_, depth-1);
            if (r != 0) {
                return r;
            }
        }
    }
    return 0;
}

XML_Node* XML_Node::findByAttr(const std::string& attr,
                               const std::string& val, int depth) const
{
    if (hasAttrib(attr) && attrib(attr) == val) {
        return const_cast<XML_Node*>(this);
    }
    if (depth > 0) {
        size_t n = nChildren();
        for (size_t i = 0; i < n; i++) {
            XML_Node* r = m_children[i]->findByAttr(attr, val, depth - 1);
            if (r != 0) {
                return r;
            }
        }
    }
    return 0;
}

const XML_Node* XML_Node::findByName(const std::string& nm, int depth) const
{
    if (name() == nm) {
        return const_cast<XML_Node*>(this);
    }
    if (depth > 0) {
        for (size_t i = 0; i < nChildren(); i++) {
            XML_Node* r = m_children[i]->findByName(nm);
            if (r != 0) {
                return r;
            }
        }
    }
    return 0;
}

XML_Node* XML_Node::findByName(const std::string& nm, int depth)
{
    if (name() == nm) {
        return this;
    }
    if (depth > 0) {
        for (size_t i = 0; i < nChildren(); i++) {
            XML_Node* r = m_children[i]->findByName(nm);
            if (r != 0) {
                return r;
            }
        }
    }
    return 0;
}

std::vector<XML_Node*> XML_Node::getChildren(const std::string& nm) const
{
    std::vector<XML_Node*> children_;
    for (size_t i = 0; i < nChildren(); i++) {
        if (caseInsensitiveEquals(child(i).name(),  nm)) {
            children_.push_back(&child(i));
        }
    }
    return children_;
}

XML_Node& XML_Node::child(const std::string& aloc) const
{
    string loc = aloc;
    while (true) {
        size_t iloc = loc.find('/');
        if (iloc != string::npos) {
            string cname = loc.substr(0,iloc);
            loc = loc.substr(iloc+1, loc.size());
            auto i = m_childindex.find(cname);
            if (i != m_childindex.end()) {
                return i->second->child(loc);
            } else {
                throw XML_NoChild(this, m_name, cname, root().m_filename,
                    lineNumber());
            }
        } else {
            auto i = m_childindex.find(loc);
            if (i != m_childindex.end()) {
                return *(i->second);
            } else {
                throw XML_NoChild(this, m_name, loc, root().m_filename,
                    lineNumber());
            }
        }
    }
}

XML_Node& XML_Node::root() const
{
    return *m_root;
}

void XML_Node::setRoot(const XML_Node& newRoot)
{
    m_root = const_cast<XML_Node*>(&newRoot);
    for (size_t i = 0; i < m_children.size(); i++) {
        m_children[i]->setRoot(newRoot);
    }
}

void XML_Node::build(const std::string& filename)
{
    ifstream fin(filename);
    if (!fin) {
        throw CanteraError("XML_Node::build",
            "Unable to open file '{}' for reading.", filename);
    }
    build(fin, filename);
}

void XML_Node::build(std::istream& f, const std::string& filename)
{
    m_filename = filename;
    XML_Reader r(f);
    XML_Node* node = this;
    bool first = true;
    while (!f.eof()) {
        map<string, string> node_attribs;
        string nm = r.readTag(node_attribs);

        if (nm == "EOF") {
            break;
        }
        if (nm == "--" && m_name == "--" && m_root == this) {
            continue;
        }
        int lnum = r.m_line;
        if (nm[nm.size() - 1] == '/') {
            string nm2 = nm.substr(0,nm.size()-1);
            if (first) {
                node->setName(nm2);
                first = false;
            } else {
                node = &node->addChild(nm2);
            }
            node->addValue("");
            node->attribs() = node_attribs;
            node->setLineNumber(lnum);
            node = node->parent();
        } else if (nm[0] != '/') {
            if (nm[0] != '!' && nm[0] != '-' && nm[0] != '?') {
                if (first) {
                    node->setName(nm);
                    first = false;
                } else {
                    node = &node->addChild(nm);
                }
                node->addValue(r.readValue());
                node->attribs() = node_attribs;
                node->setLineNumber(lnum);
            } else if (nm.substr(0,2) == "--") {
                if (nm.substr(nm.size()-2,2) == "--") {
                    node->addComment(nm.substr(2,nm.size()-4));
                }
            }
        } else {
            if (node->name() != nm.substr(1,nm.size()-1)) {
                throw XML_TagMismatch(node->name(), nm.substr(1,nm.size()-1),
                    root().m_filename, lnum);
            }
            node = node->parent();
        }
    }
}

void XML_Node::copyUnion(XML_Node* const node_dest) const
{
    node_dest->addValue(m_value);
    if (m_name == "") {
        return;
    }
    for (const auto& attr : m_attribs) {
        if (!node_dest->hasAttrib(attr.first)) {
            node_dest->addAttribute(attr.first, attr.second);
        }
    }
    const vector<XML_Node*> &vsc = node_dest->children();
    for (size_t n = 0; n < m_children.size(); n++) {
        XML_Node* sc = m_children[n];
        size_t ndc = node_dest->nChildren();
        XML_Node* dc = 0;
        if (! sc->m_iscomment) {
            for (size_t idc = 0; idc < ndc; idc++) {
                XML_Node* dcc = vsc[idc];
                if (dcc->name() == sc->name()) {
                    if (sc->hasAttrib("id") && sc->attrib("id") != dcc->attrib("id")) {
                        break;
                    }
                    if (sc->hasAttrib("name") && sc->attrib("name") != dcc->attrib("name")) {
                        break;
                    }
                    if (sc->hasAttrib("model") && sc->attrib("model") != dcc->attrib("model")) {
                        break;
                    }
                    if (sc->hasAttrib("title") && sc->attrib("title") != dcc->attrib("title")) {
                        break;
                    }
                    dc = vsc[idc];
                }
            }
        }
        if (!dc) {
            node_dest->addChild(sc->name());
            dc = vsc[ndc];
        }
        sc->copyUnion(dc);
    }
}

void XML_Node::copy(XML_Node* const node_dest) const
{
    node_dest->addValue(m_value);
    node_dest->setName(m_name);
    node_dest->setLineNumber(m_linenum);
    if (m_name == "") {
        return;
    }
    for (const auto& attr : m_attribs) {
        node_dest->addAttribute(attr.first, attr.second);
    }
    const vector<XML_Node*> &vsc = node_dest->children();

    for (size_t n = 0; n < m_children.size(); n++) {
        XML_Node* sc = m_children[n];
        size_t ndc = node_dest->nChildren();
        // Here is where we create the child node.
        node_dest->addChild(sc->name());
        XML_Node* dc = vsc[ndc];
        sc->copy(dc);
    }
}

void XML_Node::lock()
{
    m_locked = true;
    for (size_t i = 0; i < m_children.size(); i++) {
        m_children[i]->lock();
    }
}

void XML_Node::unlock()
{
    m_locked = false;
    for (size_t i = 0; i < m_children.size(); i++) {
        m_children[i]->unlock();
    }
}

void XML_Node::writeHeader(std::ostream& s)
{
    s << "<?xml version=\"1.0\"?>" << endl;
}

void XML_Node::write_int(std::ostream& s, int level, int numRecursivesAllowed) const
{
    if (m_name == "") {
        return;
    }

    string indent(level, ' ');
    if (m_iscomment) {
        // In the comment section, we test to see if there already is a space
        // beginning and ending the comment. If there already is one, we don't
        // add another one.
        s << endl << indent << "<!--";
        if (! isspace(m_value[0])) {
            s << " ";
        }
        s << m_value;
        if (! isspace(m_value[m_value.size()-1])) {
            s << " ";
        }
        s << "-->";
        return;
    }

    s << indent << "<" << m_name;
    for (const auto& attr : m_attribs) {
        s << " " << attr.first << "=\"" << attr.second << "\"";
    }
    if (m_value == "" && m_children.empty()) {
        s << "/>";
    } else {
        s << ">";

        if (m_value != "") {
            string vv = m_value;
            string::size_type ieol = vv.find('\n');
            if (ieol != string::npos) {
                while (true) {
                    ieol = vv.find('\n');
                    if (ieol != string::npos) {
                        if (ieol == 0) {
                            s << endl << indent << "  ";
                        } else {
                            size_t jf = ieol - 1;
                            for (int j = 0; j < (int) ieol; j++) {
                                if (! isspace(vv[j])) {
                                    jf = j;
                                    break;
                                }
                            }
                            s << endl << indent << "  " << vv.substr(jf,ieol-jf);
                        }
                        vv = vv.substr(ieol+1);
                    } else {
                        size_t lll = vv.size() - 1;
                        if (lll != npos) {
                            size_t jf = lll;
                            for (size_t j = 0; j < lll; j++) {
                                if (! isspace(vv[j])) {
                                    jf = j;
                                    break;
                                }
                            }
                            if (jf < lll) {
                                s << endl << indent << "  " << vv.substr(jf);
                            }
                        }
                        break;
                    }
                }
                s << endl << indent;
            } else {
                bool doSpace = true;
                bool doNewLine = false;
                size_t ll = m_value.size() - 1;
                if (ll > 25) {
                    doNewLine = true;
                }
                if (m_name == "floatArray") {
                    doNewLine = true;
                }
                if (doNewLine) {
                    doSpace = false;
                }

                if (doNewLine) {
                    s << endl << indent << "  ";
                }

                // Put spaces around a raw value field for readability
                if (doSpace && (! isspace(m_value[0]))) {
                    s << " ";
                }

                // Write out the value
                s << m_value;

                if (doSpace && (! isspace(m_value[ll]))) {
                    s << " ";
                }
                if (doNewLine) {
                    s << endl << indent;
                }
            }
        }
        if (numRecursivesAllowed > 0) {
            for (size_t i = 0; i < m_children.size(); i++) {
                s << endl;
                m_children[i]->write_int(s,level + 2, numRecursivesAllowed - 1);
            }
        }
        if (!m_children.empty()) {
            s << endl << indent;
        }
        s << "</" << m_name << ">";
    }
}

void XML_Node::write(std::ostream& s, const int level, int numRecursivesAllowed) const
{
    write_int(s, level, numRecursivesAllowed);
    s << endl;
}

XML_Node* findXMLPhase(XML_Node* root,
                       const std::string& phaseId)
{
    XML_Node* scResult = 0;
    if (!root) {
        return 0;
    }
    if (root->name() == "phase") {
        if (phaseId == "") {
            return root;
        }
        if (phaseId == root->id()) {
            return root;
        }
    }

    const vector<XML_Node*> &vsc = root->children();
    for (size_t n = 0; n < root->nChildren(); n++) {
        XML_Node* sc = vsc[n];
        if (sc->name() == "phase") {
            if (phaseId == "") {
                return sc;
            }
            if (phaseId == sc->id()) {
                return sc;
            }
        }
    }
    for (size_t n = 0; n < root->nChildren(); n++) {
        XML_Node* sc = vsc[n];
        scResult = findXMLPhase(sc, phaseId);
        if (scResult) {
            return scResult;
        }
    }
    return scResult;
}

}
