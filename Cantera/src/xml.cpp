
// simple xml functions

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <algorithm>
using namespace std;

#include "xml.h"

#define XML_INDENT 4

namespace Cantera {


    static void split(const string& src, string& file, string& id) { 
        int ipound = src.find('#');

        if (ipound >= 0) {
            id = src.substr(ipound+1,src.size());
            file = src.substr(0,ipound);
        }            
        else {
            id = "";
            file = src;
        }
    }


    ////////////////////// exceptions ////////////////////////////

    
    class XML_Error : public CanteraError {
    public:
        XML_Error(int line=0) : m_line(line) {
            m_msg = "Error in XML file";
            if (line > 0) {
                m_msg += " at line " + line;
            }
            m_msg += ".\n";
        }
        virtual ~XML_Error() {}
    protected:
        int m_line;
        string m_msg;
    };
    
    class XML_TagMismatch : public XML_Error {
    public:
        XML_TagMismatch(string opentag, string closetag, 
            int line=0) : XML_Error(line) {
            m_msg += "<" + opentag + "> paired with </" + closetag + ">.\n";
            setError("XML_TagMismatch",m_msg);
        }
        virtual ~XML_TagMismatch() {}
    };

    class XML_NoChild : public XML_Error {
    public:
        XML_NoChild(string parent, string child) {
            m_msg += "Node " + parent + " has no child named " 
                     + child + ".\n";
            setError("XML_NoChild",m_msg);
        }
        virtual ~XML_NoChild() {}
    };

    class XML_IllegalUnits : public XML_Error {
    public:
        XML_IllegalUnits(string name, string units) {
            m_msg += "Illegal units (" + units + 
                     ") specified for node " + name + ".\n";
            setError("XML_IllegalUnits",m_msg);
        }
        virtual ~XML_IllegalUnits() {}
    };



    //////////////////// XML_Reader methods ///////////////////////


    
    void XML_Reader::getchr(char& ch) {
        m_s.get(ch);
        if (ch == '\n') m_line++;
    }

    string XML_Reader::strip(const string& aline) {
        int len = aline.size();
        int i, j;
        for (i = len-1; i >= 0; i--) 
            if (aline[i] != ' ' && aline[i] != '\n') break;
        for (j = 0;  j < i; j++)
            if (aline[j] != ' ' && aline[j] != '\n') break;
        return aline.substr(j, i - j + 1);
    }

    string XML_Reader::inquotes(const string& aline) {
        int len = aline.size();
        int i, j;
        for (i = len-1; i >= 0; i--) 
            if (aline[i] == '"') break;
        for (j = 0;  j < i; j++) 
            if (aline[j] == '"') break;
        if (j == i) return "";
        else return aline.substr(j+1, i - j - 1);
    }
    

    void XML_Reader::parseTag(string line, string& name, 
			      map<string, string>& attribs) {
        int iloc;
        string attr, val;
        string s = strip(line);
        iloc = s.find(' ');
        if (iloc > 0) {
	    name = s.substr(0, iloc);
            s = strip(s.substr(iloc+1,s.size()));
            if (s[s.size()-1] == '/') name += "/";
            // get attributes
            while (1) {
                iloc = s.find('=');
                if (iloc < 0) break;
                attr = strip(s.substr(0,iloc));
                if (attr == "") break;
                s = strip(s.substr(iloc+1,s.size()));
                iloc = s.find(' ');
                if (iloc < 0) iloc = s.size();
                val = inquotes(s.substr(0,iloc));
                attribs[attr] = val;
                if (iloc < int(s.size())) 
                    s = strip(s.substr(iloc+1,s.size()));
                else
                    break;
            }
        }
        else {
            name = s;
        }
    }
    
    string XML_Reader::readTag(map<string, string>& attribs) {
        string name, tag = "";
        bool incomment = false;
        char ch  = '-';
        while (1) {
            if (m_s.eof() || (getchr(ch), ch == '<')) break;
        }
        char ch1 = ' ', ch2 = ' ';
        while (1) {
            if (m_s.eof()) { tag = "EOF"; break;}
            ch2 = ch1;
            ch1 = ch;
            getchr(ch);
            if (ch == '-') {
                if (ch1 == '-' && ch2 == '!') {
                    incomment = true;
                    tag = "-";
                }
            }
            else if (ch == '>') {
                if (incomment) {
                    if (ch1 == '-' && ch2 == '-') break;
                }
                else
                    break;
            }
            if (isprint(ch)) tag += ch;
        }
        if (incomment) {
            attribs.clear();
            return tag;
        }
        else {
            parseTag(tag, name, attribs);
            return name;
        }
    }

    string XML_Reader::readValue() {
        string tag = "";
        char ch, lastch;
        ch = '\n';
        bool front = true;
        while (1) {
            if (m_s.eof()) break;
            lastch = ch;
            getchr(ch);
            if (ch == '\n') 
                front = true;
            else if (ch != ' ') 
                front = false;
            if (ch == '<') {
                m_s.putback(ch); 
                break;
            }
            if (front && lastch == ' ' && ch == ' ') ;
            else tag += ch;
        }
        return strip(tag);
    }


    //////////////////////////  XML_Node  /////////////////////////////////


    XML_Node::XML_Node(string nm, XML_Node* p, int n) 
        : m_name(nm), m_level(0),  m_parent(p), m_nchildren(0), 
          m_n(n), m_iscomment(false) {
        if (!p) m_root = this;
        else m_root = &p->root();
    }
    

    XML_Node::~XML_Node() {
        int n = m_children.size();
        for (int i = 0; i < n; i++) {
            if (m_children[i]) {
                if (m_children[i]->parent() == this) {
                    delete m_children[i];
                    m_children[i] = 0;
                }
            }
        }
    }

    void XML_Node::addComment(string comment) {
        addChild("comment",comment);
    }

    XML_Node& XML_Node::addChild(XML_Node& node) {
        m_children.push_back(&node);
        m_nchildren = m_children.size();
        m_childindex[node.name()] = m_children.back();
        node.setRoot(root());
        return *m_children.back();
    }

    XML_Node& XML_Node::addChild(string name) { 
        int n = m_children.size();
        m_children.push_back(new XML_Node(name, this, n));
        m_nchildren = m_children.size();
        m_childindex[name] = m_children.back();
        m_children.back()->setParent(this);
        return *m_children.back();
    }

    void XML_Node::removeChild(XML_Node* node) {
        vector<XML_Node*>::iterator i;
        i = find(m_children.begin(), m_children.end(), node);
        m_children.erase(i);
        m_nchildren = m_children.size();
        m_childindex.erase(node->name());
    }

    XML_Node* XML_Node::findID(const string& id, int depth) {
        if (hasAttrib("id")) {
            if (attrib("id") == id) {
                return this;
            }
        }
        if (depth > 0) {
            XML_Node* r = 0;
            int n = nChildren();
            for (int i = 0; i < n; i++) {
                r = m_children[i]->findID(id, depth-1);
                if (r != 0) return r;
            }
        }
        return 0;
    }

    XML_Node* XML_Node::findByAttr(const string& attr, const string& val) {
        if (hasAttrib(attr)) {
            if (attrib(attr) == val) {
                return this;
            }
        }
        XML_Node* r = 0;
        int n = nChildren();
        for (int i = 0; i < n; i++) {
            r = m_children[i]->findByAttr(attr, val);
            if (r != 0) return r;
        }
        return 0;
    }

    XML_Node* XML_Node::findByName(const string& nm) {
        if (name() == nm) {
            return this;
        }
        XML_Node* r = 0;
        int n = nChildren();
        for (int i = 0; i < n; i++) {
            r = m_children[i]->findByName(nm);
            if (r != 0) return r;
        }
        return 0;
    }

    /**
     *  addChild(string name, string value):
     *    
     *    Add a child node to the current xml node, and at the
     *    same time add a value to the child
     *
     *    Resulting XML string:
     *      <name>value</name>
     *
     *    Return
     *    -------
     *     Returns a reference to the created child XML_Node object
     */
    XML_Node& XML_Node::addChild(string name, string value) {
        XML_Node& c = addChild(name);
        c.addValue(value);
        return c;
    }
    
    XML_Node& XML_Node::addChild(string name, double value, string fmt) {
        XML_Node& c = addChild(name);
        c.addValue(value,fmt);
        return c;
    }
    
    void XML_Node::addValue(string val) { 
        m_value = val;
        if (m_name == "comment") m_iscomment = true;
    }
    void XML_Node::addValue(doublereal val, string fmt) {
        char buf[30];
        sprintf(buf,fmt.c_str(),val);
        m_value = stripws(buf);
    }
    
    void XML_Node::addAttribute(string attrib, string value) {
        m_attribs[attrib] = value;
    }
    void XML_Node::addAttribute(string attrib, double value, string fmt) {
        m_attribs[attrib] = fp2str(value, fmt);
    }
    
    void XML_Node::writeHeader(ostream& s) {
        s << "<?xml version=\"1.0\"?>" << endl;
    }

    void XML_Node::build(istream& f) {
        XML_Reader r(f);
        string nm, nm2, val;
        XML_Node* node = this;
        map<string, string> attribs;
        while (!f.eof()) {
            attribs.clear();
            nm = r.readTag(attribs);
            if (nm == "EOF") break;
            int lnum = r.m_line;
            if (nm[nm.size() - 1] == '/') {
                nm2 = nm.substr(0,nm.size()-1);
                node = &node->addChild(nm2);
                node->addValue("");
                node->attribs() = attribs;
                node = node->parent();
            }
            else if (nm[0] != '/') {
                if (nm[0] != '!' && nm[0] != '-' && nm[0] != '?') {
                    node = &node->addChild(nm);
                    val = r.readValue();
                    node->addValue(val);
                    node->attribs() = attribs;
                }
                else if (nm.substr(0,2) == "--") {
                    if (nm.substr(nm.size()-2,2) == "--") {
                        node->addComment(nm.substr(2,nm.size()-4));
                    }
                }
            }
            else {
                if (node->name() != nm.substr(1,nm.size()-1)) 
                    throw XML_TagMismatch(node->name(), 
                        nm.substr(1,nm.size()-1), lnum);
                node = node->parent();
            }
        }
    } 


    void XML_Node::getChildren(string nm, 
        vector<XML_Node*>& children) const {
        int i, n = nChildren();
        for (i = 0; i < n; i++) {
            if (child(i).name() == nm) {
                children.push_back(&child(i));
            } 
        }
    }
    
    XML_Node& XML_Node::child(string loc) const {
        int iloc;
        string cname;
        map<string,XML_Node*>::const_iterator i;
        while (1) {
            iloc = loc.find('/');
            if (iloc >= 0) {
                cname = loc.substr(0,iloc);
                loc = loc.substr(iloc+1, loc.size());
                i = m_childindex.find(cname);
                //XML_Node* chld = m_childindex[cname];
                if (i != m_childindex.end()) return i->second->child(loc);
                else throw XML_NoChild(m_name, cname);
            }
            else {
                i = m_childindex.find(loc);
                if (i != m_childindex.end()) return *(i->second);
                //XML_Node* chld = m_childindex[loc];
                //if (chld) return *chld;
                else throw XML_NoChild(m_name, loc);
            }
        }
    }
    
    void XML_Node::write(ostream& s, int level) {

        if (m_name == "") return;

        m_level = level;
        string indent(level,' ');
        if (m_iscomment) {
            s << endl << indent << "<!-- " << m_value << " -->";
            return;
        }
        
        s << indent << "<" << m_name;
        map<string,string>::const_iterator b = m_attribs.begin();
        for (; b != m_attribs.end(); ++b) {
            s << " " << b->first << "=\"" << b->second << "\"";
        }
        if (m_value == "" && m_nchildren == 0) {
            s << "/>";
        }
        else {
            s << ">";
            
            if (m_value != "") {
                string vv = m_value;
                int ieol = vv.find('\n');
                if (ieol >= 0) { 
                    while (1 > 0) {
                        ieol = vv.find('\n');
                        if (ieol >= 0) {
                            s << endl << indent << vv.substr(0,ieol);
                            vv = vv.substr(ieol+1,vv.size());
                        }
                        else {
                            s << endl << indent << vv;
                            break;
                        }
                    }
                }
                else
                    s << m_value;
            }
            int i;
            for (i = 0; i < m_nchildren; i++) {
                s << endl;
                m_children[i]->write(s,level + 2);
            }
            if (m_nchildren > 0) s << endl << indent;
            s << "</" << m_name << ">";
        }
    }

    XML_Node* XML_Node::getRef() {
        if (!hasAttrib("idRef")) return this;
        XML_Node& node = *this;
        return find_XML(node["src"], &root(), node["idRef"]);
    }

    XML_Node* find_XML(string src, XML_Node* root, string id, string loc, 
        string name) {

        string file, id2;
        split(src, file, id2);
        src = file;
        if (id2 != "") id = id2;

        XML_Node *doc = 0, *r = 0;
        if (src != "") {
            doc = new XML_Node("doc");
            string spath = findInputFile(src);
            ifstream fin(spath.c_str());
            if (!fin)
                throw CanteraError("find_XML","could not open file "+src+
                    " for input.");
            doc->build(fin);
            root = 0;
        }
        else if (root) { 
            doc = root;
        }
        else {
            throw CanteraError("find_XML",
                "either root or src must be specified.");
        }

        try {
            if (id != "") 
                r = doc->findID(id);
            else if (loc != "")
                r = &doc->child(loc);
            else if (name != "") 
                r = doc->findByName(name);
            if (!r) {
                string opt = " src="+src+", loc="+loc+", id="
                             +id+", name="+name;
                throw CanteraError("find_XML", "XML element with "+opt+
                    " not found.");
            }
            return r;
        }
        catch (CanteraError) {
            
            // root was used, but element was not found. Try src.
            if (root && src != "") {
                return find_XML(src, 0, id, loc, name);
            }
            else {
                string opt = " src="+src+", loc="+loc+", id="
                             +id+", name="+name;
                throw CanteraError("find_XML", "XML element with "+opt+
                    " not found.");
                return 0;
            }
        }
    }

}


