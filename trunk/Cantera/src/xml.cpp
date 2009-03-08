
// simple xml functions

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "config.h"
#ifdef HAS_SSTREAM
#include <sstream>
#endif

#include <algorithm>
using namespace std;

#include "xml.h"

#define XML_INDENT 4

namespace Cantera {


    ////////////////////// exceptions ////////////////////////////

    
    class XML_Error : public CanteraError {
    public:
        XML_Error(int line=0) : m_line(line) {
            m_msg = "Error in XML file";
            if (line > 0) {
                m_msg += " at line " + int2str(line+1);
            }
            m_msg += ".\n";
            //setError("XML_Error",m_msg);
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
        XML_NoChild(const XML_Node* p, string parent, 
            string child, int line=0) : XML_Error(line) {
            m_msg += "           The XML Node \"" + parent + 
		"\", does not contain a required\n" +
                     "           XML child node named \"" 
                     + child + "\".\n";
#ifdef HAS_SSTREAM
            ostringstream ss(ostringstream::out);
            p->write(ss,1);
            m_msg += ss.str() + "\n";
#endif
            setError("XML_NoChild",m_msg);
        }
        virtual ~XML_NoChild() {}
    };

    class XML_IllegalUnits : public XML_Error {
    public:
        XML_IllegalUnits(string name, string units, int line=0) : XML_Error(line) {
            m_msg += "Illegal units (" + units + 
                     ") specified for node " + name + ".\n";
            setError("XML_IllegalUnits",m_msg);
        }
        virtual ~XML_IllegalUnits() {}
    };



    //////////////////// XML_Reader methods ///////////////////////


    /// Get a single character from the input stream. If the character
    /// is a new-line character, then increment the line count.
    void XML_Reader::getchr(char& ch) {
        m_s.get(ch);
        if (ch == '\n') m_line++;
    }


    /// Returns string 'aline' stripped of leading and trailing white
    /// space.
    /// @todo why is this a class method?
    string XML_Reader::strip(const string& aline) {
        int len = static_cast<int>(aline.size());
        int i, j;
        for (i = len-1; i >= 0; i--) 
            if (aline[i] != ' ' && aline[i] != '\n') break;
        for (j = 0;  j < i; j++)
            if (aline[j] != ' ' && aline[j] != '\n') break;
        return aline.substr(j, i - j + 1);
    }


    /// Looks for a substring within 'aline' enclosed in double
    /// quotes, and returns this substring (without the quotes) if
    /// found.  If not, an empty string is returned.
    /// @todo why is this a class method?
    string XML_Reader::inquotes(const string& aline) {
        int len = static_cast<int>(aline.size());
        int i, j;
        for (i = len-1; i >= 0; i--) 
            if (aline[i] == '"') break;
        for (j = 0;  j < i; j++) 
            if (aline[j] == '"') break;
        if (j == i) return "";
        else return aline.substr(j+1, i - j - 1);
    }

    /**
     * Find the first position of a character, q, in string s,
     * which is not immediately preceded by the backslash character
     * '\'
     */
    static string::size_type findUnbackslashed(string s, const char q,
					       string::size_type istart = 0) {
	string::size_type iloc, icurrent, len;
	icurrent = istart;
	len = s.size();
	while (1) {
	  iloc = s.find(q, icurrent);
	  if (iloc == string::npos || iloc == 0) {
	    return iloc;
	  }
	  char cm1 = s[iloc-1];
	  if (cm1 == '\\') {
	    if (iloc >= (len -1)) return string::npos;
	    icurrent = iloc + 1;
	  } else {
	    return iloc;
	  }
	}
    }

    /**
     *  Searches a string for the first occurrence of a valid
     *  quoted string. Quotes can start with either a single
     *  quote or a double quote, but must also end with the same
     *  type. Quotes may be commented out by preceding with a
     *  backslash character, '\\'. 
     */
    int XML_Reader::findQuotedString(const string& s, string &rstring) {
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
	if (iloc1 != string::npos) {
	  if (iloc1 < ilocStart) {
	    ilocStart = iloc1;
	    qtype = q1;
	  }
	}
	if (qtype == ' ') return 0;

	iloc1 = findUnbackslashed(s, qtype, ilocStart+1);

	if (iloc1 == string::npos) {
	  return 0;
	}
	/*
	 * Define the return string by the two endpoints.
	 * Strip the surrounding quotes as well
	 */
	rstring = s.substr(ilocStart + 1, iloc1 - 1);
	/*
	 * Return the first character position past the quotes
	 */
	return static_cast<int>(iloc1)+1;	
    }
    
    /**
     * parseTag parses XML tags, i.e., the XML elements that are
     * inbetween angle brackets.
     */
    void XML_Reader::parseTag(string tag, string& name, 
			      map<string, string>& attribs) {
        string::size_type iloc;
        string attr, val;
        string s = strip(tag);
        iloc = s.find(' ');
        if (iloc != string::npos) {
	    name = s.substr(0, iloc);
            s = strip(s.substr(iloc+1,s.size()));
            if (s[s.size()-1] == '/') {
	      name += "/";
	    }

            // get attributes
            while (1) {
                iloc = s.find('=');
                if (iloc == string::npos) break; 
                attr = strip(s.substr(0,iloc));
                if (attr == "") break;
                s = strip(s.substr(iloc+1,s.size()));
                       iloc = findQuotedString(s, val);
                attribs[attr] = val;
		if (iloc != string::npos) {
		  if (iloc < s.size()) 
		      s = strip(s.substr(iloc,s.size()));
		  else
		      break;
		}
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
        : m_name(nm), m_value(""), m_parent(p),
          m_locked(false), m_nchildren(0), 
          m_n(n), m_iscomment(false) {
        if (!p) m_root = this;
        else m_root = &p->root();
    }
    

    XML_Node::~XML_Node() {
        if (m_locked) 
            throw CanteraError("XML_Node::~XML_Node",
                "attempt to delete locked XML_Node "+name());
        int n = static_cast<int>(m_children.size());
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
        m_nchildren = static_cast<int>(m_children.size());
        m_childindex[node.name()] = m_children.back();
        node.setRoot(root());
        return *m_children.back();
    }

    XML_Node& XML_Node::addChild(string name) { 
        int n = static_cast<int>(m_children.size());
	XML_Node *xxx = new XML_Node(name, this, n);
        m_children.push_back(xxx);
        m_nchildren = static_cast<int>(m_children.size());
        m_childindex[name] = m_children.back();
        m_children.back()->setParent(this);
        return *m_children.back();
    }

    void XML_Node::removeChild(XML_Node* node) {
        vector<XML_Node*>::iterator i;
        i = find(m_children.begin(), m_children.end(), node);
        m_children.erase(i);
        m_nchildren = static_cast<int>(m_children.size());
        m_childindex.erase(node->name());
    }

    /**
     * This routine carries out a search for an XML node based
     * on both the xml element name and the attribute ID.
     * If exact matches are found for both fields, the pointer
     * to the matching XML Node is returned.
     * 
     * The ID attribute may be defaulted by setting it to "".
     * In this case the pointer to the first xml element matching the name
     * is returned.
     *
     * This algorithm does a lateral search of first generation children
     * first before diving deeper into each tree branch.
     */
    XML_Node* XML_Node::
    findNameID(const string &nameTarget, const string &idTarget) const {
	XML_Node *scResult = 0;
	XML_Node *sc;
	string idattrib = id();
	int n;
	if (name() == nameTarget) {
	  if (idTarget == "" || idTarget == idattrib) {
	    return const_cast<XML_Node*>(this);
	  }
	}
	for (n = 0; n < m_nchildren; n++) {
	  sc = m_children[n];
	  if (sc->name() == nameTarget) {
	    if (idTarget == "") return sc;
	    idattrib = sc->id();
	    if (idTarget == idattrib) return sc;
	  }
	}
	for (n = 0; n < m_nchildren; n++) {
	  sc = m_children[n];
	  scResult = sc->findNameID(nameTarget, idTarget);
	  if (scResult) return scResult;
	}
	return scResult;
    }
    

    XML_Node* XML_Node::findID(const string& id, int depth) const {
        if (hasAttrib("id")) {
            if (attrib("id") == id) {
                return const_cast<XML_Node*>(this);
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

    XML_Node* XML_Node::findByAttr(const string& attr, 
        const string& val) {
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
	    if (nm == "--" && m_name == "--" && m_root == this) {
	      continue;
	    }
            int lnum = r.m_line;
            if (nm[nm.size() - 1] == '/') {
                nm2 = nm.substr(0,nm.size()-1);
                node = &node->addChild(nm2);
                node->addValue("");
                node->attribs() = attribs;
                node->setLineNumber(lnum);
                node = node->parent();
            }
            else if (nm[0] != '/') {
                if (nm[0] != '!' && nm[0] != '-' && nm[0] != '?') {
                    node = &node->addChild(nm);
                    val = r.readValue();
                    node->addValue(val);
                    node->attribs() = attribs;
                    node->setLineNumber(lnum);
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

    void XML_Node::copyUnion(XML_Node *node_dest) {
	XML_Node *sc, *dc;
	int ndc, idc;
	node_dest->addValue(m_value);
	if (m_name == "") return;
	map<string,string>::const_iterator b = m_attribs.begin();
        for (; b != m_attribs.end(); ++b) {
	  if (! node_dest->hasAttrib(b->first)) {
	    node_dest->addAttribute(b->first, b->second);
	  }
	}
	const vector<XML_Node*> &vsc = node_dest->children();
	for (int n = 0; n < m_nchildren; n++) {
	  sc = m_children[n];
	  ndc = node_dest->nChildren();
	  dc = 0;
	  if (! sc->m_iscomment) {
	    for (idc = 0; idc < ndc; idc++) {
	      XML_Node *dcc = vsc[idc];
	      if (dcc->name() == sc->name()) {
		if (sc->hasAttrib("id")) {
		  if (sc->attrib("id") != dcc->attrib("id")) break;
		}
		if (sc->hasAttrib("name")) {
		  if (sc->attrib("name") != dcc->attrib("name")) break;
		}
		if (sc->hasAttrib("model")) {
		  if (sc->attrib("model") != dcc->attrib("model")) break;
		}
		if (sc->hasAttrib("title")) {
		  if (sc->attrib("title") != dcc->attrib("title")) break;
		}
		dc = vsc[idc];
	      }
	    }
	  }
	  if (!dc) {
	    (void) node_dest->addChild(sc->name());
	    dc = vsc[ndc];
	  }
	  sc->copyUnion(dc);
	}
    }

   void XML_Node::copy(XML_Node *node_dest) {
	XML_Node *sc, *dc;
	int ndc;
	node_dest->addValue(m_value);
	if (m_name == "") return;
	map<string,string>::const_iterator b = m_attribs.begin();
        for (; b != m_attribs.end(); ++b) {
	  node_dest->addAttribute(b->first, b->second);
	}
	const vector<XML_Node*> &vsc = node_dest->children();

	for (int n = 0; n < m_nchildren; n++) {
	  sc = m_children[n];
	  ndc = node_dest->nChildren();
	  (void) node_dest->addChild(sc->name());
	  dc = vsc[ndc];
	  sc->copy(dc);
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
        string::size_type iloc;
        string cname;
        map<string,XML_Node*>::const_iterator i;

        while (1) {
            iloc = loc.find('/');
            if (iloc != string::npos) {
                cname = loc.substr(0,iloc);
                loc = loc.substr(iloc+1, loc.size());
                i = m_childindex.find(cname);
                //XML_Node* chld = m_childindex[cname];
                if (i != m_childindex.end()) return i->second->child(loc);
                else {
                    throw XML_NoChild(this, m_name, cname, lineNumber());
                }
            }
            else {
                i = m_childindex.find(loc);
                if (i != m_childindex.end()) return *(i->second);
                //XML_Node* chld = m_childindex[loc];
                //if (chld) return *chld;
                else {
                    throw XML_NoChild(this, m_name, loc, lineNumber());
                }
            }
        }
    }
    
    /**
     * Write an XML subtree to an output stream. This is the
     * main recursive routine. It doesn't put a final endl
     * on. This is fixed up in the public method.
     */
    void XML_Node::write_int(ostream& s, int level) const {

        if (m_name == "") return;

        string indent(level, ' ');
        if (m_iscomment) {
	  /*
	   * In the comment section, we test to see if there
	   * already is a space beginning and ending the comment.
	   * If there already is one, we don't add another one.
	   */
	  s << endl << indent << "<!--";
	  if (! isspace(m_value[0])) {
	    s << " ";
	  }
	  s << m_value;
	  int ll = static_cast<int>(m_value.size()) - 1;
	  if (! isspace(m_value[ll])) {
	      s << " ";
	  }
	  s << "-->";
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
                string::size_type ieol = vv.find('\n');
                if (ieol != string::npos) { 
                    while (1 > 0) {
                        ieol = vv.find('\n');
                        if (ieol != string::npos) {
			  int jf = ieol - 1;
			  for (int j = 0; j < (int) ieol; j++) {
			    if (! isspace(vv[j])) {
			      jf = j;
			      break;
			    }
			  }
			  s << endl << indent << "  " << vv.substr(jf,ieol-jf);
			  vv = vv.substr(ieol+1);
                        }
                        else {
			  int lll = static_cast<int>(vv.size()) - 1;
			  if (lll >= 0) {
			    int jf = lll;
			    for (int j = 0; j < lll; j++) {
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
                }
                else {
		  bool doSpace = true;
		  bool doNewLine = false;
		  int ll = static_cast<int>(m_value.size()) - 1;
		  if (ll > 15) {
		    doNewLine = true;
		  }
		  if (m_name == "floatArray") {
		    doNewLine = true;
		  }
		  if (doNewLine) doSpace = false;
		
		  if (doNewLine) {
		    s << endl << indent << "  ";
		  }
		  /*
		   * Put spaces around a raw value field for readability
		   */
		  if (doSpace && (! isspace(m_value[0]))) {
		    s << " ";
		  }
		  /*
		   * Write out the value
		   */
		  s << m_value;
		
		  if (doSpace && (! isspace(m_value[ll]))) {
		    s << " ";
		  }
		  if (doNewLine) {
		    s << endl << indent;
		  }
		}
            }
            int i;
            for (i = 0; i < m_nchildren; i++) {
                s << endl;
                m_children[i]->write_int(s,level + 2);
            }
            if (m_nchildren > 0) s << endl << indent;
            s << "</" << m_name << ">";
        }
    }

    /**
     * Write an XML subtree to an output stream. This is a 
     * wrapper around the static routine write_int(). All this
     * does is add an endl on to the output stream. write_int() is
     * fine, but the last endl wasn't being written.
     * It also checks for the special name "--". If found and we
     * are at the root of the xml tree, then the block
     * is skipped and the children are processed. "--" is used
     * to denote the top of the tree.
     */
    void XML_Node::write(ostream& s, int level) const {
	if (m_name == "--" && m_root == this) {
          for (int i = 0; i < m_nchildren; i++) {
            m_children[i]->write_int(s,level);
	    s << endl;
          }
        } else {
	  write_int(s, level);
	  s << endl;
	}
    }

    /**
     * _require() 
     *    Require that the current xml node have an attribute named
     *    by the first argument, a, and that this attribute have the
     *    the string value listed in the second argument, v.
     *    If not, throw a CanteraError exception.
     */
    void XML_Node::_require(string a, string v) const {
        if (hasAttrib(a)) {
            if (attrib(a) == v) return;
        }
        string msg="XML_Node "+name()+" is required to have the value "
                   "\""+v+"\", but instead is \""+attrib(a);
        throw CanteraError("XML_Node::require",msg);
    }


        
    XML_Node * findXMLPhase(XML_Node *root, 
			    const string &idtarget) {
	XML_Node *scResult = 0;
	XML_Node *sc;
	if (!root) return 0;
	string idattrib;
	string rname = root->name();
	if (rname == "phase") {
	  if (idtarget == "") return root;
	  idattrib = root->id();
	  if (idtarget == idattrib) return root;
	  else return               0;
	}

	const vector<XML_Node*> &vsc = root->children();
        int n;
	for (n = 0; n < root->nChildren(); n++) {
	  sc = vsc[n];
	  if (sc->name() == "phase") {
	    if (idtarget == "") return sc;
	    idattrib = sc->id();
	    if (idtarget == idattrib) return sc;
	  }
	}
	for (n = 0; n < root->nChildren(); n++) {
	  sc = vsc[n];
	  if  (sc->name() != "phase") {
	    scResult = findXMLPhase(sc, idtarget);
	    if (scResult) return scResult;
	  }
	}
	return scResult;
    }

}


