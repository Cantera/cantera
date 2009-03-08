/**
 * @file xml.h
 *
 * Classes providing support for XML data files. These classes
 * implement only those aspects of XML required to read, write, and
 * manipulate CTML data files.
 */

/* $Author: hkmoffa $
 * $Revision: 1.14 $
 * $Date: 2006/07/13 20:05:11 $
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_XML_H
#define CT_XML_H

#include <string>
#include <vector>
#include <iostream>
using namespace std;

#include "ctexceptions.h"
#include "ct_defs.h"
#include "stringUtils.h"
#include <stdio.h>
#include "global.h"

#define XML_INDENT 4


namespace Cantera {


    /**
     * Class XML_Reader is designed for internal use.
     */
    class XML_Reader {
    public:
        XML_Reader(istream& input) : m_s(input), m_line(0) {}

        istream& m_s;
        int m_line;

        void getchr(char& ch);
        string strip(const string& aline);
        string inquotes(const string& aline);

	int findQuotedString(const string& aline, string &rstring);

        void parseTag(string line, string& name, map<string, string>& attribs);
        string readTag(map<string, string>& attribs);
        string readValue();
    };


    //////////////////////////  XML_Node  /////////////////////////////////

    class XML_Node {
    public:

        XML_Node(string nm = "--", XML_Node* p = 0, int n = 0);
        virtual ~XML_Node();
        void addComment(string comment);
        XML_Node& addChild(XML_Node& node);
        XML_Node& addChild(string name);
        XML_Node& addChild(string name, string value);
        XML_Node& addChild(string name, double value, string fmt="%g");
        void removeChild(XML_Node* node);
        void addValue(string val);
        void addValue(doublereal val, string fmt="%g");
        void addAttribute(string attrib, string value);
        void addAttribute(string attrib, double value, string fmt="%g");
        void writeHeader(ostream& s);
        string value() const { return m_value; }
        string value(string loc) const { return child(loc).value(); }
        void setLineNumber(int n) {m_linenum = n;}
        int lineNumber() const {return m_linenum;}
        doublereal fp_value() const { 
            return atof(m_value.c_str()); 
        }
        integer int_value() const { 
            return atoi(m_value.c_str()); 
        }
        string operator()() const { return m_value; }
        string operator()(string loc) const { return value(loc); }

	/**
	 * The operator[] is overloaded to provide a lookup capability
	 * on attributes for the current XML element.
	 *
	 * For example
	 *     xmlNode["id"] 
	 * will return the value of the attribute "id" for the current
	 * XML element. It will return the blank string if there isn't
	 * an attribute with that name.
	 */
        string operator[](string attr) const {
            return attrib(attr);
        }

        /**
         * This function searches the attibutes vector for the parameter 
         * string attribute. If a match is found, the attribute value
         * is returned as a string. If no match is found, the empty string
         * is returned.
         *
         * @param(attr) String containing the attribute to be searched for.
         */
        string attrib(string attr) const { 
            map<string,string>::const_iterator i = m_attribs.find(attr);
            if (i != m_attribs.end()) return i->second;
            return ""; 
        }

        map<string,string>& attribs() { return m_attribs; }
        XML_Node* parent() const { return m_parent; }
        XML_Node* setParent(XML_Node* p) { m_parent = p; return p; }

        bool hasChild(string ch) const {
            return (m_childindex.find(ch) != m_childindex.end());
        }
        bool hasAttrib(string a) const {
            return (m_attribs.find(a) != m_attribs.end());
        }


        string name() const { return m_name; }
        string id() const {
            if (hasAttrib("id")) return attrib("id");
            else return "";
        }
        int number() const { return m_n; }

        XML_Node& child(int n) const { return *m_children[n]; }
        vector<XML_Node*>& children()  { return m_children; }
        const vector<XML_Node*>& children() const { return m_children; }
        int nChildren() const { return m_nchildren; }

        void build(istream& f);

        void _require(string a, string v) const;

	/**
	 * This routine carries out a search for an XML node based
	 * on both the xml element name and the attribute ID.
	 * If exact matches are found for both fields, the pointer
	 * to the matching XML Node is returned.
	 * 
	 * The ID attribute may be defaulted by setting it to "".
	 * In this case the pointer to the first xml element matching the name
	 * only is returned.
	 * @internal
	 * This algorithm does a lateral search of first generation children
	 * first before diving deeper into each tree branch.
	 */
	XML_Node* findNameID(const string &nameTarget, 
			     const string &idTarget) const;

        XML_Node* findID(const string& id, int depth=100) const;
        XML_Node* findByAttr(const string& attr, const string& val);
        XML_Node* findByName(const string& nm);
        void getChildren(string name, vector<XML_Node*>& children) const;
        XML_Node& child(string loc) const;
        void write(ostream& s, int level = 0) const;
        XML_Node& root() const { return *m_root; }
        void setRoot(XML_Node& root) { m_root = &root; }
        void copyUnion(XML_Node *node_dest);
        void copy(XML_Node *node_dest);
        void lock() {m_locked = true;}
        void unlock() {m_locked = false; }

    private:
	void write_int(ostream& s, int level = 0) const;

    protected:

        string m_name;
        string m_value;
        map<string, XML_Node*> m_childindex;
        map<string, string> m_attribs;
        XML_Node* m_parent;
        XML_Node* m_root;
        bool m_locked;
        vector<XML_Node*> m_children;
        int m_nchildren;
        int m_n;
        bool m_iscomment;
        int m_linenum;
    };

  
    XML_Node * findXMLPhase(XML_Node* root, const string &id);

}

#endif

