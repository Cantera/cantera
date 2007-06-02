/**
 * @file xml.h
 *
 * Classes providing support for XML data files. These classes
 * implement only those aspects of XML required to read, write, and
 * manipulate CTML data files.
 */

/* $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_XML_H
#define CT_XML_H

#include <string>
#include <vector>
#include <iostream>


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
        XML_Reader(std::istream& input) : m_s(input), m_line(0) {}

        std::istream& m_s;
        int m_line;

        void getchr(char& ch);
        std::string strip(const std::string& aline);
        std::string inquotes(const std::string& aline);

	int findQuotedString(const std::string& aline, std::string &rstring);

        void parseTag(std::string line, std::string& name, 
            std::map<std::string, std::string>& attribs);
        std::string readTag(std::map<std::string, std::string>& attribs);
        std::string readValue();
    };


    //////////////////////////  XML_Node  /////////////////////////////////

    class XML_Node {
    public:

        XML_Node(std::string nm = "--", XML_Node* p = 0, int n = 0);
   
      XML_Node(const XML_Node &right);
      XML_Node& operator=(const XML_Node &right);

   
        virtual ~XML_Node();
        void addComment(std::string comment);
        XML_Node& addChild(XML_Node& node);
        XML_Node& addChild(std::string name);
        XML_Node& addChild(std::string name, std::string value);
        XML_Node& addChild(std::string name, double value, std::string fmt="%g");
        void removeChild(XML_Node* node);
        void addValue(std::string val);
        void addValue(doublereal val, std::string fmt="%g");
        void addAttribute(std::string attrib, std::string value);
        void addAttribute(std::string attrib, double value, std::string fmt="%g");
        void writeHeader(std::ostream& s);
        std::string value() const { return m_value; }
        std::string value(std::string loc) const { return child(loc).value(); }
        void setLineNumber(int n) {m_linenum = n;}
        int lineNumber() const {return m_linenum;}
        doublereal fp_value() const { 
            return std::atof(m_value.c_str()); 
        }
        integer int_value() const { 
            return std::atoi(m_value.c_str()); 
        }
        std::string operator()() const { return m_value; }
        std::string operator()(std::string loc) const { return value(loc); }

	/**
	 * The operator[] is overloaded to provide a lookup capability
	 * on attributes for the current XML element.
	 *
	 * For example
	 *     xmlNode["id"] 
	 * will return the value of the attribute "id" for the current
	 * XML element. It will return the blank std::string if there isn't
	 * an attribute with that name.
	 */
        std::string operator[](std::string attr) const {
            return attrib(attr);
        }

        /**
         * This function searches the attibutes vector for the parameter 
         * std::string attribute. If a match is found, the attribute value
         * is returned as a string. If no match is found, the empty string
         * is returned.
         *
         * @param(attr) Std::String containing the attribute to be searched for.
         */
        std::string attrib(std::string attr) const { 
            std::map<std::string,std::string>::const_iterator i = m_attribs.find(attr);
            if (i != m_attribs.end()) return i->second;
            return ""; 
        }

        std::map<std::string,std::string>& attribs() { return m_attribs; }
        XML_Node* parent() const { return m_parent; }
        XML_Node* setParent(XML_Node* p) { m_parent = p; return p; }

        bool hasChild(std::string ch) const {
            return (m_childindex.find(ch) != m_childindex.end());
        }
        bool hasAttrib(std::string a) const {
            return (m_attribs.find(a) != m_attribs.end());
        }

      //! Returns the name of the XML node
      /*!
       * The name is the XML node is the XML node name
       */
      std::string name() const { return m_name; }
        std::string id() const {
            if (hasAttrib("id")) return attrib("id");
            else return "";
        }
        int number() const { return m_n; }

        XML_Node& child(int n) const { return *m_children[n]; }
        std::vector<XML_Node*>& children()  { return m_children; }
        const std::vector<XML_Node*>& children() const { return m_children; }
        int nChildren() const { return m_nchildren; }

        void build(std::istream& f);

        void _require(std::string a, std::string v) const;

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
	XML_Node* findNameID(const std::string &nameTarget, 
			     const std::string &idTarget) const;

        XML_Node* findID(const std::string& id, int depth=100) const;
        XML_Node* findByAttr(const std::string& attr, const std::string& val);
        XML_Node* findByName(const std::string& nm);
        void getChildren(std::string name, std::vector<XML_Node*>& children) const;
        XML_Node& child(std::string loc) const;
        void write(std::ostream& s, int level = 0) const;
        XML_Node& root() const { return *m_root; }
        void setRoot(XML_Node& root) { m_root = &root; }
        void copyUnion(XML_Node *node_dest) const;
        void copy(XML_Node *node_dest) const;
        void lock() {m_locked = true;}
        void unlock() {m_locked = false; }

    private:
	void write_int(std::ostream& s, int level = 0) const;

    protected:

      //! XML node name of the node. 
      /*!
       *  For example, if we were in the XML_Node where
       *       <phase dim="3" id="gas">
       *       </phase>
       *  Then, this string would be equal to "phase"
       */
      std::string m_name;
        std::string m_value;
        std::map<std::string, XML_Node*> m_childindex;
        std::map<std::string, std::string> m_attribs;
        XML_Node* m_parent;
        XML_Node* m_root;
        bool m_locked;
        std::vector<XML_Node*> m_children;
        int m_nchildren;
        int m_n;
        bool m_iscomment;
        int m_linenum;
    };

  
    XML_Node * findXMLPhase(XML_Node* root, const std::string &id);

}

#endif

