/**
 * @file xml.h
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

#include "ctexceptions.h"
#include "ct_defs.h"
#include "stringUtils.h"
#include "global.h"

#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>

#define XML_INDENT 4

namespace Cantera {


  /**
   * Class XML_Reader is designed for internal use.
   *  This function reads an XML file into an XML_Node object.
   */
  class XML_Reader {
  public:

    //! Sole Constructor for the XML_Reader class
    /*!
     *   @param input   Reference to the istream object containing
     *                  the XML file
     */
    XML_Reader(std::istream& input);

    //! Input Stream containing the XML file
    std::istream& m_s;

    //! Line count
    int m_line;

    //! Read a single character from the input stream
    //! and return it
    /*!
     *  All low level reads occur through this function.
     *  The function also keeps track of the line numbers.
     *
     * @param ch   Character to be returned.
     */
    void getchr(char& ch);

    //! Returns string 'aline' stripped of leading and trailing white
    //! space.
    /*!
     *  White space is defined by the ISO C function isspace(), and
     *  includes tabs, spaces, \n. \r, \v, and \f.
     *
     * @param aline  Input line to be stripped
     *
     * @return Returns a string stripped of leading and trailing white
     *         space. 
     *
     * @todo why is this a class method?
     */
    std::string strip(const std::string& aline);

    /// Looks for a substring within 'aline' enclosed in double
    /// quotes, and returns this substring (without the quotes) if
    /// found.  If not, an empty string is returned.
    /// @todo why is this a class method?
    std::string inquotes(const std::string& aline);

    /**
     *  Searches a string for the first occurrence of a valid
     *  quoted string. Quotes can start with either a single
     *  quote or a double quote, but must also end with the same
     *  type. Quotes may be commented out by preceding with a
     *  backslash character, '\\'. 
     */
    int findQuotedString(const std::string& aline, std::string &rstring);
    
    void parseTag(std::string line, std::string& name, 
		  std::map<std::string, std::string>& attribs);
    std::string readTag(std::map<std::string, std::string>& attribs);
    std::string readValue();
  };


  //////////////////////////  XML_Node  /////////////////////////////////

  //! Class XML_Node is a tree-based representation of the contents of an XML file
  /*!
   *   Class XML_Node is a tree-based representation of the contents of an XML file.
   *
   *
   *
   *
   */
  class XML_Node {
  public:

    //! Default constructor for XML_Node, representing a tree structure
    /*!
     *  Constructor for an XML_Node, which is a node in a tree-like structure
     *  representing an XML file.
     *
     *  @param nm  Name of the node.
     *             The default name of the node is "--"
     *
     *  @param p pointer to the root for this node in the tree.
     *           The default is 0 indicating this is the top of the tree.
     */
    XML_Node(const std::string nm = "--", XML_Node * const p = 0);
   
    //! Copy constructor
    /*!
     * @param right   Object to be copied
     */
    XML_Node(const XML_Node &right);

    //! Assignment operator for XML trees
    /*!
     *  @param right    XML tree to copy
     */
    XML_Node& operator=(const XML_Node &right);

    //! Destructor for the object
    virtual ~XML_Node();

    //! Add a child node to the current node containing a comment
    /*!
     *  Child node will have the name, comment.
     *
     *  @param comment    Content of the comment
     */
    void addComment(const std::string comment);

    //! Add a child node to the current node
    /*!
     * This will add an XML_Node as a child to the current node.
     * Note, this actually adds the node. Therefore, node is changed.
     * There is no copy made of the child node.
     *
     *  @param node  Reference to a child XML_Node object
     *
     *  @return      Returns a reference to the added node
     */
    XML_Node& addChild(XML_Node& node);

    //! Add a child node to the current node with a specified name
    /*!
     * This will add an XML_Node as a child to the current node.
     * The node will be blank except for the specified name.
     *
     *  @param sname    Name of the new child
     *
     *  @return         Returns a reference to the added node
     */
    XML_Node& addChild(const std::string sname);

    //!    Add a child node to the current xml node, and at the
    //!    same time add a value to the child
    /*!
     *    Resulting XML string:
     *      <name>value</name>
     *
     *   @param   name       Name of the child XML_Node object
     *   @param   value      Value of the XML_Node - string
     *   @return  Returns a reference to the created child XML_Node object
     */
    XML_Node& addChild(const std::string name, const std::string value);

    //!    Add a child node to the current xml node, and at the
    //!    same time add a formatted value to the child
    /*!
     *  This version supplies a formatting string (printf format)
     *  to the output of the value.
     *
     *    Resulting XML string:
     *      <name>value</name>
     *
     *   @param   name       Name of the child XML_Node object
     *   @param   value      Value of the XML_Node - double.
     *   @param   fmt        Format of the output for value
     *
     *   @return  Returns a reference to the created child XML_Node object
     */
    XML_Node& addChild(const std::string name, const doublereal value, 
		       const std::string fmt="%g");

    //! Remove a child from this node's list of children
    /*!
     *  This function removes an XML_Node from the children of this node.
     *
     * @param  node  Pointer to the node to be removed. Note, this node
     *               isn't modified in any way.
     */
    void removeChild(const XML_Node * const node);

    //! Modify the value for the current node
    /*!
     * This functions fills in the m_value field of the current node
     *
     * @param val  string Value that the node will be assigned
     */
    void addValue(const std::string val);

    //! Modify the value for the current node
    /*!
     * This functions fills in the m_value field of the current node
     * with a formatted double value
     *
     * @param val  double Value that the node will be assigned
     * @param fmt  Format of the printf string conversion of the double.
     *             Default is "%g". Must be less than 63 chars
     */
    void addValue(const doublereal val, const std::string fmt="%g");

    //! Return the value of an XML node as a string
    /*!
     *  This is a simple accessor routine
     */
    std::string value() const;


    //! Overloaded parenthesis operator returns the value of the Node
    /*!
     *  @return  Returns the value of the node as a string.
     */
    std::string operator()() const;

    //!  Return the value of an XML child node as a string
    /*!
     *  @param cname  Name of the child node to the current
     *                node, for which you want the value
     */
    std::string value(const std::string cname) const;

    //!  Overloaded parenthesis operator with one augment 
    //!  returns the value of an XML child node as a string
    /*!
     *  @param cname  Name of the child node to the current
     *                node, for which you want the value
     */
    std::string operator()(std::string loc) const;

    //! Return the value of an XML node as a single double
    /*!
     *  This accesses the value string, and then tries to 
     *  interpret it as a single double value.
     */
    doublereal fp_value() const;

    //! Return the value of an XML node as a single int
    /*!
     *  This accesses the value string, and then tries to 
     *  interpret it as a single int value.
     */
    integer int_value() const;

    //! Add or modify an attribute of the current node
    /*!
     * This functions fills in the m_value field of the current node
     * with a string value
     *
     * @param attrib  String name for the attribute to be assigned
     * @param value   String value that the attribute will have
     */
    void addAttribute(const std::string & attrib, const std::string & value);

    //! Add or modify an attribute to the double, value
    /*!
     * This functions fills in the attribute field, named attrib,
     * with the double value, value. A formatting string is used.
     *
     * @param attrib  String name for the attribute to be assigned
     * @param value   double Value that the node will be assigned
     * @param fmt     Format of the printf string conversion of the double.
     *                Default is "%g". 
     */
    void addAttribute(const std::string & attrib, const double value, 
		      const std::string fmt="%g");
	
    //! The operator[] is overloaded to provide a lookup capability
    //!  on attributes for the current XML element.
    /*!
     * For example
     *     xmlNode["id"] 
     * will return the value of the attribute "id" for the current
     * XML element. It will return the blank std::string if there isn't
     * an attribute with that name.
     *
     * @param attr  attribute string to look up
     * 
     * @return  Returns a string representing the value of the attribute
     *          within the XML node. If there is no attribute
     *          with the given name, it returns the null string.
     */
    std::string operator[](const std::string & attr) const;

    //! Function returns the value of an attribute
    /*!
     * This function searches the attibutes vector for the parameter 
     * std::string attribute. If a match is found, the attribute value
     * is returned as a string. If no match is found, the empty string
     * is returned.
     *
     * @param attr  Std::String containing the attribute to be searched for.
     *
     * @return Returns  If a match is found, the attribute value
     *                  is returned as a string. If no match is found, the empty string
     *                  is returned.
     */
    std::string attrib(const std::string & attr) const;

  private:
    //! Returns a changeable value of the attributes map for the current node
    /*!
     *  Note this is a simple accessor routine. And, it is a private function.
     *  It's used in some internal copy and assignment routines
     */
    std::map<std::string,std::string>& attribs();

  public:
    //! Set the line number 
    /*!
     *  @param n   the member data m_linenum is set to n
     */
    void setLineNumber(const int n);

    //! Return the line number 
    /*!
     *  @return  returns the member data m_linenum
     */
    int lineNumber() const;

 
    //! Returns a pointer to the parent node of the current node
    XML_Node* parent() const;

    //! Sets the pointer for the parent node of the current node
    /*!
     * @param p Pointer to the parent node
     *
     * @return  Returns the pointer p
     */
    XML_Node* setParent(XML_Node * const p);

    //! Tests whether the current node has a child node with a particular name
    /*!
     * @param ch  Name of the child node to test
     *
     * @return Returns true if the child node exists, false otherwise.
     */
    bool hasChild(const std::string ch) const {
      return (m_childindex.find(ch) != m_childindex.end());
    }

    //! Tests whether the current node has an attribute with a particular name
    /*!
     * @param a  Name of the attribute to test
     *
     * @return Returns true if the attribute exists, false otherwise.
     */
    bool hasAttrib(std::string a) const {
      return (m_attribs.find(a) != m_attribs.end());
    }

    //! Returns the name of the XML node
    /*!
     * The name is the XML node is the XML node name
     */
    std::string name() const { return m_name; }

    //! Return the id attribute, if present
    /*!
     * Returns the id attribute if present. If not
     * it return the empty string
     */
    std::string id() const;

    //! Return a changeable reference to the n'th child of the current node
    /*!
     *  @param n  Number of the child to return
     */
    XML_Node& child(const int n) const ;

    //! Return an unchangeable reference to the vector of children of the current node
    /*!
     * Each of the individual XML_Node child pointers, however,
     *  is to a changeable xml node object.
     *
     *  @param n  Number of the child to return
     */
    const std::vector<XML_Node*>& children() const;

    //! return the number of children
    /*!
     *
     */
    int nChildren() const;

    //!    Require that the current xml node have an attribute named
    //!    by the first argument, a, and that this attribute have the
    //!    the string value listed in the second argument, v.
    /*!
     *   @param a  attribute name
     *   @param v  required value of the attribute
     *
     *  If the condition is not true, an exception is thrown
     */
    void _require(const std::string a, const std::string v) const;

    //! This routine carries out a recursive search for an XML node based
    //! on both the xml element name and the attribute ID.
    /*!
     * If exact matches are found for both fields, the pointer
     * to the matching XML Node is returned.
     * 
     * The ID attribute may be defaulted by setting it to "".
     * In this case the pointer to the first xml element matching the name
     * only is returned.
     *
     *  @param nameTarget  Name of the XML Node that is being searched for
     *  @param idTarget    "id" attribute of the XML Node that the routine
     *                     looks for
     *
     *  @return   Returns the pointer to the XML node that fits the criteria
     *
     * @internal
     * This algorithm does a lateral search of first generation children
     * first before diving deeper into each tree branch.
     */
    XML_Node* findNameID(const std::string &nameTarget, 
			 const std::string &idTarget) const;

    //! This routine carries out a recursive search for an XML node based
    //! on the xml element attribute ID.
    /*!
     * If exact match is found, the pointer
     * to the matching XML Node is returned. If not, 0 is returned.
     * 
     * The ID attribute may be defaulted by setting it to "".
     * In this case the pointer to the first xml element matching the name
     * only is returned.
     *
     *  @param id       "id" attribute of the XML Node that the routine
     *                  looks for
     *  @param depth    Depth of the search.
     *
     *  @return         Returns the pointer to the XML node that fits the criteria
     *
     * @internal
     * This algorithm does a lateral search of first generation children
     * first before diving deeper into each tree branch.
     */
    XML_Node* findID(const std::string& id, const int depth=100) const;

    XML_Node* findByAttr(const std::string& attr, const std::string& val) const;
    const XML_Node* findByName(const std::string& nm) const;
    XML_Node* findByName(const std::string& nm);
    void getChildren(std::string name, std::vector<XML_Node*>& children) const;

    //! Return a changeable reference to a child of the current node, 
    //! named by the argument
    /*!
     *  @param loc  Name of the child to return
     */
    XML_Node& child(std::string loc) const;

    void writeHeader(std::ostream& s);
    void write(std::ostream& s, int level = 0) const;
    XML_Node& root() const { return *m_root; }
    void setRoot(XML_Node& root) { m_root = &root; }

    //! Main routine to create an tree-like representation of an XML file
    /*!
     *   Given an input stream, this routine will read matched XML tags
     *   representing the ctml file until an EOF is read from the file. 
     *   This routine is called by the root XML_Node object.
     *
     * @param f   Input stream containing the ascii input file
     */
    void build(std::istream& f);

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
     *  Then, this string would be equal to "phase". "dim" and "id"
     *  are attributes of the XML_Node.
     */
    std::string m_name;

    //! Value of the xml node
    /*!
     *  This is the string contents of the XML node. For
     *  example. The xml node named eps:
     *
     *      <eps>
     *         valueString
     *      </eps>
     *
     *  has a m_value string containing "valueString".
     */
    std::string m_value;

    std::map<std::string, XML_Node*> m_childindex;
    std::map<std::string, std::string> m_attribs;
    XML_Node* m_parent;
    XML_Node* m_root;
    bool m_locked;
    std::vector<XML_Node*> m_children;

    //! Number of children of this node
    int m_nchildren;
    bool m_iscomment;
    int m_linenum;

   
  };

  
  XML_Node * findXMLPhase(XML_Node* root, const std::string &id);

}

#endif

