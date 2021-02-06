/**
 * @file xml.h
 * Classes providing support for XML data files. These classes
 * implement only those aspects of XML required to read, write, and
 * manipulate CTML data files.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_XML_H
#define CT_XML_H

#include "ctexceptions.h"
#include "global.h"

//@{
#define XML_INDENT 4
//@}

namespace Cantera
{
//!  Class XML_Reader reads an XML file into an XML_Node object.
/*!
 *   Class XML_Reader is designed for internal use.
 *
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
 */
class XML_Reader
{
public:
    //! Sole Constructor for the XML_Reader class
    /*!
     *  @param input   Reference to the istream object containing the XML file
     */
    XML_Reader(std::istream& input);

    //! Read a single character from the input stream and returns it
    /*!
     * All low level reads occur through this function. The function also keeps
     * track of the line numbers.
     *
     * @param ch   Character to be returned.
     */
    void getchr(char& ch);

    //!  Searches a string for the first occurrence of a valid quoted string.
    /*!
     * Quotes can start with either a single quote or a double quote, but must
     * also end with the same type. Quotes may be commented out by preceding
     * with a backslash character, '\\'.
     *
     * @param aline     This is the input string to be searched
     * @param rstring   Return value of the string that is found.
     *                  The quotes are stripped from the string.
     * @returns the integer position just after the quoted string.
     */
    int findQuotedString(const std::string& aline, std::string& rstring) const;

    //! parseTag parses XML tags, i.e., the XML elements that are in between
    //! angle brackets.
    /*!
     * @param tag            Tag to be parsed - input
     * @param name           Output string containing name of the XML
     * @param[out] attribs   map of attribute name and attribute value
     */
    void parseTag(const std::string& tag, std::string& name,
                  std::map<std::string, std::string>& attribs) const;

    //! Reads an XML tag into a string
    /*!
     * This function advances the input streams pointer
     *
     * @param attribs   map of attribute name and attribute value - output
     * @return          Output string containing name of the XML
     */
    std::string readTag(std::map<std::string, std::string>& attribs);

    //! Return the value portion of an XML element
    /*!
     * This function advances the input streams pointer
     */
    std::string readValue();

protected:
    //! Input stream containing the XML file
    std::istream& m_s;

public:
    //! Line count
    int m_line;
};

//! Class XML_Node is a tree-based representation of the contents of an XML file
/*!
 * There are routines for adding to the tree, querying and searching the tree,
 * and for writing the tree out to an output file.
 *
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
 */
class XML_Node
{
public:
    //! Constructor for XML_Node, representing a tree structure
    /*!
     * @param nm  Name of the node.
     * @param parent   Pointer to the parent for this node in the tree.
     *                 A value of 0 indicates this is the top of the tree.
     */
    explicit XML_Node(const std::string& nm="--", XML_Node* const parent=0);

    XML_Node(const XML_Node& right);
    XML_Node& operator=(const XML_Node& right);
    virtual ~XML_Node();

    //! Clear the current node and everything under it
    /*!
     * The value, attributes and children are all zeroed. The name and the
     * parent information is kept.
     */
    void clear();

    //! Merge an existing node as a child node to the current node
    /*!
     * This will merge an XML_Node as a child to the current node. Note, this
     * actually adds the node. Therefore, the current node is changed. There is
     * no copy made of the child node. The child node should not be deleted in
     * the future
     *
     * @param node  Reference to a child XML_Node object
     * @returns a reference to the added child node
     */
    XML_Node& mergeAsChild(XML_Node& node);

    // Add a child node to the current node by making a copy of an existing node
    // tree
    /*
     * This will add an XML_Node as a child to the current node. Note, this
     * actually adds the node. Therefore, node is changed. A copy is made of the
     * underlying tree
     *
     * @param node  Reference to a child XML_Node object
     * @returns a reference to the added node
     */
    XML_Node& addChild(const XML_Node& node);

    //! Add a child node to the current node with a specified name
    /*!
     * This will add an XML_Node as a child to the current node.
     * The node will be blank except for the specified name.
     *
     * @param sname    Name of the new child
     * @returns a reference to the added node
     */
    XML_Node& addChild(const std::string& sname);

    //! Add a child node to the current XML node, and at the same time add a
    //! value to the child
    /*!
     * Resulting XML string:
     *
     *      <name> value </name>
     *
     * @param   name       Name of the child XML_Node object
     * @param   value      Value of the XML_Node - string
     * @returns a reference to the created child XML_Node object
     */
    XML_Node& addChild(const std::string& name, const std::string& value);

    //! Add a child node to the current XML node, and at the same time add a
    //! formatted value to the child
    /*!
     * This version supplies a formatting string (printf format) to the output
     * of the value.
     *
     * Resulting XML string:
     *
     *      <name> value </name>
     *
     * @param   name       Name of the child XML_Node object
     * @param   value      Value of the XML_Node - double.
     * @param   fmt        Format of the output for value
     * @returns a reference to the created child XML_Node object
     */
    XML_Node& addChild(const std::string& name, const doublereal value,
                       const std::string& fmt="%g");

    //! Remove a child from this node's list of children
    /*!
     * This function removes an XML_Node from the children of this node.
     *
     * @param  node  Pointer to the node to be removed. Note, this node
     *               isn't modified in any way.
     */
    void removeChild(const XML_Node* const node);
    
    //! Add a child node to the current node containing a comment
    /*!
     *  Child node will have the name, comment.
     *
     *  @param comment    Content of the comment
     */
    void addComment(const std::string& comment);

    //! Modify the value for the current node
    /*!
     * This functions fills in the m_value field of the current node
     *
     * @param val  string Value that the node will be assigned
     */
    void addValue(const std::string& val);

    //! Modify the value for the current node
    /*!
     * This functions fills in the m_value field of the current node
     * with a formatted double value
     *
     * @param val  double Value that the node will be assigned
     * @param fmt  Format of the printf string conversion of the double.
     *             Default is "%g". Must be less than 63 chars
     */
    void addValue(const doublereal val, const std::string& fmt="%g");

    //! Return the value of an XML node as a string
    /*!
     * This is a simple accessor routine
     */
    std::string value() const;

    //! Return the value of an XML child node as a string
    /*!
     * @param cname  Name of the child node to the current node, for which you
     *               want the value
     */
    std::string value(const std::string& cname) const;

    //! The Overloaded parenthesis operator with one augment
    //! returns the value of an XML child node as a string
    /*!
     * @param cname  Name of the child node to the current node, for which you
     *               want the value
     */
    std::string operator()(const std::string& cname) const;

    //! Return the value of an XML node as a single double
    /*!
     * This accesses the value string, and then tries to interpret it as a
     * single double value.
     */
    doublereal fp_value() const;

    //! Return the value of an XML node as a single int
    /*!
     * This accesses the value string, and then tries to interpret it as a
     * single int value.
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
    void addAttribute(const std::string& attrib, const std::string& value);

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
    void addAttribute(const std::string& attrib, const doublereal value,
                      const std::string& fmt="%g");

    //! Add an integer attribute
    /*!
     * @param attrib   String name for the attribute to be assigned
     * @param value    int Value that the node will be assigned
     */
    void addAttribute(const std::string& attrib, int value);

    //! Add an unsigned integer attribute
    /*!
     * @param attrib   String name for the attribute to be assigned
     * @param value    int Value that the node will be assigned
     */
    void addAttribute(const std::string& attrib, size_t value);

    //! The operator[] is overloaded to provide a lookup capability
    //! on attributes for the current XML element.
    /*!
     * For example
     *     xmlNode["id"]
     * will return the value of the attribute "id" for the current
     * XML element. It will return the blank std::string if there isn't
     * an attribute with that name.
     *
     * @param attr  attribute string to look up
     * @returns a string representing the value of the attribute within the XML
     *          node. If there is no attribute with the given name, it returns
     *          the null string.
     */
    std::string operator[](const std::string& attr) const;

    //! Function returns the value of an attribute
    /*!
     * This function searches the attributes vector for the attribute named
     * 'attr'. If a match is found, the attribute value is returned as a
     * string. If no match is found, the empty string is returned.
     *
     * @param attr  String containing the attribute to be searched for.
     * @return  If a match is found, the attribute value is returned as a
     *          string. If no match is found, the empty string is returned.
     */
    std::string attrib(const std::string& attr) const;

private:
    //! Returns a changeable value of the attributes map for the current node
    /*!
     * Note this is a simple accessor routine. And, it is a private function.
     * It's used in some internal copy and assignment routines
     */
    std::map<std::string,std::string>& attribs();

public:
    //! Returns an unchangeable value of the attributes map for the current node
    const std::map<std::string,std::string>& attribsConst() const;

    //! Set the line number
    /*!
     *  @param n   the member data m_linenum is set to n
     */
    void setLineNumber(const int n);

    //! Return the line number
    int lineNumber() const;

    //! Returns a pointer to the parent node of the current node
    XML_Node* parent() const;

    //! Sets the pointer for the parent node of the current node
    /*!
     * @param p Pointer to the parent node
     * @returns the pointer p
     */
    XML_Node* setParent(XML_Node* const p);

    //! Tests whether the current node has a child node with a particular name
    /*!
     * @param ch  Name of the child node to test
     * @returns true if the child node exists, false otherwise.
     */
    bool hasChild(const std::string& ch) const;

    //! Tests whether the current node has an attribute with a particular name
    /*!
     * @param a  Name of the attribute to test
     * @returns true if the attribute exists, false otherwise.
     */
    bool hasAttrib(const std::string& a) const;

    //! Returns the name of the XML node
    std::string name() const {
        return m_name;
    }

    //! Sets the name of the XML node
    /*!
     * @param name_ The name of the XML node
     */
    void setName(const std::string& name_) {
        m_name = name_;
    }

    //! Return the id attribute, if present
    /*!
     * Returns the id attribute if present. If not it return the empty string
     */
    std::string id() const;

    //! Return a changeable reference to the n'th child of the current node
    /*!
     *  @param n  Number of the child to return
     */
    XML_Node& child(const size_t n) const;

    //! Return an unchangeable reference to the vector of children of the current node
    /*!
     * Each of the individual XML_Node child pointers, however, is pointing to a
     * changeable XML node object.
     */
    const std::vector<XML_Node*>& children() const;

    //! Return the number of children
    /*!
     * @param discardComments If true comments are discarded when adding up the
     *                        number of children. Defaults to false.
     */
    size_t nChildren(bool discardComments = false) const;

    //!  Boolean function indicating whether a comment
    bool isComment() const;

    //! Require that the current XML node has an attribute named by the first
    //! argument, a, and that this attribute has the string value listed in
    //! the second argument, v.
    /*!
     * @param a  attribute name
     * @param v  required value of the attribute
     *
     * If the condition is not true, an exception is thrown
     */
    void _require(const std::string& a, const std::string& v) const;

    //! This routine carries out a recursive search for an XML node based
    //! on both the XML element name and the attribute ID.
    /*!
     * If exact matches are found for both fields, the pointer
     * to the matching XML Node is returned.
     *
     * The ID attribute may be defaulted by setting it to "". In this case the
     * pointer to the first XML element matching the name only is returned.
     *
     * @param nameTarget  Name of the XML Node that is being searched for
     * @param idTarget    "id" attribute of the XML Node that the routine
     *                    looks for
     * @returns the pointer to the XML node that fits the criteria
     *
     * @internal
     * This algorithm does a lateral search of first generation children
     * first before diving deeper into each tree branch.
     */
    XML_Node* findNameID(const std::string& nameTarget,
                         const std::string& idTarget) const;

    //! This routine carries out a search for an XML node based on the XML
    //! element name, the attribute ID and an integer index.
    /*!
     * If exact matches are found for all fields, the pointer
     * to the matching XML Node is returned. The search is only carried out on
     * the current element and the child elements of the current element.
     *
     * The "id" attribute may be defaulted by setting it to "". In this case the
     * pointer to the first XML element matching the name and the Index is returned.
     *
     * @param nameTarget  Name of the XML Node that is being searched for
     * @param idTarget    "id" attribute of the XML Node that the routine
     *                    looks for
     * @param index       Integer describing the index. The index is an
     *                    attribute of the form index = "3"
     * @returns the pointer to the XML node that fits the criteria
     */
    XML_Node* findNameIDIndex(const std::string& nameTarget,
                              const std::string& idTarget, const int index) const;

    //! This routine carries out a recursive search for an XML node based
    //! on the XML element attribute "id"
    /*!
     * If exact match is found, the pointer to the matching XML Node is
     * returned. If not, 0 is returned.
     *
     * @param id       "id" attribute of the XML Node that the routine looks for
     * @param depth    Depth of the search.
     * @returns the pointer to the XML node that fits the criteria
     *
     * @internal
     * This algorithm does a lateral search of first generation children
     * first before diving deeper into each tree branch.
     */
    XML_Node* findID(const std::string& id, const int depth=100) const;

    //! This routine carries out a recursive search for an XML node based
    //! on an attribute of each XML node
    /*!
     * If exact match is found with respect to the attribute name and value of
     * the attribute, the pointer to the matching XML Node is returned. If
     * not, 0 is returned.
     *
     * @param attr     Attribute of the XML Node that the routine looks for
     * @param val      Value of the attribute
     * @param depth    Depth of the search. A value of 1 means that only the
     *                 immediate children are searched.
     * @returns the pointer to the XML node that fits the criteria
     */
    XML_Node* findByAttr(const std::string& attr, const std::string& val,
                         int depth = 100000) const;

    //! This routine carries out a recursive search for an XML node based
    //! on the name of the node.
    /*!
     * If exact match is found with respect to XML_Node name, the pointer to the
     * matching XML Node is returned. If not, 0 is returned. This is the const
     * version of the routine.
     *
     * @param nm       Name of the XML node
     * @param depth    Depth of the search. A value of 1 means that only the
     *                 immediate children are searched.
     * @returns the pointer to the XML node that fits the criteria
     */
    const XML_Node* findByName(const std::string& nm, int depth = 100000) const;

    //! This routine carries out a recursive search for an XML node based
    //! on the name of the node.
    /*!
     * If exact match is found with respect to XML_Node name, the pointer
     * to the matching XML Node is returned. If not, 0 is returned.
     * This is the non-const version of the routine.
     *
     *  @param nm       Name of the XML node
     *  @param depth    Depth of the search. A value of 1 means that only the
     *                  immediate children are searched.
     *  @returns the pointer to the XML node that fits the criteria
     */
    XML_Node* findByName(const std::string& nm, int depth = 100000);

    //! Get a vector of pointers to XML_Node containing all of the children
    //! of the current node which match the given name
    /*!
     *  @param name   Name of the XML_Node children to search for
     *  @return vector of pointers to child XML_Nodes with the matching name
     */
    std::vector<XML_Node*> getChildren(const std::string& name) const;

    //! Return a changeable reference to a child of the current node, named by
    //! the argument
    /*!
     * Note the underlying data allows for more than one XML element with the
     * same name. This routine returns the first child with the given name.
     *
     * @param loc  Name of the child to return
     */
    XML_Node& child(const std::string& loc) const;

    //! Write the header to the XML file to the specified ostream
    /*!
     * @param s   ostream to write the output to
     */
    void writeHeader(std::ostream& s);

    //! Write an XML subtree to an output stream.
    /*!
     * This is a wrapper around the static routine write_int(). All this does
     * is add an endl on to the output stream. write_int() is fine, but the
     * last endl wasn't being written.
     *
     * @param s       ostream to write to
     * @param level   Indentation level to work from
     * @param numRecursivesAllowed Number of recursive calls allowed
     */
    void write(std::ostream& s, const int level = 0, int numRecursivesAllowed = 60000) const;

    //! Return the root of the current XML_Node tree
    /*!
     * Returns a reference to the root of the current XML tree
     */
    XML_Node& root() const;

    //! Set the root XML_Node value within the current node
    /*!
     * @param root Value of the root XML_Node.
     */
    void setRoot(const XML_Node& root);

    //! Populate the XML tree from an input file
    void build(const std::string& filename);

    //! Main routine to create an tree-like representation of an XML file
    /*!
     * Given an input stream, this routine will read matched XML tags
     * representing the ctml file until an EOF is read from the file. This
     * routine is called by the root XML_Node object.
     *
     * @param f   Input stream containing the ascii input file
     * @param filename Name of the input file, used in error messages
     */
    void build(std::istream& f, const std::string& filename="[unknown]");

    //! Copy all of the information in the current XML_Node tree into the
    //! destination XML_Node tree, doing a union operation as we go
    /*!
     * Note this is a const function because the current XML_Node and its
     * children isn't altered by this operation. copyUnion() doesn't duplicate
     * existing entries in the destination XML_Node tree.
     *
     * @param node_dest  This is the XML node to receive the information
     */
    void copyUnion(XML_Node* const node_dest) const;

    //! Copy all of the information in the current XML_Node tree into the
    //! destination XML_Node tree, doing a complete copy as we go.
    /*!
     * Note this is a const function because the current XML_Node and its
     * children isn't altered by this operation.
     *
     * @param node_dest  This is the XML node to receive the information
     */
    void copy(XML_Node* const node_dest) const;

    //! Set the lock for this node and all of its children
    void lock();

    //! Unset the lock for this node and all of its children
    void unlock();

private:
    //! Write an XML subtree to an output stream.
    /*!
     * This is the main recursive routine. It doesn't put a final endl on. This
     * is fixed up in the public method. A method to only write out a limited
     * amount of the XML tree has been added.
     *
     * @param s       ostream to write to
     * @param level   Indentation level to work from
     * @param numRecursivesAllowed Number of recursive calls allowed
     */
    void write_int(std::ostream& s, int level = 0, int numRecursivesAllowed = 60000) const;

protected:
    //! XML node name of the node.
    /*!
     * For example, if we were in the XML_Node where
     *
     *      <phase dim="3" id="gas">
     *      </phase>
     *
     * Then, this string would be equal to "phase". "dim" and "id" are
     * attributes of the XML_Node.
     */
    std::string m_name;

    //! Value of the XML node
    /*!
     * This is the string contents of the XML node. For example. The XML node
     * named eps:
     *
     *      <eps>
     *         valueString
     *      </eps>
     *
     * has a m_value string containing "valueString".
     */
    std::string m_value;

    //! Name of the file from which this XML node was read. Only populated for
    //! the root node.
    std::string m_filename;

    //! Map containing an index between the node name and the
    //! pointer to the node
    /*!
     * m_childindex[node.name()] = XML_Node *pointer
     *
     * This object helps to speed up searches.
     */
    std::multimap<std::string, XML_Node*> m_childindex;

    //! Storage of attributes for a node
    /*!
     *  m_attribs[attribName] = attribValue
     */
    std::map<std::string, std::string> m_attribs;

    //! Pointer to the parent XML_Node for the current node
    /*!
     *  Note, the top node has a parent value of 0
     */
    XML_Node* m_parent;

    //! Pointer to the root XML_Node for the current node
    /*!
     *  Note, the top node has a root value equal to itself
     */
    XML_Node* m_root;

    //! Lock for this node
    /*!
     * Currently, unimplemented functionality. If locked,
     * it means you can't delete this node.
     */
    bool m_locked;

    //! Vector of pointers to child nodes
    std::vector<XML_Node*> m_children;

    //! True if the current node is a comment node
    bool m_iscomment;

    //! The member data m_linenum
    /*!
     *  Currently, unimplemented functionality
     */
    int m_linenum;
};

//! Search an XML_Node tree for a named phase XML_Node
/*!
 *  Search for a phase Node matching an id.
 *
 *  @param root         Starting XML_Node* pointer for the search
 *  @param phaseId      id of the phase to search for
 *  @returns the XML_Node pointer if the phase is found. If the phase is not
 *           found, it returns 0
 *
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
 */
XML_Node* findXMLPhase(XML_Node* root, const std::string& phaseId);

}

#endif
