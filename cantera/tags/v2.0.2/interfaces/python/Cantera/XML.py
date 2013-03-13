"""
This module provides the Python interface to C++ class XML_Node.
"""

import _cantera
import types
import tempfile
import string
import exceptions

class XML_Node:
    """A node in an XML tree."""
    def __init__(self, name="--", src="", wrap=0, root=None, preprocess=0, debug=0):
        """
        Return an instance representing a node in an XML tree.

        If 'src' is specified, then the XML tree found in file 'src' is
        constructed, and this node forms the root of the tree. The XML tree
        is saved, and a second call with the same value for 'src' will use
        the XML tree already read in, instead of reading it in again.

        If 'wrap' is greater than zero, then only a Python wrapper is
        created - no new kernel object results.
        """
        self._xml_id = 0
        self.wrap = wrap

        # create a wrapper for an existing kernel object
        if wrap > 0:
            self._xml_id = wrap

        # create an XML tree by parsing a file, and possibly
        # preprocessing it first
        elif src:
            self._xml_id = _cantera.xml_get_XML_File(src, debug)
            self.wrap = 1  # disable deleting

        # create a new empty node
        else:
            self._xml_id = _cantera.xml_new(name)

    def __del__(self):
        """Delete the node. Does nothing if this node is only a wrapper."""
        if not self.wrap:
            _cantera.xml_del(self._xml_id)

    def tag(self):
        return _cantera.xml_tag(self._xml_id)

    def id(self):
        """Return the id attribute if one exists, or else the empty string."""
        try:
            return self['id']
        except:
            return ''

    def nChildren(self):
        """Number of child elements."""
        return _cantera.xml_nChildren(self._xml_id)

    def children(self,tag=""):
        """Return a list of all child elements, or just those with a specified
        tag name.
        """
        nch = self.nChildren()
        children = []
        for n in range(nch):
            m = _cantera.xml_childbynumber(self._xml_id, n)
            ch = XML_Node(wrap = m)
            if (tag == "" or ch.tag() == tag):
                children.append(ch)
        return children

    def removeChild(self, child):
        """Remove a child and all its descendants."""
        _cantera.xml_removeChild(self._xml_id, child._xml_id)

    def addChild(self, name, value=""):
        """Add a child with tag 'name', and set its value if the value
        parameter is supplied."""
        if type(value) <> types.StringType:
            v = `value`
        else:
            v = value
        m = _cantera.xml_addChild(self._xml_id, name, v)
        return XML_Node(wrap = m)

    def hasAttrib(self, key):
        x = self.attrib(key)
        if x: return 1
        else: return 0

    def attrib(self, key):
        """Return attribute 'key', or the empty string if this attribute
        does not exist."""
        try:
            return _cantera.xml_attrib(self._xml_id, key)
        except:
            return ''

    def addAttrib(self, key, value):
        """Add attribute 'key' with value 'value'."""
        _cantera.xml_addAttrib(self._xml_id, key, value)

    def addComment(self, comment):
        """Add a comment."""
        _cantera.xml_addComment(self._xml_id, comment)

    def value(self, loc=""):
        """Return the value of this node, or, if
        the loc argument is supplied, of the node with relative
        address 'loc'."""
        if loc:
            node = self.child(loc)
            return node.value()
        else:
            return _cantera.xml_value(self._xml_id)

    def child(self, loc="", id="", name=""):
        if loc:
            m = _cantera.xml_child(self._xml_id, loc)
        elif id:
            m = _cantera.xml_findID(self._xml_id, id)
        elif name:
            m = _cantera.xml_findByName(self._xml_id, name)

        ch = XML_Node(wrap=m)
        return ch

    def __getitem__(self, key):
        """Get an attribute using the syntax node[key]"""
        return self.attrib(key)

    def __setitem__(self, key, value):
        """Set a new attribute using the syntax node[key] = value."""
        return self.addAttrib(key, value)

    def __int__(self):
        """Conversion to integer."""
        return self._xml_id

    def __call__(self, loc=''):
        """Get the value using the syntax node(loc)."""
        return self.value(loc)

    def write(self, file):
        _cantera.xml_write(self._xml_id, file)

    def __repr__(self):
        tmp = tempfile.mktemp('.xml')
        self.write(tmp)
        f = open(tmp)
        lines = f.readlines()
        f.close()
        s = ''
        for line in lines:
            s += line
        return s

def clear_XML():
    _cantera.xml_clear()


def getFloatArray(node, convert_units=0):
    sz = int(node['size'])
    return _cantera.ctml_getFloatArray(node._xml_id, convert_units, sz)
