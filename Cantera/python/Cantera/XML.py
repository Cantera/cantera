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
    def __init__(self, name="--", src="", wrap=0, root=None):
        """
        Return an instance representing a node in an XML tree.
        If 'src' is specified, then the XML tree found in file 'src' is
        constructed, and this node forms the root of the tree.
        Construct a new XML tree, with this node as the root.
        If 'wrap' is greater than zero, then a
        """
        self.wrap = wrap
        if wrap > 0:
            self._xml_id = wrap
            self._root = root
        else:
            self._xml_id = _cantera.xml_new(name)
            if src:
                _cantera.xml_build(self._xml_id, src)
            self._root = self
        
    def __del__(self):
        """Delete the node. Does nothing if this node is only a wrapper."""
        if not self.wrap:
            _cantera.xml_del(self._xml_id)

    def tag(self):
        return _cantera.xml_tag(self._xml_id)

    def id(self):
        try:
            return self['id']
        except:
            return ''

    def root(self):
        return self._root
    
    def nChildren(self):
        return _cantera.xml_nChildren(self._xml_id)

    def children(self,tag=""):
        nch = self.nChildren()
        children = []
        for n in range(nch):
            m = _cantera.xml_childbynumber(self._xml_id, n)
            ch = XML_Node(src="", wrap=m, root=self._root)
            if (tag == "" or ch.tag() == tag):
                children.append(ch)
        return children
    
    def removeChild(self, child):
        _cantera.xml_removeChild(self._xml_id, child._xml_id)            

    def addChild(self, name, value=""):
        if type(value) <> types.StringType:
            v = `value`
        else:
            v = value
        m = _cantera.xml_addChild(self._xml_id, name, v)
        return XML_Node(src="", wrap=m, root=self._root)

    def hasAttrib(self, key):
        x = self.attrib(key)
        if x: return 1
        else: return 0
        
    def attrib(self, key):
        try:
            return _cantera.xml_attrib(self._xml_id, key)
        except:
            return ''

    def addAttrib(self, key, value):
        _cantera.xml_addAttrib(self._xml_id, key, value)
        
    def value(self, loc=""):
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
            
        ch = XML_Node(src="", wrap=m, root=self._root)
        return ch

    def __getitem__(self, key):
        return self.attrib(key)

    def __setitem__(self, key, value):
        return self.addAttrib(key, value)    

    def __int__(self):
        return self._xml_id
    
    def __call__(self, loc):
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

    def getRef(self):
        if not self["idRef"]: return self
        return find_XML(src = self["src"], root = self.root(),
                        id = self["idRef"])
    

def find_XML(src = "", root = None, id = "", loc = "", name=""):
    doc = None
    r = None
    if src:
        ihash = string.find(src,'#')
        if ihash < 0:
            fname = src
        else:
            fname, idnew = string.split(src,'#')
            if idnew: id = idnew
        if fname:
            doc = XML_Node(name="doc", src=fname)
            root = None
        elif root:
            doc = root
    elif root:
        doc = root            
    else:
        raise exceptions.CanteraError("either root or src must be specified.")
    
##    try:
    if loc or id or name:
        r = doc.child(loc=loc, id=id, name=name)
    else:
        r = doc
    return r


def getFloatArray(node, convert_units=0):
    sz = int(node['size'])
    return _cantera.ctml_getFloatArray(node._xml_id, convert_units, sz)







    
        
