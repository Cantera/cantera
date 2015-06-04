function v = value(x, loc)
% VALUE  Get the value at a location in an XML_Node
% v = value(x,loc)
% value(x) returns the value of the XML element.
% value(x, loc) is shorthand for value(child(x,loc))
%
% :param x:
%     Instance of class :mat:func:`XML_Node`
% :param loc:
% :return:
%     Instance of class :mat:func:`XML_Node`
%

if nargin == 2
    c = child(x, loc);
    id = c.id;
else
    id = x.id;
end

v = ctmethods(10, 21, id);
