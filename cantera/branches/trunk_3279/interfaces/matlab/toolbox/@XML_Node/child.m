function v = child(x, loc)
% CHILD  Get the child of an XML_Node instance.
% v = child(x, loc)
% :param x:
%     Instance of class :mat:func:`XML_Node`
% :param loc:
%     String loc to search for.
% :return:
%     Instance of class :mat:func:`XML_Node`
%

id = ctmethods(10, 6, x.id, loc);
v = XML_Node('', '', id);
