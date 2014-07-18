function n = nChildren(root)
% NCHILDREN  Get the number of children of an XML_Node.
% n = nChildren(root)
% :param root:
%     Instance of class :mat:func:`XML_Node`
% :return:
%     Integer number of children of the input XML_Node
%

n = ctmethods(10, 10, root.id);
