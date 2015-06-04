function n = addChild(root, name, val)
% ADDCHILD  Add a child to the root.
% n = addChild(root, name, val)
% :param root:
%     Instance of class :mat:func:`XML_Node`
% :param name:
%     String ID of the child to be added.
% :param val:
%     String value to be added to the child.
% :return:
%     Instance of class :mat:func:`XML_Node`
%

n = ctmethods(10, 11, root.id, name, val);
