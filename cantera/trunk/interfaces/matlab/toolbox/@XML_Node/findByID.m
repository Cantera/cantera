function x = findByID(root, id)
% FINDBYID  Get an XML element given its ID.
% x = findByID(root, id)
% :param root:
%     Instance of class :mat:func:`XML_Node`
% :param id:
%     String ID of the element to search for.
% :return:
%     Instance of class :mat:func:`XML_Node`
%

index = ctmethods(10, 8, root.id, id);
x = XML_Node('', '', index);
