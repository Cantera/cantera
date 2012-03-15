function x = findByID(root, id)
% FINDBYID - Find an XML element by ID
%
index = ctmethods(10, 8, root.id, id);
x = XML_Node('','', index);
