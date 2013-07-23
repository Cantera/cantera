function x = findByName(root, name)
% FINDBYNAME - Find an XML element by name
%
index = ctmethods(10, 9, root.id, name);
x = XML_Node('','', index);
