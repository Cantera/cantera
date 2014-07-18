function x = XML_Node(name, src, wrap)
% XML_NODE  XML_Node class constructor
% x = XML_Node(name, src, wrap)
% :param name:
%     String name of the XML_Node that should be created.
% :param src:
%     String XML file name from which an instance of XML_Node
%     should be created. Reads the XML tree from the input file.
% :param wrap:
%     Specify the ID of the XML_Node.
% :return:
%     Instance of class :mat:func:`XML_Node`
%

x.id = 0;
if nargin == 3
    x.id = wrap;
elseif nargin == 2
    % read tree from a file
    x.id = ctmethods(10, 15, 0, src);
    if x.id < 0
        error(geterr);
    end
elseif nargin == 1
    x.id = ctmethods(10, 0, 0, name);
end

x = class(x, 'XML_Node');
