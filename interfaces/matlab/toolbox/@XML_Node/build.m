function x = build(x, file, pre)
% BUILD  Build an XML_Node in memory from an input file.
% x = build(x, file, pre)
% :param x:
%     Instance of class :mat:func:`XML_Node`
% :param file:
%     String input file name.
% :param pre:
%     Determine the method of building. If not specified or
%     less than zero, use XML_Node::build. Otherwise, use
%     XML_Node::get_XML_File.
% :return:
%     Instance of class :mat:func:`XML_Node`
%

if nargin < 2 || ~isa(file, 'char')
    error('Syntax error. Type "help build" for more information.')
end

if nargin == 3 && pre > 0
    ctmethods(10, 15, x.id, file)
else
    ctmethods(10, 4, x.id, file)
end
