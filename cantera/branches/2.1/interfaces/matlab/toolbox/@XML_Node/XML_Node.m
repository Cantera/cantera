function x = XML_Node(name, src, wrap)
%
% XML_Node - Cantera XML_Node class constructor
%

x.id = 0;
if nargin == 3
    x.id = wrap;
elseif nargin == 2
    % read tree from a file
    x.id = ctmethods(10,15,0,src); % newxml(name)
    if x.id < 0
        error(geterr);
    end
elseif nargin == 1
    x.id = ctmethods(10,0,0,name);
end

x = class(x,'XML_Node');
