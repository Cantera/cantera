function x = XML_Node(name, src, root, wrap)
%XML_Node Cantera XML_Node class constructor
%
x.id = 0;
x.root = 0;
if nargin == 4
   x.id = wrap;
   x.root = root;
elseif nargin > 0
   % create an empty node with name 'name'
   x.id = ctmethods(10,0,0,name); % newxml(name)
   if x.id < 0
      error(geterr);
   end
end
if nargin > 2
   x.root = root;
end

x = class(x,'XML_Node');

if nargin > 1 & nargin < 4
   build(x, src);
end
