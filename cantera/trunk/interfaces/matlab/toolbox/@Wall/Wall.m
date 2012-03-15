function x = Wall(typ)
%
if nargin == 0
    typ = 1;
end
x.index = wallmethods(0,typ);
if x.index < 0
    error(geterr);
end
x.left = -1;
x.right = -1;
x = class(x,'Wall');
