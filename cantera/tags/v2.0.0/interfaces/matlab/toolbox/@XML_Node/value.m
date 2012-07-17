function v = value(x,loc)
% VALUE - value(x) returns the value of the XML element.
%         value(x, loc) is shorthand for value(child(x,loc))
%
if nargin == 2
    c = child(x, loc);
    id = c.id;
else
    id = x.id;
end

v = ctmethods(10, 21, id);
