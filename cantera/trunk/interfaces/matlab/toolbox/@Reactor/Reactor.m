function x = Reactor(contents, typ)
% REACTOR - Create a Reactor object.
%
%    A Reactor object simulates a perfectly-stirred reactor. It has
%    a time-dependent state, and may be coupled to other reactors
%    through flow lines or through walls that may expand or
%    contract and/or conduct heat.
%
%        r1 = Reactor           % an empty reactor
%        r2 = Reactor(gas)      % a reactor containing a gas
%
%   See also: Reservoir
%
if nargin == 0
    contents = 0;
    typ = 2;
elseif nargin == 1
    typ = 2;
elseif nargin > 2
    error('too many arguments');
end

x.index = reactormethods(0,typ);
if x.index < 0
    error(geterr);
end
x.contents = contents;
x = class(x,'Reactor');

if isa(contents,'Solution')
    insert(x, contents);
end
