function x = Wall(left, right, area, k, u, q, v, kleft, kright)
%

% This is a dummy argument, it is not actually used by wall_new in ctreactor.cpp
typ = 1;

x.index = wallmethods(0, typ);

if x.index < 0
    error(geterr);
end
x.left = -1;
x.right = -1;
x = class(x, 'Wall');

if nargin >= 2
    if isa(left, 'Reactor') && isa(right, 'Reactor')
        install(x, left, right);
    else
        warning(['left and/or right were not instances of Reactor, ' ...
                 'and were not installed.'])
    end
end

if nargin >= 3
    if isnumeric(area)
        setArea(x, area);
    else
        warning('area was not a number and the area was not set')
    end
end

if nargin >= 4
    if isnumeric(k)
        setExpansionRateCoeff(x, k);
    else
        warning('k was not a number and the expansion rate coefficient was not set')
    end
end

if nargin >= 5
    if isnumeric(u)
        setHeatTransferCoeff(x, u);
    else
        warning('u was not a number and the expansion rate coefficient was not set')
    end
end

if nargin >= 6
    if isa(q, 'Func')
        setHeatFlux(x, q);
    else
        warning('q was not an instance of Func and was not set')
    end
end

if nargin >= 7
    if isa(v, 'Func')
        setVelocity(x, v)
    else
        warning('v was not an instance of Func and was not set')
    end
end

if nargin >= 8
    if ~isa(kleft, 'Kinetics')
        kleft = 0;
    end
    if nargin == 9
        if ~isa(kright, 'Kinetics')
            kright = 0;
        end
    else
        kright = 0;
    end
    if ~isa(kleft, 'Kinetics') && ~isa(kright, 'Kinetics')
        warning(['kleft and kright were not instances of class Kinetics ' ...
                 'and so were not set'])
    else
        setKinetics(x, kleft, kright);
    end
end
