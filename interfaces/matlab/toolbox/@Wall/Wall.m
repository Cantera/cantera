function x = Wall(left, right, area, k, u, q, v)
% WALL  Wall class constructor.
% x = Wall(left, right, area, k, u, q, v)
% A Wall separates two reactors, or a reactor and a reservoir. A wall has a
% finite area, may conduct heat between the two reactors on either
% side, and may move like a piston.
%
% Walls are stateless objects in Cantera, meaning that no differential
% equation is integrated to determine any wall property. Since it is the wall
% (piston) velocity that enters the energy equation, this means that it is
% the velocity, not the acceleration or displacement, that is specified.
% The wall velocity is computed from
%
% .. math:: v = K(P_{\rm left} - P_{\rm right}) + v_0(t),
%
% where :math:`K` is a non-negative constant, and :math:`v_0(t)` is a
% specified function of time. The velocity is positive if the wall is
% moving to the right.
%
% The heat flux through the wall is computed from
%
% .. math::  q = U(T_{\rm left} - T_{\rm right}) + q_0(t),
%
% where :math:`U` is the overall heat transfer coefficient for
% conduction/convection. The function
% :math:`q_0(t)` is a specified function of time. The heat flux is positive
% when heat flows from the reactor on the left to the reactor on the right.
%
% Note: all of the arguments are optional and can be activated after initial
% construction by using the various methods of the :mat:func:`Wall` class.
% Any improperly specified arguments will generate warnings; these can be ignored
% if the intention was to not use a particular argument. Thus, the velocity of
% the wall can be set by using empty strings or 0.0 for each of the arguments before
% the velocity with no harm.
%
% :param left:
%     Instance of class :mat:func:`Reactor` to be used as the bulk phase on
%     the left side of the wall. See :mat:func:`install`
% :param right:
%     Instance of class :mat:func:`Reactor` to be used as the bulk phase on
%     the right side of the wall. See :mat:func:`install`
% :param area:
%     The area of the wall in m**2. See :mat:func:`area` and :mat:func:`setArea`.
%     Defaults to 1.0 m**2 if not specified.
% :param k:
%     Expansion rate coefficient in m/(s-Pa). See :mat:func:`setExpansionRateCoeff`
%     and :mat:func:`vdot`. Defaults to 0.0 if not specified.
% :param u:
%     Heat transfer coefficient in W/(m**2-K). See :mat:func:`setHeatTransferCoeff`
%     and :mat:func:`qdot`. Defaults to 0.0 if not specified.
% :param q:
%     Heat flux in W/m**2. Must be an instance of :mat:func:`Func`. See
%     :mat:func:`setHeatFlux` and :mat:func:`qdot`. Defaults to 0.0 if not specified.
% :param v:
%     Velocity of the wall in m/s. Must be an instance of :mat:func:`Func`. See
%     :mat:func:`setVelocity` and :mat:func:`vdot`. Defaults to 0.0 if not specified.
% :return:
%     Instance of class :mat:func:`Wall`

% At the moment, only one wall type is implemented
typ = 'Wall';

x.type = char(typ);
x.index = wallmethods(0, x.type);
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
