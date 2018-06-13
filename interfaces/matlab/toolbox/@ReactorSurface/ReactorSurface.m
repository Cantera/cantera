function x = ReactorSurface(kinetics, reactor, area)
% REACTORSURFACE  ReactorSurface class constructor.
% x = ReactorSurface(kinetics, reactor, area)
%
% A surface on which heterogeneous reactions take place. The mechanism object
% (typically an instance of class :mat:func:`Interface`) must be constructed so
% that it is properly linked to the object representing the fluid in the
% reactor. The surface temperature on each side is taken to be equal to the
% temperature of the reactor.
%
% Note: all of the arguments are optional and can be activated after initial
% construction by using the various methods of the :mat:func:`ReactorSurface`
% class.
%
% :param kleft:
%     Surface reaction mechanisms for the left-facing surface. This must be an
%     instance of class :mat:func:`Kinetics`, or of a class derived from Kinetics,
%     such as :mat:func:`Interface`.
% :param reactor:
%     Instance of class :mat:func:`Reactor` to be used as the adjacent bulk
%     phase. See :mat:func:`install`
% :param area:
%     The area of the surface in m**2. See :mat:func:`area` and
%     :mat:func:`setArea`. Defaults to 1.0 m**2 if not specified.
% :return:
%     Instance of class :mat:func:`ReactorSurface`

x.index = reactorsurfacemethods(0, 0);

if x.index < 0
    error(geterr);
end
x.reactor = -1;
x = class(x, 'ReactorSurface');

if nargin >= 1
    setKinetics(x, kinetics);
end

if nargin >= 2
    if isa(reactor, 'Reactor')
        install(x, reactor);
    else
        warning(['"reactor" was not an instance of Reactor, ' ...
                 'and was not installed.'])
    end
end

if nargin >= 3
    if isnumeric(area)
        setArea(x, area);
    else
        warning('area was not a number and the area was not set')
    end
end
