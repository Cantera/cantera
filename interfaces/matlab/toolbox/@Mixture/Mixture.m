function m = Mixture(phases)
% MIXTURE  Multiphase mixture class constructor.
% m = Mixture(phases)
% Class :mat:func:`Mixture` represents mixtures of one or more phases of matter.
% To construct a mixture, supply a cell array of phases and
% mole numbers::
%
%     >> gas = Solution('gas.yaml');
%     >> graphite = Solution('graphite.yaml');
%     >> mix = Mixture({gas, 1.0; graphite, 0.1});
%
% Phases may also be added later using the addPhase method::
%
%     >> water = Solution('water.yaml');
%     >> addPhase(mix, water, 3.0);
%
% Note that the objects representing each phase compute only the
% intensive state of the phase - they do not store any information
% on the amount of this phase. Mixture objects, on the other hand,
% represent the full extensive state.
%
% Mixture objects are 'lightweight' in the sense that they do not
% store parameters needed to compute thermodynamic or kinetic
% properties of the phases. These are contained in the
% ('heavyweight') phase objects. Multiple mixture objects may be
% constructed using the same set of phase objects. Each one stores
% its own state information locally, and synchronizes the phase
% objects whenever it requires phase properties.
%
% :param phases:
%     Cell array of phases and mole numbers
% :return:
%     Instance of class :mat:func:`Mixture`
%

if nargin > 1
    error('Mixture: wrong number of arguments');
end

% create an empty mixture
m.mixindex = mixturemethods(0, 0, 0);
m.phases = phases;

m = class(m, 'Mixture');

% if phases are supplied, add them
if nargin == 1
    if ~isa(phases, 'cell')
        error('Enter phases as a cell array.');
    end

    % first column contains the phase objects, and the second column
    % the mole numbers of each phase
    [np nc] = size(phases);
    if nc ~= 2
        error('Cell array of phases should have each phase on a new row');
    end
    for n = 1:np
        addPhase(m, phases{n,1}, phases{n,2});
    end
    setTemperature(m, temperature(phases{n,1}));
    setPressure(m, pressure(phases{n,1}));
end
