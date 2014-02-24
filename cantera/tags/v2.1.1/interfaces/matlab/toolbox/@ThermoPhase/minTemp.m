function v = minTemp(p)
% MINTEMP - Minimum temperature.
%
%   The parameterizations used to represent the temperature-dependent
%   species thermodynamic properties are generally only valid in some
%   finite temperature range, which may be different for each species
%   in the phase. This method returns the lowest temperature at which
%   the parameterizations are valid for all species in the phase.
%
%   See also: maxTemp
%
v = thermo_get(p.tp_id,16);
