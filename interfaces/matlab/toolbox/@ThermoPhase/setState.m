function setState(a, job, values)
% SETSTATE  Set the state.
% a = setState(a,job,values)
% Deprecated in favor of :mat:func:`set`
%
% See also: :mat:func:`set`
%

warning('deprecated in favor of set.m')
if nargin ~= 3 || ~isa(job, 'char')
    error('Syntax error. Type "help setState" for more information.')
end

switch job
    case 'T'
        setTemperature(a, values)
    otherwise
        error(['unknown flag: ' job])
end
