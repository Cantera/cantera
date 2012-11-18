function a = setState(a,job,values)
disp('deprecated')
if nargin ~= 3 || ~isa(job,'char')
    error('Syntax error. Type "help setState" for more information.')
end

switch job
    case 'T'
        setTemperature(a,values)
    otherwise
        error(['unknown flag: ' job])
end
